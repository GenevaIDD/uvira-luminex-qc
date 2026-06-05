"""Split a multi-run xPONENT CSV into per-logical-plate single-plate CSVs.

Used for runs where one xPONENT file packs multiple plate reads — the
``Location`` field is ``N(plate_idx,well)`` and the Batch header lists
multiple RUN suffixes. Each logical plate becomes its own xPONENT-format
CSV that the rest of the pipeline can ingest unchanged, with any
single-well re-runs from other plate indices substituted into the
parent plate (e.g. a re-read of a dry-read well).

Specific to the 2026-06-03 Box1 partial-plate runs:

    plate_idx 1 → "Plate D-F"   (samples from Box1 rows D-F, plate B-D)
    plate_idx 2 → "Plate A-C"   (samples from Box1 rows A-C, plate B-D)
    plate_idx 3 → re-reads for A1 + H9 (parent: Plate A-C)
    plate_idx 5 → re-read for D4       (parent: Plate D-F)

Re-runs are matched by well: any re-run row whose well is missing from
the parent plate's main read is substituted in. The output Location
field is normalised to ``N(1,well)`` so the downstream parser is happy.

Also writes a synthesized ``*_inputfile.csv`` per logical plate so the
pipeline's authoritative Type-column classification picks up:

    A1-A10              → Standard
    A11-A12             → Background  (PBT-only blank)
    B/C/D *             → Unknown     (serum samples; Description = barcode)
    E/F/G *, H7         → Background  (PBS blank, "extra background")
    H8-H10              → Control     (Giardia NC)
    H11-H12             → Control     (NI7 / NI18)
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

LOC_RE = re.compile(r"^(\d+)\(([0-9]+),([A-H]\d+)\)$")

# Map of logical plate name → (own plate_idx, suffix used in Batch)
PLATES = {
    "PlateDF": {"own_idx": 1, "label": "PLATE_06032026_PlateDF"},
    "PlateAC": {"own_idx": 2, "label": "PLATE_06032026_PlateAC"},
}
RERUN_IDXS = (3, 4, 5)


def _parse_loc(loc: str) -> tuple[int, int, str] | None:
    m = LOC_RE.match(loc.strip())
    if not m:
        return None
    return int(m.group(1)), int(m.group(2)), m.group(3)


def _classify_well(well: str) -> tuple[str, str]:
    """Return (Type, default_description) for an inputfile.csv row."""
    row = well[0]
    col = int(well[1:])
    if row == "A" and 1 <= col <= 10:
        return "Standard", ""
    if row == "A" and col in (11, 12):
        return "Background", ""
    if row in ("B", "C", "D"):
        return "Unknown", ""  # barcode filled from data
    if row in ("E", "F", "G"):
        return "Background", ""  # PBS extra-blank
    if row == "H" and col == 7:
        return "Background", ""
    if row == "H" and col in (8, 9, 10):
        return "Control", ""  # Giardia NC
    if row == "H" and col in (11, 12):
        return "Control", ""  # NI7 / NI18
    return "Unknown", ""


def _find_blocks(rows: list[list[str]]) -> list[tuple[int, str]]:
    """Return list of (block_start_index, datatype_name)."""
    out = []
    for i, r in enumerate(rows):
        if len(r) >= 2 and r[0] == "DataType:":
            out.append((i, r[1]))
    return out


def _block_data_range(rows: list[list[str]], block_start: int) -> tuple[int, int]:
    """Inclusive [first_data_row, last_data_row] indices for a block."""
    header = block_start + 1
    first = header + 1
    last = first
    for i in range(first, len(rows)):
        r = rows[i]
        if not r or not r[0].strip():
            break
        last = i
    return first, last


def split(input_csv: Path, out_dir: Path) -> dict[str, Path]:
    """Split *input_csv* into per-plate CSVs under *out_dir*.

    Returns dict mapping logical plate name → output CSV path.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    text = input_csv.read_text(encoding="utf-8-sig")
    # Strip BOM, normalise line endings
    rows = list(csv.reader(text.splitlines()))

    blocks = _find_blocks(rows)
    if not blocks:
        raise ValueError("No DataType: blocks found in input CSV")

    # Header = everything before the first DataType: block
    header_rows = rows[: blocks[0][0]]

    # Find re-run wells: rows where plate_idx ∈ RERUN_IDXS, grouped by well
    # Take from the Median block (first one) — same wells present in all blocks
    median_start = next(i for (i, n) in blocks if n == "Median")
    first, last = _block_data_range(rows, median_start)
    rerun_wells: dict[str, int] = {}  # well → source plate_idx
    own_wells: dict[int, set[str]] = {1: set(), 2: set()}
    for i in range(first, last + 1):
        parsed = _parse_loc(rows[i][0])
        if not parsed:
            continue
        _, pidx, well = parsed
        if pidx in (1, 2):
            own_wells[pidx].add(well)
        elif pidx in RERUN_IDXS:
            rerun_wells[well] = pidx

    # Decide which re-run wells go to which parent
    rerun_assignment: dict[str, int] = {}  # well → parent plate_idx (1 or 2)
    for well, src in rerun_wells.items():
        missing_in = [p for p in (1, 2) if well not in own_wells[p]]
        if len(missing_in) == 1:
            rerun_assignment[well] = missing_in[0]
        elif len(missing_in) == 0:
            print(f"  warn: re-run well {well} (src plate {src}) "
                  f"is not missing from either parent; skipping", file=sys.stderr)
        else:
            print(f"  warn: re-run well {well} (src plate {src}) "
                  f"is missing from BOTH parents; skipping (ambiguous)",
                  file=sys.stderr)

    print(f"Re-run wells found: {sorted(rerun_wells.items())}")
    print(f"Re-run assignments (well → parent plate_idx): {rerun_assignment}")

    out_paths: dict[str, Path] = {}
    for pname, cfg in PLATES.items():
        own_idx = cfg["own_idx"]
        label = cfg["label"]

        # Build the header for this plate (override Batch line)
        out_header = []
        for r in header_rows:
            if r and r[0] == "Batch":
                out_header.append(["Batch", label])
            else:
                out_header.append(r)

        out_rows = list(out_header)

        for block_idx, (b_start, dtype) in enumerate(blocks):
            # Copy the DataType: marker + header row
            out_rows.append(["DataType:", dtype])
            header_row = rows[b_start + 1]
            out_rows.append(header_row)

            first, last = _block_data_range(rows, b_start)
            kept = []
            for i in range(first, last + 1):
                r = rows[i]
                parsed = _parse_loc(r[0])
                if not parsed:
                    continue
                seq, pidx, well = parsed
                if pidx == own_idx:
                    # Normalise Location to (1,well) so the parser doesn't care
                    new_r = list(r)
                    new_r[0] = f"{seq}(1,{well})"
                    kept.append((well, new_r))
                elif pidx in RERUN_IDXS and rerun_assignment.get(well) == own_idx:
                    new_r = list(r)
                    new_r[0] = f"{seq}(1,{well})"
                    kept.append((well, new_r))

            # Dedup by well (last-write-wins → re-runs supersede if any overlap)
            by_well: dict[str, list[str]] = {}
            for well, r in kept:
                by_well[well] = r
            # Sort by well row, then column
            def _wk(w: str) -> tuple[str, int]:
                return (w[0], int(w[1:]))
            for well in sorted(by_well, key=_wk):
                out_rows.append(by_well[well])

            # Blank separator row between blocks (matches xPONENT style)
            out_rows.append([])

        out_path = out_dir / f"PlateRunResults_{label}.csv"
        with out_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, quoting=csv.QUOTE_ALL)
            for r in out_rows:
                w.writerow(r)
        out_paths[pname] = out_path
        print(f"  wrote {out_path}  ({sum(1 for r in out_rows if r and _parse_loc(r[0])):d} data rows)")

        # Write per-plate inputfile.csv. Pulls sample names from the Median
        # block of the *output* plate so Description = sample_name for B/C/D
        # serum wells (= barcode).
        inputfile_path = out_dir / f"{label}_inputfile.csv"
        with inputfile_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, quoting=csv.QUOTE_MINIMAL)
            w.writerow(["Location", "Type", "Description"])
            # Walk own + reassigned re-run wells; same well list as the
            # Median block we just wrote.
            median_out_start = None
            for j, r in enumerate(out_rows):
                if r and r[0] == "DataType:" and len(r) > 1 and r[1] == "Median":
                    median_out_start = j
                    break
            assert median_out_start is not None
            first_o = median_out_start + 2
            for k in range(first_o, len(out_rows)):
                r = out_rows[k]
                if not r:
                    break
                loc = r[0]
                sample = r[1] if len(r) > 1 else ""
                parsed = _parse_loc(loc)
                if not parsed:
                    continue
                _, _, well = parsed
                typ, _ = _classify_well(well)
                # Description: for Unknown (B/C/D serum), use barcode from
                # sample name (e.g. "FD22124874"). For everything else,
                # leave blank so the Box xlsx merge doesn't false-match.
                desc = sample if typ == "Unknown" else ""
                w.writerow([well, typ, desc])
        print(f"  wrote {inputfile_path}")

    return out_paths


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print("usage: split_multiplate_csv.py <input.csv> [out_dir]", file=sys.stderr)
        return 2
    inp = Path(argv[1])
    out_dir = Path(argv[2]) if len(argv) > 2 else inp.parent / "split"
    split(inp, out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
