"""
# plot_vcm_vs_Ec_H.py

Aggregate multiple simulation log files, extract `(Ec_h, max(vcm))`,
and plot the trend on log-log axes.

## Dependencies

- `numpy`: numeric arrays and sorting.
- `matplotlib`: log-log plotting.
- `argparse`: command-line interface.
- `re`: header/value parsing.

## Workflow

1. Scan files in `--log_dir` filtered by `--pattern`.
2. Parse `Ec_h` from header text.
3. Parse the tabular section and compute `max(vcm)`.
4. Plot and optionally save the resulting curve.

#### Example

```bash
python3 postProcess/plot_vcm_vs_Ec_H.py --log_dir simulationCases --pattern "log*" --out vcm_vs_Ec_h.png
```
"""

import argparse
import os
import re
import sys
from typing import List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt


EC_KEY = "Ec_h"


def looks_binary(sample: bytes) -> bool:
    """
    Return `True` when a byte sample likely represents binary content.

    #### Args

    - `sample`: Leading file bytes used for a quick text/binary heuristic.

    #### Returns

    - `bool`: `True` if the sample appears non-textual.
    """
    if b"\x00" in sample:
        return True
    # count non-text-ish bytes
    text_chars = bytearray({7, 8, 9, 10, 12, 13, 27} | set(range(0x20, 0x100)))
    nontext = sum(b not in text_chars for b in sample)
    return nontext / max(1, len(sample)) > 0.30


def read_text_lines(path: str) -> Optional[List[str]]:
    """
    Read a file as text lines, returning `None` for binary or unreadable files.

    #### Args

    - `path`: Candidate log file path.

    #### Returns

    - `Optional[List[str]]`: File lines when readable text, else `None`.
    """
    try:
        with open(path, "rb") as fb:
            sample = fb.read(4096)
            if looks_binary(sample):
                return None
        with open(path, "r", errors="ignore") as f:
            return f.readlines()
    except OSError:
        return None


def extract_ec_h(lines: List[str]) -> Optional[float]:
    """
    Extract the `Ec_h` value from free-form header lines.

    #### Args

    - `lines`: Full log file lines.

    #### Returns

    - `Optional[float]`: Parsed `Ec_h` value, or `None` if unavailable.
    """
    for line in lines:
        m = re.search(rf"\b{re.escape(EC_KEY)}\s+([0-9.eE+-]+)\b", line)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                return None
    return None


def find_table_header(lines: List[str]) -> Tuple[Optional[int], Optional[List[str]]]:
    """
    Return the index and tokenized header for the first table header row.

    #### Args

    - `lines`: Full log file lines.

    #### Returns

    - `Tuple[Optional[int], Optional[List[str]]]`: Header row index and tokens.
    """
    for i, line in enumerate(lines):
        if line.strip().startswith("i"):
            headers = line.split()
            return i, headers
    return None, None


def extract_max_vcm(lines: List[str]) -> Optional[float]:
    """
    Compute the maximum `vcm` value from the parsed data table.

    #### Args

    - `lines`: Full log file lines.

    #### Returns

    - `Optional[float]`: Maximum `vcm`, or `None` if parsing fails.
    """
    header_idx, headers = find_table_header(lines)
    if headers is None or "vcm" not in headers or header_idx is None:
        return None

    vcm_col = headers.index("vcm")
    vcm_vals: List[float] = []

    for line in lines[header_idx + 1 :]:
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) <= vcm_col:
            continue
        try:
            vcm_vals.append(float(parts[vcm_col]))
        except ValueError:
            continue

    if not vcm_vals:
        return None
    return max(vcm_vals)


def iter_files(log_dir: str, pattern: str) -> List[str]:
    """
    Collect matching files from a directory using `*` wildcard matching.

    #### Args

    - `log_dir`: Directory containing log files.
    - `pattern`: Filename filter supporting `*` wildcard only.

    #### Returns

    - `List[str]`: Sorted matching file paths.
    """
    # Avoid importing glob to keep it very simple and predictable
    # Convert "log*" -> regex "^log.*$"
    pat = re.escape(pattern).replace(r"\*", ".*")
    rx = re.compile(rf"^{pat}$")

    out = []
    for name in sorted(os.listdir(log_dir)):
        full = os.path.join(log_dir, name)
        if os.path.isfile(full) and rx.match(name):
            out.append(full)
    return out


def main() -> int:
    """
    Parse logs, extract `(Ec_h, max(vcm))`, and generate the summary plot.

    #### Returns

    - `int`: Process exit code (`0` success, non-zero on errors).
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("--log_dir", type=str, required=True, help="Directory containing log files.")
    ap.add_argument("--pattern", type=str, default="*", help='Filename pattern, e.g. "log*" (supports * only).')
    ap.add_argument("--out", type=str, default=None, help="If set, save figure to this path (png/pdf/etc).")
    ap.add_argument("--include_zero_ec", action="store_true",
                    help="Include Ec_h<=0 points in the data output (still cannot be on log-log plot).")
    ap.add_argument("--show", action="store_true", help="Show the plot window.")
    args = ap.parse_args()

    if not os.path.isdir(args.log_dir):
        print(f"Error: not a directory: {args.log_dir}", file=sys.stderr)
        return 2

    paths = iter_files(args.log_dir, args.pattern)
    if not paths:
        print(f"No files matched in {args.log_dir} with pattern '{args.pattern}'", file=sys.stderr)
        return 1

    Ec_vals = []
    vmax_vals = []
    used_files = 0
    skipped = 0

    for path in paths:
        lines = read_text_lines(path)
        if lines is None:
            skipped += 1
            continue

        Ec = extract_ec_h(lines)
        vmax = extract_max_vcm(lines)

        if Ec is None or vmax is None:
            skipped += 1
            continue

        if (Ec <= 0.0) and (not args.include_zero_ec):
            # can't go on log-log, and usually not physically comparable
            skipped += 1
            continue

        Ec_vals.append(Ec)
        vmax_vals.append(vmax)
        used_files += 1

    if used_files == 0:
        print("No usable (Ec_h, max vcm) pairs found.", file=sys.stderr)
        print("Common causes:", file=sys.stderr)
        print("  - files are binary / restart files", file=sys.stderr)
        print("  - Ec_h not present in header", file=sys.stderr)
        print("  - table header not found or 'vcm' missing", file=sys.stderr)
        return 1

    Ec_vals = np.array(Ec_vals, dtype=float)
    vmax_vals = np.array(vmax_vals, dtype=float)

    # For plotting on log-log, filter non-positive Ec again
    plot_mask = Ec_vals > 0
    Ec_plot = Ec_vals[plot_mask]
    vmax_plot = vmax_vals[plot_mask]

    # Sort by Ec
    order = np.argsort(Ec_plot)
    Ec_plot = Ec_plot[order]
    vmax_plot = vmax_plot[order]

    plt.figure()
    plt.loglog(Ec_plot, vmax_plot, "o-")
    plt.xlabel("Ec_h")
    plt.ylabel("max(vcm)")
    plt.title("max(vcm) vs Ec_h")
    plt.grid(True, which="both", ls="--")

    if args.out:
        plt.savefig(args.out, dpi=300, bbox_inches="tight")
        print(f"Saved figure to: {args.out}")

    if args.show or not args.out:
        plt.show()

    print(f"Used files: {used_files}, skipped: {skipped}, total matched: {len(paths)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
