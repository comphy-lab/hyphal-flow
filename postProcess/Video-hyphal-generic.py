#!/usr/bin/env python3
"""
# Video-hyphal-generic.py

Create a mirrored-axis video from Basilisk `snapshot-*` states for
three-phase hyphal-flow runs.

## Pipeline

1. Restore each snapshot and compute drop center of mass from `f1`.
2. Sample a right-half scalar field (`vel`, `d2c`, or `trA`) on a uniform grid.
3. Optionally sample a different left-half scalar field.
4. Overlay `f1` and `f2` interface facets.
5. Write PNG frames and optionally assemble an MP4 with `ffmpeg`.

## Dependencies

- `qcc`: builds helper binaries from `getFacet-threePhase.c` and
  `getData-elastic-nonCoalescence.c`.
- `numpy`: parsing, masking, and percentile-based color scaling.
- `matplotlib`: frame rendering.
- `ffmpeg`: optional MP4 assembly (skipped with `--skip-video`).

## Default Visualization

- Left half (`r < 0`): `trA`
- Right half (`r >= 0`): `vel`
- Domain window:
  - Basilisk `x` in `0 +- 2`
  - Basilisk `y` in `y_CoM(f1) +- 4`
- Default duration: `10 s` with `fps = N_frames / duration` when `--fps` is unset.

#### Example

```bash
python3 postProcess/Video-hyphal-generic.py --case-dir simulationCases/1000
```
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
  import numpy as _np_types

  NDArray = _np_types.ndarray
  MaskedArray = _np_types.ma.MaskedArray
else:
  NDArray = Any
  MaskedArray = Any

np: Any = None
plt: Any = None
LineCollection: Any = None
WORKER_ENV_PID: int | None = None


FIELD_INDEX = {"d2c": 2, "vel": 3, "trA": 4}
FIELD_LABEL = {
  "d2c": r"$\log_{10}\!\left(\lVert\boldsymbol{\mathcal{D}}\rVert\right)$",
  "vel": r"$\lvert \mathbf{u} \rvert$",
  "trA": r"$\log_{10}\!\left(\mathrm{tr}(\mathbf{A})/3\right)$",
}

COM_HELPER_C = r"""
#include "utils.h"
#include <math.h>

scalar f1[], f2[];

int main (int argc, char const *argv[]) {
  if (argc < 2) {
    fprintf (stderr, "Usage: %s snapshot\n", argv[0]);
    return 1;
  }

  if (!restore (file = argv[1])) {
    fprintf (stderr, "Failed to restore snapshot: %s\n", argv[1]);
    return 2;
  }

  double wt = 0., xcm = 0., ycm = 0.;
  foreach (reduction(+:wt) reduction(+:xcm) reduction(+:ycm)) {
    double w = f1[]*sq(Delta);
    wt += w;
    xcm += x*w;
    ycm += y*w;
  }

  if (wt <= 0.) {
    printf ("nan nan\n");
    return 0;
  }

  xcm /= wt;
  ycm /= wt;
  printf ("%.16g %.16g\n", xcm, ycm);
  return 0;
}
"""


def parse_args() -> argparse.Namespace:
  """
  Parse command-line options for frame extraction and video assembly.

  #### Returns

  - `argparse.Namespace`: Parsed command-line arguments.
  """
  parser = argparse.ArgumentParser(
    description="Create a generic hyphal-flow video with f1/f2 interfaces.",
  )
  parser.add_argument(
    "--case-dir",
    type=Path,
    default=Path.cwd(),
    help="Case directory containing `intermediate/snapshot-*` (default: cwd).",
  )
  parser.add_argument(
    "--snap-glob",
    default="intermediate/snapshot-*",
    help="Snapshot glob pattern relative to `case-dir`.",
  )
  parser.add_argument(
    "--field",
    choices=tuple(FIELD_INDEX.keys()),
    default="vel",
    help="Scalar field from getData-elastic-nonCoalescence.",
  )
  parser.add_argument(
    "--ny",
    type=int,
    default=400,
    help="Number of grid points along y for sampled scalar field.",
  )
  parser.add_argument(
    "--cpus",
    "--CPUs",
    dest="cpus",
    type=int,
    default=4,
    help="Number of worker processes for snapshot rendering (default: 4).",
  )
  parser.add_argument(
    "--fps",
    type=float,
    default=None,
    help="Output FPS. If omitted, computed from --duration.",
  )
  parser.add_argument(
    "--duration",
    type=float,
    default=10.0,
    help="Target output duration in seconds when --fps is not set (default: 10).",
  )
  parser.add_argument(
    "-o",
    "--output",
    default="video-hyphal-generic.mp4",
    help="Output MP4 path (relative paths resolved inside case-dir).",
  )
  parser.add_argument(
    "--frames-dir",
    type=Path,
    help="Directory for PNG frames (default: case-dir/Video).",
  )
  parser.add_argument(
    "--clean-frames",
    dest="clean_frames",
    action="store_true",
    default=True,
    help="Delete existing PNG files in frames-dir before rendering (default).",
  )
  parser.add_argument(
    "--no-clean-frames",
    dest="clean_frames",
    action="store_false",
    help="Keep existing PNG files in frames-dir.",
  )
  parser.add_argument(
    "--skip-video",
    action="store_true",
    help="Only write PNG frames and skip ffmpeg MP4 assembly.",
  )
  parser.add_argument(
    "--ffmpeg",
    default="ffmpeg",
    help="ffmpeg executable name/path.",
  )
  parser.add_argument("--max-frames", type=int, help="Render only first N frames.")
  parser.add_argument("--start-time", type=float, help="Skip snapshots with t < this.")
  parser.add_argument("--end-time", type=float, help="Skip snapshots with t > this.")
  parser.add_argument(
    "--vmin",
    type=float,
    help="Right-field color scale minimum (default: vel->0, d2c->-3, trA->-0.5).",
  )
  parser.add_argument(
    "--vmax",
    type=float,
    help="Right-field color scale maximum (default: vel->0.1, d2c->0, trA->0.5).",
  )
  parser.add_argument("--cmap", default=None, help="Right-field colormap (default: field-based auto).")
  parser.add_argument(
    "--left-field",
    choices=tuple(FIELD_INDEX.keys()),
    default="trA",
    help="Field on the left half (r < 0). Default: trA.",
  )
  parser.add_argument(
    "--no-left-field",
    action="store_true",
    help="Disable left-half field rendering.",
  )
  parser.add_argument(
    "--left-D2",
    dest="left_d2",
    action="store_true",
    help="Use D2 field on the left half (equivalent to --left-field d2c).",
  )
  parser.add_argument(
    "--left-tra",
    action="store_true",
    help="Use tr(A) field on the left half (equivalent to --left-field trA).",
  )
  parser.add_argument(
    "--left-vmin",
    type=float,
    help="Left-field color scale minimum (default: vel->0, d2c->-3, trA->-0.5).",
  )
  parser.add_argument(
    "--left-vmax",
    type=float,
    help="Left-field color scale maximum (default: vel->0.1, d2c->0, trA->0.5).",
  )
  parser.add_argument("--left-cmap", default=None, help="Left-field colormap (default: field-based auto).")

  # Requested defaults: x in 0 ± 2 and y in CoM ± 4.
  parser.add_argument("--x-center", type=float, default=0.0)
  parser.add_argument("--x-half-width", type=float, default=2.0)
  parser.add_argument("--y-half-width", type=float, default=4.0)
  return parser.parse_args()


def ensure_python_dependencies() -> None:
  """
  Import `numpy`/`matplotlib` lazily and configure plotting defaults.

  #### Raises

  - `RuntimeError`: Required plotting dependencies are not installed.
  """
  global np, plt, LineCollection

  missing = []
  _matplotlib: Any | None = None
  try:
    import numpy as _np
  except ModuleNotFoundError:
    missing.append("numpy")
    _np = None

  try:
    import matplotlib as _matplotlib
    _matplotlib.use("Agg")
    import matplotlib.pyplot as _plt
    from matplotlib.collections import LineCollection as _LineCollection
  except ModuleNotFoundError:
    missing.append("matplotlib")
    _plt = None
    _LineCollection = None

  if missing:
    missing_str = ", ".join(sorted(set(missing)))
    raise RuntimeError(
      f"Missing required Python packages: {missing_str}. "
      "Install them (e.g. `pip install numpy matplotlib`)."
    )

  # Publication-style defaults with LaTeX typography.
  assert _matplotlib is not None
  _matplotlib.rcParams["font.family"] = "serif"
  _matplotlib.rcParams["font.serif"] = ["Computer Modern Roman", "DejaVu Serif"]
  _matplotlib.rcParams["mathtext.fontset"] = "cm"
  _matplotlib.rcParams["text.usetex"] = True
  _matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
  _matplotlib.rcParams["axes.linewidth"] = 2.5

  np = _np
  plt = _plt
  LineCollection = _LineCollection


def ensure_plotting_runtime() -> None:
  """
  Ensure plotting globals are available in the current process.
  """
  if np is None or plt is None or LineCollection is None:
    ensure_python_dependencies()


def configure_worker_environment(cache_root: Path | None) -> None:
  """
  Configure process-local matplotlib and LaTeX cache directories.
  """
  global WORKER_ENV_PID

  if cache_root is None:
    return

  pid = os.getpid()
  if WORKER_ENV_PID == pid:
    return

  worker_root = cache_root / f"worker-{pid}"
  mpl_dir = worker_root / "mplconfig"
  tex_var_dir = worker_root / "texmf-var"
  tex_cfg_dir = worker_root / "texmf-config"
  for path in (mpl_dir, tex_var_dir, tex_cfg_dir):
    path.mkdir(parents=True, exist_ok=True)

  os.environ["MPLCONFIGDIR"] = str(mpl_dir)
  os.environ["TEXMFVAR"] = str(tex_var_dir)
  os.environ["TEXMFCONFIG"] = str(tex_cfg_dir)
  os.environ.setdefault("OMP_NUM_THREADS", "1")

  WORKER_ENV_PID = pid


def run_capture(cmd: list[str], cwd: Path | None = None) -> str:
  """
  Run a subprocess and return combined `stdout` + `stderr` text.

  #### Args

  - `cmd`: Command and arguments passed to `subprocess.run`.
  - `cwd`: Working directory used for command execution.

  #### Returns

  - `str`: Concatenated standard output and error output.

  #### Raises

  - `subprocess.CalledProcessError`: The command exits with a non-zero code.
  """
  result = subprocess.run(
    cmd,
    cwd=cwd,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
    check=True,
  )
  return (result.stdout or "") + (result.stderr or "")


def snapshot_argument(snapshot: Path, case_dir: Path) -> str:
  """
  Return a case-relative snapshot path when possible.
  """
  try:
    return str(snapshot.relative_to(case_dir))
  except ValueError:
    return str(snapshot)


def snapshot_time(path: Path) -> float:
  """
  Extract time from a filename of form `snapshot-<time>`.
  """
  name = path.name
  if "snapshot-" not in name:
    return math.inf
  raw = name.split("snapshot-", 1)[1]
  try:
    return float(raw)
  except ValueError:
    return math.inf


def list_snapshots(case_dir: Path, pattern: str) -> list[Path]:
  """
  Collect and sort snapshot files by simulation time.
  """
  snapshots = sorted(case_dir.glob(pattern), key=snapshot_time)
  return [p for p in snapshots if p.is_file()]


def compile_get_helper(source: Path, output: Path) -> None:
  """
  Compile one get* helper with `qcc`.
  """
  cmd = [
    "qcc",
    "-O2",
    "-Wall",
    "-disable-dimensions",
    source.name,
    "-o",
    str(output),
    "-lm",
  ]
  subprocess.run(cmd, check=True, cwd=source.parent)


def precompile_get_helpers(script_dir: Path, build_dir: Path) -> tuple[Path, Path, Path]:
  """
  Pre-processing step: compile get* helper binaries before rendering.

  #### Returns

  - `tuple[Path, Path, Path]`: `(facet_bin, data_bin, com_bin)`.
  """
  if shutil.which("qcc") is None:
    raise RuntimeError("qcc not found in PATH.")

  facet_src = script_dir / "getFacet-threePhase.c"
  data_src = script_dir / "getData-elastic-nonCoalescence.c"
  com_src = build_dir / "getDropCom-f1.c"

  if not facet_src.exists() or not data_src.exists():
    raise FileNotFoundError(
      "Required files not found in postProcess/: "
      "getFacet-threePhase.c and/or getData-elastic-nonCoalescence.c"
    )

  com_src.write_text(COM_HELPER_C, encoding="utf-8")

  facet_bin = build_dir / "getFacet-threePhase"
  data_bin = build_dir / "getData-elastic-nonCoalescence"
  com_bin = build_dir / "getDropCom-f1"

  compile_get_helper(facet_src, facet_bin)
  compile_get_helper(data_src, data_bin)
  compile_get_helper(com_src, com_bin)

  return facet_bin, data_bin, com_bin


def parse_facet_segments(raw: str) -> NDArray:
  """
  Parse `output_facets` text into a `N x 2 x 2` segment array.
  """
  points: list[list[float]] = []
  for line in raw.splitlines():
    line = line.strip()
    if not line:
      continue
    vals = line.split()
    if len(vals) < 2:
      continue
    try:
      points.append([float(vals[0]), float(vals[1])])
    except ValueError:
      continue

  if len(points) < 2:
    return np.empty((0, 2, 2), dtype=float)

  usable = len(points) - (len(points) % 2)
  arr = np.asarray(points[:usable], dtype=float)
  return arr.reshape(-1, 2, 2)


def get_facets(snapshot: Path, facet_bin: Path, case_dir: Path, include_f1: bool) -> NDArray:
  """
  Extract interface segments for `f1` (`include_f1=True`) or `f2`.
  """
  mode = "true" if include_f1 else "false"
  raw = run_capture(
    [str(facet_bin), snapshot_argument(snapshot, case_dir), mode],
    cwd=case_dir,
  )
  return parse_facet_segments(raw)


def get_drop_com(snapshot: Path, com_bin: Path, case_dir: Path) -> tuple[float, float]:
  """
  Compute `f1` center of mass `(xcm, ycm)` for one snapshot.
  """
  raw = run_capture(
    [str(com_bin), snapshot_argument(snapshot, case_dir)],
    cwd=case_dir,
  ).strip().splitlines()

  if not raw:
    return float("nan"), float("nan")
  vals = raw[-1].split()
  if len(vals) < 2:
    return float("nan"), float("nan")
  try:
    return float(vals[0]), float(vals[1])
  except ValueError:
    return float("nan"), float("nan")


def get_field_grid(
  snapshot: Path,
  data_bin: Path,
  case_dir: Path,
  field_key: str,
  xmin: float,
  ymin: float,
  xmax: float,
  ymax: float,
  ny: int,
) -> tuple[NDArray, NDArray, MaskedArray]:
  """
  Sample one derived field on a uniform grid inside `[xmin, xmax] x [ymin, ymax]`.
  """
  raw = run_capture(
    [
      str(data_bin),
      snapshot_argument(snapshot, case_dir),
      f"{xmin:.16g}",
      f"{ymin:.16g}",
      f"{xmax:.16g}",
      f"{ymax:.16g}",
      str(ny),
      "0",
      "0",
      "0",
    ],
    cwd=case_dir,
  )

  rows = []
  for line in raw.splitlines():
    line = line.strip()
    if not line:
      continue
    vals = line.split()
    if len(vals) < 5:
      continue
    try:
      rows.append([float(vals[0]), float(vals[1]), float(vals[2]), float(vals[3]), float(vals[4])])
    except ValueError:
      continue

  if not rows:
    raise RuntimeError(f"No field data parsed for snapshot: {snapshot}")

  arr = np.asarray(rows, dtype=float)
  x = arr[:, 0]
  y = arr[:, 1]
  field = arr[:, FIELD_INDEX[field_key]]

  x_unique = np.unique(x)
  y_unique = np.unique(y)
  ix = np.searchsorted(x_unique, x)
  iy = np.searchsorted(y_unique, y)

  grid = np.full((len(y_unique), len(x_unique)), np.nan, dtype=float)
  grid[iy, ix] = field

  invalid = (~np.isfinite(grid)) | (np.abs(grid) > 1e20)
  return x_unique, y_unique, np.ma.array(grid, mask=invalid)


def grid_extent(x: NDArray, y: NDArray) -> list[float]:
  """
  Convert center-coordinates to image extent bounds.
  """
  dx = float(np.median(np.diff(x))) if len(x) > 1 else 1.0
  dy = float(np.median(np.diff(y))) if len(y) > 1 else 1.0
  return [x[0] - 0.5 * dx, x[-1] + 0.5 * dx, y[0] - 0.5 * dy, y[-1] + 0.5 * dy]


def map_segments_xy_to_rz(segments: NDArray) -> NDArray:
  """
  Map Basilisk `(x, y)` segments to plotting `(r, z)` with `r = y`, `z = x`.
  """
  if len(segments) == 0:
    return segments
  return segments[..., [1, 0]]


def mirror_segments_about_r0(segments_rz: NDArray) -> NDArray:
  """
  Mirror `(r, z)` segments about `r = 0` and return combined segments.
  """
  if len(segments_rz) == 0:
    return segments_rz
  mirrored = segments_rz.copy()
  mirrored[..., 0] *= -1.0
  return np.concatenate([mirrored, segments_rz], axis=0)


def mirror_field_xy_to_rz(field_xy: MaskedArray, r_pos: NDArray) -> tuple[NDArray, MaskedArray]:
  """
  Build a full `(-r, +r)` field in `(r, z)` coordinates from positive-`r` data.
  """
  field_pos = np.ma.array(field_xy.T, copy=False)  # (nz, nr_pos)
  r_positive = np.asarray(r_pos, dtype=float)
  r_negative = -r_positive[::-1]
  field_negative = field_pos[:, ::-1]
  r_full = np.concatenate([r_negative, r_positive])
  field_full = np.ma.concatenate([field_negative, field_pos], axis=1)
  return r_full, field_full


def render_frame(
  frame_path: Path,
  t: float,
  xcm: float,
  ycm: float,
  x: NDArray,
  y: NDArray,
  field: MaskedArray,
  left_field: MaskedArray | None,
  left_field_key: str | None,
  f1_segments: NDArray,
  f2_segments: NDArray,
  args: argparse.Namespace,
  vmin: float | None,
  vmax: float | None,
  left_vmin: float | None,
  left_vmax: float | None,
) -> None:
  """
  Render one PNG frame with field background and `f1/f2` interface overlays.
  """
  fig, ax = plt.subplots(figsize=(10.5, 5.5), dpi=180)

  # Mapping: Basilisk x -> plotted z (vertical), Basilisk y -> plotted r (horizontal).
  # We mirror about r=0 so the full domain is visible.
  r_full, field_rz = mirror_field_xy_to_rz(field, y)
  extent_rz = grid_extent(r_full, x)

  left_image = None
  if left_field is not None:
    _, left_rz = mirror_field_xy_to_rz(left_field, y)
    left_plot = np.ma.array(left_rz, copy=True)
    right_plot = np.ma.array(field_rz, copy=True)
    left_mask = r_full >= 0.0
    right_mask = r_full < 0.0
    left_plot[:, left_mask] = np.ma.masked
    right_plot[:, right_mask] = np.ma.masked
    left_image = ax.imshow(
      left_plot,
      origin="lower",
      extent=extent_rz,
      cmap=args.left_cmap,
      vmin=left_vmin,
      vmax=left_vmax,
      aspect="equal",
      interpolation="nearest",
    )
    display_field = right_plot
  else:
    display_field = field_rz

  image = ax.imshow(
    display_field,
    origin="lower",
    extent=extent_rz,
    cmap=args.cmap,
    vmin=vmin,
    vmax=vmax,
    aspect="equal",
    interpolation="nearest",
  )

  f1_rz = mirror_segments_about_r0(map_segments_xy_to_rz(f1_segments))
  f2_rz = mirror_segments_about_r0(map_segments_xy_to_rz(f2_segments))

  if len(f2_rz):
    # High-contrast two-layer styling so f2 remains visible on both light/dark backgrounds.
    ax.add_collection(LineCollection(f2_rz, colors="black", linewidths=3.4, alpha=0.95))
    ax.add_collection(LineCollection(f2_rz, colors="#00E5FF", linewidths=1.9, alpha=0.98))
  if len(f1_rz):
    ax.add_collection(LineCollection(f1_rz, colors="black", linewidths=2.4, alpha=1.0))

  if np.isfinite(xcm):
    # Requested: hardcode r_CoM to zero.
    ax.plot(0.0, xcm, "o", ms=7.0, mfc="gold", mec="black", mew=1.2)

  # Window is defined in Basilisk coordinates and then mapped:
  # x in [x_center ± x_half_width] -> z-axis
  # r is centered at 0 and shown as [-y_half_width, +y_half_width].
  xmin = args.x_center - args.x_half_width
  xmax = args.x_center + args.x_half_width
  ax.set_xlim(-args.y_half_width, args.y_half_width)
  ax.set_ylim(xmin, xmax)
  ax.set_aspect("equal")

  # Requested: no ticks/labels on the coordinate axes.
  ax.set_xlabel("")
  ax.set_ylabel("")
  ax.set_xticks([])
  ax.set_yticks([])
  ax.tick_params(
    axis="both",
    which="both",
    bottom=False,
    top=False,
    left=False,
    right=False,
    labelbottom=False,
    labelleft=False,
  )
  for spine in ax.spines.values():
    spine.set_linewidth(2.5)

  ax.set_title(
    rf"$t={t:0.4f},\ \mathrm{{CoM}}(f_1):\ z={xcm:0.4f}$",
    fontsize=32,
    pad=16,
  )

  fig.tight_layout()

  # Slim vertical colorbars on both sides of the main axes.
  l, b, w, h = ax.get_position().bounds
  cbar_w = 0.018
  cbar_gap = 0.012

  # Right colorbar for the right-half field.
  cax_right = fig.add_axes([l + w + cbar_gap, b, cbar_w, h])
  cbar = fig.colorbar(image, cax=cax_right, orientation="vertical")
  cbar.set_label(FIELD_LABEL[args.field], fontsize=20, labelpad=8)
  cbar.ax.tick_params(labelsize=16, width=1.4, length=5, direction="out")
  cbar.outline.set_linewidth(1.4)

  # Left colorbar for the left-half field (if enabled).
  if left_image is not None and left_field_key is not None:
    cax_left = fig.add_axes([l - cbar_gap - cbar_w, b, cbar_w, h])
    cbar_left = fig.colorbar(left_image, cax=cax_left, orientation="vertical")
    cbar_left.set_label(FIELD_LABEL[left_field_key], fontsize=20, labelpad=8)
    cbar_left.ax.tick_params(labelsize=16, width=1.4, length=5, direction="out")
    cbar_left.outline.set_linewidth(1.4)
    cbar_left.ax.yaxis.set_ticks_position("left")
    cbar_left.ax.yaxis.set_label_position("left")

  fig.savefig(frame_path, bbox_inches="tight")
  plt.close(fig)


def format_fps(fps: float) -> str:
  """
  Format FPS for ffmpeg CLI arguments.
  """
  return f"{fps:.6f}".rstrip("0").rstrip(".")


def default_limits_for_field(field_key: str) -> tuple[float | None, float | None]:
  """
  Return fixed default color limits for known fields.

  - `d2c`: [-3, 0]
  - `vel`: [0, 0.1]
  - `trA`: [-0.5, 0.5]
  """
  if field_key == "d2c":
    return -3.0, 0.0
  if field_key == "vel":
    return 0.0, 0.1
  if field_key == "trA":
    return -0.1, 0.1
  return None, None


def default_cmap_for_field(field_key: str) -> str:
  """
  Return a default colormap per field.
  """
  if field_key == "d2c":
    return "hot_r"
  if field_key == "trA":
    return "bwr"
  return "viridis"


def render_single_snapshot(
  idx: int,
  snapshot: Path,
  case_dir: Path,
  frames_dir: Path,
  facet_bin: Path,
  data_bin: Path,
  com_bin: Path,
  args: argparse.Namespace,
  use_vmin: float | None,
  use_vmax: float | None,
  use_left_vmin: float | None,
  use_left_vmax: float | None,
  worker_cache_root: Path | None,
) -> tuple[int, Path]:
  """
  Render one frame for a snapshot and return `(index, frame_path)`.
  """
  configure_worker_environment(worker_cache_root)
  ensure_plotting_runtime()

  t = snapshot_time(snapshot)
  xcm, ycm = get_drop_com(snapshot, com_bin, case_dir)

  xmin = args.x_center - args.x_half_width
  xmax = args.x_center + args.x_half_width
  ymin = 0.0
  ymax = args.y_half_width

  x, y, field = get_field_grid(snapshot, data_bin, case_dir, args.field, xmin, ymin, xmax, ymax, args.ny)
  left_field = None
  if args.left_field is not None:
    x_left, y_left, left_field = get_field_grid(
      snapshot, data_bin, case_dir, args.left_field, xmin, ymin, xmax, ymax, args.ny
    )
    if len(x_left) != len(x) or len(y_left) != len(y):
      raise RuntimeError("Left/right field grids are inconsistent.")

  f1_segments = get_facets(snapshot, facet_bin, case_dir, include_f1=True)
  f2_segments = get_facets(snapshot, facet_bin, case_dir, include_f1=False)

  frame_path = frames_dir / f"frame_{idx:06d}.png"
  render_frame(
    frame_path=frame_path,
    t=t,
    xcm=xcm,
    ycm=ycm,
    x=x,
    y=y,
    field=field,
    left_field=left_field,
    left_field_key=args.left_field,
    f1_segments=f1_segments,
    f2_segments=f2_segments,
    args=args,
    vmin=use_vmin,
    vmax=use_vmax,
    left_vmin=use_left_vmin,
    left_vmax=use_left_vmax,
  )
  return idx, frame_path


def render_snapshots(
  snapshots: list[Path],
  case_dir: Path,
  frames_dir: Path,
  facet_bin: Path,
  data_bin: Path,
  com_bin: Path,
  args: argparse.Namespace,
  use_vmin: float | None,
  use_vmax: float | None,
  use_left_vmin: float | None,
  use_left_vmax: float | None,
  worker_cache_root: Path | None,
) -> None:
  """
  Render all snapshots, batching work in chunks of `args.cpus`.
  """
  tasks = list(enumerate(snapshots))
  total = len(tasks)

  if args.cpus <= 1:
    for idx, snapshot in tasks:
      _, frame_path = render_single_snapshot(
        idx=idx,
        snapshot=snapshot,
        case_dir=case_dir,
        frames_dir=frames_dir,
        facet_bin=facet_bin,
        data_bin=data_bin,
        com_bin=com_bin,
        args=args,
        use_vmin=use_vmin,
        use_vmax=use_vmax,
        use_left_vmin=use_left_vmin,
        use_left_vmax=use_left_vmax,
        worker_cache_root=None,
      )
      print(f"[{idx + 1}/{total}] wrote {frame_path}", file=sys.stderr)
    return

  with ProcessPoolExecutor(max_workers=args.cpus) as executor:
    for start in range(0, total, args.cpus):
      batch = tasks[start:start + args.cpus]
      futures = [
        executor.submit(
          render_single_snapshot,
          idx,
          snapshot,
          case_dir,
          frames_dir,
          facet_bin,
          data_bin,
          com_bin,
          args,
          use_vmin,
          use_vmax,
          use_left_vmin,
          use_left_vmax,
          worker_cache_root,
        )
        for idx, snapshot in batch
      ]
      batch_results = [future.result() for future in futures]
      for idx, frame_path in sorted(batch_results, key=lambda item: item[0]):
        print(f"[{idx + 1}/{total}] wrote {frame_path}", file=sys.stderr)


def main() -> int:
  """
  Execute snapshot discovery, frame rendering, and optional MP4 assembly.

  #### Returns

  - `int`: Process exit code (`0` for success, non-zero for errors).
  """
  args = parse_args()
  if args.left_tra:
    args.left_field = "trA"
  if args.left_d2:
    args.left_field = "d2c"
  if args.no_left_field:
    args.left_field = None

  # Resolve automatic colormap defaults after field-selection flags.
  if args.cmap is None:
    args.cmap = default_cmap_for_field(args.field)
  if args.left_field is not None and args.left_cmap is None:
    args.left_cmap = default_cmap_for_field(args.left_field)

  try:
    ensure_python_dependencies()
  except RuntimeError as exc:
    print(f"Error: {exc}", file=sys.stderr)
    return 1

  case_dir = args.case_dir.resolve()
  script_dir = Path(__file__).resolve().parent

  snapshots = list_snapshots(case_dir, args.snap_glob)
  if args.start_time is not None:
    snapshots = [s for s in snapshots if snapshot_time(s) >= args.start_time]
  if args.end_time is not None:
    snapshots = [s for s in snapshots if snapshot_time(s) <= args.end_time]
  if args.max_frames is not None:
    snapshots = snapshots[: args.max_frames]

  if not snapshots:
    print(f"No snapshots found with pattern '{args.snap_glob}' in {case_dir}", file=sys.stderr)
    return 1
  if args.ny <= 2:
    print("--ny must be > 2", file=sys.stderr)
    return 1
  if args.duration <= 0:
    print("--duration must be > 0", file=sys.stderr)
    return 1
  if args.y_half_width <= 0:
    print("--y-half-width must be > 0", file=sys.stderr)
    return 1
  if args.cpus <= 0:
    print("--cpus must be > 0", file=sys.stderr)
    return 1
  if not args.skip_video and shutil.which(args.ffmpeg) is None:
    print(f"ffmpeg executable not found: {args.ffmpeg}", file=sys.stderr)
    return 1

  out_path = Path(args.output)
  if not out_path.is_absolute():
    out_path = case_dir / out_path
  out_path.parent.mkdir(parents=True, exist_ok=True)

  if args.frames_dir is None:
    frames_dir = case_dir / "Video"
  else:
    frames_dir = args.frames_dir if args.frames_dir.is_absolute() else (case_dir / args.frames_dir)
  frames_dir.mkdir(parents=True, exist_ok=True)

  if args.clean_frames:
    for old_png in frames_dir.glob("*.png"):
      old_png.unlink()

  temp_objects: list[tempfile.TemporaryDirectory[str]] = []
  temp_build = tempfile.TemporaryDirectory(prefix="video-hyphal-tools-", dir=case_dir)
  temp_objects.append(temp_build)
  build_dir = Path(temp_build.name)
  worker_cache_root: Path | None = None
  if args.cpus > 1:
    temp_worker_cache = tempfile.TemporaryDirectory(prefix="video-hyphal-worker-cache-", dir=case_dir)
    temp_objects.append(temp_worker_cache)
    worker_cache_root = Path(temp_worker_cache.name)

  try:
    print("Pre-processing: compiling get* helpers...", file=sys.stderr)
    facet_bin, data_bin, com_bin = precompile_get_helpers(script_dir, build_dir)

    first = snapshots[0]
    xmin0 = args.x_center - args.x_half_width
    xmax0 = args.x_center + args.x_half_width
    # Physical extraction in Basilisk y is only non-negative for this setup.
    # We mirror about r=0 at plotting stage.
    ymin0 = 0.0
    ymax0 = args.y_half_width
    _, _, field0 = get_field_grid(first, data_bin, case_dir, args.field, xmin0, ymin0, xmax0, ymax0, args.ny)

    fixed_vmin, fixed_vmax = default_limits_for_field(args.field)
    need_auto_right = (
      (args.vmin is None and fixed_vmin is None) or
      (args.vmax is None and fixed_vmax is None)
    )
    auto_vmin, auto_vmax = None, None
    if need_auto_right:
      valid0 = field0.compressed()
      if valid0.size:
        auto_vmin = float(np.percentile(valid0, 2.0))
        auto_vmax = float(np.percentile(valid0, 98.0))
        if not np.isfinite(auto_vmin) or not np.isfinite(auto_vmax) or auto_vmin == auto_vmax:
          auto_vmin = float(np.nanmin(valid0))
          auto_vmax = float(np.nanmax(valid0))

    use_vmin = args.vmin if args.vmin is not None else (fixed_vmin if fixed_vmin is not None else auto_vmin)
    use_vmax = args.vmax if args.vmax is not None else (fixed_vmax if fixed_vmax is not None else auto_vmax)

    use_left_vmin, use_left_vmax = args.left_vmin, args.left_vmax
    if args.left_field is not None:
      fixed_left_vmin, fixed_left_vmax = default_limits_for_field(args.left_field)
      need_auto_left = (
        (use_left_vmin is None and fixed_left_vmin is None) or
        (use_left_vmax is None and fixed_left_vmax is None)
      )
      auto_left_vmin, auto_left_vmax = None, None
      if need_auto_left:
        _, _, field0_left = get_field_grid(
          first, data_bin, case_dir, args.left_field, xmin0, ymin0, xmax0, ymax0, args.ny
        )
        valid0_left = field0_left.compressed()
        if valid0_left.size:
          auto_left_vmin = float(np.percentile(valid0_left, 2.0))
          auto_left_vmax = float(np.percentile(valid0_left, 98.0))
          if (not np.isfinite(auto_left_vmin) or not np.isfinite(auto_left_vmax)
              or auto_left_vmin == auto_left_vmax):
            auto_left_vmin = float(np.nanmin(valid0_left))
            auto_left_vmax = float(np.nanmax(valid0_left))

      if use_left_vmin is None:
        use_left_vmin = fixed_left_vmin if fixed_left_vmin is not None else auto_left_vmin
      if use_left_vmax is None:
        use_left_vmax = fixed_left_vmax if fixed_left_vmax is not None else auto_left_vmax

    render_snapshots(
      snapshots=snapshots,
      case_dir=case_dir,
      frames_dir=frames_dir,
      facet_bin=facet_bin,
      data_bin=data_bin,
      com_bin=com_bin,
      args=args,
      use_vmin=use_vmin,
      use_vmax=use_vmax,
      use_left_vmin=use_left_vmin,
      use_left_vmax=use_left_vmax,
      worker_cache_root=worker_cache_root,
    )

    if args.skip_video:
      print(f"Frames written to: {frames_dir}", file=sys.stderr)
      return 0

    fps = float(args.fps) if args.fps is not None else (len(snapshots) / args.duration)
    fps = max(1e-6, fps)
    fps_str = format_fps(fps)

    # Matches requested ffmpeg style:
    # ffmpeg -framerate <fps> -pattern_type glob -i 'Video/*.png'
    #   -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -c:v libx264 -r <fps> -pix_fmt yuv420p out.mp4
    cmd = [
      args.ffmpeg,
      "-y",
      "-framerate",
      fps_str,
      "-pattern_type",
      "glob",
      "-i",
      str(frames_dir / "*.png"),
      "-vf",
      "pad=ceil(iw/2)*2:ceil(ih/2)*2",
      "-c:v",
      "libx264",
      "-r",
      fps_str,
      "-pix_fmt",
      "yuv420p",
      str(out_path),
    ]
    subprocess.run(cmd, check=True)
    print(
      f"Wrote video: {out_path} | fps={fps_str} | frames={len(snapshots)} | duration~{args.duration}s",
      file=sys.stderr,
    )
    return 0

  except subprocess.CalledProcessError as exc:
    print(f"Command failed with exit code {exc.returncode}: {exc.cmd}", file=sys.stderr)
    return 2
  except Exception as exc:  # noqa: BLE001
    print(f"Error: {exc}", file=sys.stderr)
    return 2
  finally:
    for obj in temp_objects:
      obj.cleanup()


if __name__ == "__main__":
  raise SystemExit(main())
