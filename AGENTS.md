# AGENTS.md

This file is the authoritative project guidance for `hyphal-flow`.

## Project Layout (CoMPhy Standard)

- `src-local/`: project-specific Basilisk headers and the C/shell parameter parsers.
- `simulationCases/`: the `hypha.c` entry point and generated `simulationCases/<CaseNo>/` output.
- `postProcess/`: extraction, plotting, and snapshot video.
- `legacy/`: superseded Cartesian capillary solvers and the previous in-tree VE headers. Do not compile them.
- Repository root: `runSimulation.sh`, `runParameterSweep.sh`, `*.params`, and site batch wrappers.

## Solver

There is one production case: `simulationCases/hypha.c`.

- Axisymmetric tube, three non-coalescing phases (drop `f1`, film `f2`, cytoplasm).
- Log-conformation stress from MultiRheoFlow
  `log-conform-viscoelastic-scalar-2D.h`.
- Per-phase `G` and `lambda` are interpolated in
  `three-phase-nonCoalescing-VE.h`.
- `DriveMode=bond` (periodic body force) or `DriveMode=pressure` (open ends).

Do not revive `hypha-capillary.c` or `Hypha-capillary-length.c` as extra
entry points. They mixed Cartesian geometry with inconsistent boundary
conditions; the useful knobs now live in `hypha.c` / `default.params`.

## Build and Run

```bash
bash runSimulation.sh default.params
bash runSimulation.sh default.params --mpi --CPUs 8
bash runParameterSweep.sh sweep.params --dry-run
bash runParameterSweep.sh sweep.params
```

- Runners source `src-local/parse_params.sh`. Do not add a third parser.
- `.project_config` is optional. `qcc` must be on `PATH`.
- Never hardcode a machine-local `qcc` path.
- Compile against CoMPhy Basilisk `v2026-08-30` (`basilisk-C` ref-locked installer).

## Parameters

- C: `params_init_from_argv` in `hypha.c`, typed reads in `src-local/params.h`.
- Shell: `parse_param_file`, `get_param`, `set_param_in_file` in `src-local/parse_params.sh`.
- Base file: `default.params`. Sweeps: `sweep.params`, `sweep-128.params`.

## Repository Rules

- Keep headers in `src-local/` and the live solver in `simulationCases/`.
- Do not delete case output directories unless explicitly requested.
- Do not commit binaries, restart files, logs, `basilisk/`, `.comphy-basilisk/`, or `.docker_mode`.
- `CLAUDE.md` is local (`@AGENTS.md`) and gitignored. `.coderabbit.yaml` is committed.
- Keep `README.md` synchronized with actual paths. Repository trees in the README must use a plain fenced block (no language tag) and `path - description` lines.

## Agent File Policy

- `AGENTS.md` is committed and maintained.
- `CLAUDE.md` is a local pointer file and is gitignored.
