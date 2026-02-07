# AGENTS.md

This file is the authoritative project guidance for `hyphal-flow`.

## Project Layout (CoMPhy Standard)

- `src-local/`: project-specific Basilisk headers and model extensions.
- `simulationCases/`: simulation entry points (`.c`) and case output folders.
- `postProcess/`: post-processing and plotting utilities.
- Repository root: run scripts and parameter files (`*.params`, `*.sbatch`).

## Build and Run

- Compile single case (OpenMP):
  ```bash
  qcc -Wall -O2 -fopenmp -I$(PWD)/src-local simulationCases/hypha.c -o hypha -lm
  ```
- Compile MPI case:
  ```bash
  CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions -I$(PWD)/src-local simulationCases/hypha.c -o hyphaMPI -lm
  ```
- Run parameter sweep on Hamilton:
  ```bash
  sbatch runSweepHamilton.sbatch
  ```
  - Uses `sweep.params` by default for sweep values (`SWEEP_Ec_h`).
  - You can provide a custom sweep file:
    ```bash
    sbatch runSweepHamilton.sbatch my-sweep.params
    ```
  - Uses `default.params` (or `BASE_CONFIG` from the sweep file) as base config.
  - Writes outputs under `simulationCases/<CaseNo>/`.

## Repository Rules

- Keep simulation sources in `simulationCases/`; keep headers in `src-local/`.
- Do not delete or rewrite existing case output directories unless explicitly requested.
- Do not commit generated binaries, restart files, or logs.
- Do not commit local Basilisk checkouts (`basilisk/`) or qcc-lsp workspace state (`.comphy-basilisk/`).
- Keep `README.md` synchronized with actual repository paths and scripts.

## Agent File Policy

- `AGENTS.md` is committed and maintained.
- `CLAUDE.md` is a local pointer file and is gitignored.
