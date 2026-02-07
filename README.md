# hyphal-flow

Basilisk C simulations for drop transport in a single fungal hypha branch with
three-phase, non-coalescing viscoelastic modeling.

## Repository Layout

```text
hyphal-flow/
├── src-local/              # Project-specific Basilisk headers
├── simulationCases/        # Main simulation entry points and case outputs
├── postProcess/            # Analysis and plotting utilities
├── default.params          # Base parameters for sweeps
├── sweep.params            # Sweep metadata/config values
├── sweep-128.params        # 128-case sweep configuration
├── runSweepHamilton.sbatch # Hamilton batch sweep script
├── runSweepHamilton-serial-128.sbatch # Serial compile + parallel launch script
└── AGENTS.md               # Authoritative project instructions
```

## Requirements

- Basilisk `qcc`
- C compiler (GCC/Clang)
- `mpicc` and MPI runtime for MPI builds
- Python 3 (for plotting scripts and sweep value generation)

## Quick Start

### 1. Compile a simulation case

OpenMP build:

```bash
qcc -Wall -O2 -fopenmp -I$(PWD)/src-local simulationCases/hypha.c -o hypha -lm
```

MPI build:

```bash
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions \
  -I$(PWD)/src-local simulationCases/hypha.c -o hyphaMPI -lm
```

### 2. Run locally

```bash
./hypha
```

or with MPI:

```bash
mpirun -np 8 ./hyphaMPI
```

### 3. Run a cluster sweep

The Hamilton job script performs an `Ec_h` sweep defined in the selected sweep
config file and writes results to `simulationCases/<CaseNo>/`.

```bash
sbatch runSweepHamilton.sbatch
```

By default it reads `sweep.params`. You can pass a different sweep file:

```bash
sbatch runSweepHamilton.sbatch my-sweep.params
```

For the serial executable workflow (compile without MPI, launch all cases in
parallel), use:

```bash
sbatch runSweepHamilton-serial-128.sbatch
sbatch runSweepHamilton-serial-128.sbatch my-sweep-128.params
```

## Parameter Files

- `default.params`: base parameter values used by the sweep script.
- `sweep.params`: default sweep input for `runSweepHamilton.sbatch` (contains
  `CASE_START`, optional `CASE_END`, and `SWEEP_Ec_h=...` list).
- `sweep-128.params`: default sweep input for
  `runSweepHamilton-serial-128.sbatch` (128 `Ec_h` values, `CASE_START=1506`).

## Post-Processing

`postProcess/` contains extraction and plotting utilities. Typical usage is to
run plotting scripts from inside a case directory containing `log` and/or
`hypha-def-log`, for example:

```bash
python3 ../postProcess/plot_vcm_vs_time.py
python3 ../postProcess/plot_hypha_width_vs_time.py
```

For multi-case analysis (after collecting logs into one directory):

```bash
python3 postProcess/plot_vcm_vs_Ec_H.py --log_dir path/to/logs --pattern "log*"
```

## Notes

- `simulationCases/hypha.c` reads `Ec_h` from a parameter interface
  (`param_double(...)`), so runtime parameter handling must be available in the
  configured Basilisk toolchain.
- Generated binaries, logs, restart files, local Basilisk checkouts, and
  `.comphy-basilisk/` are intentionally not committed.

## License

GPL-3.0. See `LICENSE`.
