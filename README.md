# Hyphal Flow

[![Basilisk](https://img.shields.io/badge/Basilisk-C-blue)](http://basilisk.fr)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

Axisymmetric Basilisk simulations of a viscoelastic drop transported
through a liquid-filled, deformable hypha.

## Physical model

The canonical case uses three material regions:

- `f1`: a finite-relaxation Oldroyd-B drop;
- `f2`: a soft, incompressible finite-strain neo-Hookean
  Kelvin--Voigt-type outer solid;
- `1 - f1 - f2`: a finite-relaxation Oldroyd-B carrier liquid.

The Newtonian `Oh_*` term is the solvent viscosity. For a relaxing phase,
`Ec_*` is the modulus and `De_*` the relaxation time, giving an implied
polymer viscosity `Ec_* De_*`. The outer solid uses a very large
`De_solid`, so its conformation stress does not relax on the simulated
timescale; `Oh_solid` supplies the viscous dashpot.

Drop and solid interfaces use separate VoF fields. This prevents them from
becoming one tracer, but does not by itself guarantee a resolved liquid film.
Every run therefore logs VoF overlap, phase volumes and phase-weighted stress
maxima.

The constitutive update is the MultiRheoFlow tensor log-conform solver.
One regularised conformation tensor is used, with a phase-weighted modulus
and a modulus-weighted relaxation rate. That is a practical
diffuse-interface approximation, not three independent material-history
tensors.

## Basilisk (required)

This repository targets CoMPhy Basilisk `v2026-08-30`. First-time install
(or reinstall):

```bash
curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/v2026-08-30/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=v2026-08-30 --hard
```

Subsequent runs (reuse the same checkout if the ref matches):

```bash
curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/v2026-08-30/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=v2026-08-30
```

When a newer stable release is available, replace `v2026-08-30` in both
the script URL and `--ref` with the same
[release tag](https://github.com/comphy-lab/basilisk-C/releases).

The runner also accepts a current `qcc` on `PATH`, or `$BASILISK/qcc`, and
detects the current `set_prolongation()` API versus the older scalar-`dirty`
API.

## Repository structure

```
├── src-local/ - project headers and runtime parameter API
│   ├── parse_params.h - low-level key=value loader
│   ├── params.h - typed accessors (`param_int`, `param_double`, ...)
│   ├── parse_params.sh - shared shell helper for runners
│   ├── log-conform-viscoelastic.h - MultiRheoFlow tensor log-conform
│   ├── three-phase-rheology.h - normalised three-material properties
│   └── reduced-three-phase-nonCoalescing.h - reduced gravity on f1 and f2
├── simulationCases/ - solver entry point and generated case directories
│   └── hyphal-flow.c - canonical axisymmetric three-rheology case
├── postProcess/ - extraction, plotting and parallel video
│   ├── getData-elastic-nonCoalescence.c - uniform-grid field sampler
│   ├── getFacet-threePhase.c - interface facet extraction
│   ├── Video-hyphal-generic.py - parallel snapshot video
│   ├── plot_vcm_vs_time.py - drop centre-of-mass velocity
│   ├── plot_vcm_vs_time-2.py - alternate CoM velocity plot
│   ├── plot_hypha_width_vs_time.py - film height from hypha-def-log
│   ├── plot_vcm_vs_Ec_H.py - multi-case Ec comparison
│   └── requirements.txt - Python dependencies for plotting and video
├── .github/ - docs CI, issue templates and generated HTML
├── runSimulation.sh - single-case compile and run
├── runParameterSweep.sh - Cartesian product of SWEEP_* keys
├── default.params - base runtime parameters
├── sweep.params - default parameter sweep
├── runSingleCaseHamilton.sbatch - site batch wrapper for one MPI case
├── runSweepHamilton.sbatch - site batch wrapper for an MPI sweep
├── Makefile - portable compile of hyphal-flow.c
├── AGENTS.md - repository operating manual
└── LICENSE - GPL-3.0
```

## Running the code

On a shared-HPC login node, use the repository batch wrappers for
compilation and execution. Login nodes are for editing, status and
submission only.

### Scripts (recommended)

```bash
bash runSimulation.sh default.params --dry-run
bash runSimulation.sh default.params
bash runSimulation.sh default.params --mpi --CPUs 8
bash runParameterSweep.sh sweep.params --dry-run
bash runParameterSweep.sh sweep.params
```

Each case runs in `simulationCases/<CaseNo>/` unless `--output-root` is set.
The binary is invoked with `case.params` as `argv[1]`.

`SWEEP_*` keys form a Cartesian product. `CASE_START` and `CASE_END` must
match the generated count exactly.

### Manual compilation

These commands are for a local workstation or an allocated compute node.

```bash
qcc -O2 -Wall -disable-dimensions -I$(pwd)/src-local simulationCases/hyphal-flow.c -o hyphal-flow -lm
./hyphal-flow case.params
```

MPI:

```bash
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions -I$(pwd)/src-local simulationCases/hyphal-flow.c -o hyphal-flowMPI -lm
mpirun -np 8 ./hyphal-flowMPI case.params
```

## Managed compute

Simulation output does not belong in the source checkout. The Hamilton
wrappers require a reserved external `RUN_ROOT` and the site-provided
Basilisk/MPI environment:

```bash
sbatch --export=ALL,RUN_ROOT=/reserved/run/root runSingleCaseHamilton.sbatch default.params
```

They do not install Basilisk or duplicate parameter parsing inside the job.

## Post-processing and documentation

`postProcess/Video-hyphal-generic.py` compiles extraction helpers once,
accepts `--cpus`/`--CPUs`, isolates worker caches and writes deterministic
frame names. Use `--max-frames` and `--skip-video` for reduced checks.

Build the documentation site from the repository root:

```bash
.github/scripts/build.sh --force-rebuild
```

The generated site is written to `.github/docs/` and deployed through GitHub
Pages.

## Licence

GPL-3.0. See [LICENSE](LICENSE).
