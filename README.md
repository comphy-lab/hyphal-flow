# hyphal-flow

[![Basilisk](https://img.shields.io/badge/Basilisk-C-blue)](http://basilisk.fr)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

Axisymmetric Basilisk simulations of drop transport in a single fungal
hypha branch, with three non-coalescing phases and log-conformation
viscoelasticity.

## Overview

A drop moves through a thin liquid film that lines a hyphal tube filled
with cytoplasm. The drop and film occupy separate volume-fraction
fields and do not coalesce. Polymeric stress is evolved with the
MultiRheoFlow log-conformation solver, so each phase can be Newtonian,
Oldroyd-B, or Kelvin--Voigt according to its modulus and relaxation
time.

## Physics

### Problem

The tube is treated as axisymmetric. The drop is an ellipse of
semi-axes $(1.5, 1)$ in capillary units. The film/cytoplasm interface
is a tanh profile whose far-field radius is set by `hf`. Streamwise
driving is either a Bond body force on drop and film (periodic tube)
or an imposed pressure difference between open ends.

### Governing equations

Incompressible momentum with solvent viscosity $\mu_s$ and polymeric
stress $\mathbf{T}$:

$$\nabla \cdot \mathbf{u} = 0$$

$$\rho\left(\partial_t\mathbf{u} + \mathbf{u}\cdot\nabla\mathbf{u}\right)
= -\nabla p + \nabla\cdot(2\mu_s\mathbf{D}) + \nabla\cdot\mathbf{T}
+ \rho\mathbf{a}$$

The conformation tensor $\mathbf{A}$ is integrated in log form
$\Psi = \log\mathbf{A}$. For Oldroyd-B,

$$\mathbf{T} = G_p(\mathbf{A} - \mathbf{I})$$

with relaxation time $\lambda$. A very large $\lambda$ recovers
Kelvin--Voigt elasticity.

### Dimensionless groups

| Parameter | Definition | Role |
|-----------|------------|------|
| $Oh$ | $\mu/\sqrt{\rho\sigma R}$ | solvent viscosity of each phase |
| $Ec$ | $G_p R/\sigma$ | elasto-capillary number |
| $De$ | $\lambda/\sqrt{\rho R^3/\sigma}$ | Deborah number |
| $Bo$ | body-force strength in capillary units | Bond driving |
| $hf$ | far-field film radius / drop radius | constriction geometry |

Capillary scaling uses drop radius $R$ and interfacial tension
$\sigma$. Phase labels are drop (`d`), film (`h`/`f`), and cytoplasm
(`c`).

## Basilisk (Required)

This repository targets CoMPhy Basilisk `v2026-08-30`. First-time
install (or reinstall):

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

## Repository Structure

```
├── src-local/ - Project headers and runtime parameter API
│   ├── parse_params.h - Low-level key=value loader
│   ├── params.h - Typed accessors (`param_int`, `param_double`, ...)
│   ├── parse_params.sh - Shared shell helper for runners
│   ├── log-conform-viscoelastic-scalar-2D.h - MultiRheoFlow log-conform (2D/axi)
│   ├── log-conform-viscoelastic-scalar-3D.h - MultiRheoFlow log-conform (3D)
│   ├── log-conform-viscoelastic.h - MultiRheoFlow tensor log-conform
│   ├── three-phase-nonCoalescing-VE.h - Three-phase coupling to Gp and lambda
│   └── reduced-three-phase-nonCoalescing.h - Reduced gravity on f1 and f2
├── simulationCases/ - Solver entry point and generated case directories
│   ├── hypha.c - Unified axisymmetric hypha solver
│   ├── Makefile - Optional qcc include of Basilisk Makefile.defs
│   └── runCodesInParallel.sh - Compatibility wrapper; use runSimulation.sh
├── postProcess/ - Extraction, plots, and snapshot video
│   ├── getData-elastic-nonCoalescence.c - Uniform-grid field sampler
│   ├── getFacet-threePhase.c - Interface facet extraction
│   ├── Video-hyphal-generic.py - Parallel snapshot video
│   ├── plot_vcm_vs_time.py - Drop centre-of-mass velocity
│   ├── plot_vcm_vs_time-2.py - Alternate CoM velocity plot
│   ├── plot_hypha_width_vs_time.py - Film height from hypha-def-log
│   └── plot_vcm_vs_Ec_H.py - Multi-case Ec_h comparison
├── legacy/ - Superseded Cartesian capillary solvers and old VE headers
├── .github/ - Docs CI, issue templates, and generated HTML
├── runSimulation.sh - Single-case compile and run
├── runParameterSweep.sh - Cartesian product of SWEEP_* keys
├── default.params - Base runtime parameters
├── sweep.params - Default Ec_h sweep
├── sweep-128.params - Serial-128 sweep list
├── runSingleCaseHamilton.sbatch - Site batch wrapper for one MPI case
├── runSweepHamilton.sbatch - Site batch wrapper for an MPI sweep
├── runSweepHamilton-serial-128.sbatch - Site batch wrapper for serial-compiled sweep
├── AGENTS.md - Repository operating manual
└── LICENSE - GPL-3.0
```

## Running the Code

On a shared-HPC login node, use the repository batch wrappers for
compilation and execution. Login nodes are for editing, status, and
submission only.

### Scripts (recommended)

```bash
bash runSimulation.sh default.params
bash runSimulation.sh default.params --mpi --CPUs 8
bash runParameterSweep.sh sweep.params --dry-run
bash runParameterSweep.sh sweep.params
```

Each case runs in `simulationCases/<CaseNo>/` with `case.params` as
`argv[1]`.

### Manual compilation

These commands are for a local workstation or an allocated compute node.

**Serial / OpenMP:**

```bash
qcc -O2 -Wall -disable-dimensions -I$(pwd)/src-local simulationCases/hypha.c -o hypha -lm
./hypha case.params
```

**MPI:**

```bash
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions -I$(pwd)/src-local simulationCases/hypha.c -o hyphaMPI -lm
mpirun -np 8 ./hyphaMPI case.params
```

### Parameter file

| Key | Meaning | Default |
|-----|---------|---------|
| `CaseNo` | Output directory `simulationCases/<CaseNo>/` | `1000` |
| `MAXlevel` | Maximum wavelet level | `12` (`13` in `default.params`) |
| `tmax` | End time | `200` |
| `DriveMode` | `bond` or `pressure` | `bond` |
| `Bond` | Streamwise body force (`bond` mode) | `1` |
| `PL`, `PR` | End pressures (`pressure` mode) | `1`, `0` |
| `Ohd`, `Ohf`, `Ohc` | Solvent Ohnesorge numbers | `1`, `1`, `0.01` |
| `Ec_d`, `Ec_h`, `Ec_c` | Elasto-capillary numbers | `0`, `0`, `0` |
| `De_d`, `De_h`, `De_c` | Deborah numbers | `0`, `1e30`, `0` |
| `hf` / `hr` | Far-field film radius | `0.90` |
| `Ldomain` | Domain length | `80` |

A large `De_h` is a Kelvin--Voigt film. Set `NanDetector=1` to abort on
the first non-finite field value.

## Post-Processing

Run plotting scripts from a case directory that contains `log` and/or
`hypha-def-log`:

```bash
python3 ../postProcess/plot_vcm_vs_time.py
python3 ../postProcess/plot_hypha_width_vs_time.py
```

Snapshot video (compiles helpers once, then frames in parallel):

```bash
python3 postProcess/Video-hyphal-generic.py --case-dir simulationCases/1000 --cpus 4
```

## Authors

- **Vatsal Sanjay** (Durham University), [vatsal.sanjay@comphy-lab.org](mailto:vatsal.sanjay@comphy-lab.org)

Log-conformation headers are imported from
[MultiRheoFlow](https://github.com/comphy-lab/MultiRheoFlow).

## License

GNU General Public License v3.0. See [LICENSE](LICENSE).
