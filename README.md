# Hyphal Flow

[Hyphal Flow](https://github.com/comphy-lab/hyphal-flow) is a CoMPhy Lab
Basilisk project for an axisymmetric viscoelastic drop transported through a
liquid-filled, deformable hypha.

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

The current implementation uses one regularised conformation tensor with a
phase-weighted modulus and modulus-weighted relaxation rate. It is a practical
diffuse-interface approximation, not three independent material-history
tensors.

## Requirements

- a current [Basilisk](https://github.com/comphy-lab/basilisk-C) `qcc` on
  `PATH`, or an executable `$BASILISK/qcc`;
- a C compiler and `libm`;
- `mpicc` and `mpirun` for MPI runs;
- Bash, `awk` and `shasum` for the runners;
- Python 3, NumPy and Matplotlib for post-processing; ffmpeg for video output.

## Quick start

Inspect a clean-clone run without writing output:

```bash
bash runSimulation.sh default.params --dry-run
```

Compile and run one case into an explicit output root:

```bash
bash runSimulation.sh default.params --output-root /path/to/run-root
```

The runner creates `<output-root>/<CaseNo>/`, copies `case.params`, records
source and parameter hashes, compiles `simulationCases/hyphal-flow.c`, and
executes the binary with `case.params` as `argv[1]`. It refuses a non-empty
case directory unless `--resume` is supplied and the recorded hashes match.

Parameter sweeps use the same case path:

```bash
bash runParameterSweep.sh sweep.params --dry-run
bash runParameterSweep.sh sweep.params --output-root /path/to/run-root
```

`SWEEP_*` keys form a Cartesian product. `CASE_START` and `CASE_END` must match
the generated count exactly.

## CoMPhy smoke campaign

The reduced smoke inputs exercise a Newtonian control, the outer-solid branch,
the two liquid-rheology branches, and the complete requested model. The driver
also exercises a two-case sweep, two-rank MPI execution and an input-matched
restart:

```bash
bash runSmokeTests.sh --output-root /path/to/smoke-root
```

These are compile-and-run smoke checks. They establish that the selected code
paths advance, terminate and emit finite diagnostics; they are not a
convergence study, numerical verification or physical validation.

## Repository structure

```
├── simulationCases/ - supported simulation entry points
│   └── hyphal-flow.c - canonical axisymmetric three-rheology case
├── src-local/ - project-specific Basilisk and parameter modules
│   ├── parse_params.h - low-level C key-value parser
│   ├── params.h - typed runtime parameter accessors
│   ├── parse_params.sh - shared safe shell parser
│   ├── three-phase-rheology.h - normalised three-material properties
│   ├── log-conformation-rheology.h - shared constitutive update
│   └── reduced-three-phase-nonCoalescing.h - phase body-force potentials
├── smoke/ - reduced execution-contract parameter files
├── postProcess/ - extraction, plotting and parallel video tools
├── runSimulation.sh - single-case compile, run and resume contract
├── runParameterSweep.sh - Cartesian sweep expansion and execution
├── runSmokeTests.sh - bounded CoMPhy smoke campaign
├── runSingleCaseHamilton.sbatch - thin MPI scheduler wrapper
├── runSweepHamilton.sbatch - thin MPI sweep scheduler wrapper
├── default.params - reference full-rheology configuration
├── sweep.params - reference solid-modulus sweep
└── Makefile - portable canonical compile target
```

Historical cases and headers are retained for provenance but are excluded from
the supported runner and smoke claims.

## Managed compute

Simulation output does not belong in the source checkout. The Hamilton
wrappers require a reserved external `RUN_ROOT` and the site-provided
Basilisk/MPI environment:

```bash
sbatch --export=ALL,RUN_ROOT=/reserved/run/root runSingleCaseHamilton.sbatch default.params
```

They do not install Basilisk or duplicate parameter parsing inside the job.

## Post-processing and documentation

`postProcess/Video-hyphal-generic.py` compiles extraction helpers once, accepts
`--cpus`/`--CPUs`, isolates worker caches and writes deterministic frame names.
Use `--max-frames` and `--skip-video` for reduced checks.

Build the documentation site from the repository root:

```bash
.github/scripts/build.sh --force-rebuild
```

The generated site is written to `.github/docs/` and deployed through GitHub
Pages.

## Licence

GPL-3.0. See `LICENSE`.
