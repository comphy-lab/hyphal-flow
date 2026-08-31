# AGENTS.md

This is the operating manual for `comphy-lab/hyphal-flow`.

## Scientific contract

The supported model is `simulationCases/hyphal-flow.c`:

- `f1`: finite-relaxation Oldroyd-B drop;
- `f2`: incompressible finite-strain neo-Hookean Kelvin--Voigt-type outer
  solid (`De_solid` in the practical no-relaxation limit);
- `1 - f1 - f2`: finite-relaxation Oldroyd-B carrier liquid.

`Oh_*` is the Newtonian/solvent viscosity. For each relaxing phase the
implied polymer viscosity is `Ec_* De_*`. The current solver uses one
regularised conformation tensor with phase-weighted modulus and relaxation
rate; it does not provide three independent interface tensors.

Separate VoF tracers keep the drop and solid interfaces distinct. That alone
does not prove a positive carrier-film thickness: inspect the overlap and
phase-volume diagnostics in each run log.

## Supported layout

- `simulationCases/hyphal-flow.c`: only supported simulation entry point.
- `src-local/`: runtime parameter, phase-property, conformation and forcing
  headers.
- `smoke/`: reduced execution inputs. They establish that code paths compile,
  advance, terminate and emit finite diagnostics; they are not verification or
  validation cases.
- `postProcess/`: extraction and visualisation tools.
- `simulationCases/legacy/`, `src-local/legacy/`, `scripts/legacy/`:
  provenance only. Do not route production or smoke runs through them.

## Build and run

Use the root runner so `case.params`, source hashes and resume checks stay
consistent:

```bash
bash runSimulation.sh default.params --dry-run
bash runSimulation.sh default.params --output-root /path/to/run-root
```

`--resume` is allowed only when the stored source closure and parameter hash
match. Never overwrite or silently reuse a non-empty case directory.

Sweep expansion is owned by one runner:

```bash
bash runParameterSweep.sh sweep.params --dry-run
bash runParameterSweep.sh sweep.params --output-root /path/to/run-root
```

For managed compute, output belongs in the reserved run root, not in this
checkout. The Hamilton wrappers require `RUN_ROOT` and a site-provided
Basilisk/MPI environment; they never install software inside a job.

## Smoke evidence

The CoMPhy smoke driver runs four reduced physical cases, a two-case sweep,
the full case with two MPI ranks and an input-matched restart:

```bash
bash runSmokeTests.sh --output-root /path/to/smoke-root
```

Record the exact source commit, parameters and output root with any receipt.
A successful smoke establishes execution only. Numerical verification needs an
independent exact/manufactured comparator; physical validation needs
independent data.

## Source and documentation discipline

- Ground Basilisk API and numerical changes in current upstream source before
  editing solver headers.
- Keep `README.md` and generated docs aligned with actual paths and commands.
- Build docs with `.github/scripts/build.sh --force-rebuild` and inspect the
  generated homepage/tree before committing.
- `AGENTS.md` is committed. Local `CLAUDE.md` contains only `@AGENTS.md` and is
  gitignored.
- Never commit local Basilisk trees, `.docker_mode`, `.comphy-basilisk`, build
  output, case output, restart files or raw simulation data.

## Research documentation boundary

<!-- documentation-boundary-v1 -->

Treat `README.md`, generated docs, figures, reports and collaborator-facing
material as public candidates. Use precise, evidence-proportionate prose and
keep provisional values, failed hypotheses, runtime paths, job details and
debugging lore on internal project/run surfaces until explicitly approved for
a named target.

Every figure task, including exploratory and diagnostic output, must follow
the CoMPhy `publication-plots` workflow. Parallel frame rendering must use
process-safe mathtext rather than `text.usetex=True`.
