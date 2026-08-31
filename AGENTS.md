# AGENTS.md

This is the operating manual for `comphy-lab/hyphal-flow`.

## Scientific contract

The supported model is `simulationCases/hyphal-flow.c`:

- `f1`: finite-relaxation Oldroyd-B drop;
- `f2`: incompressible finite-strain neo-Hookean Kelvin--Voigt-type outer
  solid (`De_solid` in the practical no-relaxation limit);
- `1 - f1 - f2`: finite-relaxation Oldroyd-B carrier liquid.

`Oh_*` is the Newtonian/solvent viscosity. For each relaxing phase the
implied polymer viscosity is `Ec_* De_*`. Conformation stress comes from
MultiRheoFlow `log-conform-viscoelastic.h`, included before
`three-phase-rheology.h`. The solver uses one regularised conformation
tensor with phase-weighted modulus and modulus-weighted relaxation rate; it
does not provide three independent interface tensors.

Separate VoF tracers keep the drop and solid interfaces distinct. That
alone does not prove a positive carrier-film thickness: inspect the overlap
and phase-volume diagnostics in each run log.

## Layout (CoMPhy standard)

- `src-local/`: project headers and the C/shell parameter parsers.
- `simulationCases/`: `hyphal-flow.c` and generated `simulationCases/<CaseNo>/`.
- `postProcess/`: extraction, plotting and snapshot video.
- Root: `runSimulation.sh`, `runParameterSweep.sh`, `default.params`,
  `sweep.params`, site batch wrappers, and `Makefile`.

Reduced compile/run checks use `runSimulation.sh default.params`
(optionally `--compile-only`) with an external `--output-root`.

## Build and run

```bash
bash runSimulation.sh default.params --dry-run
bash runSimulation.sh default.params
bash runSimulation.sh default.params --mpi --CPUs 8
bash runParameterSweep.sh sweep.params --dry-run
bash runParameterSweep.sh sweep.params
```

`--resume` is allowed only when the stored source closure and parameter hash
match. Never overwrite or silently reuse a non-empty case directory.

For managed compute, output belongs in the reserved run root, not in this
checkout. The Hamilton wrappers require `RUN_ROOT` and a site-provided
Basilisk/MPI environment; they never install software inside a job.

## Parameters

- C: `params_init_from_argv` in `hyphal-flow.c`, typed reads in `src-local/params.h`.
- Shell: one helper, `src-local/parse_params.sh`.
- Base file: `default.params`. Sweep file: `sweep.params`.

## Source and documentation discipline

- Ground Basilisk API and numerical changes in current upstream source before
  editing solver headers.
- Keep `README.md` and generated docs aligned with actual paths and commands.
- Repository trees in the README use a plain fenced block (no language tag)
  and `path - description` lines.
- Build docs with `.github/scripts/build.sh --force-rebuild` and inspect the
  generated homepage tree before committing.
- `AGENTS.md` is committed. Local `CLAUDE.md` contains only `@AGENTS.md` and
  is gitignored. `.coderabbit.yaml` is committed.
- Never commit local Basilisk trees, `.docker_mode`, `.comphy-basilisk`,
  build output, case output, restart files or raw simulation data.

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
