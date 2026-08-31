# Legacy solvers and headers

These files are kept for provenance. They are not compiled by
`runSimulation.sh`.

The production entry point is `simulationCases/hypha.c`, which is
axisymmetric and parameterised. The two capillary files were Cartesian
experiments: one mixed periodic and pressure boundary conditions, the
other added open-end pressure driving and extra abort diagnostics.

`src-local/` here holds the previous in-tree log-conformation and
three-phase headers. Viscoelastic stress is now taken from MultiRheoFlow
and coupled through `src-local/three-phase-nonCoalescing-VE.h`.
