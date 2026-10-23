# Mechanical constraint validation on geomech

The clean line still checks direction counts and a coordinate-free rotational
heuristic locally. It also skips this check for separate mechanics grids. A
locally unsupported rank can belong to a globally supported body, while nodes
that pass direction counts can still leave rigid-body rotation unconstrained.

This port gathers connectivity and displacement masks by global mechanics-grid
vertex identity, checks the rank of the six rigid-body modes for each connected
body, and broadcasts the result. Deficient bodies fail before the solve; a
separate supported body cannot hide an unsupported component. It replaces the
existing heuristic in the VEM initialization path and includes separate grids.

## Validation on 23 October 2026

The native suite compiled and linked against clean geomech common/grid and
passed sequentially on one, two and four MPI ranks. It covers free bodies,
fixed-point/collinear supports, distributed supports, duplicate overlap masks,
disconnected bodies, translation/scaling and collective bad-index rejection.
The instantiated simulator translation unit passed a syntax check. Compilation
was sequential with an 8 GiB process ceiling; MPI checks used a 4 GiB ceiling.
Full simulator linking, nine regressions, negative input simulations and fresh
SIMPLE/WELLREF comparisons remain required before merging.

Six whole-body modes are a necessary check, not proof of stiffness definiteness,
conditioning or absence of point/line-hinge mechanisms. Large-grid root memory,
scheduled BC changes, TPSA and conflicting prescribed values are outside scope.
