# Completion identity on geomech

A WSEED names a parent cell. Its selected completion may instead belong to an
LGR and have a different grid-local index. Store both identities: retain the
parent for activation and use (completion global index, LGR number) for property
lookup and mapping into MPI-local completion lists.

At hnil/geomech 7324e1c, localSeedCells indexes the filtered local list with a
full-schedule ordinal, and property/filter-cake lookup uses the WSEED parent
index without a grid qualifier. These paths can miss the selected completion or
use another grid's identically numbered completion. Existing grid-qualified
lookup in opm-common does not resolve the missing identity in Fracture itself.

The port requires identity at every Fracture::init call, keeps current OpmLog
logging, and leaves fracture physics and the one-fracture-per-WSEED policy intact.
It does not import the old experimental grid history.

## Validation on 23 October 2026

The native identity test compiled and linked against clean geomech common/grid,
and passed with mpiexec -n 1. Cases cover colliding indices in different grids,
parent/completion differences, reordered local completions and absent/nonlocal
completions. Fracture.cpp, FractureModel.cpp and the instantiated simulator
translation unit passed syntax checks. Compilation was sequential and capped
at 8 GiB per process. Full simulator linking, nine regression runs and fresh
SIMPLE/WELLREF comparisons remain required before merging.

The old line's 440-to-zero warning result is historical evidence, not a measured
result for this port. Filter-cake concentration/zero-thickness behavior, nested
LGR and schedule changes remain separate investigations.
