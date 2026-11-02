# Complete owner-vertex stencils for parallel VEM

VEM stiffness and body-force assembly uses every local cell incident to a vertex. A vertex owner must therefore have its complete incident-cell star. Face-neighbor overlap, even with multiple layers, can omit incident cells at a partition junction and change an owned mechanics row.

Geomechanics now defaults to the existing `add-corners=true` and `add-corners-by-vertex=true` dependency options. A collective validator counts unique incident cells through face vertices, including conforming face/edge refinement vertices. It compares the owner's locally available count with the global count accumulated from unique cell owners. NNC faces without geometric vertices do not contribute. An incomplete stencil is rejected collectively before timestepping, and checked again when setting up the parallel VEM solver. Cutting timesteps cannot repair this grid configuration.

The test assembles the production stiffness and nonuniform body-force operators on a 4x4x4 grid. It compares every owned row with a serial reference in physical-coordinate DOF numbering. On four ranks a fixed partition junction with face-only overlap must produce mismatched rows and be rejected; vertex-complete corner overlap must reproduce the serial rows. Both startup and solver checks must agree. One-rank assembly must also reproduce the reference.

```sh
cmake --build build --target test_mechanical_vertex_stencils --parallel 1
ctest --test-dir build -R '^Unit.MechanicalVertexStencils\.' --output-on-failure -j1
```

The supported dependency branches are `hnil/opm-simulators:geomech` and `hnil/opm-grid:geomech`, which already expose vertex-corner selection. The regression covers a Cartesian partition junction; it is not a complete nonconforming-grid, separate-mechanics-grid or arbitrary-rank test matrix. Serial and no-MPI execution have a dedicated CTest registration; no-MPI compilation must be verified separately. No fracture constitutive, propagation, storage, acceptance-band or coupling-policy changes are part of this patch.
