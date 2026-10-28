# Aperture storage across mesh changes

The aperture-storage option remains opt-in. Its old state is now a geometric
checkpoint separate from the indexed state retained for rollback. For current
triangle K_j and checkpoint triangles K_i^n, the old cell volume is

    V_old,j = sum_i area(K_j intersect K_i^n) * w_i^n
    storage_j = (area(K_j) * w_j - V_old,j) / dt.

This preserves the integrated old volume when the new domain covers the old
one. Newly added area contains zero old fluid. If an experimental mesh change
removes part of the old domain, only the overlapping old volume is represented;
the transfer does not silently redistribute removed fluid elsewhere. Any physical
handling of a shrinking domain is a separate model decision.

The checkpoint is immutable during growth and retries. A new accepted-step
checkpoint replaces it in `moveForwardInTime`; rebuilding a trial mesh does not
modify it. The original aperture, pressure and mesh checkpoints still serve
`resetFracture`. Transfers always start from the checkpoint, avoiding cumulative
coarsen/refine interpolation diffusion. The current triangle geometry is cached;
it is checked again even if the cell count is unchanged, covering reordering and
boundary-node movement.

The transfer projects coplanar triangles into an orthonormal fracture plane,
clips convex triangles, and integrates the piecewise-constant old aperture.
It is independent of RegularTrimesh levels, cell references and mapper indices.
Noncoplanar transfers are rejected rather than assigned a projected volume.
Source bounding boxes sorted along the first plane axis reject disjoint pairs.
A geometry change has worst-case quadratic search cost, matching the dense
mechanics matrix's asymptotic order; unchanged geometry reuses cached volumes.

`test_fracture_storage` checks analytic identity, permutation, refinement,
coarsening, partial overlap, genuinely new area, checkpoint replacement,
rollback/retry, rotated and translated geometry, invalid-plane rejection and
zero aperture. It also builds actual RegularTrimesh/Dune grids, including
multiresolution and smoothed growth, and checks integrated and per-cell volumes.
The production coarsening fixture demonstrates that the old indexed formula
misrepresents old volume. MPI1/MPI4 run these local fracture tests independently
on each rank; they do not constitute a distributed fracture-storage test.
The retry unit checks the geometric snapshot, rather than forcing a rejected
simulator timestep through `resetFracture`. Full
one/four-rank WELLREF simulations validate the integrated path separately.

This change does not alter the selected approximate-Newton policy or add the
missing area/dt exact-Jacobian contribution. It also does not establish a full
well/leakoff/accepted-step conservation audit; that needs source/time diagnostics.
