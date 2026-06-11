# Pure mechanics discretization: stability choices and high-aspect-ratio robustness

This note summarizes a review (June 2026) of the linear elasticity discretizations in
this module — the first-order virtual element method (VEM) in `vem.cpp` and the Q1
finite element path — with focus on why **all stability variants degrade on cells
with high aspect ratio**, and what was changed about it.

## Structure of the discretizations

The VEM element stiffness follows Gain et al. (2014), DOI:10.1016/j.cma.2014.05.005:

```
Ke = Wc^T D Wc * |E|   +   (I-P)^T S (I-P)
     `consistency'         `stabilization'
```

The consistency term (`EWcDWct` in `final_assembly`) is *exact for linear
displacement fields by construction, on any polyhedron, at any aspect ratio*.
All aspect-ratio trouble lives in the stabilization matrix `S`
(`getStabilityMatrix`), which must give the non-polynomial part of the element
kernel a stiffness comparable to the genuine deformation modes.

The FEM path (Q1, trilinear hexahedra, full 2x2x2 Gauss integration; see
`vemutils.cpp` and `elasticity_solver_impl.hpp`/opm-upscaling) is a correct standard
implementation — its aspect-ratio problems are inherent to the element, not bugs.

## Why each variant degrades for high aspect ratios

On a cell with edge lengths `(hx, hy, hz)` the diagonal of the consistency matrix
scales like `E*|E|/h_d^2` for a dof in direction `d`: the genuine deformation modes
span a stiffness range of `(max h / min h)^2`.  Any *isotropic* stabilization must
either under-stabilize the stiff directions (spurious, nearly-zero-energy modes →
hourglassing, solver breakdown) or over-stabilize the soft ones (locking →
overstiff response).

- **SIMPLE** (Gain 2014): `alpha = |E| tr(D) / tr(N^T N)`, scalar x identity.
  `tr(N^T N)` is dominated by the *long* axis → alpha too small → spurious soft
  modes in the short direction.
- **HARMONIC** (Andersen 2017): `alpha ~ |E| tr(D) * sum(1/diag(N^T N))`, dominated
  by the *short* axis → larger alpha, fails later but is still one isotropic number.
- **D_RECIPE**: per-dof diagonal of the consistency matrix — the right,
  naturally anisotropic idea — but the floor `max(tr(D)/9 * cbrt(|E|), .)` is
  isotropic and based on the *geometric-mean* length.  On a high-aspect cell that
  floor clobbers the soft-direction entries and re-introduces locking.  (It also
  used `cbrt` even in 2D.)
- **EXPERIMENTAL**: constant, geometry-independent — not robust by design.
- **Q1 FEM, full integration**: classic *shear locking* in bending (parasitic shear
  energy grows ~ aspect^2) and *volumetric locking* as `nu -> 0.5`.

## What was added (June 2026)

New `StabilityChoice` values in `vem.hpp`:

- **`ANISO_DIAG` (8)** — recommended default for general (polyhedral) grids.
  D_RECIPE's consistency diagonal, but with an *anisotropic* floor
  `tr(D)/9 * |E| / h_d^2` per coordinate direction `d`, where `h_d` are
  per-direction length scales from the second moments of the element's corner
  cloud (`compute_axis_lengths`; equals the edge lengths for a box).  Reduces to
  D_RECIPE behavior on isotropic cells; tracks the aspect ratio automatically.
- **`ANISO_HARMONIC` (9)** — purely geometric anisotropic diagonal
  `tr(D)/9 * |E| / h_d^2` (no consistency diagonal); the anisotropic
  generalization of HARMONIC.
- **`FEM` (7, now actually implemented end-to-end)** — standard Q1 elements on every
  cell where they are valid (topological hexahedra with positive Jacobian at all
  quadrature points, in 2D quadrilaterals); VEM with ANISO_DIAG stabilization on
  all other cells.  Works both through the Dune/CpGrid path (`vemutils.cpp`,
  which additionally restricts to "nice" corner-point cells) and through the pure
  array-based API (`assemble_stiffness_matrix_3D/2D` with hex detection in
  `hex_corner_ordering`).

  The hex classification in `hex_corner_ordering` requires: 8 corners, 6 quad
  faces, 12 distinct edges *each shared by exactly two faces* (closed surface),
  *every vertex of degree 3*, a consistent bottom/top pairing, and positive
  trilinear Jacobians at all 8 Gauss points *and* all 8 corners.  The
  combinatorial part is watertight: the cube is the unique 3-regular
  quadrangulation of the sphere with 6 faces, so no non-hex face complex can
  pass (the counting conditions alone — 8 corners/6 quads/12 edges — would not
  exclude non-manifold complexes from corrupt corner-point cells).  The
  geometric part (Jacobian sign at finitely many points) is the standard
  practical criterion; it covers everything the 2x2x2 quadrature evaluates, but
  does not prove global injectivity of the trilinear map — strongly twisted
  cells that pass it would also be mis-integrated by any standard Q1 code.
  Cells failing any check fall back to VEM.

  The standalone Q1 element (`stiffness_matrix_fem_hex_3D`) is verified to
  reproduce the legacy opm-upscaling `Opm::Elasticity` Q1 element stiffness
  entry-by-entry (rel. diff ~1e-15, ~1e-12 at cell aspect 1000;
  `Unit.PatchRecovery` test).  The differences are plumbing, not discretization:
  the standalone version needs no Dune grid (array API, used by the unit tests),
  derives the corner ordering from face topology, uses |detJ| (robust to
  inverted orientation), and adds the optional B-bar modification; the
  Dune/CpGrid assembly path continues to use the opm-upscaling element (with the
  same B-bar two-pass added for FEM_BBAR).
- **`FEM_BBAR` (10)** — as FEM, but with the B-bar / mean-dilatation modification
  (Hughes 1980): the dilatational part of the strain-displacement matrix is
  replaced by its element average.  Cures *volumetric* locking (important for
  nu -> 0.5, e.g. undrained response); it does **not** cure bending/shear locking
  (see "further options" below).

Other changes:

- The VEM fallback for non-hex cells under FEM choices now uses ANISO_DIAG
  (was D_RECIPE).
- **Diagonal (Jacobi) scaling** of the assembled system, `A -> S A S` with
  `S = diag(1/sqrt(A_ii))`: JSON option `"mech_diagonal_scaling": true` under
  `fractureparam` (plumbed `eclgeomechmodel.hh -> VemElasticitySolver`), or
  `vem::diagonal_scale_system()` for the array-based API.  High-aspect cells make
  the diagonal magnitudes span orders of magnitude; equilibration helps the
  iterative solver and the AMG setup independently of which element/stability is
  used.  The solution is recovered as `x = S y` (handled inside
  `VemElasticitySolver::solve()`).

## Test

`examples/test_aspect_ratio.cpp` (ctest target `Test.AspectRatio`) is standalone
(array-based API only, dense in-file solver) and runs:

1. 3D + 2D patch tests on stretched grids (hard failure if a linear field is not
   reproduced) — also with diagonal scaling on.
2. Single-hex eigenvalue check at aspect 1/10/100: exactly 6 rigid-body zero modes,
   no spurious modes (hard failure for the new variants, warning for legacy ones).
3. Cantilever bending sweep at fixed element count and growing element aspect
   ratio, compared to the Timoshenko solution (diagnostic table: locking shows as
   ratios falling below 1).
4. Volumetric locking sweep: the same cantilever at nu=0.4999 across aspect
   ratios.  Hard checks: FEM_BBAR must be softer than FEM at every aspect ratio
   and must essentially cure the locking at aspect 1.

Measured (tip deflection / Timoshenko, June 2026): at nu=0.4999 **all VEM
variants collapse to ~0.03 already at aspect 1** — not only does the consistency
D stiffen, the stabilization scalings are all proportional to `tr(D)`, which is
dominated by lambda as nu -> 0.5, so the kernel (bending/hourglass) modes become
nearly rigid.  Plain Q1 gives 0.105 at aspect 1 (classic volumetric locking);
FEM_BBAR gives 0.915 (cured), degrading with aspect only through the remaining
shear locking.  Consequence: for near-incompressible response (undrained
poroelasticity, volumetric thermal loads) FEM_BBAR is currently the only good
choice.  The natural VEM fix — not yet implemented — is to scale the
stabilization with the *shear* modulus only (e.g. replace `tr(D)` by `~18 mu`,
or use the deviatoric part of D), since the VEM consistency term applies the
volumetric constraint only cell-wise (one constraint per cell, B-bar-like) and
is therefore not the locking bottleneck.

Note also (mesh refinement at fixed aspect 1, cantilever L=5): all variants
converge to 1 (e.g. SIMPLE 0.68 -> 0.94, HARMONIC 0.27 -> 0.74, D_RECIPE
0.46 -> 0.87, FEM 0.87 -> 0.98 for 10x2x2 -> 30x6x6).  The large spread on the
coarse mesh is pre-asymptotic: with 1-2 linear elements through the thickness
most of the bending energy lives in the projection kernel, so the deflection is
controlled directly by the magnitude of S (or by parasitic shear for Q1) — the
variants differ in that magnitude by factors of a few, which is exactly the
spread seen at aspect 1.

## Load regularization for under-resolved fronts (June 2026)

Sharp, cell-wise constant temperature/pressure fields (e.g. an early-time cold-water
plume) produce eigenstrain jumps whose stress extrema are grid dependent — and at
the *edges* of the cooled region the elastic stress is log-singular, so the
pointwise minimum stress does not converge under refinement at all.  The
`"smooth_force"` option (one-ring averaging) has a grid-dependent smoothing radius
and therefore does not remove this dependence.

New option `"smooth_force_length"` (meters, under `fractureparam`): smooths the
mechanics potential force with `vem::diffuseCellVector`, a conservative TPFA
Helmholtz filter `(V + l^2 L) u = V u0` whose kernel decays exponentially with the
given *physical* length.  Choose the length from physics, e.g. the conductive front
width `sqrt(thermal diffusivity * time)`, so the regularized stress extrema are
grid independent and physically meaningful.  Takes precedence over `smooth_force`.
In parallel the filter falls back to explicit sub-steps with per-step
owner-to-all communication (exact, but capped at 1000 steps with a warning).

For *evaluating* stress, two patch-recovery options are wired in (both verified
linear-exact also on aspect-ratio-1000 cells, `Unit.PatchRecovery` test):

- `"patch_recovery_stress": "true"` (under `fractureparam`): replaces the raw
  cell-wise stress (`linstress_`, hence everything downstream of
  `problem.stress()`, including output and fracture sampling) by its patch
  recovery (`vem::patchRecoveryCells6`).
- `"reservoir": { "interpolate_stress": "true" }`: the fracture model evaluates
  the reservoir stress at each fracture cell centroid by *trilinear interpolation
  of patch-recovered nodal stress* (`vem::patchRecoveryNodal6` + CpGrid
  `geometry().local()`), giving a traction that is continuous along the fracture
  surface instead of jumping at each reservoir cell boundary (relevant when
  reservoir cells are large compared to the fracture mesh).  Falls back to the
  containing-cell value on degenerate cells where the trilinear inversion fails.

Even with recovery, a criterion based on pointwise interface-cell minima remains
non-convergent under refinement at a sharp front edge — combine with
`smooth_force_length` or the analytic plateau `E alpha dT/(1-nu)` inside the
cooled zone.

## Recommended configuration (high-aspect grids, thermal injection)

```json
"fractureparam": {
    "vem_stability_choice": "9",          // ANISO_HARMONIC: flat across aspect ratios
    "mech_diagonal_scaling": "true",      // equilibrates the system for AMG
    "smooth_force_length": "<front width in m, e.g. sqrt(kappa*t)>",
    "patch_recovery_stress": "true",
    "reservoir": { "interpolate_stress": "true" }
}
```

Reasoning: the dominant deformation regime for an injection-cooled layer is
laterally constrained (oedometric), where both VEM and Q1 are consistent — the
element choice matters less than (a) not over-stabilizing the thin-cell direction
(HARMONIC and D_RECIPE lock; ANISO_HARMONIC does not) and (b) regularizing the
under-resolved front load and the sampled stress (the three options above).
Pure VEM with ANISO_HARMONIC treats every cell — including degenerate
corner-point cells — uniformly, avoiding the FEM/VEM stabilization interface that
the combined FEM choice introduces exactly where pinch-outs cluster.  Use
FEM/FEM_BBAR (7/10) as a cross-check on clean Cartesian decks, and FEM_BBAR if
materials with nu >~ 0.45 appear.  At nu ~ 0.3, volumetric locking is not a
concern in the sequential coupling (the mech element always sees drained moduli).

## Further options (not implemented)

## Further options (not implemented)

- **EAS / incompatible modes (Wilson Q6, Taylor QM6)** for the Q1 path: the
  standard cure for *bending/shear* locking of elongated hexes (~9 internal dofs,
  condensed per element).  This is the highest-payoff next step for logically
  Cartesian grids; B-bar alone leaves the shear locking visible in the bending
  table of the test.
- **Principal-axis length scales**: `compute_axis_lengths` uses coordinate-axis
  second moments, which is right for (logically) Cartesian cells; for strongly
  sheared/rotated cells use the eigenvectors of the full inertia tensor and apply
  the floor in the principal frame.
- **Spectral stabilization**: build S from the eigenstructure of the consistency
  matrix so every kernel mode gets stiffness comparable to the smallest genuine
  deformation mode — the most robust general-mesh option, at the cost of a small
  eigensolve per element.
- **Stabilization-free VEM** (Berrone/Borio/Beirao da Veiga line of work): enlarge
  the strain projection space so no artificial S is needed at all.
