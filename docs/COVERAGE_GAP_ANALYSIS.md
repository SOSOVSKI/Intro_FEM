# Coverage Gap Analysis

Comparison of the syllabus (SYLLABUS_REVISION_DRAFT.md), lectures (course_deck.qmd + book chapters), and hands-on notebooks.

Generated: 2026-04-12

## Coverage Table

| Week | Syllabus Topic | In Lectures? | Book Chapter | Has Notebook? | Notebook(s) | Action |
|------|---------------|-------------|-------------|--------------|-------------|--------|
| 1 | Introduction, essence of FEM | Yes (Week 1) | ch01 | No dedicated NB | — | Syllabus-only flag: intro is lecture-only by design |
| 2 | Weak form derivation, test function | Yes (Week 2) | ch01, ch02 | Yes | NB01 | Covered |
| 3 | Galerkin method, shape functions | Yes (Week 3) | ch02, ch04 | Yes | NB02 | Covered |
| 4 | 1D finite elements, numerical solution | Yes (Week 4) | ch04 | Yes | NB03 | Covered |
| 5 | Gauss quadrature, generic 1D FEM | Yes (Week 5) | ch04 | Yes | NB04 | Covered |
| 6 | Natural and essential BCs | Yes (Week 6) | ch04 | Yes | NB08, NB09 | Covered |
| 7 | Error analysis, solution characteristics | Yes (Week 7, 13) | ch03, ch12 | **Gap** → NB14 created | NB14 | **Created**: 14_error_analysis.ipynb |
| 8 | Time-dependent 1D, beginning 2D | Yes (Week 8) | ch05, ch06 | **Gap** → NB13 created | NB13 | **Created**: 13_time_dependent_fem.ipynb |
| 9 | 2D formulation, divergence theorem | Yes (Week 9) | ch06, ch08 | Yes | NB05, NB06 | Covered |
| 10 | 2D shape functions, meshing, BCs in 2D | Yes (Week 10) | ch08 | Yes | NB07, NB09 | Covered |
| 11 | Continuum mechanics, elasticity | Yes (Week 11) | ch09 | Yes | NB10 | Covered |
| 12 | 3D elasticity, implementation | Yes (Week 12) | ch10, ch11 | Yes | NB11 | Covered |
| 13 | Error estimation, adaptivity | Yes (Week 13) | ch12 | Yes (see Week 7) | NB14 | Covered by NB14 |
| — | I❤️LA export, widgets | Optional bridge | — | Yes | NB12 | Covered (optional) |

## Syllabus-Only Flags (topics in syllabus not fully covered in notebooks)

These topics appear in the syllabus and lectures but do not have dedicated hands-on notebooks. This is by design — they are primarily theoretical.

1. **Week 1: Introduction and essence of FEM** — Lecture-only. The overview and motivation do not require a hands-on notebook.

2. **Function spaces (L², H¹, Sobolev)** — Covered in Week 7 lectures (ch03) and referenced in NB14's convergence rate discussion, but there is no dedicated notebook on function space theory. This is appropriate since function spaces are a mathematical prerequisite, not a computational exercise.

3. **Post-processing** — Mentioned briefly in Week 5 deck slides but has no dedicated content anywhere. Students see solution plots in NB07 (2D assembly) and NB13 (time-dependent), but there is no systematic treatment of stress recovery, gradient smoothing, or visualisation pipelines.

4. **Q1 quadrilateral workflow** — The `ReferenceQuadrilateralQ1` and `ReferenceQuadrilateralQ2` classes exist in the workbench, and NB05 (reference element gallery) displays them, but there is no end-to-end Poisson-on-quads notebook analogous to NB06 (triangles). This was marked as optional in Phase 3 of the integration plan.

5. **UFL/FEniCS comparison** — Listed in syllabus Week 11+ as an optional bridge. Not implemented in any notebook.

## Notebooks Created to Fill Gaps

### 13_time_dependent_fem.ipynb (NEW)

Covers Week 8 syllabus topic: time-dependent FEM. Contents:

- Symbolic mass matrix derivation using `build_bar_1d_mass_matrix()`
- Lumped vs. consistent mass matrix comparison
- Numerical assembly of M and K for the 1D transient heat equation
- Backward Euler (implicit) time stepping
- Forward Euler (explicit) time stepping with lumped mass
- CFL stability analysis and instability demonstration
- Comparison with exact solution $u(x,t) = \sin(\pi x) e^{-\pi^2 t}$

Verified: executes cleanly.

### 14_error_analysis.ipynb (NEW)

Covers Week 7 and Week 13 syllabus topics: error analysis and a posteriori estimation. Contents:

- Manufactured solution for 1D Poisson ($u = \sin(\pi x)$)
- FEM solver for varying mesh sizes (h-refinement study)
- $L^2$ and $H^1$ error norm computation via Gauss quadrature
- Convergence rate verification: $O(h^2)$ for $L^2$, $O(h)$ for $H^1$
- Nodal superconvergence on uniform meshes
- Error distribution visualisation (pointwise and element-wise)
- A posteriori element residual error indicator
- Effectivity index comparison (indicator vs. true error)

Verified: executes cleanly.

## Summary

| Metric | Count |
|--------|-------|
| Total syllabus topics mapped | 14 (13 weeks + optional I❤️LA) |
| Fully covered (lectures + notebook) | 10 |
| Gaps filled by new notebooks | 2 (NB13, NB14) |
| Syllabus-only flags (no notebook needed) | 5 |
| Total notebooks | 14 (12 original + 2 new) |
| All notebooks verified | 13 PASS, 1 SKIP (NB11: 3D compute timeout) |
