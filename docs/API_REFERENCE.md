# API Reference

Complete reference of all public functions and classes in Symbolic FEM Workbench.

## symbols.py

### Types

**Domain1D** = `tuple[sp.Symbol, sp.Expr, sp.Expr]`
- A 1D domain tuple: (coordinate_symbol, lower_bound, upper_bound)

**Domain2D** = `tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], tuple[sp.Symbol, sp.Expr, sp.Expr]]`
- A 2D domain as two Domain1D tuples, one per coordinate

### Classes

**ScalarField(name: str, expression: sp.Expr)**
- A thin wrapper around a SymPy scalar field expression. Provides transparent access to the expression's attributes.

### Functions

**make_field_1d(name: str, x: sp.Symbol) -> sp.Expr**
- Create a symbolic 1D scalar field as a SymPy Function application (undefined function). Returns Function(name)(x).

**make_field_2d(name: str, x: sp.Symbol, y: sp.Symbol) -> sp.Expr**
- Create a symbolic 2D scalar field as a SymPy Function application. Returns Function(name)(x, y).

---

## forms.py

### Classes

**DomainIntegral(integrand: sp.Expr, domain: Domain1D)**
- A domain integral ∫ integrand dx over a 1D domain. Call `.as_integral()` to get a SymPy Integral object.

**BoundaryContribution(expr: sp.Expr, x: sp.Symbol, location: sp.Expr, label: str | None = None)**
- A boundary term from integration by parts. Call `.evaluate()` to substitute the coordinate and simplify.

**WeightedResidual(domain_integrals: tuple[DomainIntegral, ...], boundary_terms: tuple[BoundaryContribution, ...] = ())**
- The complete weighted residual before splitting into bilinear/linear forms. Call `.as_expression()` to combine into a single SymPy expression.

**WeakForm(bilinear: sp.Expr, linear: sp.Expr, trial: sp.Expr, test: sp.Expr)**
- The split weak form with bilinear part a(u,v) and linear part F(v).

---

## transforms.py

**weighted_residual(lhs: sp.Expr, rhs: sp.Expr, test: sp.Expr, domain: Domain1D) -> WeightedResidual**
- Construct the weighted residual from a strong-form equation (lhs = rhs). Multiplies residual by test function and integrates.

**integrate_divergence_1d(flux: sp.Expr, test: sp.Expr, x: sp.Symbol, domain: Domain1D) -> Tuple[DomainIntegral, Tuple[BoundaryContribution, BoundaryContribution]]**
- Apply integration by parts to a flux-divergence term in 1D. Returns (domain_integral, (boundary_left, boundary_right)).

**drop_dirichlet_boundary(boundary_term: BoundaryContribution, test: sp.Expr) -> sp.Expr**
- Enforce a homogeneous Dirichlet condition by setting the test function to zero at a boundary.

**apply_neumann_flux(boundary_term: BoundaryContribution, flux_expr: sp.Expr, prescribed_flux: sp.Expr) -> sp.Expr**
- Apply a prescribed Neumann (natural) boundary condition by substituting the known flux value.

**split_linear_weak_form(expr: sp.Expr, trial: sp.Expr, test: sp.Expr) -> WeakForm**
- Split a weak-form expression into bilinear part a(u,v) and linear part F(v).

**grad_2d(expr: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> sp.Matrix**
- Return the 2D gradient as a column vector [du/dx, du/dy]^T.

**pullback_gradient_2d(gradient_ref: sp.Matrix, jacobian: sp.Matrix) -> sp.Matrix**
- Map a reference-space gradient into physical coordinates using J^{-T}.

**substitute_field(expr: sp.Expr, field: sp.Expr, replacement: sp.Expr, *coordinates: sp.Symbol) -> sp.Expr**
- Replace a symbolic field with a finite-element expansion, handling derivatives of all orders.

**substitute_fe(expr: sp.Expr, replacements: Mapping[sp.Expr, sp.Expr], *coordinates: sp.Symbol) -> sp.Expr**
- Apply multiple field substitutions (convenience wrapper around substitute_field).

**gateaux_derivative(form: sp.Expr, trial_var: sp.Expr, increment_var: sp.Expr) -> sp.Expr**
- Compute the Gateaux (directional) derivative for linearization (used in Newton-Raphson).

---

## fe_spaces.py

### Classes

**LinearElement1D(x: sp.Symbol, x0: sp.Expr, x1: sp.Expr)**
- A two-node linear (P1) finite element on interval [x0, x1]. Provides shape functions N₁, N₂ and element length.

### Functions

**local_trial_expansion(shape_functions, dofs) -> sp.Expr**
- Build the local trial field expansion u_h = Σ Nᵢ dᵢ.

**local_test_expansion(shape_functions, dofs) -> sp.Expr**
- Build the local test field expansion v_h = Σ Nᵢ wᵢ (Galerkin method).

---

## reference.py

### Classes

**ReferenceIntervalP1(xi: sp.Symbol)**
- Reference linear interval element on [0, 1]. Provides shape functions N₁ = 1 - ξ, N₂ = ξ and constant derivatives.

**ReferenceTriangleP1(xi: sp.Symbol, eta: sp.Symbol)**
- Reference P1 triangle with vertices at (0,0), (1,0), (0,1). Barycentric shape functions N₁ = 1-ξ-η, N₂ = ξ, N₃ = η.

**ReferenceQuadrilateralQ1(xi: sp.Symbol, eta: sp.Symbol)**
- Reference bilinear quadrilateral on [-1,1]². Four nodes at corners, bilinear shape functions.

**ReferenceQuadrilateralQ2(xi: sp.Symbol, eta: sp.Symbol)**
- Reference biquadratic quadrilateral on [-1,1]². Nine nodes (corners + edges + center), tensor-product Lagrange shape functions.

**ReferenceHexahedronQ1(xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol)**
- Reference trilinear hexahedron on [-1,1]³. Eight corner nodes.

**ReferenceTetrahedronP1(xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol)**
- Reference linear tetrahedron with vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1). Barycentric shape functions.

**ReferenceTetrahedronP2(xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol)**
- Reference quadratic tetrahedron with 10 nodes. Second-order barycentric shape functions.

**AffineTriangleMap2D(xi: sp.Symbol, eta: sp.Symbol, x1: sp.Expr, y1: sp.Expr, x2: sp.Expr, y2: sp.Expr, x3: sp.Expr, y3: sp.Expr)**
- Affine mapping from reference triangle to physical triangle. Provides Jacobian, determinant, inverse-transpose, and coordinate mappings.

---

## quadrature.py

**integrate_reference_triangle_exact(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol) -> sp.Expr**
- Exact symbolic integration over the reference P1 triangle [0,1]×[0,1-ξ].

**integrate_reference_quadrilateral_exact(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol) -> sp.Expr**
- Exact symbolic integration over the reference quadrilateral [-1,1]².

**integrate_reference_tetra_exact(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol) -> sp.Expr**
- Exact symbolic integration over the reference P1 tetrahedron.

**integrate_reference_hexahedron_exact(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol) -> sp.Expr**
- Exact symbolic integration over the reference hexahedron [-1,1]³.

**triangle_one_point_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol) -> sp.Expr**
- One-point quadrature on the reference triangle (centroid rule). Exact for linear polynomials.

**triangle_three_point_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol) -> sp.Expr**
- Three-point quadrature on the reference triangle. Exact for quadratic polynomials.

**triangle_six_point_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol) -> sp.Expr**
- Six-point quadrature on the reference triangle. Exact for quartic polynomials.

**quadrilateral_gauss_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol, order: int = 2) -> sp.Expr**
- Gauss-Legendre tensor-product quadrature on the reference quad [-1,1]². Supports orders 1, 2, 3.

**hexahedron_gauss_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol, order: int = 2) -> sp.Expr**
- Gauss-Legendre tensor-product quadrature on the reference hex [-1,1]³.

**tetrahedron_one_point_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol) -> sp.Expr**
- One-point quadrature on the reference tetrahedron (centroid rule).

**tetrahedron_four_point_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol) -> sp.Expr**
- Four-point quadrature on the reference tetrahedron. Exact for quadratic polynomials.

---

## extract.py

**extract_coefficient_matrix(expr: sp.Expr, trial_dofs, test_dofs) -> sp.Matrix**
- Extract the local element stiffness matrix from a discretized bilinear form. Returns an (n_test, n_trial) Matrix.

**extract_coefficient_vector(expr: sp.Expr, test_dofs) -> sp.Matrix**
- Extract the local element load vector from a discretized linear form. Returns a column Matrix of length n_test.

---

## validate.py

**ensure_same_variable(expected: sp.Symbol, actual: sp.Symbol) -> None**
- Raise ValueError if expected and actual coordinate symbols don't match.

**field_dependency(expr: sp.Expr, field: sp.Expr) -> bool**
- Check whether an expression depends on a given symbolic field, including through derivatives.

**split_terms(expr: sp.Expr) -> list[sp.Expr]**
- Expand an expression and return its additive terms as a list.

---

## workflow.py

**build_bar_1d_local_problem() -> dict[str, sp.Expr | sp.Matrix]**
- Build the complete 1D bar/Poisson local element problem. Returns Ke (2×2), fe (2×1), and all intermediate expressions.

**build_bar_1d_mass_matrix() -> dict[str, sp.Expr | sp.Matrix]**
- Build the consistent mass matrix for a 1D bar element. Returns Me (2×2) and symbols.

**build_poisson_triangle_p1_local_problem() -> dict[str, sp.Expr | sp.Matrix]**
- Build a teaching-first local P1 triangle problem for scalar Poisson/diffusion. Returns Ke (3×3), fe (3×1), and all intermediate expressions.

---

## assembly.py

**assemble_dense_matrix(global_matrix: np.ndarray, local_matrix: np.ndarray, connectivity: Sequence[int]) -> None**
- Add a local element matrix into a dense global matrix in place.

**assemble_dense_vector(global_vector: np.ndarray, local_vector: np.ndarray, connectivity: Sequence[int]) -> None**
- Add a local element vector into a dense global vector in place.

**apply_dirichlet_by_reduction(global_matrix: np.ndarray, global_vector: np.ndarray, constrained_dofs: Sequence[int], prescribed_values: Sequence[float] | None = None) -> tuple[np.ndarray, np.ndarray, list[int]]**
- Return the reduced linear system after enforcing Dirichlet values.

**expand_reduced_solution(free_solution: Sequence[float], ndofs: int, free_dofs: Sequence[int], constrained_dofs: Sequence[int], prescribed_values: Sequence[float] | None = None) -> np.ndarray**
- Expand a reduced solution back into the full nodal vector.

---

## viz.py

**plot_triangle_shape_functions(reference_element, xi_sym, eta_sym, ax=None)**
- Plot 2D contour plots of triangle shape functions on the reference element. Returns matplotlib Figure.

**plot_affine_mapping(geom, xi_sym, eta_sym, vertex_subs, ax=None)**
- Visualize the affine mapping from reference to physical triangle. Returns matplotlib Figure.

**plot_mesh_solution(nodes, elements, u, ax=None, title="FEM Solution")**
- Plot a scalar FEM solution on a triangular mesh. Returns matplotlib Figure.
