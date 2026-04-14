"""Small workflow helpers for course-ready examples."""

from __future__ import annotations

import sympy as sp

from .extract import extract_coefficient_matrix, extract_coefficient_vector
from .fe_spaces import LinearElement1D, local_test_expansion, local_trial_expansion
from .quadrature import integrate_reference_triangle_exact
from .reference import AffineTriangleMap2D, ReferenceTriangleP1, ReferenceTetrahedronP1
from .symbols import make_field_1d, make_field_2d
from .transforms import (
    apply_neumann_flux,
    drop_dirichlet_boundary,
    grad_2d,
    integrate_divergence_1d,
    pullback_gradient_2d,
    split_linear_weak_form,
    substitute_fe,
)


def build_bar_1d_local_problem() -> dict[str, sp.Expr | sp.Matrix]:
    """Build the complete 1D bar/Poisson local element problem symbolically. Derives the weak form from the strong form d/dx(EA·du/dx) + q = 0, applies integration by parts, enforces homogeneous Dirichlet at x=0 and Neumann traction P at x=L, substitutes P1 shape functions, and extracts the 2×2 element stiffness matrix and 2×1 load vector. Returns a dict with keys: x, L, E, A, q, P (symbols); u, v (fields); weak_bilinear, weak_linear (forms); Ke (2×2 Matrix), fe (2×1 Matrix)."""
    x, L = sp.symbols("x L", positive=True)
    E, A = sp.symbols("E A", positive=True)
    q, P = sp.symbols("q P", real=True)

    u = make_field_1d("u", x)
    v = make_field_1d("v", x)

    domain = (x, 0, L)
    flux = E * A * sp.diff(u, x)
    domain_term, (left_boundary, right_boundary) = integrate_divergence_1d(flux, v, x, domain)
    left_bc = drop_dirichlet_boundary(left_boundary, v)
    right_bc = apply_neumann_flux(right_boundary, flux_expr=flux, prescribed_flux=P)

    weak_expr = domain_term.as_integral() - sp.Integral(q * v, domain) + left_bc + right_bc
    weak_form = split_linear_weak_form(weak_expr, u, v)

    element = LinearElement1D(x=x, x0=0, x1=L)
    N = element.shape_functions
    d0, d1 = sp.symbols("d0 d1", real=True)
    w0, w1 = sp.symbols("w0 w1", real=True)
    uh = local_trial_expansion(N, [d0, d1])
    vh = local_test_expansion(N, [w0, w1])

    a_disc = substitute_fe(weak_form.bilinear, {u: uh, v: vh}, x)
    l_disc = substitute_fe(weak_form.linear, {v: vh}, x)

    Ke = extract_coefficient_matrix(a_disc, [d0, d1], [w0, w1])
    fe = extract_coefficient_vector(l_disc, [w0, w1])

    return {
        "x": x,
        "L": L,
        "E": E,
        "A": A,
        "q": q,
        "P": P,
        "u": u,
        "v": v,
        "weak_bilinear": weak_form.bilinear,
        "weak_linear": weak_form.linear,
        "Ke": sp.simplify(Ke),
        "fe": sp.simplify(fe),
    }


def build_poisson_triangle_p1_local_problem() -> dict[str, sp.Expr | sp.Matrix]:
    """Build a teaching-first local P1 triangle problem for scalar Poisson/diffusion."""
    x, y = sp.symbols("x y", real=True)
    u = make_field_2d("u", x, y)
    v = make_field_2d("v", x, y)

    xi, eta = sp.symbols("xi eta", real=True)
    x1, y1, x2, y2, x3, y3 = sp.symbols("x1 y1 x2 y2 x3 y3", real=True)
    f = sp.symbols("f", real=True)

    ref = ReferenceTriangleP1(xi=xi, eta=eta)
    geom = AffineTriangleMap2D(xi=xi, eta=eta, x1=x1, y1=y1, x2=x2, y2=y2, x3=x3, y3=y3)

    d0, d1, d2 = sp.symbols("d0 d1 d2", real=True)
    w0, w1, w2 = sp.symbols("w0 w1 w2", real=True)

    uh = local_trial_expansion(ref.shape_functions, [d0, d1, d2])
    vh = local_test_expansion(ref.shape_functions, [w0, w1, w2])

    grad_uh_ref = grad_2d(uh, xi, eta)
    grad_vh_ref = grad_2d(vh, xi, eta)
    grad_uh_phys = pullback_gradient_2d(grad_uh_ref, geom.jacobian)
    grad_vh_phys = pullback_gradient_2d(grad_vh_ref, geom.jacobian)

    bilinear_integrand = sp.simplify((grad_uh_phys.T * grad_vh_phys)[0] * geom.det_jacobian)
    linear_integrand = sp.simplify(f * vh * geom.det_jacobian)

    bilinear_expr = integrate_reference_triangle_exact(bilinear_integrand, xi, eta)
    linear_expr = integrate_reference_triangle_exact(linear_integrand, xi, eta)

    Ke = extract_coefficient_matrix(bilinear_expr, [d0, d1, d2], [w0, w1, w2])
    fe = extract_coefficient_vector(linear_expr, [w0, w1, w2])

    unit_subs = {x1: 0, y1: 0, x2: 1, y2: 0, x3: 0, y3: 1}
    Ke_unit = sp.simplify(Ke.subs(unit_subs))
    fe_unit = sp.simplify(fe.subs(unit_subs))

    return {
        "x": x,
        "y": y,
        "u": u,
        "v": v,
        "xi": xi,
        "eta": eta,
        "f": f,
        "reference_triangle": ref,
        "geometry": geom,
        "bilinear_integrand_reference": bilinear_integrand,
        "linear_integrand_reference": linear_integrand,
        "Ke": sp.simplify(Ke),
        "fe": sp.simplify(fe),
        "Ke_unit_right_triangle": Ke_unit,
        "fe_unit_right_triangle": fe_unit,
    }


def build_bar_1d_mass_matrix() -> dict[str, sp.Expr | sp.Matrix]:
    """Build the consistent mass matrix for a 1D bar element.

    Derives M_e from the L2 inner product ∫ ρA·u·v dx on a two-node
    linear element. Used in time-dependent problems: M·a_ddot + K·a = F.

    Returns a dict with keys: x, L, rho, A (symbols); Me (2×2 Matrix).
    """
    x, L = sp.symbols("x L", positive=True)
    rho, A = sp.symbols("rho A", positive=True)

    element = LinearElement1D(x=x, x0=0, x1=L)
    N = element.shape_functions

    d0, d1 = sp.symbols("d0 d1", real=True)
    w0, w1 = sp.symbols("w0 w1", real=True)

    uh = local_trial_expansion(N, [d0, d1])
    vh = local_test_expansion(N, [w0, w1])

    mass_integrand = rho * A * uh * vh
    mass_integral = sp.integrate(mass_integrand, (x, 0, L))
    mass_expanded = sp.expand(mass_integral)

    Me = extract_coefficient_matrix(mass_expanded, [d0, d1], [w0, w1])

    return {
        "x": x,
        "L": L,
        "rho": rho,
        "A": A,
        "Me": sp.simplify(Me),
    }


def build_elasticity_triangle_p1_2d(
    formulation: str = "plane_stress",
) -> dict[str, sp.Expr | sp.Matrix]:
    """Build local element matrices for 2D linear elasticity on a P1 triangle.

    Demonstrates the full pipeline: constitutive matrix D, gradient pullback,
    B-matrix construction, and K_e = t A B^T D B.

    Parameters
    ----------
    formulation : str
        Either "plane_stress" or "plane_strain".

    Returns a dict with: E_sym, nu (material); D (3×3 constitutive); B (3×6
    strain-displacement); Ke (6×6 element stiffness); geometry; etc.
    """
    from .elasticity import (
        plane_stress_D, plane_strain_D, B_matrix_triangle_2d,
        element_stiffness_BtDB,
    )

    E_sym, nu = sp.symbols("E nu", positive=True)
    t_sym = sp.Symbol("t", positive=True)  # thickness

    if formulation == "plane_stress":
        D = plane_stress_D(E_sym, nu)
    elif formulation == "plane_strain":
        D = plane_strain_D(E_sym, nu)
    else:
        raise ValueError(f"Unknown formulation: {formulation!r}")

    xi, eta = sp.symbols("xi eta", real=True)
    x1, y1, x2, y2, x3, y3 = sp.symbols("x1 y1 x2 y2 x3 y3", real=True)

    ref = ReferenceTriangleP1(xi=xi, eta=eta)
    geom = AffineTriangleMap2D(xi=xi, eta=eta, x1=x1, y1=y1,
                                x2=x2, y2=y2, x3=x3, y3=y3)

    # Physical-space gradients via J^{-T}
    JinvT = geom.inverse_transpose_jacobian
    dN_dx = []
    dN_dy = []
    for grad_ref in ref.shape_gradients_reference:
        grad_phys = sp.simplify(JinvT * grad_ref)
        dN_dx.append(grad_phys[0])
        dN_dy.append(grad_phys[1])

    B = B_matrix_triangle_2d(dN_dx, dN_dy)
    Ke = element_stiffness_BtDB(B, D, geom.area, thickness=t_sym)

    # Specialise for unit right triangle
    unit_subs = {x1: 0, y1: 0, x2: 1, y2: 0, x3: 0, y3: 1}
    Ke_unit = sp.simplify(Ke.subs(unit_subs))

    return {
        "E": E_sym,
        "nu": nu,
        "t": t_sym,
        "D": D,
        "formulation": formulation,
        "B": sp.simplify(B),
        "Ke": sp.simplify(Ke),
        "Ke_unit_right_triangle": Ke_unit,
        "reference_triangle": ref,
        "geometry": geom,
        "xi": xi, "eta": eta,
        "dN_dx": dN_dx, "dN_dy": dN_dy,
    }


def build_elasticity_tetra_p1_3d() -> dict[str, sp.Expr | sp.Matrix]:
    """Build local element matrices for 3D linear elasticity on a P1 tetrahedron.

    Derives the full pipeline: 6×6 constitutive matrix D, gradient pullback
    via 3×3 Jacobian, B-matrix (6×12), and K_e = V B^T D B (12×12).

    Returns a dict with: E, nu (material); D (6×6); B (6×12);
    Ke (12×12); vertex symbols; etc.
    """
    from .elasticity import isotropic_3d_D, B_matrix_tetra_3d, element_stiffness_BtDB

    E_sym, nu = sp.symbols("E nu", positive=True)

    xi, eta, zeta = sp.symbols("xi eta zeta", real=True)
    ref = ReferenceTetrahedronP1(xi=xi, eta=eta, zeta=zeta)

    # Generic tetrahedron vertices
    x1, y1, z1 = sp.symbols("x1 y1 z1", real=True)
    x2, y2, z2 = sp.symbols("x2 y2 z2", real=True)
    x3, y3, z3 = sp.symbols("x3 y3 z3", real=True)
    x4, y4, z4 = sp.symbols("x4 y4 z4", real=True)

    # Affine map: x = x1 + (x2-x1)ξ + (x3-x1)η + (x4-x1)ζ  (similarly for y, z)
    J = sp.Matrix([
        [x2 - x1, x3 - x1, x4 - x1],
        [y2 - y1, y3 - y1, y4 - y1],
        [z2 - z1, z3 - z1, z4 - z1],
    ])
    det_J = sp.simplify(J.det())
    volume = sp.simplify(sp.Abs(det_J) / 6)
    JinvT = sp.simplify(J.inv().T)

    # Physical-space gradients
    dN_dx, dN_dy, dN_dz = [], [], []
    for grad_ref in ref.shape_gradients_reference:
        grad_phys = sp.simplify(JinvT * grad_ref)
        dN_dx.append(grad_phys[0])
        dN_dy.append(grad_phys[1])
        dN_dz.append(grad_phys[2])

    D = isotropic_3d_D(E_sym, nu)
    B = B_matrix_tetra_3d(dN_dx, dN_dy, dN_dz)
    Ke = element_stiffness_BtDB(B, D, volume)

    # Specialise for unit right tetrahedron
    unit_subs = {
        x1: 0, y1: 0, z1: 0,
        x2: 1, y2: 0, z2: 0,
        x3: 0, y3: 1, z3: 0,
        x4: 0, y4: 0, z4: 1,
    }
    Ke_unit = sp.simplify(Ke.subs(unit_subs))

    return {
        "E": E_sym,
        "nu": nu,
        "D": D,
        "B": sp.simplify(B),
        "Ke": sp.simplify(Ke),
        "Ke_unit_right_tetra": Ke_unit,
        "jacobian": J,
        "det_jacobian": det_J,
        "volume": volume,
        "reference_tetra": ref,
        "xi": xi, "eta": eta, "zeta": zeta,
        "x1": x1, "y1": y1, "z1": z1,
        "x2": x2, "y2": y2, "z2": z2,
        "x3": x3, "y3": y3, "z3": z3,
        "x4": x4, "y4": y4, "z4": z4,
        "dN_dx": dN_dx, "dN_dy": dN_dy, "dN_dz": dN_dz,
    }
