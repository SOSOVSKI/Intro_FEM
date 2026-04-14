from __future__ import annotations

import sympy as sp

from symbolic_fem_workbench.quadrature import (
    integrate_reference_quadrilateral_exact,
    integrate_reference_tetra_exact,
    integrate_reference_hexahedron_exact,
    quadrilateral_gauss_rule,
    tetrahedron_four_point_rule,
    hexahedron_gauss_rule,
)
from symbolic_fem_workbench.reference import (
    ReferenceQuadrilateralQ1,
    ReferenceQuadrilateralQ2,
    ReferenceTetrahedronP1,
    ReferenceTetrahedronP2,
    ReferenceHexahedronQ1,
)



def _check_partition_of_unity(element) -> None:
    assert sp.simplify(sum(element.shape_functions) - 1) == 0



def _check_gradient_sum_zero(element, dim: int) -> None:
    grad_sum = sum(element.shape_gradients_reference, sp.zeros(dim, 1))
    assert sp.simplify(grad_sum) == sp.zeros(dim, 1)



def _check_nodal_kronecker_property(element) -> None:
    coords = []
    if hasattr(element, "xi"):
        coords.append(element.xi)
    if hasattr(element, "eta"):
        coords.append(element.eta)
    if hasattr(element, "zeta"):
        coords.append(element.zeta)

    for i, N in enumerate(element.shape_functions):
        for j, node in enumerate(element.nodes):
            subs = {coord: value for coord, value in zip(coords, node)}
            expected = 1 if i == j else 0
            assert sp.simplify(N.subs(subs) - expected) == 0



def test_quadrilateral_q1_properties() -> None:
    xi, eta = sp.symbols("xi eta")
    element = ReferenceQuadrilateralQ1(xi, eta)
    _check_partition_of_unity(element)
    _check_gradient_sum_zero(element, 2)
    _check_nodal_kronecker_property(element)



def test_quadrilateral_q2_properties() -> None:
    xi, eta = sp.symbols("xi eta")
    element = ReferenceQuadrilateralQ2(xi, eta)
    _check_partition_of_unity(element)
    _check_gradient_sum_zero(element, 2)
    _check_nodal_kronecker_property(element)



def test_tetrahedron_p1_properties() -> None:
    xi, eta, zeta = sp.symbols("xi eta zeta")
    element = ReferenceTetrahedronP1(xi, eta, zeta)
    _check_partition_of_unity(element)
    _check_gradient_sum_zero(element, 3)
    _check_nodal_kronecker_property(element)



def test_tetrahedron_p2_properties() -> None:
    xi, eta, zeta = sp.symbols("xi eta zeta")
    element = ReferenceTetrahedronP2(xi, eta, zeta)
    _check_partition_of_unity(element)
    _check_gradient_sum_zero(element, 3)
    _check_nodal_kronecker_property(element)



def test_hexahedron_q1_properties() -> None:
    xi, eta, zeta = sp.symbols("xi eta zeta")
    element = ReferenceHexahedronQ1(xi, eta, zeta)
    _check_partition_of_unity(element)
    _check_gradient_sum_zero(element, 3)
    _check_nodal_kronecker_property(element)



def test_quadrilateral_gauss_rule_matches_exact_for_polynomial() -> None:
    xi, eta = sp.symbols("xi eta")
    expr = 1 + xi**2 + eta**2 + xi**2 * eta**2 + xi * eta
    exact = integrate_reference_quadrilateral_exact(expr, xi, eta)
    approx = quadrilateral_gauss_rule(expr, xi, eta, order=2)
    assert sp.simplify(exact - approx) == 0



def test_tetrahedron_four_point_rule_matches_exact_for_quadratic() -> None:
    xi, eta, zeta = sp.symbols("xi eta zeta")
    expr = 1 + xi + eta + zeta + xi * eta + zeta**2
    exact = integrate_reference_tetra_exact(expr, xi, eta, zeta)
    approx = tetrahedron_four_point_rule(expr, xi, eta, zeta)
    assert sp.simplify(exact - approx) == 0



def test_hexahedron_gauss_rule_matches_exact_for_polynomial() -> None:
    xi, eta, zeta = sp.symbols("xi eta zeta")
    expr = 1 + xi**2 + eta * zeta + xi**2 * eta**2 * zeta**2
    exact = integrate_reference_hexahedron_exact(expr, xi, eta, zeta)
    approx = hexahedron_gauss_rule(expr, xi, eta, zeta, order=2)
    assert sp.simplify(exact - approx) == 0
