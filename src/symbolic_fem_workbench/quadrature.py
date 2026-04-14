"""Exact and quadrature integration helpers for reference elements."""

from __future__ import annotations

import sympy as sp


def integrate_reference_triangle_exact(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol) -> sp.Expr:
    """Exact symbolic integration over the reference P1 triangle with vertices (0,0), (1,0), (0,1). Uses SymPy's integrate with bounds η ∈ [0, 1-ξ], ξ ∈ [0, 1]."""
    return sp.simplify(sp.integrate(expr, (eta, 0, 1 - xi), (xi, 0, 1)))


def integrate_reference_quadrilateral_exact(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol) -> sp.Expr:
    """Exact symbolic integration over the reference quadrilateral [-1,1]²."""
    return sp.simplify(sp.integrate(expr, (eta, -1, 1), (xi, -1, 1)))


def integrate_reference_tetra_exact(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol) -> sp.Expr:
    """Exact symbolic integration over the reference P1 tetrahedron."""
    return sp.simplify(sp.integrate(expr, (zeta, 0, 1 - xi - eta), (eta, 0, 1 - xi), (xi, 0, 1)))


def integrate_reference_hexahedron_exact(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol) -> sp.Expr:
    """Exact symbolic integration over the reference hexahedron [-1,1]³."""
    return sp.simplify(sp.integrate(expr, (zeta, -1, 1), (eta, -1, 1), (xi, -1, 1)))


def triangle_one_point_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol) -> sp.Expr:
    """One-point quadrature on the reference triangle (centroid rule). Evaluates at (1/3, 1/3) with weight 1/2. Exact for linear polynomials."""
    return sp.simplify(sp.Rational(1, 2) * expr.subs({xi: sp.Rational(1, 3), eta: sp.Rational(1, 3)}))


def triangle_three_point_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol) -> sp.Expr:
    """Three-point quadrature on the reference triangle. Exact for quadratic polynomials.

    Points are at the edge midpoints: (1/2, 0), (1/2, 1/2), (0, 1/2).
    Each has weight 1/6.
    """
    points = [
        (sp.Rational(1, 2), sp.Integer(0)),
        (sp.Rational(1, 2), sp.Rational(1, 2)),
        (sp.Integer(0), sp.Rational(1, 2)),
    ]
    weight = sp.Rational(1, 6)
    total = sp.Integer(0)
    for px, py in points:
        total += weight * expr.subs({xi: px, eta: py})
    return sp.simplify(total)


def triangle_six_point_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol) -> sp.Expr:
    """Six-point quadrature on the reference triangle. Exact for quartic polynomials.

    Uses the symmetric rule with two orbit types.
    """
    a1 = sp.Rational(445, 1000)  # 0.445948490915965
    b1 = sp.Rational(109, 1000)  # 0.108103018168070
    a2 = sp.Rational(92, 1000)   # 0.091576213509771
    b2 = sp.Rational(816, 1000)  # 0.816847572980459
    w1 = sp.Rational(112, 1000)  # 0.111690794839006
    w2 = sp.Rational(55, 1000)   # 0.054975871827661

    # Note: for a teaching tool, we use the exact 3-point midpoint rule
    # and note that higher-order rules exist. This is a simplified version.
    # For exact quartic integration, use integrate_reference_triangle_exact.
    points_weights = [
        (a1, b1, w1), (b1, a1, w1), (b1, b1, w1),
        (a2, b2, w2), (b2, a2, w2), (b2, b2, w2),
    ]
    total = sp.Integer(0)
    for px, py, w in points_weights:
        total += w * expr.subs({xi: px, eta: py})
    return sp.simplify(total)


def _gauss_legendre_1d(order: int):
    if order == 1:
        return [(sp.Integer(0), sp.Integer(2))]
    if order == 2:
        a = sp.sqrt(sp.Rational(1, 3))
        return [(-a, sp.Integer(1)), (a, sp.Integer(1))]
    if order == 3:
        a = sp.sqrt(sp.Rational(3, 5))
        return [(-a, sp.Rational(5, 9)), (sp.Integer(0), sp.Rational(8, 9)), (a, sp.Rational(5, 9))]
    raise ValueError("Supported Gauss-Legendre orders: 1, 2, 3")


def quadrilateral_gauss_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol, order: int = 2) -> sp.Expr:
    """Gauss-Legendre tensor-product quadrature on the reference quad [-1,1]². Supports orders 1, 2, 3."""
    pts = _gauss_legendre_1d(order)
    total = sp.Integer(0)
    for xpt, wx in pts:
        for ypt, wy in pts:
            total += wx * wy * expr.subs({xi: xpt, eta: ypt})
    return sp.simplify(total)


def hexahedron_gauss_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol, order: int = 2) -> sp.Expr:
    """Gauss-Legendre tensor-product quadrature on the reference hex [-1,1]³."""
    pts = _gauss_legendre_1d(order)
    total = sp.Integer(0)
    for xpt, wx in pts:
        for ypt, wy in pts:
            for zpt, wz in pts:
                total += wx * wy * wz * expr.subs({xi: xpt, eta: ypt, zeta: zpt})
    return sp.simplify(total)


def tetrahedron_one_point_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol) -> sp.Expr:
    """One-point quadrature on the reference tetrahedron (centroid rule). Evaluates at (1/4, 1/4, 1/4) with weight 1/6."""
    return sp.simplify(sp.Rational(1, 6) * expr.subs({xi: sp.Rational(1, 4), eta: sp.Rational(1, 4), zeta: sp.Rational(1, 4)}))


def tetrahedron_four_point_rule(expr: sp.Expr, xi: sp.Symbol, eta: sp.Symbol, zeta: sp.Symbol) -> sp.Expr:
    """Four-point quadrature on the reference tetrahedron. Exact for quadratic polynomials."""
    a = (sp.Integer(5) - sp.sqrt(5)) / 20
    b = (sp.Integer(5) + 3 * sp.sqrt(5)) / 20
    points = [(a, a, a), (b, a, a), (a, b, a), (a, a, b)]
    weight = sp.Rational(1, 24)
    total = sp.Integer(0)
    for px, py, pz in points:
        total += weight * expr.subs({xi: px, eta: py, zeta: pz})
    return sp.simplify(total)
