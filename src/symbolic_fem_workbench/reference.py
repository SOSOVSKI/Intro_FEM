"""Reference-element and affine-mapping helpers for teaching examples."""

from __future__ import annotations

from dataclasses import dataclass
import sympy as sp


def _lagrange_1d_m1(s: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    return (sp.simplify((1 - s) / 2), sp.simplify((1 + s) / 2))


def _lagrange_1d_m1_p1_q2(s: sp.Symbol) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    return (
        sp.simplify(s * (s - 1) / 2),
        sp.simplify(1 - s**2),
        sp.simplify(s * (s + 1) / 2),
    )


@dataclass(frozen=True)
class ReferenceIntervalP1:
    """Reference linear interval element on [0, 1]. Two nodes at endpoints.

    Provides shape functions N₁ = 1 - ξ, N₂ = ξ and their constant derivatives.
    This is the 1D counterpart to ReferenceTriangleP1 for consistency.
    """
    xi: sp.Symbol

    @property
    def nodes(self):
        return ((0,), (1,))

    @property
    def shape_functions(self):
        return (1 - self.xi, self.xi)

    @property
    def shape_gradients_reference(self):
        return (sp.Matrix([[-1]]), sp.Matrix([[1]]))


@dataclass(frozen=True)
class ReferenceTriangleP1:
    """Reference P1 triangle with vertices at (0,0), (1,0), (0,1). Provides barycentric shape functions N₁ = 1-ξ-η, N₂ = ξ, N₃ = η and their constant gradients."""
    xi: sp.Symbol
    eta: sp.Symbol

    @property
    def nodes(self):
        return ((0, 0), (1, 0), (0, 1))

    @property
    def shape_functions(self):
        return (1 - self.xi - self.eta, self.xi, self.eta)

    @property
    def shape_gradients_reference(self):
        return (
            sp.Matrix([[-1], [-1]]),
            sp.Matrix([[1], [0]]),
            sp.Matrix([[0], [1]]),
        )

    @property
    def exact_domain(self):
        return ((self.eta, 0, 1 - self.xi), (self.xi, 0, 1))


@dataclass(frozen=True)
class ReferenceQuadrilateralQ1:
    """Reference bilinear quadrilateral on [-1,1]². Four nodes at corners, bilinear shape functions."""
    xi: sp.Symbol
    eta: sp.Symbol

    @property
    def nodes(self):
        return ((-1, -1), (1, -1), (1, 1), (-1, 1))

    @property
    def shape_functions(self):
        xi, eta = self.xi, self.eta
        return (
            sp.simplify((1 - xi) * (1 - eta) / 4),
            sp.simplify((1 + xi) * (1 - eta) / 4),
            sp.simplify((1 + xi) * (1 + eta) / 4),
            sp.simplify((1 - xi) * (1 + eta) / 4),
        )

    @property
    def shape_gradients_reference(self):
        return tuple(sp.Matrix([[sp.diff(N, self.xi)], [sp.diff(N, self.eta)]]) for N in self.shape_functions)


@dataclass(frozen=True)
class ReferenceQuadrilateralQ2:
    """Reference biquadratic quadrilateral on [-1,1]². Nine nodes (corners + edges + center), tensor-product Lagrange shape functions."""
    xi: sp.Symbol
    eta: sp.Symbol

    @property
    def nodes(self):
        return ((-1, -1), (0, -1), (1, -1), (1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (0, 0))

    @property
    def shape_functions(self):
        lxi = _lagrange_1d_m1_p1_q2(self.xi)
        leta = _lagrange_1d_m1_p1_q2(self.eta)
        return (
            sp.simplify(lxi[0] * leta[0]),
            sp.simplify(lxi[1] * leta[0]),
            sp.simplify(lxi[2] * leta[0]),
            sp.simplify(lxi[2] * leta[1]),
            sp.simplify(lxi[2] * leta[2]),
            sp.simplify(lxi[1] * leta[2]),
            sp.simplify(lxi[0] * leta[2]),
            sp.simplify(lxi[0] * leta[1]),
            sp.simplify(lxi[1] * leta[1]),
        )

    @property
    def shape_gradients_reference(self):
        return tuple(sp.Matrix([[sp.diff(N, self.xi)], [sp.diff(N, self.eta)]]) for N in self.shape_functions)


@dataclass(frozen=True)
class ReferenceHexahedronQ1:
    """Reference trilinear hexahedron on [-1,1]³. Eight corner nodes, trilinear shape functions."""
    xi: sp.Symbol
    eta: sp.Symbol
    zeta: sp.Symbol

    @property
    def nodes(self):
        return (
            (-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
            (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1),
        )

    @property
    def shape_functions(self):
        lxi = _lagrange_1d_m1(self.xi)
        leta = _lagrange_1d_m1(self.eta)
        lz = _lagrange_1d_m1(self.zeta)
        return (
            sp.simplify(lxi[0] * leta[0] * lz[0]),
            sp.simplify(lxi[1] * leta[0] * lz[0]),
            sp.simplify(lxi[1] * leta[1] * lz[0]),
            sp.simplify(lxi[0] * leta[1] * lz[0]),
            sp.simplify(lxi[0] * leta[0] * lz[1]),
            sp.simplify(lxi[1] * leta[0] * lz[1]),
            sp.simplify(lxi[1] * leta[1] * lz[1]),
            sp.simplify(lxi[0] * leta[1] * lz[1]),
        )

    @property
    def shape_gradients_reference(self):
        return tuple(
            sp.Matrix([[sp.diff(N, self.xi)], [sp.diff(N, self.eta)], [sp.diff(N, self.zeta)]])
            for N in self.shape_functions
        )


@dataclass(frozen=True)
class ReferenceTetrahedronP1:
    """Reference linear tetrahedron with vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1). Barycentric shape functions."""
    xi: sp.Symbol
    eta: sp.Symbol
    zeta: sp.Symbol

    @property
    def nodes(self):
        return ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))

    @property
    def shape_functions(self):
        L1 = 1 - self.xi - self.eta - self.zeta
        L2 = self.xi
        L3 = self.eta
        L4 = self.zeta
        return (L1, L2, L3, L4)

    @property
    def shape_gradients_reference(self):
        return tuple(
            sp.Matrix([[sp.diff(N, self.xi)], [sp.diff(N, self.eta)], [sp.diff(N, self.zeta)]])
            for N in self.shape_functions
        )


@dataclass(frozen=True)
class ReferenceTetrahedronP2:
    """Reference quadratic tetrahedron with 10 nodes. Second-order barycentric shape functions."""
    xi: sp.Symbol
    eta: sp.Symbol
    zeta: sp.Symbol

    @property
    def nodes(self):
        return (
            (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
            (sp.Rational(1, 2), 0, 0), (sp.Rational(1, 2), sp.Rational(1, 2), 0),
            (0, sp.Rational(1, 2), 0), (0, 0, sp.Rational(1, 2)),
            (sp.Rational(1, 2), 0, sp.Rational(1, 2)), (0, sp.Rational(1, 2), sp.Rational(1, 2)),
        )

    @property
    def shape_functions(self):
        L1 = 1 - self.xi - self.eta - self.zeta
        L2 = self.xi
        L3 = self.eta
        L4 = self.zeta
        return (
            sp.simplify(L1 * (2 * L1 - 1)),
            sp.simplify(L2 * (2 * L2 - 1)),
            sp.simplify(L3 * (2 * L3 - 1)),
            sp.simplify(L4 * (2 * L4 - 1)),
            sp.simplify(4 * L1 * L2),
            sp.simplify(4 * L2 * L3),
            sp.simplify(4 * L1 * L3),
            sp.simplify(4 * L1 * L4),
            sp.simplify(4 * L2 * L4),
            sp.simplify(4 * L3 * L4),
        )

    @property
    def shape_gradients_reference(self):
        return tuple(
            sp.Matrix([[sp.diff(N, self.xi)], [sp.diff(N, self.eta)], [sp.diff(N, self.zeta)]])
            for N in self.shape_functions
        )


@dataclass(frozen=True)
class AffineTriangleMap2D:
    """Affine mapping from the reference triangle to a physical triangle defined by vertices (x₁,y₁), (x₂,y₂), (x₃,y₃). Provides the Jacobian matrix J, its determinant, the inverse transpose J^{-T} for gradient pullback, and the coordinate mapping formulas."""
    xi: sp.Symbol
    eta: sp.Symbol
    x1: sp.Expr
    y1: sp.Expr
    x2: sp.Expr
    y2: sp.Expr
    x3: sp.Expr
    y3: sp.Expr

    @property
    def jacobian(self) -> sp.Matrix:
        return sp.Matrix([[self.x2 - self.x1, self.x3 - self.x1], [self.y2 - self.y1, self.y3 - self.y1]])

    @property
    def det_jacobian(self) -> sp.Expr:
        return sp.simplify(self.jacobian.det())

    @property
    def area(self) -> sp.Expr:
        return sp.simplify(self.det_jacobian / 2)

    @property
    def inverse_transpose_jacobian(self) -> sp.Matrix:
        return sp.simplify(self.jacobian.inv().T)

    @property
    def x_map(self) -> sp.Expr:
        return sp.expand(self.x1 + (self.x2 - self.x1) * self.xi + (self.x3 - self.x1) * self.eta)

    @property
    def y_map(self) -> sp.Expr:
        return sp.expand(self.y1 + (self.y2 - self.y1) * self.xi + (self.y3 - self.y1) * self.eta)
