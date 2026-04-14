"""Structured form containers used by the symbolic FEM workflow."""

from __future__ import annotations

from dataclasses import dataclass
import sympy as sp

from .symbols import Domain1D


@dataclass(frozen=True)
class DomainIntegral:
    """A domain integral ∫ integrand dx over a 1D domain. Call .as_integral() to get a SymPy Integral object."""
    integrand: sp.Expr
    domain: Domain1D

    def as_integral(self) -> sp.Integral:
        return sp.Integral(self.integrand, self.domain)


@dataclass(frozen=True)
class BoundaryContribution:
    """A boundary term from integration by parts. Stores the expression, coordinate, evaluation point, and optional label. Call .evaluate() to substitute the coordinate and simplify."""
    expr: sp.Expr
    x: sp.Symbol
    location: sp.Expr
    label: str | None = None

    def evaluate(self) -> sp.Expr:
        return sp.simplify(self.expr.subs(self.x, self.location))


@dataclass(frozen=True)
class WeightedResidual:
    """The complete weighted residual before splitting into bilinear/linear forms. Contains domain integrals and boundary contributions. Call .as_expression() to combine into a single SymPy expression."""
    domain_integrals: tuple[DomainIntegral, ...]
    boundary_terms: tuple[BoundaryContribution, ...] = ()

    def as_expression(self) -> sp.Expr:
        return sum((di.as_integral() for di in self.domain_integrals), sp.Integer(0)) + sum(
            (bt.evaluate() for bt in self.boundary_terms), sp.Integer(0)
        )


@dataclass(frozen=True)
class WeakForm:
    """The split weak form with bilinear part a(u,v) and linear part F(v). The bilinear form contains terms with both trial and test functions. The linear form is the right-hand side."""
    bilinear: sp.Expr
    linear: sp.Expr
    trial: sp.Expr
    test: sp.Expr
