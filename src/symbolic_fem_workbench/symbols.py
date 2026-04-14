"""Symbol and domain helpers for teaching-first symbolic FEM workflows."""

from __future__ import annotations

from dataclasses import dataclass
import sympy as sp

Domain1D = tuple[sp.Symbol, sp.Expr, sp.Expr]
Domain2D = tuple[tuple[sp.Symbol, sp.Expr, sp.Expr], tuple[sp.Symbol, sp.Expr, sp.Expr]]


@dataclass(frozen=True)
class ScalarField:
    """A thin wrapper around a SymPy scalar field expression."""

    name: str
    expression: sp.Expr

    def __sympy__(self) -> sp.Expr:
        return self.expression

    def __getattr__(self, item):
        return getattr(self.expression, item)

    def __repr__(self) -> str:
        return f"ScalarField(name={self.name!r}, expression={self.expression!r})"


def make_field_1d(name: str, x: sp.Symbol) -> sp.Expr:
    """Create a symbolic 1D scalar field as a SymPy Function application. Returns Function(name)(x), which SymPy treats as an undefined function of x, suitable for differentiation and integration."""
    return sp.Function(name)(x)


def make_field_2d(name: str, x: sp.Symbol, y: sp.Symbol) -> sp.Expr:
    """Create a symbolic 2D scalar field as a SymPy Function application. Returns Function(name)(x, y)."""
    return sp.Function(name)(x, y)
