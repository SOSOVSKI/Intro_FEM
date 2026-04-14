"""Finite-element interpolation helpers."""

from __future__ import annotations

from dataclasses import dataclass
import sympy as sp


@dataclass(frozen=True)
class LinearElement1D:
    """A two-node linear (P1) finite element on an interval [x0, x1]. Provides symbolic shape functions N₁ = (x₁ - x) / L and N₂ = (x - x₀) / L where L = x₁ - x₀."""
    x: sp.Symbol
    x0: sp.Expr
    x1: sp.Expr

    @property
    def length(self) -> sp.Expr:
        return sp.simplify(self.x1 - self.x0)

    @property
    def shape_functions(self) -> tuple[sp.Expr, sp.Expr]:
        L = self.length
        return (
            sp.simplify((self.x1 - self.x) / L),
            sp.simplify((self.x - self.x0) / L),
        )


def local_trial_expansion(shape_functions, dofs) -> sp.Expr:
    """Build the local trial field expansion u_h = Σ Nᵢ dᵢ. Parameters: shape_functions - tuple of shape function expressions, dofs - list of DOF symbols. Returns the expanded polynomial."""
    return sp.expand(sum(N * d for N, d in zip(shape_functions, dofs)))


def local_test_expansion(shape_functions, dofs) -> sp.Expr:
    """Build the local test field expansion v_h = Σ Nᵢ wᵢ. Identical to local_trial_expansion (Galerkin method uses the same basis for trial and test)."""
    return sp.expand(sum(N * d for N, d in zip(shape_functions, dofs)))
