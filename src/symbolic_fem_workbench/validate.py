"""Small validation helpers used across the symbolic FEM workflow."""

from __future__ import annotations

import sympy as sp


def ensure_same_variable(expected: sp.Symbol, actual: sp.Symbol) -> None:
    """Raise ValueError if expected and actual coordinate symbols don't match."""
    if expected != actual:
        raise ValueError(f"Expected variable {expected}, got {actual}")


def field_dependency(expr: sp.Expr, field: sp.Expr) -> bool:
    """Check whether an expression depends on a given symbolic field, including through derivatives. Traverses the expression tree looking for the field or any Derivative involving it."""
    if expr.has(field):
        return True
    field_func = getattr(field, "func", None)
    for node in sp.preorder_traversal(expr):
        if isinstance(node, sp.Derivative):
            if node.expr == field:
                return True
            if field_func is not None and getattr(node.expr, "func", None) == field_func:
                return True
    return False


def split_terms(expr: sp.Expr) -> list[sp.Expr]:
    """Expand an expression and return its additive terms as a list."""
    expanded = sp.expand(expr)
    return list(expanded.as_ordered_terms())
