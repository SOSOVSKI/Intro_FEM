"""Very small backend wrapper around the I❤️LA printer.

This module intentionally stops at source generation. Compilation to NumPy is left to the
external I❤️LA toolchain so the course can keep that boundary explicit.
"""

from __future__ import annotations

import sympy as sp

from ..printers.iheartla_printer import iheartla_matrix_definition, iheartla_scalar_definition


def emit_iheartla_scalar(name: str, expr: sp.Expr, declarations: dict[str, str] | None = None) -> str:
    return iheartla_scalar_definition(name, expr, declarations)


def emit_iheartla_matrix(name: str, matrix: sp.Matrix, declarations: dict[str, str] | None = None) -> str:
    return iheartla_matrix_definition(name, matrix, declarations)
