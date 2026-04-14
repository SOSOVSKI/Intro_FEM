"""Minimal I❤️LA printer helpers for course examples."""

from __future__ import annotations

import sympy as sp


def _expr_to_iheartla(expr: sp.Expr) -> str:
    return sp.sstr(sp.simplify(expr)).replace("**", "^")


def iheartla_scalar_definition(name: str, expr: sp.Expr, declarations: dict[str, str] | None = None) -> str:
    lines = [f"{name} = {_expr_to_iheartla(expr)}"]
    if declarations:
        where_parts = [f"{sym} in {domain}" for sym, domain in declarations.items()]
        lines.append("where " + ", ".join(where_parts))
    return "\n".join(lines)


def iheartla_matrix_definition(name: str, matrix: sp.Matrix, declarations: dict[str, str] | None = None) -> str:
    rows = []
    for i in range(matrix.rows):
        row = " & ".join(_expr_to_iheartla(matrix[i, j]) for j in range(matrix.cols))
        rows.append(row)
    body = " \\\n".join(rows)
    lines = [f"{name} = \\begin{{bmatrix}} {body} \\end{{bmatrix}}"]
    if declarations:
        where_parts = [f"{sym} in {domain}" for sym, domain in declarations.items()]
        lines.append("where " + ", ".join(where_parts))
    return "\n".join(lines)
