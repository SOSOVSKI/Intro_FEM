"""Notebook UI helpers for symbolic FEM exploration.

Provides tabbed displays and graceful fallbacks when optional visualization
backends (ipywidgets, pandas, matplotlib) are unavailable.
"""

from __future__ import annotations

import io
from dataclasses import dataclass
from typing import Any

import sympy as sp


@dataclass(frozen=True)
class NumericViewOptions:
    """Numeric rendering options for matrix/vector previews."""

    truncate: bool = True
    digits: int = 4


def symbolic_text(expr: Any, *, pretty: bool = True) -> str:
    """Return symbolic text using ``sympy.pretty`` or plain ``repr``."""
    if pretty:
        try:
            return sp.pretty(expr)
        except Exception:
            return repr(expr)
    return repr(expr)


def matrix_preview_text(matrix_like: Any, *, truncate: bool = True, digits: int = 4) -> str:
    """Return textual matrix/vector preview with optional truncation."""
    mat = sp.Matrix(matrix_like)
    if truncate:
        out = mat.evalf(digits)
    else:
        out = mat
    return sp.pretty(out)


def matrix_full_precision_csv_bytes(matrix_like: Any) -> bytes:
    """Export full-precision matrix/vector CSV bytes for reproducibility."""
    mat = sp.Matrix(matrix_like)
    buf = io.StringIO()
    rows, cols = mat.shape
    for i in range(rows):
        buf.write(",".join(str(mat[i, j]) for j in range(cols)))
        buf.write("\n")
    return buf.getvalue().encode("utf-8")


def visualization_hint(problem: dict[str, Any]) -> str:
    """Build a best-effort textual visualization hint from reference/mesh data."""
    ref = problem.get("reference_triangle") or problem.get("reference_tetrahedron")
    geom = problem.get("geometry")
    hints: list[str] = []
    if ref is not None:
        hints.append(f"Reference element: {type(ref).__name__}")
    if geom is not None:
        hints.append("Geometry map available (affine).")
    if not hints:
        hints.append("No reference/mesh objects found; showing textual fallback only.")
    return "\n".join(hints)


def build_main_panel(problem: dict[str, Any]) -> Any:
    """Create main panel with tabs: Math View, Table View, Codegen/Export.

    Falls back to a plain dict if ipywidgets is unavailable.
    """
    math_content = {
        "symbolic": symbolic_text(problem, pretty=False),
        "hint": visualization_hint(problem),
    }
    table_content = {
        "Ke": matrix_preview_text(problem["Ke"]) if "Ke" in problem else "N/A",
        "fe": matrix_preview_text(problem["fe"]) if "fe" in problem else "N/A",
    }
    export_content = {
        "Ke_csv": matrix_full_precision_csv_bytes(problem["Ke"]) if "Ke" in problem else b"",
        "fe_csv": matrix_full_precision_csv_bytes(problem["fe"]) if "fe" in problem else b"",
    }

    try:
        import ipywidgets as widgets

        math = widgets.Accordion(children=[widgets.Textarea(value=math_content["symbolic"], layout=widgets.Layout(width="100%", height="240px")),
                                           widgets.HTML(value=f"<pre>{math_content['hint']}</pre>")])
        math.set_title(0, "Symbolic form")
        math.set_title(1, "Reference/mesh hints")

        table = widgets.Accordion(children=[widgets.Textarea(value=table_content["Ke"], layout=widgets.Layout(width="100%", height="180px")),
                                            widgets.Textarea(value=table_content["fe"], layout=widgets.Layout(width="100%", height="120px"))])
        table.set_title(0, "Assembled matrix preview")
        table.set_title(1, "Assembled vector preview")

        export = widgets.VBox([
            widgets.HTML("<b>Full-precision CSV ready for download/export from bytes payload.</b>"),
            widgets.HTML(f"<pre>Ke bytes: {len(export_content['Ke_csv'])}\nfe bytes: {len(export_content['fe_csv'])}</pre>"),
        ])

        tabs = widgets.Tab(children=[math, table, export])
        tabs.set_title(0, "Math View")
        tabs.set_title(1, "Table View")
        tabs.set_title(2, "Codegen/Export")
        return tabs
    except Exception:
        return {
            "Math View": math_content,
            "Table View": table_content,
            "Codegen/Export": export_content,
        }
