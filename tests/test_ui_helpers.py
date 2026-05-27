import sympy as sp

from symbolic_fem_workbench.ui import (
    build_main_panel,
    matrix_full_precision_csv_bytes,
    matrix_preview_text,
    symbolic_text,
)


def test_symbolic_text_pretty_and_repr():
    x = sp.symbols("x")
    expr = sp.sin(x) + x**2
    assert "sin" in symbolic_text(expr, pretty=False)
    assert "2" in symbolic_text(expr, pretty=True)


def test_matrix_preview_truncation_toggle():
    m = sp.Matrix([[sp.Rational(1, 3), sp.pi]])
    trunc = matrix_preview_text(m, truncate=True, digits=3)
    full = matrix_preview_text(m, truncate=False)
    assert "0.333" in trunc or "0.3333" in trunc
    assert "π" in full or "pi" in full


def test_full_precision_csv_export():
    m = sp.Matrix([[sp.Rational(1, 3), sp.sqrt(2)]])
    csv_bytes = matrix_full_precision_csv_bytes(m)
    text = csv_bytes.decode("utf-8")
    assert "1/3" in text
    assert "sqrt(2)" in text


def test_build_main_panel_fallback_or_widget():
    payload = {"Ke": sp.eye(2), "fe": sp.Matrix([1, 2])}
    panel = build_main_panel(payload)
    if isinstance(panel, dict):
        assert "Math View" in panel
        assert "Table View" in panel
        assert "Codegen/Export" in panel
    else:
        # ipywidgets Tab-like object
        assert hasattr(panel, "children")
