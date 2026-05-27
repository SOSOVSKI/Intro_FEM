"""Small validation helpers used across the symbolic FEM workflow."""

from __future__ import annotations

import sympy as sp
import traceback


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


def validate_symbolic_inputs(
    *,
    trial_dofs: list[sp.Symbol],
    test_dofs: list[sp.Symbol],
    coefficients: dict[str, object] | None = None,
    required_coefficients: tuple[str, ...] = (),
    selected_boundaries: tuple[str, ...] = (),
) -> None:
    """Validate common UI inputs before symbolic compute steps.

    Checks:
    - trial/test DOFs are non-empty and dimension-compatible
    - required coefficients exist and are not None
    - at least one boundary is selected (when boundary-sensitive workflows run)
    """
    if not trial_dofs or not test_dofs:
        raise ValueError("Trial/test DOF lists must be non-empty.")
    if len(trial_dofs) != len(test_dofs):
        raise ValueError(
            "Trial/test DOF order mismatch: "
            f"{len(trial_dofs)} trial vs {len(test_dofs)} test."
        )

    coefficient_map = coefficients or {}
    missing = [name for name in required_coefficients if coefficient_map.get(name) is None]
    if missing:
        raise ValueError(f"Missing required coefficient(s): {', '.join(missing)}")

    if required_coefficients and not selected_boundaries:
        raise ValueError("At least one boundary selection is required.")


def run_compute_step(step_name: str, compute_fn):
    """Run a symbolic compute step with readable error + traceback details."""
    try:
        return {"ok": True, "result": compute_fn(), "error": None, "traceback": None}
    except Exception as exc:  # pragma: no cover - defensive wrapper
        return {
            "ok": False,
            "result": None,
            "error": f"{step_name} failed: {exc}",
            "traceback": traceback.format_exc(),
        }


def sanity_checks_panel_data(
    *,
    trial_dofs: list[sp.Symbol],
    test_dofs: list[sp.Symbol],
    coefficients: dict[str, object] | None,
    required_coefficients: tuple[str, ...],
    selected_boundaries: tuple[str, ...],
) -> dict[str, object]:
    """Build lightweight sanity-check data suitable for a small UI panel."""
    issues: list[str] = []
    try:
        validate_symbolic_inputs(
            trial_dofs=trial_dofs,
            test_dofs=test_dofs,
            coefficients=coefficients,
            required_coefficients=required_coefficients,
            selected_boundaries=selected_boundaries,
        )
    except ValueError as exc:
        issues.append(str(exc))

    return {
        "ready": not issues,
        "issues": issues,
        "trial_dof_count": len(trial_dofs),
        "test_dof_count": len(test_dofs),
        "selected_boundary_count": len(selected_boundaries),
    }


def downstream_enabled(session_state: dict[str, object], required_keys: tuple[str, ...]) -> bool:
    """Return True only when all prerequisite artifacts exist in session state."""
    return all(key in session_state and session_state[key] is not None for key in required_keys)
