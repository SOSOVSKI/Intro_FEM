from __future__ import annotations

import sympy as sp

from symbolic_fem_workbench.workflow import build_poisson_triangle_p1_local_problem


def test_triangle_p1_unit_right_triangle() -> None:
    data = build_poisson_triangle_p1_local_problem()
    f = data["f"]

    Ke_expected = sp.Matrix(
        [
            [1, sp.Rational(-1, 2), sp.Rational(-1, 2)],
            [sp.Rational(-1, 2), sp.Rational(1, 2), 0],
            [sp.Rational(-1, 2), 0, sp.Rational(1, 2)],
        ]
    )
    fe_expected = sp.Matrix([sp.Rational(1, 6) * f, sp.Rational(1, 6) * f, sp.Rational(1, 6) * f])

    assert sp.simplify(data["Ke_unit_right_triangle"] - Ke_expected) == sp.zeros(3, 3)
    assert sp.simplify(data["fe_unit_right_triangle"] - fe_expected) == sp.zeros(3, 1)
