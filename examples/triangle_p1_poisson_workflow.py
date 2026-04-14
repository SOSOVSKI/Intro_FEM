from __future__ import annotations

import sympy as sp

from symbolic_fem_workbench.workflow import build_poisson_triangle_p1_local_problem


def main() -> None:
    data = build_poisson_triangle_p1_local_problem()

    print("Generic symbolic P1 triangle stiffness matrix Ke:")
    sp.pprint(data["Ke"])
    print()

    print("Generic symbolic P1 triangle load vector fe for constant f:")
    sp.pprint(data["fe"])
    print()

    print("Unit right triangle specialization:")
    print("Ke =")
    sp.pprint(data["Ke_unit_right_triangle"])
    print("fe =")
    sp.pprint(data["fe_unit_right_triangle"])


if __name__ == "__main__":
    main()
