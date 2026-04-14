from __future__ import annotations

import numpy as np
import sympy as sp

from symbolic_fem_workbench.assembly import (
    apply_dirichlet_by_reduction,
    assemble_dense_matrix,
    assemble_dense_vector,
    expand_reduced_solution,
)
from symbolic_fem_workbench.workflow import build_poisson_triangle_p1_local_problem


def test_square_four_triangle_manual_assembly_center_value() -> None:
    data = build_poisson_triangle_p1_local_problem()
    x1, y1 = data["geometry"].x1, data["geometry"].y1
    x2, y2 = data["geometry"].x2, data["geometry"].y2
    x3, y3 = data["geometry"].x3, data["geometry"].y3
    f = data["f"]

    ke_fn = sp.lambdify((x1, y1, x2, y2, x3, y3), data["Ke"], "numpy")
    fe_fn = sp.lambdify((x1, y1, x2, y2, x3, y3, f), data["fe"], "numpy")

    nodes = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [0.5, 0.5],
        ],
        dtype=float,
    )
    elements = [(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)]

    K = np.zeros((5, 5), dtype=float)
    F = np.zeros(5, dtype=float)

    for conn in elements:
        coords = nodes[list(conn)]
        x1v, y1v = coords[0]
        x2v, y2v = coords[1]
        x3v, y3v = coords[2]
        K_local = np.asarray(ke_fn(x1v, y1v, x2v, y2v, x3v, y3v), dtype=float)
        F_local = np.asarray(fe_fn(x1v, y1v, x2v, y2v, x3v, y3v, 1.0), dtype=float).reshape(-1)
        assemble_dense_matrix(K, K_local, conn)
        assemble_dense_vector(F, F_local, conn)

    K_red, F_red, free = apply_dirichlet_by_reduction(K, F, constrained_dofs=[0, 1, 2, 3])
    u_free = np.linalg.solve(K_red, F_red)
    u = expand_reduced_solution(u_free, ndofs=5, free_dofs=free, constrained_dofs=[0, 1, 2, 3])

    assert np.allclose(u[:4], 0.0)
    assert np.isclose(u[4], 1.0 / 12.0)
