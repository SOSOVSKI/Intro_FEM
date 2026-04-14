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


def main() -> None:
    data = build_poisson_triangle_p1_local_problem()

    # Turn the symbolic local tensors into callable numerical kernels.
    x1, y1 = data["geometry"].x1, data["geometry"].y1
    x2, y2 = data["geometry"].x2, data["geometry"].y2
    x3, y3 = data["geometry"].x3, data["geometry"].y3
    f = data["f"]

    ke_fn = sp.lambdify((x1, y1, x2, y2, x3, y3), data["Ke"], "numpy")
    fe_fn = sp.lambdify((x1, y1, x2, y2, x3, y3, f), data["fe"], "numpy")

    # Tiny mesh: unit square split into four triangles around a center node.
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
    elements = [
        (0, 1, 4),
        (1, 2, 4),
        (2, 3, 4),
        (3, 0, 4),
    ]

    ndofs = len(nodes)
    K_global = np.zeros((ndofs, ndofs), dtype=float)
    F_global = np.zeros(ndofs, dtype=float)

    source_value = 1.0

    for elem_id, conn in enumerate(elements):
        coords = nodes[list(conn)]
        x1v, y1v = coords[0]
        x2v, y2v = coords[1]
        x3v, y3v = coords[2]

        K_local = np.asarray(ke_fn(x1v, y1v, x2v, y2v, x3v, y3v), dtype=float)
        F_local = np.asarray(fe_fn(x1v, y1v, x2v, y2v, x3v, y3v, source_value), dtype=float).reshape(-1)

        print(f"Element {elem_id}, connectivity {conn}")
        print("K_local =")
        print(K_local)
        print("F_local =")
        print(F_local)
        print()

        assemble_dense_matrix(K_global, K_local, conn)
        assemble_dense_vector(F_global, F_local, conn)

    print("Assembled global matrix K:")
    print(K_global)
    print()
    print("Assembled global vector F:")
    print(F_global)
    print()

    # Homogeneous Dirichlet conditions on the outer boundary nodes.
    constrained = [0, 1, 2, 3]
    K_red, F_red, free = apply_dirichlet_by_reduction(K_global, F_global, constrained)

    print("Reduced system K_ff u_f = F_f")
    print("K_ff =")
    print(K_red)
    print("F_f =")
    print(F_red)
    print()

    u_free = np.linalg.solve(K_red, F_red)
    u = expand_reduced_solution(u_free, ndofs=ndofs, free_dofs=free, constrained_dofs=constrained)

    print("Solved nodal values u =")
    print(u)
    print()
    print("Center-node value u[4] =", u[4])
    print("Expected for this tiny mesh with f = 1 and zero boundary values: 1/12 =", 1.0 / 12.0)


if __name__ == "__main__":
    main()
