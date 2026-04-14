"""Small reference-element library demo.

This script prints a compact summary of the currently supported reference
families. It is meant for inspection in class or in a notebook, not for heroic
numerical work.
"""

from __future__ import annotations

import sympy as sp

from symbolic_fem_workbench.reference import (
    ReferenceTriangleP1,
    ReferenceQuadrilateralQ1,
    ReferenceQuadrilateralQ2,
    ReferenceTetrahedronP1,
    ReferenceTetrahedronP2,
    ReferenceHexahedronQ1,
)



def summarize_element(name: str, element) -> None:
    print(f"\n{name}")
    print("-" * len(name))
    print(f"number of nodes: {len(element.nodes)}")
    print("nodes:")
    for i, node in enumerate(element.nodes, start=1):
        print(f"  {i}: {node}")
    print("shape functions:")
    for i, N in enumerate(element.shape_functions, start=1):
        print(f"  N{i} = {sp.simplify(N)}")


if __name__ == "__main__":
    xi, eta, zeta = sp.symbols("xi eta zeta")

    summarize_element("Triangle P1", ReferenceTriangleP1(xi, eta))
    summarize_element("Quadrilateral Q1", ReferenceQuadrilateralQ1(xi, eta))
    summarize_element("Quadrilateral Q2", ReferenceQuadrilateralQ2(xi, eta))
    summarize_element("Tetrahedron P1 / Tet4", ReferenceTetrahedronP1(xi, eta, zeta))
    summarize_element("Tetrahedron P2 / Tet10", ReferenceTetrahedronP2(xi, eta, zeta))
    summarize_element("Hexahedron Q1 / Hex8", ReferenceHexahedronQ1(xi, eta, zeta))
