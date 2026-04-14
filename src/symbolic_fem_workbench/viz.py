"""Visualization helpers for teaching FEM concepts.

Provides matplotlib-based functions for plotting shape functions,
reference-to-physical mappings, and FEM solution fields.
"""

from __future__ import annotations

import numpy as np
import sympy as sp


def plot_triangle_shape_functions(reference_element, xi_sym, eta_sym, ax=None):
    """Plot 2D contour plots of triangle shape functions on the reference element.

    Parameters
    ----------
    reference_element
        A reference triangle object with .shape_functions property.
    xi_sym, eta_sym
        SymPy symbols for the reference coordinates.
    ax
        Optional matplotlib axes (creates new figure if None).

    Returns
    -------
    fig : matplotlib Figure
    """
    import matplotlib.pyplot as plt

    n_funcs = len(reference_element.shape_functions)
    if ax is None:
        fig, axes = plt.subplots(1, n_funcs, figsize=(4 * n_funcs, 4))
        if n_funcs == 1:
            axes = [axes]
    else:
        fig = ax[0].get_figure() if hasattr(ax, '__len__') else ax.get_figure()
        axes = ax if hasattr(ax, '__len__') else [ax]

    xi_vals = np.linspace(0, 1, 60)
    eta_vals = np.linspace(0, 1, 60)
    XI, ETA = np.meshgrid(xi_vals, eta_vals)
    mask = XI + ETA <= 1.0

    for idx, (a, N_sym) in enumerate(zip(axes, reference_element.shape_functions)):
        N_fn = sp.lambdify((xi_sym, eta_sym), N_sym, "numpy")
        Z = N_fn(XI, ETA)
        Z = np.where(mask, Z, np.nan)
        c = a.contourf(XI, ETA, Z, levels=20, cmap="viridis")
        a.set_title(f"$N_{{{idx + 1}}}$", fontsize=14)
        a.set_xlabel(r"$\xi$")
        a.set_ylabel(r"$\eta$")
        a.set_aspect("equal")
        plt.colorbar(c, ax=a, shrink=0.7)

    plt.tight_layout()
    return fig


def plot_affine_mapping(geom, xi_sym, eta_sym, vertex_subs, ax=None):
    """Visualize the affine mapping from reference to physical triangle.

    Parameters
    ----------
    geom : AffineTriangleMap2D
    xi_sym, eta_sym : sp.Symbol
    vertex_subs : dict
        Substitution dict mapping symbolic vertices to numerical values,
        e.g. {x1: 0, y1: 0, x2: 2, y2: 0, x3: 1, y3: 1.5}.
    ax : optional matplotlib axes pair (ax_ref, ax_phys)

    Returns
    -------
    fig : matplotlib Figure
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, (ax_ref, ax_phys) = plt.subplots(1, 2, figsize=(10, 4))
    else:
        ax_ref, ax_phys = ax
        fig = ax_ref.get_figure()

    # Reference triangle
    ref_pts = np.array([[0, 0], [1, 0], [0, 1], [0, 0]])
    ax_ref.fill(ref_pts[:3, 0], ref_pts[:3, 1], alpha=0.3, color="steelblue")
    ax_ref.plot(ref_pts[:, 0], ref_pts[:, 1], "b-o", markersize=8)
    ax_ref.set_title("Reference Triangle")
    ax_ref.set_xlabel(r"$\xi$")
    ax_ref.set_ylabel(r"$\eta$")
    ax_ref.set_aspect("equal")

    # Physical triangle
    x_map_fn = sp.lambdify((xi_sym, eta_sym), geom.x_map.subs(vertex_subs), "numpy")
    y_map_fn = sp.lambdify((xi_sym, eta_sym), geom.y_map.subs(vertex_subs), "numpy")

    ref_corners_xi = [0, 1, 0]
    ref_corners_eta = [0, 0, 1]
    phys_x = [float(x_map_fn(xi_c, eta_c)) for xi_c, eta_c in zip(ref_corners_xi, ref_corners_eta)]
    phys_y = [float(y_map_fn(xi_c, eta_c)) for xi_c, eta_c in zip(ref_corners_xi, ref_corners_eta)]

    phys_pts = np.array(list(zip(phys_x, phys_y)) + [(phys_x[0], phys_y[0])])
    ax_phys.fill(phys_pts[:3, 0], phys_pts[:3, 1], alpha=0.3, color="coral")
    ax_phys.plot(phys_pts[:, 0], phys_pts[:, 1], "r-o", markersize=8)
    ax_phys.set_title("Physical Triangle")
    ax_phys.set_xlabel("x")
    ax_phys.set_ylabel("y")
    ax_phys.set_aspect("equal")

    plt.tight_layout()
    return fig


def plot_mesh_solution(nodes, elements, u, ax=None, title="FEM Solution"):
    """Plot a scalar FEM solution on a triangular mesh.

    Parameters
    ----------
    nodes : ndarray of shape (n_nodes, 2)
    elements : list of 3-tuples (connectivity)
    u : ndarray of shape (n_nodes,) — nodal values
    ax : optional matplotlib axes
    title : str

    Returns
    -------
    fig : matplotlib Figure
    """
    import matplotlib.pyplot as plt
    from matplotlib.tri import Triangulation

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    else:
        fig = ax.get_figure()

    tri_conn = np.array(elements)
    triang = Triangulation(nodes[:, 0], nodes[:, 1], tri_conn)

    tpc = ax.tripcolor(triang, u, shading="flat", cmap="viridis")
    ax.triplot(triang, "k-", linewidth=0.5)
    for i, (x, y) in enumerate(nodes):
        ax.annotate(f"{i}", (x, y), fontsize=9,
                    xytext=(3, 3), textcoords="offset points")
    plt.colorbar(tpc, ax=ax, label="u")
    ax.set_title(title)
    ax.set_aspect("equal")
    return fig
