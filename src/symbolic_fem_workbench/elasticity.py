"""Elasticity helpers for 2D and 3D vector FEM problems.

Provides constitutive matrices (D), strain-displacement matrices (B),
and element stiffness computation for plane stress, plane strain,
and full 3D isotropic linear elasticity using Voigt notation.
"""

from __future__ import annotations

import sympy as sp


# ---------------------------------------------------------------------------
# 2D constitutive matrices
# ---------------------------------------------------------------------------

def plane_stress_D(E: sp.Expr, nu: sp.Expr) -> sp.Matrix:
    """Constitutive matrix for plane stress (thin plates).

    Voigt ordering: [σ_xx, σ_yy, τ_xy]^T = D [ε_xx, ε_yy, γ_xy]^T

    D = E/(1-ν²) [[1, ν, 0], [ν, 1, 0], [0, 0, (1-ν)/2]]
    """
    coeff = E / (1 - nu**2)
    return coeff * sp.Matrix([
        [1, nu, 0],
        [nu, 1, 0],
        [0, 0, (1 - nu) / 2],
    ])


def plane_strain_D(E: sp.Expr, nu: sp.Expr) -> sp.Matrix:
    """Constitutive matrix for plane strain (thick/constrained bodies).

    Voigt ordering: [σ_xx, σ_yy, τ_xy]^T = D [ε_xx, ε_yy, γ_xy]^T

    D = E/((1+ν)(1-2ν)) [[1-ν, ν, 0], [ν, 1-ν, 0], [0, 0, (1-2ν)/2]]
    """
    coeff = E / ((1 + nu) * (1 - 2 * nu))
    return coeff * sp.Matrix([
        [1 - nu, nu, 0],
        [nu, 1 - nu, 0],
        [0, 0, (1 - 2 * nu) / 2],
    ])


# ---------------------------------------------------------------------------
# 3D constitutive matrix
# ---------------------------------------------------------------------------

def isotropic_3d_D(E: sp.Expr, nu: sp.Expr) -> sp.Matrix:
    """Full 3D isotropic constitutive matrix (6×6).

    Voigt ordering: [σ_xx, σ_yy, σ_zz, τ_yz, τ_xz, τ_xy]^T
                  = D [ε_xx, ε_yy, ε_zz, γ_yz, γ_xz, γ_xy]^T
    """
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    return sp.Matrix([
        [lam + 2*mu, lam, lam, 0, 0, 0],
        [lam, lam + 2*mu, lam, 0, 0, 0],
        [lam, lam, lam + 2*mu, 0, 0, 0],
        [0, 0, 0, mu, 0, 0],
        [0, 0, 0, 0, mu, 0],
        [0, 0, 0, 0, 0, mu],
    ])


# ---------------------------------------------------------------------------
# 2D B-matrix for triangles
# ---------------------------------------------------------------------------

def B_matrix_triangle_2d(
    dN_dx: list[sp.Expr],
    dN_dy: list[sp.Expr],
) -> sp.Matrix:
    """Build the 2D strain-displacement matrix B for a triangle element.

    For n_nodes nodes with physical-space shape function derivatives,
    the B matrix maps nodal displacements [u1, v1, u2, v2, ...]^T
    to Voigt strains [ε_xx, ε_yy, γ_xy]^T.

    B is 3 × (2 * n_nodes):
        B = [[dN1/dx, 0, dN2/dx, 0, ...],
             [0, dN1/dy, 0, dN2/dy, ...],
             [dN1/dy, dN1/dx, dN2/dy, dN2/dx, ...]]

    Parameters
    ----------
    dN_dx, dN_dy : list of sp.Expr
        Physical-space partial derivatives of shape functions.
    """
    n = len(dN_dx)
    B = sp.zeros(3, 2 * n)
    for i in range(n):
        B[0, 2*i] = dN_dx[i]
        B[1, 2*i + 1] = dN_dy[i]
        B[2, 2*i] = dN_dy[i]
        B[2, 2*i + 1] = dN_dx[i]
    return B


# ---------------------------------------------------------------------------
# 3D B-matrix for tetrahedra
# ---------------------------------------------------------------------------

def B_matrix_tetra_3d(
    dN_dx: list[sp.Expr],
    dN_dy: list[sp.Expr],
    dN_dz: list[sp.Expr],
) -> sp.Matrix:
    """Build the 3D strain-displacement matrix B for a tetrahedral element.

    Maps nodal displacements [u1, v1, w1, u2, v2, w2, ...]^T
    to Voigt strains [ε_xx, ε_yy, ε_zz, γ_yz, γ_xz, γ_xy]^T.

    B is 6 × (3 * n_nodes).
    """
    n = len(dN_dx)
    B = sp.zeros(6, 3 * n)
    for i in range(n):
        col = 3 * i
        B[0, col] = dN_dx[i]       # ε_xx
        B[1, col + 1] = dN_dy[i]   # ε_yy
        B[2, col + 2] = dN_dz[i]   # ε_zz
        B[3, col + 1] = dN_dz[i]   # γ_yz
        B[3, col + 2] = dN_dy[i]
        B[4, col] = dN_dz[i]       # γ_xz
        B[4, col + 2] = dN_dx[i]
        B[5, col] = dN_dy[i]       # γ_xy
        B[5, col + 1] = dN_dx[i]
    return B


# ---------------------------------------------------------------------------
# Element stiffness from B and D
# ---------------------------------------------------------------------------

def element_stiffness_BtDB(
    B: sp.Matrix,
    D: sp.Matrix,
    volume_or_area: sp.Expr,
    thickness: sp.Expr | None = None,
) -> sp.Matrix:
    """Compute element stiffness K_e = (t) * |Ω_e| * B^T D B.

    For constant-strain elements (P1 triangle, P1 tet) where B and D are
    constant over the element, the integral reduces to a simple product.

    Parameters
    ----------
    B : sp.Matrix
        Strain-displacement matrix.
    D : sp.Matrix
        Constitutive matrix.
    volume_or_area : sp.Expr
        Element area (2D) or volume (3D).
    thickness : sp.Expr or None
        Plate thickness for 2D problems (defaults to 1 if omitted).
    """
    t = thickness if thickness is not None else sp.Integer(1)
    return sp.simplify(t * volume_or_area * B.T * D * B)
