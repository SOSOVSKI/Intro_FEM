"""Tests for the elasticity module."""
from __future__ import annotations
import sympy as sp
from symbolic_fem_workbench.elasticity import (
    plane_stress_D, plane_strain_D, isotropic_3d_D,
    B_matrix_triangle_2d, B_matrix_tetra_3d,
    element_stiffness_BtDB,
)
from symbolic_fem_workbench import build_elasticity_triangle_p1_2d


def test_plane_stress_D_symmetry():
    E, nu = sp.symbols("E nu", positive=True)
    D = plane_stress_D(E, nu)
    assert D == D.T


def test_plane_strain_D_symmetry():
    E, nu = sp.symbols("E nu", positive=True)
    D = plane_strain_D(E, nu)
    assert D == D.T


def test_3d_D_symmetry():
    E, nu = sp.symbols("E nu", positive=True)
    D = isotropic_3d_D(E, nu)
    assert sp.simplify(D - D.T) == sp.zeros(6, 6)


def test_B_matrix_triangle_shape():
    dN_dx = [sp.Rational(-1, 1), sp.Rational(1, 1), sp.Integer(0)]
    dN_dy = [sp.Rational(-1, 1), sp.Integer(0), sp.Rational(1, 1)]
    B = B_matrix_triangle_2d(dN_dx, dN_dy)
    assert B.shape == (3, 6)


def test_B_matrix_tetra_shape():
    dN = [sp.Integer(i) for i in range(4)]
    B = B_matrix_tetra_3d(dN, dN, dN)
    assert B.shape == (6, 12)


def test_Ke_elasticity_symmetry():
    data = build_elasticity_triangle_p1_2d("plane_stress")
    Ke = data["Ke_unit_right_triangle"]
    assert sp.simplify(Ke - Ke.T) == sp.zeros(6, 6)


def test_Ke_elasticity_plane_strain():
    data = build_elasticity_triangle_p1_2d("plane_strain")
    Ke = data["Ke_unit_right_triangle"]
    assert Ke.shape == (6, 6)
    assert sp.simplify(Ke - Ke.T) == sp.zeros(6, 6)
