"""Tests for 3D elasticity."""
from __future__ import annotations
import sympy as sp
from symbolic_fem_workbench.elasticity import (
    isotropic_3d_D, B_matrix_tetra_3d,
)
from symbolic_fem_workbench.reference import ReferenceTetrahedronP1


def test_3d_D_positive_diagonal():
    """Test that the 3D constitutive matrix has positive diagonal."""
    E, nu = sp.symbols("E nu", positive=True)
    D = isotropic_3d_D(E, nu)
    # For physically valid parameters, diagonal should be positive
    # Test with nu = 1/4 (Poisson's ratio)
    D_num = D.subs({E: 1, nu: sp.Rational(1, 4)})
    for i in range(6):
        assert D_num[i, i] > 0, f"Diagonal element {i} is not positive"


def test_3d_D_symmetry():
    """Test that the 3D constitutive matrix is symmetric."""
    E, nu = sp.symbols("E nu", positive=True)
    D = isotropic_3d_D(E, nu)
    assert sp.simplify(D - D.T) == sp.zeros(6, 6), "D matrix is not symmetric"


def test_3d_B_matrix_tetra_shape():
    """Test that B matrix has correct shape (6×12) for tetrahedron."""
    dN = [sp.Integer(i) for i in range(4)]
    B = B_matrix_tetra_3d(dN, dN, dN)
    assert B.shape == (6, 12), f"Expected B shape (6, 12), got {B.shape}"


def test_reference_tetra_p1_nodes():
    """Test that reference tetrahedron has correct node coordinates."""
    xi, eta, zeta = sp.symbols("xi eta zeta", real=True)
    ref = ReferenceTetrahedronP1(xi=xi, eta=eta, zeta=zeta)
    expected = ((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
    assert ref.nodes == expected, f"Reference tet nodes mismatch: {ref.nodes}"


def test_reference_tetra_p1_shape_functions():
    """Test that reference tetrahedron shape functions sum to 1."""
    xi, eta, zeta = sp.symbols("xi eta zeta", real=True)
    ref = ReferenceTetrahedronP1(xi=xi, eta=eta, zeta=zeta)
    N = ref.shape_functions
    assert len(N) == 4, f"Expected 4 shape functions, got {len(N)}"

    # Sum should be 1
    N_sum = sp.simplify(sum(N))
    assert N_sum == 1, f"Shape functions do not sum to 1: {N_sum}"


def test_reference_tetra_p1_gradients():
    """Test that reference tetrahedron has 4 gradient vectors of shape (3,1)."""
    xi, eta, zeta = sp.symbols("xi eta zeta", real=True)
    ref = ReferenceTetrahedronP1(xi=xi, eta=eta, zeta=zeta)
    grads = ref.shape_gradients_reference
    assert len(grads) == 4, f"Expected 4 gradients, got {len(grads)}"
    for i, grad in enumerate(grads):
        assert grad.shape == (3, 1), f"Gradient {i} has wrong shape: {grad.shape}"


def test_3d_D_size():
    """Test that 3D D matrix is 6×6."""
    E, nu = sp.symbols("E nu", positive=True)
    D = isotropic_3d_D(E, nu)
    assert D.shape == (6, 6), f"Expected 6×6 D, got {D.shape}"
