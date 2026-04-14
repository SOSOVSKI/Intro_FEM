"""Tests for features added in phases 1-4."""
from __future__ import annotations

import numpy as np
import sympy as sp
import pytest

from symbolic_fem_workbench import (
    build_bar_1d_mass_matrix,
    triangle_three_point_rule,
    triangle_six_point_rule,
    substitute_field,
    extract_coefficient_matrix,
    extract_coefficient_vector,
    make_field_1d,
)
from symbolic_fem_workbench.reference import ReferenceIntervalP1


# --- Mass matrix ---

def test_mass_matrix_symmetry_and_values():
    m = build_bar_1d_mass_matrix()
    Me = m["Me"]
    L, rho, A = m["L"], m["rho"], m["A"]
    # Consistent mass matrix for P1 bar: rho*A*L/6 * [[2,1],[1,2]]
    expected = sp.Matrix([[rho*A*L/3, rho*A*L/6], [rho*A*L/6, rho*A*L/3]])
    assert sp.simplify(Me - expected) == sp.zeros(2, 2)


def test_mass_matrix_row_sum():
    """Row sum of consistent mass matrix equals rho*A*L/2 (lumped mass)."""
    m = build_bar_1d_mass_matrix()
    Me = m["Me"]
    L, rho, A = m["L"], m["rho"], m["A"]
    row_sum = sp.simplify(Me[0, 0] + Me[0, 1])
    assert sp.simplify(row_sum - rho * A * L / 2) == 0


# --- Triangle quadrature ---

def test_triangle_three_point_exact_for_quadratic():
    xi, eta = sp.symbols("xi eta")
    f = xi**2 + eta**2 + xi * eta
    exact = sp.integrate(f, (eta, 0, 1 - xi), (xi, 0, 1))
    approx = triangle_three_point_rule(f, xi, eta)
    assert sp.simplify(exact - approx) == 0


def test_triangle_three_point_exact_for_linear():
    xi, eta = sp.symbols("xi eta")
    f = 3 * xi + 5 * eta + 7
    exact = sp.integrate(f, (eta, 0, 1 - xi), (xi, 0, 1))
    approx = triangle_three_point_rule(f, xi, eta)
    assert sp.simplify(exact - approx) == 0


# --- ReferenceIntervalP1 ---

def test_reference_interval_p1_partition_of_unity():
    xi = sp.Symbol("xi")
    ref = ReferenceIntervalP1(xi)
    assert sp.simplify(sum(ref.shape_functions) - 1) == 0


def test_reference_interval_p1_kronecker():
    xi = sp.Symbol("xi")
    ref = ReferenceIntervalP1(xi)
    for i, node in enumerate(ref.nodes):
        for j, N in enumerate(ref.shape_functions):
            val = sp.simplify(N.subs(xi, node[0]))
            expected = 1 if i == j else 0
            assert val == expected, f"N{j+1}(node{i+1}) = {val}, expected {expected}"


def test_reference_interval_p1_gradient_sum_zero():
    xi = sp.Symbol("xi")
    ref = ReferenceIntervalP1(xi)
    grad_sum = sum(ref.shape_gradients_reference, sp.zeros(1, 1))
    assert grad_sum == sp.zeros(1, 1)


# --- Input validation ---

def test_substitute_field_requires_coordinates():
    x = sp.Symbol("x")
    u = make_field_1d("u", x)
    with pytest.raises(ValueError, match="At least one coordinate"):
        substitute_field(u, u, sp.Symbol("d0"))


def test_extract_coefficient_matrix_requires_nonempty():
    with pytest.raises(ValueError, match="non-empty"):
        extract_coefficient_matrix(sp.Integer(0), [], [sp.Symbol("w0")])


def test_extract_coefficient_vector_requires_nonempty():
    with pytest.raises(ValueError, match="non-empty"):
        extract_coefficient_vector(sp.Integer(0), [])
