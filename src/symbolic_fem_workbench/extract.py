"""Helpers for extracting local element matrices and vectors from symbolic forms."""

from __future__ import annotations

import sympy as sp


def extract_coefficient_matrix(expr: sp.Expr, trial_dofs, test_dofs) -> sp.Matrix:
    """Extract the local element stiffness matrix from a discretized bilinear form. After Galerkin substitution, the bilinear form is a polynomial in trial DOFs (d₀, d₁, ...) and test DOFs (w₀, w₁, ...). This function extracts the coefficient of each (wᵢ, dⱼ) pair to build the element matrix K_e[i,j]. Parameters: expr - the discretized bilinear form expression, trial_dofs - list of trial DOF symbols, test_dofs - list of test DOF symbols. Returns a SymPy Matrix of shape (n_test, n_trial)."""
    if not trial_dofs or not test_dofs:
        raise ValueError("trial_dofs and test_dofs must be non-empty lists of SymPy symbols.")
    expanded = sp.expand(expr)
    rows = []
    for wt in test_dofs:
        row = []
        for du in trial_dofs:
            row.append(sp.simplify(expanded.coeff(wt).coeff(du)))
        rows.append(row)
    return sp.Matrix(rows)


def extract_coefficient_vector(expr: sp.Expr, test_dofs) -> sp.Matrix:
    """Extract the local element load vector from a discretized linear form. After Galerkin substitution, the linear form is a polynomial in test DOFs (w₀, w₁, ...). This function extracts the coefficient of each wᵢ to build the element load vector f_e[i]. Parameters: expr - the discretized linear form expression, test_dofs - list of test DOF symbols. Returns a SymPy column Matrix of length n_test."""
    if not test_dofs:
        raise ValueError("test_dofs must be a non-empty list of SymPy symbols.")
    expanded = sp.expand(expr)
    return sp.Matrix([sp.simplify(expanded.coeff(wt)) for wt in test_dofs])
