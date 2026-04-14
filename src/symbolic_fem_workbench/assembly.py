"""Minimal manual-assembly helpers for teaching examples.

These helpers are intentionally small and explicit. They support classroom examples
without trying to hide topology or boundary-condition reasoning.
"""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np


def assemble_dense_matrix(global_matrix: np.ndarray, local_matrix: np.ndarray, connectivity: Sequence[int]) -> None:
    """Add a local element matrix into a dense global matrix in place."""
    for a_local, a_global in enumerate(connectivity):
        for b_local, b_global in enumerate(connectivity):
            global_matrix[a_global, b_global] += local_matrix[a_local, b_local]



def assemble_dense_vector(global_vector: np.ndarray, local_vector: np.ndarray, connectivity: Sequence[int]) -> None:
    """Add a local element vector into a dense global vector in place."""
    local_vector = np.asarray(local_vector).reshape(-1)
    for a_local, a_global in enumerate(connectivity):
        global_vector[a_global] += local_vector[a_local]



def apply_dirichlet_by_reduction(
    global_matrix: np.ndarray,
    global_vector: np.ndarray,
    constrained_dofs: Sequence[int],
    prescribed_values: Sequence[float] | None = None,
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """Return the reduced linear system after enforcing Dirichlet values.

    Parameters
    ----------
    global_matrix, global_vector
        Full linear system ``K u = F``.
    constrained_dofs
        Global dof indices where the solution is prescribed.
    prescribed_values
        Values at the constrained dofs. Defaults to homogeneous conditions.
    """
    K = np.asarray(global_matrix, dtype=float)
    F = np.asarray(global_vector, dtype=float).reshape(-1)

    constrained = list(constrained_dofs)
    if prescribed_values is None:
        prescribed_values = [0.0 for _ in constrained]
    if len(prescribed_values) != len(constrained):
        raise ValueError("prescribed_values must have the same length as constrained_dofs")

    all_dofs = list(range(K.shape[0]))
    free = [i for i in all_dofs if i not in constrained]

    u_c = np.asarray(prescribed_values, dtype=float).reshape(-1)
    K_ff = K[np.ix_(free, free)]
    K_fc = K[np.ix_(free, constrained)]
    F_f = F[free] - K_fc @ u_c
    return K_ff, F_f, free



def expand_reduced_solution(
    free_solution: Sequence[float],
    ndofs: int,
    free_dofs: Sequence[int],
    constrained_dofs: Sequence[int],
    prescribed_values: Sequence[float] | None = None,
) -> np.ndarray:
    """Expand a reduced solution back into the full nodal vector."""
    full = np.zeros(ndofs, dtype=float)
    full[list(free_dofs)] = np.asarray(free_solution, dtype=float).reshape(-1)
    if prescribed_values is None:
        prescribed_values = [0.0 for _ in constrained_dofs]
    full[list(constrained_dofs)] = np.asarray(prescribed_values, dtype=float).reshape(-1)
    return full


def apply_dirichlet_by_row_substitution(
    global_matrix: np.ndarray,
    global_vector: np.ndarray,
    constrained_dofs: Sequence[int],
    prescribed_values: Sequence[float] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Enforce Dirichlet BCs by modifying the global system in-place (row substitution).

    This is the textbook approach that most students learn first:
    1. For each constrained DOF i with prescribed value g_i:
       - Zero out row i of K
       - Set K[i, i] = 1
       - Subtract K[:, i] * g_i from the RHS for all other rows
       - Zero out column i of K
       - Set F[i] = g_i

    After modification, the system K u = F can be solved directly and
    u[i] = g_i for every constrained DOF.

    Parameters
    ----------
    global_matrix, global_vector
        The assembled system K u = F.  **Modified in-place** and also returned.
    constrained_dofs
        Indices of DOFs with prescribed Dirichlet values.
    prescribed_values
        The prescribed values.  Defaults to zeros (homogeneous).

    Returns
    -------
    K_mod, F_mod : np.ndarray
        The modified system (same arrays, returned for convenience).
    """
    K = np.asarray(global_matrix, dtype=float)
    F = np.asarray(global_vector, dtype=float).reshape(-1)
    constrained = list(constrained_dofs)
    if prescribed_values is None:
        prescribed_values = [0.0] * len(constrained)

    for loc, dof in enumerate(constrained):
        g = prescribed_values[loc]
        # Subtract column contribution from RHS
        F -= K[:, dof] * g
        # Zero row and column
        K[dof, :] = 0.0
        K[:, dof] = 0.0
        # Set diagonal to 1
        K[dof, dof] = 1.0
        # Set RHS to prescribed value
        F[dof] = g

    return K, F


def apply_dirichlet_by_lifting(
    global_matrix: np.ndarray,
    global_vector: np.ndarray,
    constrained_dofs: Sequence[int],
    prescribed_values: Sequence[float] | None = None,
) -> tuple[np.ndarray, np.ndarray, list[int], np.ndarray]:
    """Enforce Dirichlet BCs using the lifting approach.

    The lifting approach decomposes the solution as u = u_0 + u_D, where:
    - u_D is a "lift" vector that satisfies the Dirichlet conditions
      (u_D[i] = g_i for constrained DOFs, zero elsewhere)
    - u_0 ∈ V_0 is the homogeneous part (u_0[i] = 0 for constrained DOFs)

    Substituting into a(u, v) = F(v):
        a(u_0 + u_D, v) = F(v)
        a(u_0, v) = F(v) - a(u_D, v)

    In matrix form this becomes:
        K_ff u_0f = F_f - K_fc g_c

    This is mathematically equivalent to the reduction approach but shows
    more clearly how the Dirichlet data enters as a correction to the RHS.

    The approach follows the treatment in Dokken's FEniCS tutorial:
    the lift u_D lives in the full space, the modified RHS is F - K u_D,
    and the constrained rows/columns are then eliminated.

    Parameters
    ----------
    global_matrix, global_vector
        Assembled system K u = F.  Not modified.
    constrained_dofs
        Indices of Dirichlet DOFs.
    prescribed_values
        Prescribed values at the Dirichlet DOFs.  Defaults to zeros.

    Returns
    -------
    K_ff : np.ndarray
        Reduced stiffness (free-free block).
    F_modified : np.ndarray
        Modified load vector = F_f - K_fc g_c.
    free_dofs : list[int]
    u_lift : np.ndarray
        The full lift vector (zero except at constrained DOFs).
    """
    K = np.asarray(global_matrix, dtype=float)
    F = np.asarray(global_vector, dtype=float).reshape(-1)
    ndofs = K.shape[0]

    constrained = list(constrained_dofs)
    if prescribed_values is None:
        prescribed_values = [0.0] * len(constrained)

    # Build the lift vector
    u_lift = np.zeros(ndofs, dtype=float)
    for loc, dof in enumerate(constrained):
        u_lift[dof] = prescribed_values[loc]

    # Compute modified RHS:  F_modified = F - K @ u_lift
    F_modified = F - K @ u_lift

    # Extract free DOFs
    free = [i for i in range(ndofs) if i not in constrained]
    K_ff = K[np.ix_(free, free)]
    F_f = F_modified[free]

    return K_ff, F_f, free, u_lift


__all__ = [
    "assemble_dense_matrix",
    "assemble_dense_vector",
    "apply_dirichlet_by_reduction",
    "expand_reduced_solution",
    "apply_dirichlet_by_row_substitution",
    "apply_dirichlet_by_lifting",
]
