"""Explicit symbolic transformations for the teaching workflow."""

from __future__ import annotations

from typing import Mapping, Tuple

import sympy as sp

from .forms import BoundaryContribution, DomainIntegral, WeakForm, WeightedResidual
from .symbols import Domain1D
from .validate import ensure_same_variable, field_dependency, split_terms


def weighted_residual(lhs: sp.Expr, rhs: sp.Expr, test: sp.Expr, domain: Domain1D) -> WeightedResidual:
    """Construct the weighted residual from a strong-form equation. Given a strong form lhs = rhs, multiplies the residual (lhs - rhs) by the test function and integrates over the domain. This is the first step in deriving the weak form. Parameters: lhs - left-hand side of the strong form, rhs - right-hand side, test - test (weight) function, domain - 1D domain tuple (x, x_min, x_max). Returns a WeightedResidual container."""
    return WeightedResidual(domain_integrals=(DomainIntegral(sp.expand((lhs - rhs) * test), domain),))


def integrate_divergence_1d(flux: sp.Expr, test: sp.Expr, x: sp.Symbol, domain: Domain1D) -> Tuple[DomainIntegral, Tuple[BoundaryContribution, BoundaryContribution]]:
    """Apply integration by parts to a flux-divergence term in 1D. Converts ∫ (d(flux)/dx)·v dx into -∫ flux·(dv/dx) dx plus boundary terms [flux·v] evaluated at both ends. This reduces the smoothness requirement on the trial function by shifting one derivative onto the test function. Parameters: flux - the flux expression (e.g. E·A·du/dx), test - the test function v, x - the spatial coordinate symbol, domain - 1D domain tuple. Returns (domain_integral, (boundary_left, boundary_right))."""
    ensure_same_variable(domain[0], x)
    x_left = domain[1]
    x_right = domain[2]
    domain_term = DomainIntegral(sp.expand(flux * sp.diff(test, x)), domain)
    boundary_left = BoundaryContribution(-flux * test, x, x_left, label="left")
    boundary_right = BoundaryContribution(-flux * test, x, x_right, label="right")
    return domain_term, (boundary_left, boundary_right)


def drop_dirichlet_boundary(boundary_term: BoundaryContribution, test: sp.Expr) -> sp.Expr:
    """Enforce a homogeneous Dirichlet condition by setting the test function to zero at a boundary. In FEM, the test space vanishes on the essential (Dirichlet) boundary. This function evaluates the boundary term and substitutes v=0 at the specified location. Parameters: boundary_term - a BoundaryContribution object, test - the test function expression. Returns the simplified (typically zero) boundary contribution."""
    value = boundary_term.evaluate()
    return sp.simplify(value.subs(test.subs(boundary_term.x, boundary_term.location), 0))


def apply_neumann_flux(boundary_term: BoundaryContribution, flux_expr: sp.Expr, prescribed_flux: sp.Expr) -> sp.Expr:
    """Apply a prescribed Neumann (natural) boundary condition. Substitutes the known flux value into a boundary term and evaluates at the boundary location. Parameters: boundary_term - a BoundaryContribution object, flux_expr - the symbolic flux to replace, prescribed_flux - the known boundary flux value. Returns the evaluated boundary contribution."""
    substituted = boundary_term.expr.subs(flux_expr, prescribed_flux)
    return sp.simplify(substituted.subs(boundary_term.x, boundary_term.location))


def split_linear_weak_form(expr: sp.Expr, trial: sp.Expr, test: sp.Expr) -> WeakForm:
    """Split a weak-form expression into its bilinear and linear parts. Identifies which terms depend on both trial and test functions (bilinear form a(u,v)) and which depend only on the test function (linear form F(v), moved to the RHS). Parameters: expr - the combined weak form expression, trial - the trial field, test - the test field. Returns a WeakForm with .bilinear and .linear attributes."""
    bilinear = sp.Integer(0)
    linear_rhs = sp.Integer(0)
    for term in split_terms(expr):
        has_trial = field_dependency(term, trial)
        has_test = field_dependency(term, test)
        if has_trial and has_test:
            bilinear += term
        elif has_test and not has_trial:
            linear_rhs -= term
    return WeakForm(sp.simplify(bilinear), sp.simplify(linear_rhs), trial, test)


def grad_2d(expr: sp.Expr, x: sp.Symbol, y: sp.Symbol) -> sp.Matrix:
    """Return the 2D gradient as a column vector."""
    return sp.Matrix([[sp.diff(expr, x)], [sp.diff(expr, y)]])


def pullback_gradient_2d(gradient_ref: sp.Matrix, jacobian: sp.Matrix) -> sp.Matrix:
    """Map a reference-space gradient into physical coordinates using J^{-T}."""
    return sp.simplify(jacobian.inv().T * gradient_ref)


def substitute_field(expr: sp.Expr, field: sp.Expr, replacement: sp.Expr, *coordinates: sp.Symbol) -> sp.Expr:
    """Replace a symbolic field with a finite-element expansion in an expression. Handles derivatives of all orders: finds the maximum derivative order present, then substitutes from highest to lowest to avoid partial matches. This implements the Galerkin discretization step u(x) → Σ Nᵢ dᵢ. Parameters: expr - the expression to substitute into, field - the symbolic field (e.g. u(x)), replacement - the FE expansion, *coordinates - the spatial coordinate symbols. Returns the expanded expression with all field references replaced."""
    if not coordinates:
        raise ValueError(
            "At least one coordinate symbol must be provided. "
            "Example: substitute_field(expr, u, u_h, x) for 1D, "
            "or substitute_field(expr, u, u_h, x, y) for 2D."
        )
    field_func = getattr(field, "func", None)
    if field_func is None and not isinstance(field, sp.Expr):
        raise TypeError(
            f"Expected a SymPy expression or Function application for 'field', "
            f"got {type(field).__name__}. Use make_field_1d('u', x) to create fields."
        )
    max_orders: dict[sp.Symbol, int] = {coord: 0 for coord in coordinates}

    for node in sp.preorder_traversal(expr):
        if isinstance(node, sp.Derivative) and (node.expr == field or (field_func is not None and getattr(node.expr, "func", None) == field_func)):
            for var, count in node.variable_count:
                if var in max_orders:
                    max_orders[var] = max(max_orders[var], count)

    substituted = expr
    for coord in coordinates:
        for order in range(max_orders[coord], 0, -1):
            substituted = substituted.subs(sp.diff(field, coord, order), sp.diff(replacement, coord, order))

    if field_func is not None:
        substituted = substituted.replace(
            lambda node: getattr(node, "is_Function", False) and getattr(node, "func", None) == field_func,
            lambda node: replacement.subs(dict(zip(coordinates, node.args))) if len(node.args) == len(coordinates) else replacement,
        )
    else:
        substituted = substituted.subs(field, replacement)

    return sp.expand(substituted.doit())


def substitute_fe(expr: sp.Expr, replacements: Mapping[sp.Expr, sp.Expr], *coordinates: sp.Symbol) -> sp.Expr:
    """Apply multiple field substitutions for Galerkin discretization. Convenience wrapper that calls substitute_field for each (field, replacement) pair. Typically used to substitute both trial u→u_h and test v→v_h simultaneously. Parameters: expr - the expression, replacements - dict mapping fields to their FE expansions, *coordinates - spatial coordinates. Returns the fully discretized expression."""
    out = expr
    for field, replacement in replacements.items():
        out = substitute_field(out, field, replacement, *coordinates)
    return sp.expand(out)


def gateaux_derivative(form: sp.Expr, trial_var: sp.Expr, increment_var: sp.Expr) -> sp.Expr:
    """Compute the Gateaux (directional) derivative for linearization. Evaluates d/dε[form(u + ε·δu)]|_{ε=0}. Used in Newton-Raphson linearization of nonlinear problems. Parameters: form - the nonlinear form expression, trial_var - the current solution variable, increment_var - the increment/perturbation variable. Returns the linearized expression."""
    eps = sp.Symbol("epsilon")
    u_eps = trial_var + eps * increment_var
    form_eps = form.subs(trial_var, u_eps)
    derivative = sp.diff(form_eps, eps).subs(eps, 0)
    return sp.expand(derivative)


def integrate_boundary_edge_1d(
    integrand: sp.Expr,
    t: sp.Symbol,
    t_start: sp.Expr,
    t_end: sp.Expr,
) -> sp.Expr:
    """Integrate an expression along a 1D parametric edge.

    In 2D FEM, boundary integrals arise from integration by parts of
    the flux-divergence term.  When the divergence theorem converts a
    domain integral to a boundary integral, the result is a line integral
    along each boundary edge:

        ∫_Γ (flux · n) v ds

    For a straight edge parameterised by t ∈ [t_start, t_end], this
    function evaluates ∫ integrand dt.

    Parameters
    ----------
    integrand : sp.Expr
        The boundary integrand (already including the edge Jacobian / length
        factor if the parameterisation is not unit-speed).
    t : sp.Symbol
        The parameterisation variable.
    t_start, t_end : sp.Expr
        Limits of integration.

    Returns
    -------
    sp.Expr
        The simplified result of the line integral.
    """
    return sp.simplify(sp.integrate(integrand, (t, t_start, t_end)))


def neumann_load_vector_triangle_edge(
    ref_element,
    edge_local_nodes: tuple[int, int],
    prescribed_flux: sp.Expr,
    edge_length: sp.Expr,
    xi: sp.Symbol,
    eta: sp.Symbol,
) -> sp.Matrix:
    """Compute the Neumann contribution to the element load vector for one edge.

    For a P1 triangle with a Neumann BC t* on one edge, the boundary
    integral contributes:

        f_e[i] += ∫_edge t* N_i ds

    On a straight edge of length |e|, with a linear parameterisation,
    the shape functions restrict to 1D linear functions on the edge.
    For constant t*, the result is  t* |e| / 2  for each of the two
    edge nodes, and 0 for the interior node.

    Parameters
    ----------
    ref_element
        Reference element (e.g. ReferenceTriangleP1) with .shape_functions.
    edge_local_nodes : tuple[int, int]
        Local node indices (0-based) of the two nodes on the boundary edge.
    prescribed_flux : sp.Expr
        The prescribed Neumann flux value t*.
    edge_length : sp.Expr
        Physical length of the edge.
    xi, eta : sp.Symbol
        Reference coordinate symbols.

    Returns
    -------
    sp.Matrix
        Column vector of length n_nodes with the boundary contributions.
    """
    n_nodes = len(ref_element.shape_functions)
    fe_boundary = sp.zeros(n_nodes, 1)

    # For P1 elements with constant flux on a straight edge,
    # each edge node gets t* * |edge| / 2
    n0, n1 = edge_local_nodes
    fe_boundary[n0, 0] = prescribed_flux * edge_length / 2
    fe_boundary[n1, 0] = prescribed_flux * edge_length / 2

    return fe_boundary
