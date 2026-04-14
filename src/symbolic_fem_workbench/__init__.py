"""symbolic_fem_workbench

A teaching-first symbolic FEM helper package.
"""

from .symbols import Domain1D, Domain2D, make_field_1d, make_field_2d
from .forms import DomainIntegral, BoundaryContribution, WeightedResidual, WeakForm
from .fe_spaces import LinearElement1D, local_trial_expansion, local_test_expansion
from .reference import ReferenceIntervalP1, ReferenceTriangleP1, AffineTriangleMap2D, ReferenceTetrahedronP1
from .quadrature import (
    integrate_reference_triangle_exact,
    triangle_one_point_rule,
    triangle_three_point_rule,
    triangle_six_point_rule,
)
from .transforms import (
    weighted_residual,
    integrate_divergence_1d,
    drop_dirichlet_boundary,
    apply_neumann_flux,
    split_linear_weak_form,
    grad_2d,
    pullback_gradient_2d,
    substitute_field,
    substitute_fe,
    gateaux_derivative,
    integrate_boundary_edge_1d,
    neumann_load_vector_triangle_edge,
)
from .extract import extract_coefficient_matrix, extract_coefficient_vector
from .workflow import (
    build_bar_1d_local_problem,
    build_bar_1d_mass_matrix,
    build_poisson_triangle_p1_local_problem,
    build_elasticity_tetra_p1_3d,
)
from .assembly import (
    assemble_dense_matrix,
    assemble_dense_vector,
    apply_dirichlet_by_reduction,
    expand_reduced_solution,
    apply_dirichlet_by_row_substitution,
    apply_dirichlet_by_lifting,
)
from .elasticity import (
    plane_stress_D,
    plane_strain_D,
    isotropic_3d_D,
    B_matrix_triangle_2d,
    B_matrix_tetra_3d,
    element_stiffness_BtDB,
)
from .workflow import build_elasticity_triangle_p1_2d
from . import viz

__all__ = [
    "Domain1D",
    "Domain2D",
    "make_field_1d",
    "make_field_2d",
    "DomainIntegral",
    "BoundaryContribution",
    "WeightedResidual",
    "WeakForm",
    "LinearElement1D",
    "ReferenceIntervalP1",
    "ReferenceTriangleP1",
    "ReferenceTetrahedronP1",
    "AffineTriangleMap2D",
    "local_trial_expansion",
    "local_test_expansion",
    "integrate_reference_triangle_exact",
    "triangle_one_point_rule",
    "triangle_three_point_rule",
    "triangle_six_point_rule",
    "weighted_residual",
    "integrate_divergence_1d",
    "drop_dirichlet_boundary",
    "apply_neumann_flux",
    "split_linear_weak_form",
    "grad_2d",
    "pullback_gradient_2d",
    "substitute_field",
    "substitute_fe",
    "gateaux_derivative",
    "integrate_boundary_edge_1d",
    "neumann_load_vector_triangle_edge",
    "extract_coefficient_matrix",
    "extract_coefficient_vector",
    "build_bar_1d_local_problem",
    "build_bar_1d_mass_matrix",
    "build_poisson_triangle_p1_local_problem",
    "build_elasticity_tetra_p1_3d",
    "assemble_dense_matrix",
    "assemble_dense_vector",
    "apply_dirichlet_by_reduction",
    "expand_reduced_solution",
    "apply_dirichlet_by_row_substitution",
    "apply_dirichlet_by_lifting",
    "plane_stress_D",
    "plane_strain_D",
    "isotropic_3d_D",
    "B_matrix_triangle_2d",
    "B_matrix_tetra_3d",
    "element_stiffness_BtDB",
    "build_elasticity_triangle_p1_2d",
    "viz",
]
