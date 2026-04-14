#%%
# 2D Timoshenko Beam Example in FEniCSx
#

#
import dolfinx
from dolfinx import fem, mesh, io
from mpi4py import MPI
import ufl
import numpy as np
import pyvista
import basix

# --- 1. Problem Parameters ---
# Geometry
L = 10.0  # Length of the beam
n_elem = 10 # Number of beam elements
h = 0.5   # Height of the beam cross-section
b = 0.2   # Width of the beam cross-section (out-of-plane)

# Material Properties (Steel)
E = 70e9 # Young's modulus (Pa)
nu = 0.3   # Poisson's ratio
G = E / (2 * (1 + nu)) # Shear modulus (μ)

# Load
q_val = -4000.0 # Uniformly distributed load in y-direction (N/m)

# Cross-section properties
S = b * h         # Area
I = b * h**3 / 12 # Second moment of area (for bending in x-y plane)
kappa = 5.0 / 6.0   # Shear correction factor for rectangular cross-section

# --- 2. Mesh Generation ---
# Create a 1D mesh embedded in a 2D space (gdim=2, tdim=1)
domain = mesh.create_interval(MPI.COMM_WORLD, n_elem, [0, L])

# --- 3. Function Space Definition ---
# We use a mixed formulation with two independent fields:
# - Displacement `u`: A 2D vector field.
# - Rotation `theta`: A scalar field.
Ue = basix.ufl.element("P", domain.basix_cell(), 1, shape=(2,))
Te = basix.ufl.element("P", domain.basix_cell(), 1)
#%%
W = fem.functionspace(domain, basix.ufl.mixed_element([Ue, Te]))

# create subspaces for displacement and rotation
V_u, _ = W.sub(0).collapse()
V_theta, _ = W.sub(1).collapse()
#%%
# --- 4. Variational Formulation ---
# Define trial and test functions for the mixed system
(du, dtheta) = ufl.TrialFunctions(W)
(v_u, v_theta) = ufl.TestFunctions(W)

# Define generalized strains for the 2D Timoshenko beam
def get_strains(u, theta):
    # Axial strain (stretching)
    delta = u[0].dx(0)
    # Bending curvature
    chi = theta.dx(0)
    # Shear strain
    gamma = u[1].dx(0) - theta
    return ufl.as_vector([delta, chi, gamma])
#%%
# Get strains for trial and test functions
strains_trial = get_strains(du, dtheta)
strains_test = get_strains(v_u, v_theta)

# Define constitutive law (diagonal matrix)
# [Normal Force, Bending Moment, Shear Force]
C = ufl.diag(ufl.as_vector([E * S, E * I, kappa * G * S]))

# Define bilinear form (internal virtual work)
# a = integral( (C * strain_trial) . strain_test ) * dx

dx_shear = ufl.dx(metadata={"quadrature_degree": 0}) # reduced integration
a = ( E * S * strains_trial[0] * strains_test[0] * ufl.dx + # Axial part
      E * I * strains_trial[1] * strains_test[1] * ufl.dx + # Bending part
      kappa * G * S * strains_trial[2] * strains_test[2] * dx_shear ) # Shear part (reduced integration)
#%%
# Define linear form (external virtual work)
q = fem.Constant(domain, dolfinx.default_scalar_type(q_val))
l_form = q * v_u[1] * ufl.dx #  load in the y-direction

# --- 5. Boundary Conditions ---
# We model a cantilever beam, so it's clamped (fixed) at x=0.
# This means displacement (u_x, u_y) and rotation (theta) are all zero.
def clamped_boundary(x):
    return np.isclose(x[0], 0)

# Locate DOFs for displacement and rotation at the clamped end
clamped_dofs_u = fem.locate_dofs_geometrical((W.sub(0), V_u), clamped_boundary)
clamped_dofs_theta = fem.locate_dofs_geometrical((W.sub(1), V_theta), clamped_boundary)

# Create zero-value functions for the Dirichlet BCs
u0 = fem.Function(V_u, name="u_clamped")
u0.x.array[:] = 0.0
theta0 = fem.Function(V_theta, name="theta_clamped")
theta0.x.array[:] = 0.0
#%%
# Apply the boundary conditions
bc_u = fem.dirichletbc(u0, clamped_dofs_u, W.sub(0))
bc_theta = fem.dirichletbc(theta0, clamped_dofs_theta, W.sub(1))
bcs = [bc_u, bc_theta]

# --- 6. Solve the Linear System ---
# w_sol will hold the solution for both displacement and rotation
w_sol = fem.Function(W, name="solution")
#%%
# Set up and solve the linear problem
from dolfinx.fem.petsc import LinearProblem
problem = LinearProblem(a, l_form, bcs=bcs, u=w_sol,
                                  petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
problem.solve()
#%%
# --- 7. Post-processing  ---
# Extract displacement and rotation from the solution
u_sol = w_sol.sub(0).collapse()
theta_sol = w_sol.sub(1).collapse()

# Analytical solution for Timoshenko beam tip deflection under uniform load `q`
w_analytical_T = (q_val * L**4) / (8 * E * I) + (q_val * L**2) / (2 * kappa * G * S)
print(f"Timoshenko analytical tip deflection: {w_analytical_T:.6e}")
'''
For the Euler-Bernoulli beam, you are expected to find the analytical solution yourself

w_analytical_E = ...
'''
print(f"Euler Bernoulli analytical tip deflection:")
# Get the computed deflection at the tip (x=L)
print("\n--- Solution Results ---")
print(f"Maximum vertical displacement: {np.max(np.abs(u_sol.x.array[1::2])):.6f} m")
print(f"Tip displacement (x-direction): {u_sol.x.array[-2]:.6e} m")
print(f"Tip displacement (y-direction): {u_sol.x.array[-1]:.6e} m")
#%%

# --- 8. Visualization with PyVista ---
from dolfinx import plot
# The 1D mesh needs to be converted to a 3D object for pyvista plotting
gdim = domain.geometry.dim
topology, cell_types, geometry = plot.vtk_mesh(domain, gdim)
grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
#%%
# Add the computed displacement vector to the grid
# We need to create a 3D vector from our 2D displacement solution for pyvista
u_plot = np.zeros((geometry.shape[0], 3))
value_dim = V_u.dofmap.bs # Block size of the vector function space, which is 2 since u_sol has 2 components
u_plot[:, :value_dim] = u_sol.x.array.reshape(-1, value_dim)
grid["Displacement"] = u_plot

# Warp the mesh by the displacement vector to show the deformed shape
factor = 20.  # Factor to exaggerate the deformation for visualization
warped = grid.warp_by_vector("Displacement", factor=factor) 
# Create plotter
plotter = pyvista.Plotter()
plotter.add_text("Deformed Beam Structure (Deformation exaggerated)", font_size=15)
plotter.add_mesh(grid, style='wireframe', color='gray', label='Undeformed')
plotter.add_mesh(warped, show_edges=True, line_width=5, color='red', label='Deformed')
plotter.add_legend()
plotter.view_xy() # View in the x-y plane
plotter.show()

# %%
