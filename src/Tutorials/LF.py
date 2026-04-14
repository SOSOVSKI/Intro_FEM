#%%
from dolfinx import mesh
import ufl
import numpy as np
from mpi4py import MPI 
# Define the geometry of the beam
Length = 20.0  # Length of the beam
Height = 1  # Height of the beam
Nx = 20  # Number of elements in the x-direction
Ny = 4  # Number of elements in the y-direction
#create the rectangular mesh
domain = mesh.create_rectangle(MPI.COMM_WORLD,
                                     [np.array([0, 0]), np.array([Length, Height])],
                                     [Nx, Ny],
                                     cell_type=mesh.CellType.quadrilateral)
gdim = domain.geometry.dim # dimension 
#%%
from dolfinx import fem

# Define function space
Vdegree = 1  # Degree of the polynomial
V = fem.functionspace(domain, ("P", Vdegree, (gdim,)))
# V is a vector function space with gdim=2
# Vdegree=1 means linear elements (P1)

#Define the test and trial functions
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
#%%
# Define Dirichlet boundary condition (fixed left edge)
def left_boundary(x):
    return np.isclose(x[0], 0)
#locate the dofs on the left boundary
fixed_side = fem.locate_dofs_geometrical(V, left_boundary)
# Define Dirichlet boundary condition
bc = [fem.dirichletbc(np.zeros(gdim), fixed_side, V)]

#%%
from dolfinx import default_scalar_type
# Define body force (gravity)
rho = fem.Constant(domain, 2700.) # density of the material 
g = fem.Constant(domain, 9.81) # gravitational acceleration
# Define body force vector
f = fem.Constant(domain, default_scalar_type((0, -rho * g)))
#%%
def bottom_boundary(x):
    return np.isclose(x[1], 0.0)
# Mark the boundaries with unique tags
fdim = gdim - 1  # The boundary dimension for 2D problems
bottom_facet = mesh.locate_entities_boundary(domain, fdim, bottom_boundary)
# Create a MeshTag: we tag the bottom boundary with a unique tag
bottom_tag = np.full_like(bottom_facet, 1)
facet_tag = mesh.meshtags(domain, fdim, bottom_facet, bottom_tag)
# Define a function for the tracitons
Traction_Function = fem.Function(V,name="Traction_Function")
# Define the traction function
tractions = lambda x: np.vstack((np.zeros_like(x[0], dtype=default_scalar_type), -0. * x[0]))
# Assign the traction function to the Traction_Function
Traction_Function.interpolate(tractions)
# Define facet normal vector
n_vector = ufl.FacetNormal(domain)
#%%
ds = ufl.Measure("ds", domain=domain, subdomain_data=facet_tag)
# Define material properties
E = fem.Constant(domain, 70e9)  # Young's modulus in MPa
nu = fem.Constant(domain, 0.3)  # Poisson's ratio
# Define Lame's parameters
lame = E * nu / ((1 + nu) * (1 - 2 * nu))  # Lame's first parameter
mu = E / (2 * (1 + nu))  # Lame's second parameter - shear modulus

def strain(v):
    """Define the strain tensor"""
    return ufl.sym(ufl.grad(v)) # the symmetric gradient

def stress(u):
     return lame * ufl.tr(strain(u))*ufl.Identity(gdim) + 2 * mu * strain(u)  # the stress tensor
#%%
# Define the bilinear form a(u,v)
a = ufl.inner(stress(u), strain(v)) * ufl.dx # ufl.dx is the integration measure over the domain
# Define the linear form L(v)
# L = ufl.dot(f, v) * ufl.dx 
# if you wish to apply the traction boundary condition you can add it to the linear form as follows:
L_with_tractions = L = ufl.dot(f, v) * ufl.dx  + ufl.dot(Traction_Function,v) * ds(1)  # dS(1) is the integration measure over the bottom boundary

#%%
from dolfinx.fem.petsc import LinearProblem
problem = LinearProblem(a, L, bcs=bc, petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
u_solution = problem.solve()
#%%
import pyvista as pv
#set off-screen rendering for Jupyter notebooks
pv.OFF_SCREEN = True
from dolfinx import plot
pv.set_jupyter_backend('trame')
# Create a PyVista plotter
plotter = pv.Plotter()
# Create a PyVista mesh from the FEniCSx mesh
topology, cell_types, geometry = plot.vtk_mesh(V)
grid = pv.UnstructuredGrid(topology, cell_types, geometry)
#Since Pyvista expects a 3D vector field, we need to reshape the solution
u_2d = u_solution.x.array.reshape(geometry.shape[0], gdim)
u_3d = np.zeros((geometry.shape[0], 3))
u_3d[:, :gdim] = u_2d  # Copy the x and y components
# Add the displacements ato the grid
grid["u"] = u_3d
# plot the undeformed mesh
undeformed = plotter.add_mesh(grid, style="wireframe", color="k")
# plot the deformed mesh
#create a warped mesh, and scale the displacements by a factor of your choice
factor = 10.
warped = grid.warp_by_vector("u",factor=factor)
#add the deformed mesh to the plotter
deformed = plotter.add_mesh(warped, show_edges=True)
plotter.show_axes()
#%%
#plot the two meshes. 
if not pv.OFF_SCREEN:
    plotter.show()
else:
    disp_figure = plotter.screenshot("displacements.png")
#%%
# Compute the stress from the solution
stress_solution = stress(u_solution)
#transfer to deviatoric stress
deviatoric_stress = stress_solution - ufl.tr(stress_solution) / gdim * ufl.Identity(gdim)
# Compute the Von Mises stress
von_mises_stress = ufl.sqrt(3.0 / 2.0 * ufl.inner(deviatoric_stress, deviatoric_stress))
# Create a FunctionSpace for the stress 
stress_space = fem.functionspace(domain, ("DG", 0))
# Create an expression for the Von Mises stress
Mises_expression = fem.Expression(von_mises_stress, stress_space.element.interpolation_points())
VM_stress = fem.Function(stress_space)
VM_stress.interpolate(Mises_expression)
#%%
warped.cell_data["Von_Mises_Stress"] = VM_stress.x.array
warped.set_active_scalars("Von_Mises_Stress")
plotter = pv.Plotter()
plotter.add_mesh(warped)
plotter.show_axes()
if not pv.OFF_SCREEN:
    plotter.show()
else:
    VM_Stress_figure = plotter.screenshot(f"Von_Mises_Stress.png")
# %%
