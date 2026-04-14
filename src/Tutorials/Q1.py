#%%
import matplotlib as mpl
import matplotlib.pyplot as plt
import pyvista
import ufl
import numpy as np

from petsc4py import PETSc
from mpi4py import MPI
#%%
from dolfinx import fem, mesh, io, plot,geometry
from dolfinx.fem.petsc import assemble_vector, assemble_matrix, create_vector, apply_lifting, set_bc
#%%
# Define temporal parameters
t = 0  # Start time
T = 5000.0  # Final time
num_steps = 100
dt = T / num_steps  # time step size

# Define mesh
nx, ny = 50, 50
domain = mesh.create_rectangle(MPI.COMM_WORLD,[np.array([0, 0]), np.array([2, 2])],
                               [nx, ny], mesh.CellType.quadrilateral)
#%%
# #convert to trapezoid 

# points = domain.geometry.x.reshape(-1, 3)
# x = points[:, 0]
# y = points[:, 1]
# points[:, 0] = x +(y / 4.0) * (1-x)
# tdim = domain.topology.dim
# domain.topology.create_connectivity(tdim, 0)
# cell_to_vertex_map = domain.topology.connectivity(tdim, 0)
#%%
# Define physical parameters for a material like steel
rho = fem.Constant(domain, PETSc.ScalarType(2700))  # Density (kg/m^3)
c = fem.Constant(domain, PETSc.ScalarType(500))   # Specific heat capacity (J/(kg*K))
k = fem.Constant(domain, PETSc.ScalarType(75))    # Thermal conductivity (W/(m*K))

V = fem.functionspace(domain, ("P", 1))
pyvista.set_jupyter_backend('client')
topology, cells, geom = plot.vtk_mesh(V)
grid = pyvista.UnstructuredGrid(topology, cells, geom)
plotter = pyvista.Plotter(window_size=(600, 600))
renderer = plotter.add_mesh(grid, show_edges=True)
plotter.view_xy()
plotter.show()

#%%
# Create initial condition
def initial_condition_smooth(x, a=10,b=5):
    """Initial condition function."""
    return 500*np.exp(-a * (x[0]-1)**2 -b*(x[1]-1)**2)

def initial_condition_sharp(x):
    '''
    top hat initial condition: return 500 if 0.5<x[0]<1.5 and 0.5<x[1]<1.5, else return 0
    '''
    return np.where((0.75 < x[0]) & (x[0] < 1.25) & (0.5 < x[1]) & (x[1] < 1.5), 500, 0)

initial_condition = initial_condition_sharp
u_n = fem.Function(V)
u_n.name = "u_n"
u_n.interpolate(initial_condition)

# Create boundary condition
fdim = domain.topology.dim - 1
boundary_facets = mesh.locate_entities_boundary(
    domain, fdim, lambda x: np.full(x.shape[1], True, dtype=bool))
bc = fem.dirichletbc(PETSc.ScalarType(0), fem.locate_dofs_topological(V, fdim, boundary_facets), V)
#%%
xdmf = io.XDMFFile(domain.comm, "diffusion.xdmf", "w")
xdmf.write_mesh(domain)

# Define solution variable, and interpolate initial solution for visualization in Paraview
uh = fem.Function(V)
uh.name = "uh"
uh.interpolate(initial_condition)
xdmf.write_function(uh, t)
#%%
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
f = fem.Constant(domain, PETSc.ScalarType(0))
# Define the bilinear and linear forms
# Equation: rho*c*u*v*dx + dt*k*nabla_u*nabla_v*dx = (rho*c*u_n + dt*f)*v*dx
a = rho * c * u * v * ufl.dx + dt * k * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = (rho * c * u_n + dt * f) * v * ufl.dx
#%%
bilinear_form = fem.form(a)
linear_form = fem.form(L)
#%%
A = assemble_matrix(bilinear_form, bcs=[bc])
A.assemble()
#to print the matrix use:
# print(A.convert("dense").getDenseArray())
b = create_vector(linear_form)
#%%
solver = PETSc.KSP().create(domain.comm)
solver.setOperators(A)
solver.setType(PETSc.KSP.Type.PREONLY)
solver.getPC().setType(PETSc.PC.Type.LU)
#%%
pyvista.set_jupyter_backend('trame')

grid = pyvista.UnstructuredGrid(*plot.vtk_mesh(V))

plotter = pyvista.Plotter()
plotter.open_gif("u_time.gif", fps=10)

grid.point_data["uh"] = uh.x.array
warped = grid.warp_by_scalar("uh", factor=0.01)

viridis = mpl.colormaps.get_cmap("viridis").resampled(25)
sargs = dict(title_font_size=25, label_font_size=20, fmt="%.2e", color="black",
             position_x=0.1, position_y=0.8, width=0.8, height=0.1)

renderer = plotter.add_mesh(warped, show_edges=True, lighting=False,
                            cmap=viridis, scalar_bar_args=sargs,
                            clim=[0, max(uh.x.array)])
#%%
# Setup for plotting temperature profile at x=1

probe_points = np.zeros((ny + 1, 3))
probe_points[:, 0] = 1
probe_points[:, 1] = np.linspace(0, 2, ny + 1)

# Create a bounding box tree for efficient point location
bb_tree = geometry.bb_tree(domain, domain.topology.dim)

# Find cells containing the probe points
cells = []
points_on_proc = []

for i, point in enumerate(probe_points):
    cell_candidates = geometry.compute_collisions_points(bb_tree, point.reshape(1, -1))
    if len(cell_candidates.links(0)) > 0:
        # Take the first colliding cell
        cell = cell_candidates.links(0)[0]
        cells.append([cell])
        points_on_proc.append(point)
    else:
        # Point not found on this process 
        pass

# Convert to numpy arrays for evaluation
if len(points_on_proc) > 0:
    points_on_proc = np.array(points_on_proc)
    cells = np.array([c[0] for c in cells], dtype=np.int32)
else:
    points_on_proc = np.array([]).reshape(0, 3)
    cells = np.array([], dtype=np.int32)

print(f"Found {len(points_on_proc)} probe points on this process")
y_coords_plot = points_on_proc[:, 1]
#%%
Ts=[]
for i in range(num_steps):
    t += dt

    # Update the right hand side reusing the initial vector
    with b.localForm() as loc_b:
        loc_b.set(0)
    assemble_vector(b, linear_form)

    # Apply Dirichlet boundary condition to the vector
    apply_lifting(b, [bilinear_form], [[bc]])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    set_bc(b, [bc])

    # Solve linear problem
    solver.solve(b, uh.x.petsc_vec)
    uh.x.scatter_forward()

    # Update solution at previous time step (u_n)
    u_n.x.array[:] = uh.x.array

    # Write solution to file
    xdmf.write_function(uh, t)
    # Update plot
    new_warped = grid.warp_by_scalar("uh", factor=0.01)
    warped.points[:, :] = new_warped.points
    warped.point_data["uh"][:] = uh.x.array
    plotter.write_frame()

    # Evaluate and plot temperature profile every 10 steps
    if (i + 1) % 10 == 0:
        print(f"Plotting profile at t = {t:.2f}")
        
        if len(points_on_proc) > 0:
            # Evaluate solution at the valid probe points
            temp_at_x1 = uh.eval(points_on_proc, cells)
            
            Ts.append(temp_at_x1)

plotter.close()
xdmf.close()

# Plot the temperature profile at x=1
Ts = np.array(Ts)  # Shape: (num_timesteps, 51)
plt.figure(figsize=(10, 6))

# Plot each time step as a separate line
for i, temp_profile in enumerate(Ts):
    plt.plot(y_coords_plot, temp_profile, marker='o', linestyle='-', 
             markersize=4, label=f'Time step {i}')

plt.title("Temperature Profile at x=0 over Time")
plt.xlabel("y-coordinate")
plt.ylabel("Temperature (K)")
plt.grid()
plt.legend()
plt.show()
#%%
