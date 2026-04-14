
# Example of implementing a simple 1D FEM solver for a specific problem
def solve_specific_1d_problem(n_elements=10):
    """
    Solve a specific 1D problem using FEM
    Example: -d^2y/dx^2 = 1 with y(0) = y(1) = 0
    """
    Hprint("\n" + "=" * 70)
    Hprint("EXAMPLE: SOLVING A SPECIFIC 1D PROBLEM")
    Hprint("-d²y/dx² = 1 with y(0) = y(1) = 0")
    Hprint("=" * 70)
    
    Hprint("Step 11: Set up a concrete problem and discretize the domain:")
    Hprint("EXPLANATION: We'll solve the Poisson equation -d²y/dx² = 1 with homogeneous")
    Hprint("Dirichlet boundary conditions y(0) = y(1) = 0. This has the exact solution")
    Hprint("y(x) = 0.5x(1-x), which represents the deflection of a simply supported beam")
    Hprint("under uniform load. We divide the domain [0,1] into n_elements equal elements.")
    
    # Number of nodes (one more than elements)
    n_nodes = n_elements + 1
    
    # Element length
    h = 1.0 / n_elements
    
    Hprint(f"\nStep 12: Initialize the global system (n_elements={n_elements}):")
    Hprint("EXPLANATION: We create a global stiffness matrix K and load vector F. For our")
    Hprint("1D problem with piecewise linear elements, K will be a tridiagonal matrix of")
    Hprint(f"size {n_nodes}×{n_nodes}, and F will be a vector of size {n_nodes}.")
    
    # Initialize stiffness matrix and load vector
    K = np.zeros((n_nodes, n_nodes))
    F = np.zeros(n_nodes)
    
    Hprint("\nStep 13: Assemble the system element by element:")
    Hprint("EXPLANATION: We loop through each element and compute its contribution to the")
    Hprint("global system. For each element, we compute a 2×2 element stiffness matrix and")
    Hprint("a 2×1 element load vector. These are then assembled into the global K and F at")
    Hprint("the appropriate locations based on the element's global node numbers.")
    
    # Assemble the system
    for e in range(n_elements):
        # Element indices
        i = e
        j = e + 1
        
        # Element stiffness matrix (for -d²y/dx²)
        # For linear elements, this is [[1, -1], [-1, 1]]/h
        # Derived from integrating phi_i' * phi_j' over the element
        k_e = np.array([[1, -1], [-1, 1]]) / h
        
        # Element load vector (for right side = 1)
        # For linear elements and constant source term = 1
        # This is [h/2, h/2], derived from integrating phi_i over the element
        f_e = np.array([h/2, h/2])
        
        # Add to global system (assembly process)
        K[i:i+2, i:i+2] += k_e
        F[i:i+2] += f_e
    
    Hprint("\nStep 14: Apply boundary conditions:")
    Hprint("EXPLANATION: We need to enforce the Dirichlet boundary conditions y(0) = y(1) = 0.")
    Hprint("A straightforward approach is to simply remove the corresponding rows and columns")
    Hprint("from the system. Since the first and last nodes are on the boundary, we remove the")
    Hprint("first and last rows/columns from K and the first and last entries from F.")
    
    # Apply boundary conditions (y(0) = y(1) = 0)
    # Remove first and last rows/columns from system
    K_reduced = K[1:-1, 1:-1]
    F_reduced = F[1:-1]
    
    Hprint("\nStep 15: Solve the linear system:")
    Hprint("EXPLANATION: We now solve the reduced linear system K_reduced·a = F_reduced for")
    Hprint("the unknown coefficients a. These coefficients represent the solution values at")
    Hprint(f"the interior nodes. For our problem, this is a system of {n_nodes-2} equations.")
    
    # Solve the system
    a = np.linalg.solve(K_reduced, F_reduced)
    
    Hprint("\nStep 16: Reconstruct the full solution:")
    Hprint("EXPLANATION: We create a full solution vector including the boundary nodes.")
    Hprint("The values at the boundary nodes are set to 0 (per our boundary conditions),")
    Hprint("and the values at the interior nodes are set to the computed coefficients.")
    
    # Full solution vector (including boundary nodes)
    y = np.zeros(n_nodes)
    y[1:-1] = a
    
    Hprint("\nStep 17: Compare with the exact solution:")
    Hprint("EXPLANATION: For this problem, we know the exact solution is y(x) = 0.5x(1-x).")
    Hprint("We compare our FEM approximation with this exact solution to assess the accuracy.")
    Hprint("We expect the FEM solution to converge to the exact solution as we increase the")
    Hprint("number of elements.")
    
    # Exact solution for comparison
    x_exact = np.linspace(0, 1, 100)
    y_exact = 0.5 * x_exact * (1 - x_exact)
    
    # FEM solution for plotting
    x_fem = np.linspace(0, 1, n_nodes)
    
    # Plot results
    plt.figure(figsize=(10, 6))
    plt.plot(x_exact, y_exact, 'b-', label='Exact solution: y = 0.5x(1-x)')
    plt.plot(x_fem, y, 'ro-', label=f'FEM solution ({n_elements} elements)')
    plt.title('1D FEM Solution: -d²y/dx² = 1, y(0) = y(1) = 0')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.legend()
    plt.show()
    
    return x_fem, y, x_exact, y_exact
#%%

def solve_specific_1d_problem_animated(n_elements=10):
    """
    Solves -d^2y/dx^2 = 1 with y(0)=y(1)=0 using n_elements.
    Returns the matplotlib figure object.
    """
    # --- Step 11: Setup (Calculations only) ---
    n_nodes = n_elements + 1
    h = 1.0 / n_elements

    # --- Step 12: Initialize ---
    K = np.zeros((n_nodes, n_nodes))
    F = np.zeros(n_nodes)

    # --- Step 13: Assemble ---
    for e in range(n_elements):
        i = e
        j = e + 1
        k_e = np.array([[1, -1], [-1, 1]]) / h
        f_e = np.array([h/2, h/2])
        K[i:i+2, i:i+2] += k_e
        F[i:i+2] += f_e

    # --- Step 14: Apply BCs ---
    K_reduced = K[1:-1, 1:-1]
    F_reduced = F[1:-1]

    # --- Step 15: Solve ---
    
    if K_reduced.shape[0] > 0:
         try:
            a = np.linalg.solve(K_reduced, F_reduced)
         except np.linalg.LinAlgError:
            print(f"Warning: Singular matrix for n_elements={n_elements}. Cannot solve.")
            a = np.zeros(K_reduced.shape[0]) 
    else:
        a = np.array([]) # No interior nodes to solve for


    # --- Step 16: Reconstruct Solution ---
    y = np.zeros(n_nodes)
    if len(a) > 0:
      y[1:-1] = a

    
    x_exact = np.linspace(0, 1, 100)
    y_exact = 0.5 * x_exact * (1 - x_exact)
    x_fem = np.linspace(0, 1, n_nodes)

    # Create the figure object
    fig = plt.figure(figsize=(8, 5)) 
    ax = fig.add_subplot(1, 1, 1) 

    ax.plot(x_exact, y_exact, 'b-', label='Exact solution: y = 0.5x(1-x)')
    ax.plot(x_fem, y, 'ro-', label=f'FEM solution ({n_elements} elements)')
    ax.set_title(f'1D FEM Solution (n_elements = {n_elements})') 
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.grid(True)
    ax.legend()

    plt.close(fig) 

    return fig 

#%%
from ipywidgets import interact, IntSlider
@interact(n_elements=IntSlider(min=1, max=50, step=1, value=2, description='# Elements:'))
def interactive_solver_plot(n_elements):

    fig = solve_specific_1d_problem_animated(n_elements)

    return fig

# %%
