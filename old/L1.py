#%%
"""
Introduction to the Finite Element Method
1D Setup Implementation using SymPy
"""

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, Function, Derivative, Integral, diff, integrate, Matrix, Symbol,latex
from sympy.printing import pprint
from IPython.display import Math, display, HTML
#%%
# Define a global function to pretty print math equations using LaTeX
def Lprint(x):
    display(Math(latex(x)))

def Hprint(x):
    display(HTML(f"<span style='font-size: 30px;'>{x}</span>"))
#%%
"""
Introduction to the Finite Element Method
1D Setup Implementation using SymPy
"""


# Define symbolic variables and functions
x, L = symbols('x L', real=True, positive=True)
A = Function('A')(x)  # Coefficient function A(x)
B = Function('B')(x)  # Coefficient function B(x)
C = Function('C')(x)  # Coefficient function C(x)
F = Function('F')(x)  # Right hand side function F(x)

y = Function('y')(x)  # Solution function
v = Function('v')(x)  # Test function

# Define the differential operator (strong form)
def differential_operator(y):
    """Define the strong form of the 1D differential equation"""
    return A * diff(y, x, 2) + B * diff(y, x) + C * y

# Display the strong form of the equation
def display_strong_form():
    """Display the strong form of the equation"""
    Hprint("Strong form of the equation:")
    eq = sp.Eq(differential_operator(y), F)
    Lprint(eq)
    Hprint("\nWith boundary conditions:")
    Lprint(sp.Eq(y.subs(x, 0), 0))
    Lprint(sp.Eq(y.subs(x, L), 0))
    Hprint("where x ∈ [0, L]")

# Derivation of the weak form
def derive_weak_form():
    """Derive the weak form of the equation step by step with detailed explanations"""
    Hprint("Step 1: Multiply the strong form by a test function v(x):")
    Hprint("EXPLANATION: We start with the strong form of our differential equation and multiply")
    Hprint("both sides by a test function v(x). This is the first step in converting to the weak")
    Hprint("form. The test function v(x) will eventually vanish at the boundaries, helping us")
    Hprint("incorporate boundary conditions naturally.")
    eq1 = sp.Eq(v * differential_operator(y), v * F)
    Lprint(eq1)
    
    Hprint("\nStep 2: Integrate over the domain [0, L]:")
    Hprint("EXPLANATION: We integrate both sides of the equation over the entire domain [0, L].")
    Hprint("This transforms our pointwise equation into an integral equation that will hold over")
    Hprint("the whole domain. This step moves us from a local formulation to a global one, which")
    Hprint("is a key aspect of the finite element method.")
    lhs = sp.Integral(v * differential_operator(y), (x, 0, L))
    rhs = sp.Integral(v * F, (x, 0, L))
    eq2 = sp.Eq(lhs, rhs)
    Lprint(eq2)
    
    Hprint("\nStep 3: Expand the left-hand side:")
    Hprint("EXPLANATION: We expand the left-hand side to separate the terms with different")
    Hprint("derivatives of y. This makes it easier to apply integration by parts to the")
    Hprint("second-derivative term. Breaking down the differential operator allows us to")
    Hprint("handle each term appropriately based on the order of derivatives.")
    lhs_expanded = sp.Integral(v * A * diff(y, x, 2), (x, 0, L)) + \
                  sp.Integral(v * B * diff(y, x), (x, 0, L)) + \
                  sp.Integral(v * C * y, (x, 0, L))
    eq3 = sp.Eq(lhs_expanded, rhs)
    Lprint(eq3)
    
    # Integration by parts for the second derivative term
    Hprint("\nStep 4: Integration by parts for the term with second derivative:")
    Hprint("EXPLANATION: This is the critical step that reduces the order of derivatives in our")
    Hprint("equation. We apply integration by parts to the term containing the second derivative")
    Hprint("of y. The formula for integration by parts is:")
    Hprint("  ∫u·dv = [u·v] - ∫v·du")
    Hprint("For our case:")
    Hprint("  u = v(x)·A(x) and dv = d²y/dx² dx")
    Hprint("  du = d(v(x)·A(x))/dx dx and v = dy/dx")
    Hprint("Using the product rule:")
    Hprint("  d/dx(v·A·dy/dx) = v·A·d²y/dx² + (v·A)'·dy/dx")
    Hprint("Rearranging:")
    Hprint("  v·A·d²y/dx² = d/dx(v·A·dy/dx) - (v·A)'·dy/dx")
    Hprint("This substitution allows us to replace the second derivative with first derivatives,")
    Hprint("reducing the continuity requirements on our approximation functions.")
    
    # Apply integration by parts
    ibp_term = sp.Integral(sp.diff(v * A * diff(y, x), x), (x, 0, L)) - \
              sp.Integral(sp.diff(v * A, x) * diff(y, x), (x, 0, L))
    
    # Evaluate the boundary term from integration by parts
    boundary_term = (v * A * diff(y, x)).subs(x, L) - (v * A * diff(y, x)).subs(x, 0)
    
    print("After integration by parts:")
    Lprint(sp.Eq(sp.Integral(v * A * diff(y, x, 2), (x, 0, L)), 
                 boundary_term - sp.Integral(sp.diff(v * A, x) * diff(y, x), (x, 0, L))))
    
    # Apply boundary conditions: v(0) = v(L) = 0
    Hprint("\nStep 5: Apply boundary conditions v(0) = v(L) = 0:")
    Hprint("EXPLANATION: We choose our test functions v(x) to vanish at the domain boundaries")
    Hprint("(v(0) = v(L) = 0). This is a crucial step that eliminates the boundary terms from")
    Hprint("integration by parts. By carefully selecting test functions that satisfy these")
    Hprint("conditions, we naturally incorporate the essential (Dirichlet) boundary conditions")
    Hprint("into our formulation. This is one of the elegant aspects of the finite element method.")
    Hprint("The boundary term [v·A·dy/dx]₀ᴸ vanishes because v(0) = v(L) = 0.")
    
    # Final weak form
    Hprint("\nStep 6: Assemble the final weak form:")
    Hprint("EXPLANATION: Now we collect all terms to form the complete weak formulation. This")
    Hprint("form requires lower continuity requirements on our solution (only first derivatives")
    Hprint("appear). The weak form is equivalent to the strong form for sufficiently smooth")
    Hprint("solutions, but it allows us to find approximate solutions with less smoothness.")
    Hprint("This is the foundation for the finite element approximation where we'll represent")
    Hprint("the solution as a combination of piecewise polynomial functions.")
    weak_lhs = -sp.Integral(sp.diff(v * A, x) * diff(y, x), (x, 0, L)) + \
               sp.Integral(v * B * diff(y, x), (x, 0, L)) + \
               sp.Integral(v * C * y, (x, 0, L))
    
    weak_eq = sp.Eq(weak_lhs, rhs)
    Lprint(weak_eq)
    
    return weak_eq

# Define shape and test functions for FEM discretization
def define_shape_functions(N):
    """Define shape functions for discretization with N elements"""
    Hprint(f"Step 7: Define shape functions for {N} elements:")
    Hprint("EXPLANATION: Now we discretize the problem by representing both the solution y(x)")
    Hprint("and test function v(x) as linear combinations of shape functions. This transforms our")
    Hprint("continuous problem into a discrete one with a finite number of unknowns (the coefficients).")
    Hprint("The shape functions are typically chosen to be simple piecewise polynomials (often linear)")
    Hprint("that are non-zero only over a small portion of the domain .")
    Hprint("This 'local support' property leads to sparse matrices in the final system.")
    
    # Symbol for shape functions
    phi = [Function(f'phi_{i}')(x) for i in range(1, N+1)]
    
    Hprint("\nStep 8: Express solution and test functions using shape functions:")
    Hprint("EXPLANATION: We approximate both the solution y(x) and test function v(x) using the")
    Hprint("same set of shape functions. The solution is written as a sum of shape functions with")
    Hprint("unknown coefficients a_j. Similarly, the test function is written as a sum with coefficients b_i.")
    Hprint("This is the Galerkin approach, where we use the same functions for both approximation and testing.")
    
    # Assume solution as linear combination of shape functions
    a = [Symbol(f'a_{j}') for j in range(1, N+1)]
    y_approx = sum(a[j-1] * phi[j-1] for j in range(1, N+1))
    Lprint(sp.Eq(y, y_approx))
    
    # Test function as linear combination of shape functions
    b = [Symbol(f'b_{i}') for i in range(1, N+1)]
    v_approx = sum(b[i-1] * phi[i-1] for i in range(1, N+1))
    Lprint(sp.Eq(v, v_approx))
    
    Hprint("\nWhere:")
    Hprint("a_j are unknown constants we need to solve for")
    Hprint("b_i are arbitrary constants for the test function")
    Hprint("phi_i(x), phi_j(x) are known shape functions (typically piecewise linear functions)")
    Hprint("\nNOTE: In the finite element method, we typically choose shape functions that are")
    Hprint("equal to 1 at their corresponding node and 0 at all other nodes. This makes it easy")
    Hprint("to enforce boundary conditions and interpret the coefficients a_j as the solution values")
    Hprint("at the nodes.")
    
    return phi, a, b, y_approx, v_approx

# Assemble the FEM system
def assemble_fem_system(weak_form, phi, a, b, y_approx, v_approx, N):
    """Assemble the FEM system using the weak form and shape functions"""
    Hprint("Step 9: Substitute discretized functions into the weak form:")
    Hprint("EXPLANATION: We now substitute our approximations for y(x) and v(x) into the weak form.")
    Hprint("This transforms the continuous weak form into a discrete algebraic system of equations.")
    Hprint("Since the test function coefficients b_i are arbitrary, we can collect terms and set")
    Hprint("up a system of equations for the unknown coefficients a_j.")
    
    # Replace y and v with their approximations in the weak form
    # This is a simplified assembly process - in practice, you'd loop through elements
    
    Hprint("\nStep 10: Assemble the system matrix and right-hand side vector:")
    Hprint("EXPLANATION: When we substitute the discretized functions and collect terms, we")
    Hprint("obtain a linear system of equations K·a = F, where:")
    Hprint("  - K is the stiffness matrix (or system matrix)")
    Hprint("  - a is the vector of unknown coefficients")
    Hprint("  - F is the load vector (right-hand side)")
    Hprint("Each entry K[i,j] is computed by evaluating an integral that involves the shape functions")
    Hprint("phi_i and phi_j and their derivatives. Similarly, each entry F[i] comes from an integral")
    Hprint("involving phi_i and the force function F(x).")
    Hprint("\nFor a simplified case where A(x) = B(x) = C(x) = 1, the entries would be:")
    Hprint("  K[i,j] = ∫(phi_i'·phi_j')dx + ∫(phi_i·phi_j')dx + ∫(phi_i·phi_j)dx")
    Hprint("  F[i] = ∫(phi_i·F(x))dx")
    Hprint("Where phi_i' represents the derivative of phi_i with respect to x.")
    
    # Example for constant coefficients
    A_val = 1
    B_val = 1 
    C_val = 1
    
    # Define the stiffness matrix K and load vector F
    K = sp.zeros(N, N)
    F_vec = sp.zeros(N, 1)
    
    Hprint("\nStiffness matrix K and load vector F will have the form:")
    Lprint(sp.Matrix(K))
    Lprint(sp.Matrix(F_vec))
    
    Hprint("\nNOTE: In practice, the assembly process is typically done by looping over elements,")
    Hprint("computing the element contributions, and adding them to the global system. This")
    Hprint("element-by-element approach is more efficient and aligns with the local support")
    Hprint("property of shape functions.")
    
    return K, F_vec

# Main function to demonstrate the FEM process
def finite_element_method_1d_demo():
    """Demonstrate the 1D FEM process"""
    Hprint("=" * 70)
    Hprint("INTRODUCTION TO THE FINITE ELEMENT METHOD")
    Hprint("Lecture #1: Introduction + 1D Setup")
    Hprint("=" * 70)
    
    # Display the strong form
    display_strong_form()
    
    Hprint("\n" + "=" * 70)
    Hprint("DERIVING THE WEAK FORM")
    Hprint("=" * 70)
    weak_form = derive_weak_form()
    
    Hprint("\n" + "=" * 70)
    Hprint("DISCRETIZATION USING SHAPE FUNCTIONS")
    Hprint("=" * 70)
    N = 5  # Number of elements for example
    phi, a, b, y_approx, v_approx = define_shape_functions(N)
    
    Hprint("\n" + "=" * 70)
    Hprint("ASSEMBLING THE FEM SYSTEM")
    Hprint("=" * 70)
    K, F_vec = assemble_fem_system(weak_form, phi, a, b, y_approx, v_approx, N)
    
    Hprint("\n" + "=" * 70)
    Hprint("SOLVING THE FEM SYSTEM")
    Hprint("=" * 70)
    Hprint("The final step would be solving the linear system:")
    Lprint(sp.Eq(sp.Matrix(K) * sp.Matrix([a[i] for i in range(N)]), sp.Matrix(F_vec)))
    
    Hprint("\nAfter solving for the unknown coefficients a_j, we can reconstruct the solution y(x)")
    
# Run the demonstration
if __name__ == "__main__":
    finite_element_method_1d_demo()

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

# Uncomment to run the specific example
#%%
def show_sol(n):
    x_fem, y, x_exact, y_exact =solve_specific_1d_problem(20)
# %%
