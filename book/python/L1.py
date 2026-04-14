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
from IPython.display import Math, display, HTML, Markdown
import re as _re
#%%
# Define a global function to pretty print math equations using LaTeX
def Lprint(x):
    display(Math(latex(x)))

def Hprint(x):
    """Display text as Markdown (works in both HTML and PDF output)."""
    text = str(x)
    # Convert common HTML tags to Markdown equivalents
    text = _re.sub(r'<b>(.*?)</b>', r'**\1**', text, flags=_re.DOTALL)
    text = text.replace('<br>', '\n').replace('<br />', '\n').replace('<br/>', '\n')
    display(Markdown(text))

def Lspace():
    pass  # spacing handled naturally by Markdown paragraph breaks
#%%


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

def derive_weak_form():
    """Derive the weak form of the equation step by step with explanations"""
    Hprint("<b>Step 1: Multiply the strong form by a test function v(x):</b>")
    Lspace()
    Hprint("EXPLANATION: We start with the strong form of our differential equation and multiply both sides by a test function v(x)")
    Lspace()
    Hprint("This is the first step in converting to the weak form.")
    Lspace()
    Hprint(" The test function v(x) will eventually vanish at the boundaries, helping us incorporate boundary conditions naturally.")
    Lspace()
    eq1 = sp.Eq(v * differential_operator(y), v * F)
    Lprint(eq1)

    Lspace()
    Hprint("<b>Step 2: Integrate over the domain [0, L]:</b>")
    Lspace()
    Hprint("EXPLANATION:<br> We integrate both sides of the equation over the entire domain [0, L].")
    Hprint("<br> This transforms our pointwise equation into an integral equation that will hold over the whole domain.")
    Hprint(" This step moves us from a local formulation to a global one, which is a key aspect of the finite element method.")
    Lspace()
    lhs = sp.Integral(v * differential_operator(y), (x, 0, L))
    rhs = sp.Integral(v * F, (x, 0, L))
    eq2 = sp.Eq(lhs, rhs)
    Lprint(eq2)

    Hprint("<b>Step 3: Expand the left-hand side:</b>")
    Lspace()
    Hprint("EXPLANATION: ")
    Lspace()
    Hprint("We expand the left-hand side to separate the terms with different derivatives of y. <br> This makes it easier to apply integration by parts to the second-derivative term.")
    Lspace()
    Hprint("Breaking down the differential operator allows us to handle each term appropriately based on the order of derivatives.")
    Lspace()    
    lhs_expanded = sp.Integral(v * A * diff(y, x, 2), (x, 0, L)) + \
                  sp.Integral(v * B * diff(y, x), (x, 0, L)) + \
                  sp.Integral(v * C * y, (x, 0, L))
    eq3 = sp.Eq(lhs_expanded, rhs)
    Lprint(eq3)
 
    Hprint("<b>Step 4: Integration by parts for the term with second derivative:</b>")
    Lspace()
    Hprint("EXPLANATION:<br> ")
    Hprint("This is the critical step that reduces the order of derivatives in our equation. <br> We apply integration by parts to the term containing the second derivative of y.")
    Hprint(" The formula for integration by parts is:")
    Lspace()
    Hprint("  ∫u·dv = [u·v] - ∫v·du")
    Hprint("For our case:<br>")
    Hprint("  u = v(x)·A(x) and dv = d²y/dx² dx")
    Lspace()
    Hprint("  du = d(v(x)·A(x))/dx dx and v = dy/dx")
    Lspace()
    Hprint("Using the product rule:")
    Lspace()
    Hprint("  d/dx(v·A·dy/dx) = v·A·d²y/dx² + (v·A)'·dy/dx")
    Hprint("Rearranging:")
    Lspace()
    Hprint("  v·A·d²y/dx² = d/dx(v·A·dy/dx) - (v·A)'·dy/dx")
    Lspace()
    Hprint("This substitution allows us to replace the second derivative with first derivatives, reducing the continuity requirements on our approximation functions.")
    Lspace()
    
    # # Apply integration by parts
    # ibp_term = sp.Integral(sp.diff(v * A * diff(y, x), x), (x, 0, L)) - \
    #           sp.Integral(sp.diff(v * A, x) * diff(y, x), (x, 0, L))
    # Lprint(ibp_term)
    # Evaluate the boundary term from integration by parts
    boundary_term = (v * A * diff(y, x)).subs(x, L) - (v * A * diff(y, x)).subs(x, 0)
    
    Hprint("After integration by parts:")
    Lspace()
    Lprint(sp.Eq(sp.Integral(v * A * diff(y, x, 2), (x, 0, L)), 
                 boundary_term - sp.Integral(sp.diff(v * A, x) * diff(y, x), (x, 0, L))))

    Hprint("<b>Step 5: Apply boundary conditions v(0) = v(L) = 0:</b>")
    Lspace()
    Hprint("EXPLANATION: ")
    Lspace()
    Hprint("We choose our test functions v(x) to vanish at the domain boundaries <br>(v(0) = v(L) = 0). ")
    Hprint("This is a crucial step that eliminates the boundary terms from integration by parts. <br> By carefully selecting test functions that satisfy these conditions, we naturally incorporate the essential (Dirichlet) boundary conditions into our formulation.")
    Hprint(" This is one of the elegant aspects of the finite element method.")
    Lspace()
    Hprint("The boundary term [v·A·dy/dx]₀ᴸ vanishes because v(0) = v(L) = 0.")

    # Final weak form
    Hprint("<b>Step 6: Assemble the final weak form:</b>")
    Lspace()
    Hprint("EXPLANATION: <br> Now we collect all terms to form the complete weak formulation.")
    Hprint("This form requires lower continuity requirements on our solution (only first derivatives appear).")
    Lspace()
    Hprint(" The weak form is equivalent to the strong form for sufficiently smooth solutions, but it allows us to find approximate solutions with less smoothness.")
    Lspace()
    Hprint("This is the foundation for the finite element approximation where we'll represent the solution as a combination of piecewise polynomial functions.")
    Lspace()
    weak_lhs = -sp.Integral(sp.diff(v * A, x) * diff(y, x), (x, 0, L)) + \
               sp.Integral(v * B * diff(y, x), (x, 0, L)) + \
               sp.Integral(v * C * y, (x, 0, L))
    
    weak_eq = sp.Eq(weak_lhs, rhs)
    Lprint(weak_eq)
    
    return weak_eq

# Define shape and test functions for FEM discretization
def define_shape_functions(N):
    """Define shape functions for discretization with N elements"""
    Hprint(f"<b>Step 7: Define shape functions for N elements:</b><br>")
    Hprint("EXPLANATION: <br>")
    Hprint("Now we discretize the problem by representing both the solution y(x) and test function v(x) as linear combinations of shape functions. ")
    Lspace()
    Hprint("This transforms our continuous problem into a discrete one with a finite number of unknowns (the coefficients).")
    Lspace()
    Hprint("The shape functions are typically chosen to be simple piecewise polynomials (often linear) that are non-zero only over a small portion of the domain .")
    Lspace()
    Hprint("This 'local support' property often leads to sparse matrices in the final system.")
    Lspace()
    # Symbol for shape functions
    phi = [Function(f'phi_{i}')(x) for i in range(1, N+1)]
    Lprint(phi)
    
    Hprint("<b>Step 8: Express solution and test functions using shape functions:</b>")
    Lspace()
    Hprint("EXPLANATION:<br> ")
    Hprint("We approximate both the solution y(x) and test function v(x) using the same set of shape functions. ")
    Lspace()
    Hprint("The solution is written as a sum of shape functions with unknown coefficients a_j. <br> Similarly, the test function is written as a sum with coefficients b_i.")
    Lspace()
    Hprint("This is the Galerkin approach, where we use the same functions for both approximation and testing.")
    Lspace()
    # Assume solution as linear combination of shape functions
    a = [Symbol(f'a_{j}') for j in range(1, N+1)]
    y_approx = sum(a[j-1] * phi[j-1] for j in range(1, N+1))
    Lprint(sp.Eq(y, y_approx))
    
    # Test function as linear combination of shape functions
    b = [Symbol(f'b_{i}') for i in range(1, N+1)]
    v_approx = sum(b[i-1] * phi[i-1] for i in range(1, N+1))
    Lprint(sp.Eq(v, v_approx))
    
    Hprint("\nWhere:<br>")
    Hprint("- aⱼ are unknown constants we need to solve for ")
    Lspace()
    Hprint("- bᵢ are arbitrary constants for the test function")
    Lspace()
    Hprint("- φᵢ(x), φⱼ(x) are known shape functions (typically piecewise linear functions)")
    Lspace()
    Hprint("NOTE: In the finite element method, we often choose shape functions that are equal to 1 at their corresponding node and 0 at all other nodes.<br>")
    Hprint(" This makes it easy to enforce boundary conditions and interpret the coefficients aⱼ as the solution values at the nodes.")
    
    return phi, a, b, y_approx, v_approx

# Assemble the FEM system
def assemble_fem_system(weak_form, phi, a, b, y_approx, v_approx, N):
    """Assemble the FEM system using the weak form and shape functions"""
    Hprint("<b> Step 9: Substitute discretized functions into the weak form:</b>")
    Lspace()
    Hprint("EXPLANATION:<br> We now substitute our approximations for y(x) and v(x) into the weak form.")
    Lspace()
    Hprint("This transforms the continuous weak form into a discrete algebraic system of equations.")
    Lspace()
    Hprint("Since the test function coefficients bᵢ are arbitrary, we can collect terms and set up a system of equations for the unknown coefficients aⱼ.")
    Lspace()
    # This is a simplified assembly process - in practice, you'd loop through elements
    
    Hprint("<b>Step 10: Assemble the system matrix and right-hand side vector:</b>")
    Hprint("EXPLANATION: <br>")
    Hprint("When we substitute the discretized functions and collect terms, we obtain a linear system of equations K·a = F, where:")
    Lspace()
    Hprint("  - K is the stiffness matrix (or system matrix)")
    Lspace()
    Hprint("  - a is the vector of unknown coefficients")
    Lspace()
    Hprint("  - F is the load vector (right-hand side)")
    Lspace()
    Hprint("Each entry K[i,j] is computed by evaluating an integral that involves the shape functions phi_i and phi_j and their derivatives.")
    Lspace()
    Hprint(" Similarly, each entry F[i] comes from an integral involving phi_i and the force function F(x).")
    Lspace()
    Hprint("For a simplified case where A(x) = B(x) = C(x) = 1, the entries would be:")
    Lspace()
    display(Math(r"K_{ij} = \int_0^L \frac{d\phi_i}{dx} \frac{d\phi_j}{dx} dx + \int_0^L \phi_i(x) \frac{d\phi_j}{dx} dx + \int_0^L \phi_i(x) \phi_j(x) dx"))
    
    Lspace()
    display(Math(r"F_i = \int_0^L \phi_i(x) F(x) dx"))
    
    Lspace()
    
    
    # Example for constant coefficients
    A_val = 1
    B_val = 1 
    C_val = 1
    
    # Define the stiffness matrix K and load vector F
    K = sp.zeros(N, N)
    F_vec = sp.zeros(N, 1)
    Lspace()
    
    Hprint("NOTE: In practice, the assembly process is typically done by looping over elements,computing the element contributions, and adding them to the global system.")
    Hprint(" This element-by-element approach is more efficient and aligns with the local support property of shape functions.")
    
    return K, F_vec


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
    
if __name__ == "__main__":
    finite_element_method_1d_demo()
#%%