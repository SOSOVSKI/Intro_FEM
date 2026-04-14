#%%
"""
Finite Element Method - Lecture 2
Deriving the Algebraic System & Galerkin's Method
Based on SymPy
"""

import sympy as sp
from sympy import symbols, Function, Derivative, Integral, diff, integrate, Matrix, Symbol, Idx, IndexedBase, Sum, Eq, Q, Wild
from sympy.printing import pprint, latex
from IPython.display import Math, display, HTML, Markdown

#%% --- Helper Functions ---

def print_styled(text, size="1.1em"):
    """Display text as Markdown (works in both HTML and PDF output)."""
    display(Markdown(str(text)))

# Function to pretty print math equations using LaTeX
def Lprint(expr):
    """Displays SymPy expressions using LaTeX in Jupyter/Quarto."""
    display(Math(latex(expr)))

#%% --- Symbolic Definitions ---

print_styled("<h2>Symbolic Definitions</h2>")

# Independent variable and domain length
x, L = symbols('x L', real=True, positive=True)

# Indices for sums and matrix notation
i = Symbol('i', integer=True, positive=True)
j = Symbol('j', integer=True, positive=True)
N = Symbol('N', integer=True, positive=True) # Number of basis functions / nodes (approx)

# Unknown solution function y(x) and test function v(x)
y = Function('y')(x)
v = Function('v')(x)

# Coefficients in the differential equation (ASSUMED CONSTANT for this lecture)
A, B, C = symbols('A B C', real=True) # No longer functions of x

# Right-hand side function F(x)
F = Function('F')(x)

# Shape functions (basis functions) phi_i(x) and phi_j(x)
# Using IndexedBase for phi_i, phi_j representation
phi = Function('phi') # Represents the set of shape functions phi_k(x)

# Unknown coefficients a_j and arbitrary coefficients b_i
a = IndexedBase('a') # Represents the vector of unknown coefficients a_j
b = IndexedBase('b') # Represents the vector of arbitrary coefficients b_i


#%% --- Recall Weak Form (from Lecture 1) ---
print_styled("<h2>1. Recall Weak Form (from Lecture 1)</h2>")
print_styled("Assuming constant A, B, C and Dirichlet BCs v(0)=v(L)=0:")

# Define the terms symbolically. Note: In lecture 1, (v*A)' was used.
# If A is constant, (v*A)' = A*dv/dx.
weak_term1 = -Integral(A * diff(v, x) * diff(y, x), (x, 0, L))
weak_term2 = Integral(v * B * diff(y, x), (x, 0, L))
weak_term3 = Integral(v * C * y, (x, 0, L))
rhs_term = Integral(v * F, (x, 0, L))

weak_form_eq = Eq(weak_term1 + weak_term2 + weak_term3, rhs_term)
print_styled("The weak form is: Find y(x) (such that y(0)=y(L)=0) for which the following holds for all valid test functions v(x) (with v(0)=v(L)=0):")
Lprint(weak_form_eq)

#%% --- Discretization using Shape Functions ---
print_styled("<h2>2. Substitute Discretized Approximations</h2>")
print_styled("We approximate the solution y(x) and test function v(x) using shape functions phi_k(x):")

# Define the approximations using Sum
y_approx = Sum(a[j] * phi(j,x), (j, 1, N))
v_approx = Sum(b[i] * phi(i,x), (i, 1, N)) # Using index i for v

print_styled("Approximate solution y(x):")
Lprint(Eq(y, y_approx, evaluate=False))
print_styled("Test function v(x):")
Lprint(Eq(v, v_approx, evaluate=False))

print_styled("Derivatives:")
dy_dx_approx = Derivative(y_approx, x) # Keep derivative symbolic for now
# Or symbolically differentiate term-by-term:
# dy_dx_approx = Sum(a[j] * diff(phi(j,x), x), (j, 1, N))
dv_dx_approx = Derivative(v_approx, x)
# dv_dx_approx = Sum(b[i] * diff(phi(i,x), x), (i, 1, N))

Lprint(Eq(Derivative(y,x), dy_dx_approx, evaluate=False))
Lprint(Eq(Derivative(v,x), dv_dx_approx, evaluate=False))

print_styled("Now, substitute these sums into the weak form.")
print_styled("Substituting into the first term: -Integral(A * dv/dx * dy/dx, (x, 0, L))")
# Replace v and y with their sum representations. Be careful with indices.
# SymPy doesn't easily simplify nested sums within integrals symbolically in the general case.
# We demonstrate the conceptual substitution.
term1_subst = -Integral(A * Sum(b[i] * diff(phi(i,x), x), (i, 1, N)) * Sum(a[j] * diff(phi(j,x), x), (j, 1, N)), (x, 0, L))
print_styled("Term 1 becomes:")
Lprint(term1_subst)

print_styled("Substituting into the second term: Integral(B * v * dy/dx, (x, 0, L))")
term2_subst = Integral(B * Sum(b[i] * phi(i,x), (i, 1, N)) * Sum(a[j] * diff(phi(j,x), x), (j, 1, N)), (x, 0, L))
print_styled("Term 2 becomes:")
Lprint(term2_subst)

print_styled("Substituting into the third term: Integral(C * v * y, (x, 0, L))")
term3_subst = Integral(C * Sum(b[i] * phi(i,x), (i, 1, N)) * Sum(a[j] * phi(j,x), (j, 1, N)), (x, 0, L))
print_styled("Term 3 becomes:")
Lprint(term3_subst)

print_styled("Substituting into the right-hand side term: Integral(v * F, (x, 0, L))")
rhs_subst = Integral(Sum(b[i] * phi(i,x), (i, 1, N)) * F, (x, 0, L))
print_styled("RHS term becomes:")
Lprint(rhs_subst)

#%% --- Rearranging the Summation and Integrals ---
print_styled("<h2>3. Rearrange Sums and Integrals</h2>")
print_styled("Assuming sufficient continuity, we can swap the order of integration and summation.")
print_styled("We factor out the arbitrary coefficients b_i from each term.")

print_styled("Consider Term 1:")
print_styled("$-A \\int_0^L \\left( \\sum_{i=1}^N b_i \\frac{d\\phi_i}{dx} \\right) \\left( \\sum_{j=1}^N a_j \\frac{d\\phi_j}{dx} \\right) dx = \\sum_{i=1}^N b_i \\left[ -A \\sum_{j=1}^N \\left( \\int_0^L \\frac{d\\phi_i}{dx} \\frac{d\\phi_j}{dx} dx \\right) a_j \\right]$")
# Symbolic representation of the inner part:
term1_inner = Sum( -A * Integral(diff(phi(i,x), x) * diff(phi(j,x), x), (x, 0, L)) * a[j], (j, 1, N))
# Full term rearranged
term1_rearranged = Sum(b[i] * term1_inner, (i, 1, N))
# Lprint(Eq(term1_subst, term1_rearranged)) # SymPy might struggle to prove this equality automatically

print_styled("Term 2:")
print_styled("$B \\int_0^L \\left( \\sum_{i=1}^N b_i \\phi_i \\right) \\left( \\sum_{j=1}^N a_j \\frac{d\\phi_j}{dx} \\right) dx = \\sum_{i=1}^N b_i \\left[ B \\sum_{j=1}^N \\left( \\int_0^L \\phi_i \\frac{d\\phi_j}{dx} dx \\right) a_j \\right]$")
term2_inner = Sum( B * Integral(phi(i,x) * diff(phi(j,x), x), (x, 0, L)) * a[j], (j, 1, N))
term2_rearranged = Sum(b[i] * term2_inner, (i, 1, N))

print_styled("Term 3:")
print_styled("$C \\int_0^L \\left( \\sum_{i=1}^N b_i \\phi_i \\right) \\left( \\sum_{j=1}^N a_j \\phi_j \\right) dx = \\sum_{i=1}^N b_i \\left[ C \\sum_{j=1}^N \\left( \\int_0^L \\phi_i \\phi_j dx \\right) a_j \\right]$")
term3_inner = Sum( C * Integral(phi(i,x) * phi(j,x), (x, 0, L)) * a[j], (j, 1, N))
term3_rearranged = Sum(b[i] * term3_inner, (i, 1, N))

print_styled("RHS:")
print_styled("$\\int_0^L \\left( \\sum_{i=1}^N b_i \\phi_i \\right) F(x) dx = \\sum_{i=1}^N b_i \\left( \\int_0^L \\phi_i F(x) dx \\right)$")
rhs_inner = Integral(phi(i,x) * F, (x, 0, L))
rhs_rearranged = Sum(b[i] * rhs_inner, (i, 1, N))


#%% --- Obtaining N Equations ---
print_styled("<h2>4. Obtain N Equations from Arbitrary b_i</h2>")

print_styled("Combining the rearranged terms:")
final_sum_eq = Eq(term1_rearranged + term2_rearranged + term3_rearranged, rhs_rearranged)
# Lprint(final_sum_eq) # This looks messy, let's represent it conceptually

print_styled("$\\sum_{i=1}^N b_i \\left( [\\text{Term 1 Inner}]_i + [\\text{Term 2 Inner}]_i + [\\text{Term 3 Inner}]_i \\right) = \\sum_{i=1}^N b_i [\\text{RHS Inner}]_i$")
print_styled("Rearranging everything to one side:")
print_styled("$\\sum_{i=1}^N b_i \\left( [\\text{Term 1 Inner}]_i + [\\text{Term 2 Inner}]_i + [\\text{Term 3 Inner}]_i - [\\text{RHS Inner}]_i \\right) = 0$")

print_styled("Since this equation must hold for *any* choice of coefficients $b_i$, the term within the parenthesis must be zero for *each* value of $i$ from 1 to N.")
print_styled("This gives us N distinct equations:")
equation_for_i = Eq(term1_inner + term2_inner + term3_inner, rhs_inner)
print_styled("For each $i = 1, 2, ..., N$:")
Lprint(equation_for_i)

#%% --- Matrix Form (System of Linear Equations) ---
print_styled("<h2>5. Assemble the Matrix System K a = F</h2>")
print_styled("The $i$-th equation involves a sum over $j$ of terms multiplying the unknown coefficients $a_j$.")
print_styled("Let's rewrite the $i$-th equation to isolate the sum over $a_j$:")

print_styled("$\\sum_{j=1}^N \\left[ -A \\int_0^L \\frac{d\\phi_i}{dx} \\frac{d\\phi_j}{dx} dx + B \\int_0^L \\phi_i \\frac{d\\phi_j}{dx} dx + C \\int_0^L \\phi_i \\phi_j dx \\right] a_j = \\int_0^L \\phi_i F(x) dx$")

print_styled("This is a system of N linear equations for the N unknowns $a_1, a_2, ..., a_N$. We can write it in matrix form:")
print_styled("$ K a = F $")

print_styled("Where:")
print_styled(" - $a$ is the column vector of unknown coefficients $[a_1, a_2, ..., a_N]^T$.")
print_styled(" - $F$ is the column vector (load vector) with entries $F_i$.")
print_styled(" - $K$ is the $N \\times N$ matrix (stiffness matrix) with entries $K_{ij}$.")

# Define K_ij and F_i symbolically
K_ij = -A * Integral(diff(phi(i,x), x) * diff(phi(j,x), x), (x, 0, L)) + \
        B * Integral(phi(i,x) * diff(phi(j,x), x), (x, 0, L)) + \
        C * Integral(phi(i,x) * phi(j,x), (x, 0, L))

F_i = Integral(phi(i,x) * F, (x, 0, L))

print_styled("The components are defined as:")
Lprint(Eq(Symbol('K_{ij}'), K_ij))
Lprint(Eq(Symbol('F_i'), F_i))

# Represent the matrix equation symbolically (for small N, e.g., N=3, if desired)
# N_example = 3
# K_matrix = Matrix(N_example, N_example, lambda r, c: K_ij.subs({i: r + 1, j: c + 1}))
# a_vector = Matrix(N_example, 1, lambda r, c: a[r + 1])
# F_vector = Matrix(N_example, 1, lambda r, c: F_i.subs({i: r + 1}))
# Lprint(Eq(K_matrix * a_vector, F_vector, evaluate=False)) # Display K a = F

print_styled("Now we have transformed the original differential equation into a system of algebraic equations $Ka = F$, which can be solved using linear algebra.")


#%% --- Boundary Conditions ---
print_styled("<h2>6. Note on Boundary Conditions</h2>")
print_styled("We assumed $y(0)=0$ and $y(L)=0$. How are these applied to the matrix system?")
print_styled("If using nodal basis functions where $\phi_k(x_{node}) = \delta_{knode}$, then $y(x_1) = y(0) = \sum a_j \phi_j(x_1) = a_1$. Similarly $y(x_N) = y(L) = a_N$.")
print_styled("So, for $y(0)=y(L)=0$, we know $a_1 = 0$ and $a_N = 0$.")
print_styled("These known values must be incorporated when solving $Ka=F$. Common methods include:")
print_styled("  - Modifying the K matrix and F vector (e.g., setting rows/columns to enforce the condition).")
print_styled("  - Partitioning the system into known and unknown degrees of freedom.")
print_styled("Handling boundary conditions precisely is a key part of FEM implementation and will be detailed later.")


#%% --- Introduction to Weighted Residuals and Galerkin ---
print_styled("<h2>7. Weighted Residual Method and Galerkin's Method</h2>")

print_styled("Let's revisit the core idea from a different perspective.")
print_styled("Our original differential equation is $A(y) = F$, where $A$ is the differential operator.")
print_styled("$A(y) = A_{coeff} \\frac{d^2y}{dx^2} + B_{coeff} \\frac{dy}{dx} + C_{coeff} y$") # Use different symbols A_coeff etc to avoid clash

A_op = lambda func: A * diff(func, x, 2) + B * diff(func, x) + C * func
strong_form_op = Eq(A_op(y), F)
# Lprint(strong_form_op) # Display operator form

print_styled("Our approximate solution $y^N(x) = \\sum_{j=1}^N a_j \\phi_j(x)$ likely does *not* perfectly satisfy the equation.")
print_styled("The difference is the **residual**, $r^N(x)$:")
r_N = Symbol('r^N')
# Need to use y_approx here for the definition
residual_def = Eq(r_N, A_op(y_approx) - F, evaluate=False)
Lprint(residual_def)

print_styled("The goal of FEM is to make this residual 'small' in some sense.")
print_styled("The **Method of Weighted Residuals** requires the residual to be orthogonal to a set of *weighting functions*, $w_i(x)$:")
# w = IndexedBase('w')
w = Function('w') # Represents the set of weighting functions w_i(x)
weighted_residual_eq = Eq(Integral(r_N * w(i,x), (x, 0, L)), 0)
print_styled("$\\int_0^L r^N(x) w_i(x) dx = 0$, for $i = 1, 2, ..., N$")
# Lprint(weighted_residual_eq)

print_styled("Different choices for $w_i$ lead to different methods (e.g., Collocation, Subdomain, Least Squares).")
print_styled("The **Galerkin Method** makes a specific choice: use the *same* functions for weighting as for the basis/shape functions.")
print_styled("Choose $w_i(x) = \\phi_i(x)$:")
galerkin_eq = Eq(Integral(r_N * phi(i,x), (x, 0, L)), 0)
print_styled("Galerkin's Method: $\\int_0^L r^N(x) \\phi_i(x) dx = 0$, for $i = 1, 2, ..., N$")
Lprint(galerkin_eq)

print_styled("Substituting the definition of $r^N$:")
galerkin_expanded = Eq(Integral((A_op(y_approx) - F) * phi(i,x), (x, 0, L)), 0)
Lprint(galerkin_expanded)

print_styled("$\\int_0^L \\left( A \\frac{d^2y^N}{dx^2} + B \\frac{dy^N}{dx} + C y^N - F \\right) \\phi_i(x) dx = 0$")

print_styled("If we perform integration by parts on the second derivative term (similar to deriving the weak form), we recover the *exact same set of equations* as we derived earlier by substituting into the weak form!")
print_styled("$\\int_0^L \\left( -A \\frac{d\\phi_i}{dx} \\frac{dy^N}{dx} + B \\phi_i \\frac{dy^N}{dx} + C \\phi_i y^N - \\phi_i F \\right) dx + [A \\phi_i \\frac{dy^N}{dx}]_0^L = 0$")
print_styled("The boundary term vanishes if $\phi_i(0) = \phi_i(L) = 0$ (as required for test functions corresponding to essential BCs).")
print_styled("Substituting $y^N = \\sum a_j \\phi_j$ leads back to the $Ka=F$ system.")

print_styled("Therefore, the Galerkin method provides a powerful theoretical justification for the procedure we followed: it minimizes the residual in a way that makes it orthogonal to the space spanned by the basis functions themselves.")

#%%