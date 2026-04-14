# Getting Started with Symbolic FEM Workbench

## What is This Package?

**Symbolic FEM Workbench** is a teaching-first finite element method package. It is **not a numerical solver**. Instead, it helps students learn FEM by building finite-element formulations *symbolically* using SymPy.

Use this package to:
- Derive weak forms from strong-form PDEs (strong form → weighted residual → weak form)
- Apply integration by parts and enforce boundary conditions step-by-step
- Discretize weak forms with finite-element shape functions
- Extract local element matrices (stiffness K, mass M, load vectors f) as SymPy expressions
- Visualize reference elements and their mappings to physical coordinates

This is not for solving large problems—it's for understanding the mathematics of FEM.

## Installation

1. Make sure you have Python 3.9+ installed.
2. Install the package and its dependencies:
   ```bash
   cd /path/to/IntroFem
   uv sync
   ```

This creates a virtual environment with all required packages (SymPy, NumPy, Jupyter, pytest, etc.).

## Running Your First Example

After installation, run the 1D bar/Poisson example:

```bash
python examples/bar_1d_workflow.py
```

You should see the element stiffness matrix and load vector printed to the console:

```
Ke =
[E*A/L    -E*A/L  ]
[-E*A/L    E*A/L  ]
fe =
[P*L/2 - q*L/2]
[P*L/2 - q*L/2]
```

## Walkthrough: bar_1d_workflow.py

Let's break down what happens in `examples/bar_1d_workflow.py`:

```python
from symbolic_fem_workbench.workflow import build_bar_1d_local_problem

problem = build_bar_1d_local_problem()
```

This single function call does everything:

1. **Define symbols and fields** (inside the function):
   - Spatial coordinate: `x`
   - Domain length: `L`
   - Material properties: `E` (Young's modulus), `A` (cross-section area)
   - Loads: `q` (distributed), `P` (point load at right end)
   - Fields: `u(x)` (trial), `v(x)` (test)

2. **Build the strong form**:
   - Strong form: `d/dx(EA·du/dx) + q = 0` on `[0, L]`
   - Boundary conditions: `u(0) = 0` (Dirichlet), `EA·du/dx|_{x=L} = P` (Neumann)

3. **Apply integration by parts**:
   - Multiply residual by test function `v`
   - Integrate: `∫ (d(EA·du/dx)/dx)·v dx + ∫ q·v dx = 0`
   - Use integration by parts: `∫ EA·du/dx·dv/dx dx - [EA·du/dx·v]_boundary + ∫ q·v dx = 0`

4. **Enforce boundary conditions**:
   - At `x=0`: `v(0) = 0` (test space vanishes on Dirichlet boundary)
   - At `x=L`: boundary term becomes `P·v(L)`

5. **Discretize with P1 shape functions**:
   - Shape functions on `[0, L]`: `N₁ = (L-x)/L`, `N₂ = x/L`
   - Trial expansion: `u_h = N₁·d₀ + N₂·d₁`
   - Test expansion: `v_h = N₁·w₀ + N₂·w₁`

6. **Extract Ke and fe**:
   - Substitute expansions into the bilinear form `a(u_h, v_h)`
   - Extract coefficient of each `wᵢ·dⱼ` pair → element stiffness matrix Ke
   - Substitute into linear form `F(v_h)` (the RHS)
   - Extract coefficient of each `wᵢ` → element load vector fe

The result is a dict with:
- `problem["Ke"]`: the 2×2 element stiffness matrix
- `problem["fe"]`: the 2×1 element load vector
- All intermediate symbols and forms for inspection

## Using Jupyter Notebooks

Interactive tutorials are in `notebooks/`:

```bash
jupyter notebook notebooks/01_strong_to_weak_form.ipynb
```

These notebooks walk through the FEM derivation step-by-step with visualizations and explanations.

## Package Structure

- **`symbols.py`**: Create symbolic fields (`make_field_1d`, `make_field_2d`)
- **`forms.py`**: Containers for domain integrals, boundary terms, weak forms
- **`transforms.py`**: Core operations (integration by parts, field substitution, weak-form splitting)
- **`fe_spaces.py`**: Shape functions and element expansions (`LinearElement1D`, `local_trial_expansion`)
- **`reference.py`**: Reference-element classes (triangles, quads, tetrahedra) and affine mappings
- **`quadrature.py`**: Exact and quadrature integration rules
- **`extract.py`**: Extract matrices and vectors from symbolic expressions
- **`validate.py`**: Input validation helpers
- **`workflow.py`**: High-level functions (`build_bar_1d_local_problem`, `build_bar_1d_mass_matrix`)
- **`assembly.py`**: Manual assembly helpers for teaching examples
- **`viz.py`**: Visualization helpers (shape functions, mappings, solutions)

## Common Pitfalls

### 1. SymPy simplify() is slow
Calling `sp.simplify()` on complex expressions can take minutes. Use `sp.expand()` instead for most operations. Simplification is useful only when the expression is already small.

### 2. Create fields with make_field_1d, not raw symbols
**Wrong:**
```python
u = sp.symbols('u')  # This is a scalar symbol, not a field
```

**Correct:**
```python
from symbolic_fem_workbench import make_field_1d
u = make_field_1d('u', x)  # This is u(x), undefined function application
```

The difference matters for differentiation and substitution.

### 3. Use expand() before extracting coefficients
Before calling `extract_coefficient_matrix()` or `extract_coefficient_vector()`, always expand your expression:

```python
a_disc = sp.expand(a_disc)  # Expand first
Ke = extract_coefficient_matrix(a_disc, trial_dofs, test_dofs)
```

Unexpanded expressions may have nested structure that `coeff()` cannot navigate.

## Next Steps

1. Read `API_REFERENCE.md` for detailed function signatures
2. Explore `examples/` for more case studies
3. Work through the Jupyter notebooks
4. Modify existing examples and observe how the matrices change

Happy learning!
