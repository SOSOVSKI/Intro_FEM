# Integration Plan: symbolic_fem_workbench into IntroFEM Lectures

## 1. Project Survey Summary

### 1.1 Lecture Content Inventory

The course has **five Quarto/reveal.js lecture decks** plus supporting Python modules:

| Lecture File | Topic | Slides (approx) | Key FEM Content |
|---|---|---|---|
| `Lecture1.qmd` | Introduction to FEM (1D) | ~10 | Strong form, weak form derivation, shape functions, assembly, 1D solver demo |
| `Lecture2.qmd` | Algebraic system & Galerkin | ~15 | Discretization with indexed basis functions, K·a = F derivation, boundary conditions, Galerkin method |
| `L3.qmd` | 1D Elasticity & Implementation | ~30 | Bar element, constitutive relations, reference element mapping, Gaussian quadrature, element assembly, post-processing |
| `Lectures.qmd` | Time-Dependent & 2D/3D | ~50 | Mass matrix, theta-method, 2D heat equation, divergence theorem, triangular/rectangular elements, error estimation, adaptivity |
| `Function_Spaces.qmd` | Mathematical Foundations | ~15 | Hilbert spaces, L^2, H^1, Sobolev spaces, weak derivatives |

Supporting Python modules: `old/L1.py` (symbolic display + 1D solver), `old/L2.py` (Galerkin derivation helpers), `src/Lectures/L3.py` (FEM diagrams, shape function visualizations).

All lectures are **Quarto reveal.js** presentations rendered to HTML, using **KaTeX** for math and embedded **Python code cells** (some interactive via ipywidgets).

### 1.2 Workbench Capabilities

The `symbolic_fem_workbench` package provides:

**Core pipeline (transforms.py, forms.py, extract.py):**
- Weighted residual construction from strong form
- Integration by parts for divergence terms
- Dirichlet/Neumann boundary condition handling
- Bilinear/linear form splitting
- Field substitution with FE expansions (Galerkin discretization)
- Local stiffness matrix and load vector extraction by coefficient matching

**Element library (fe_spaces.py, reference.py):**
- 1D: `LinearElement1D` (P1 interval)
- 2D: `ReferenceTriangleP1`, `ReferenceQuadrilateralQ1`, `ReferenceQuadrilateralQ2`
- 3D: `ReferenceTetrahedronP1/P2`, `ReferenceHexahedronQ1`
- `AffineTriangleMap2D` with Jacobian, determinant, and J^{-T} for gradient pullback

**Integration (quadrature.py):**
- Exact symbolic integration on triangle, quad, tet, hex
- Gauss-Legendre rules for quads/hexes (order 1-3)
- One-point and four-point rules for tet

**Assembly (assembly.py):**
- Dense matrix/vector assembly with connectivity
- Dirichlet BC by reduction (condensation)
- Solution expansion from reduced system

**Code generation (printers/, codegen/):**
- I-heart-LA notation export for local tensors

**Pre-built workflows (workflow.py):**
- `build_bar_1d_local_problem()` -- complete 1D bar/Poisson pipeline
- `build_poisson_triangle_p1_local_problem()` -- complete 2D P1 triangle pipeline

### 1.3 Examples and Tests

| File | What it demonstrates | Teaching value |
|---|---|---|
| `bar_1d_workflow.py` | Minimal 1D pipeline (12 lines) | Excellent starter |
| `triangle_p1_poisson_workflow.py` | 2D Poisson symbolic derivation | Good for 2D introduction |
| `manual_assembly_square_4tri.py` | Full pipeline: symbolic local tensors -> lambdify -> numerical assembly -> Dirichlet -> solve | Excellent end-to-end demo |
| `reference_element_library_demo.py` | Reference element shape function inspection (6 element types) | Good reference material |
| `test_triangle_p1.py` | Verifies Ke, fe against hand-calculated values for unit right triangle | Formula validation |
| `test_reference_elements_extended.py` | Partition of unity, gradient sum, Kronecker delta, quadrature accuracy | Axiom verification |
| `test_manual_assembly_square_4tri.py` | End-to-end solve with u(center) = 1/12 check | Pipeline correctness |


---

## 2. Lecture-to-Workbench Mapping

### 2.1 Direct Mappings

| Lecture | Week | Workbench Capability | Integration Opportunity |
|---|---|---|---|
| **Lecture 1** (Intro, 1D weak form) | 1-2 | `weighted_residual()`, `integrate_divergence_1d()` | Live derivation of weak form from strong form; show each step symbolically |
| **Lecture 2** (Galerkin, algebraic system) | 2-3 | `split_linear_weak_form()`, `substitute_fe()`, `extract_coefficient_matrix()` | Derive K·a = F from weak form; show matrix entries emerge from bilinear form |
| **L3** (1D elements, quadrature) | 3-5 | `LinearElement1D`, `build_bar_1d_local_problem()`, `quadrature.py` | Element stiffness derivation, exact vs. Gauss integration comparison |
| **Lectures.qmd** (2D/3D formulation) | 8-10 | `ReferenceTriangleP1`, `AffineTriangleMap2D`, `pullback_gradient_2d()`, `build_poisson_triangle_p1_local_problem()` | Reference element mapping, gradient pullback, local 2D tensor extraction |
| **Lectures.qmd** (time-dependent) | 7-8 | Gap: no mass matrix workflow yet | Opportunity: add `build_bar_1d_mass_matrix()` |
| **Lectures.qmd** (error estimation) | 10+ | Gap: no error estimation in workbench | Low priority -- stays as lecture-only theory |
| **Function_Spaces** | 5-6 | Not directly related | No workbench integration needed |

### 2.2 Coverage Gaps

These lecture topics have no corresponding workbench support:

1. **Mass matrix / time-dependent problems** -- Lectures.qmd covers theta-method and M*a_ddot + K*a = R, but no workbench workflow for mass matrix extraction.
2. **Vector problems (2D elasticity)** -- Lectures.qmd discusses plane stress/strain, but the workbench only handles scalar equations.
3. **Higher-order elements in workflows** -- Reference elements for Q1, Q2, P2 exist but have no corresponding workflow.py pipelines.
4. **Error estimation and adaptivity** -- Purely theoretical in lectures; no symbolic counterpart needed.
5. **Boundary integrals in 2D** -- No edge integration helpers for Neumann/Robin BCs on triangle edges.


---

## 3. Proposed Interactive Examples and Demos

### 3.1 New Jupyter Notebook Demos (for discussion sections)

#### Demo 1: "Strong Form to Weak Form" (Week 2)
- **Workbench functions**: `weighted_residual()`, `integrate_divergence_1d()`, `drop_dirichlet_boundary()`, `apply_neumann_flux()`
- **What students see**: Step-by-step symbolic transformation with each intermediate expression printed and explained.
- **Interactive element**: Students modify the strong form (change coefficients, add terms) and re-derive.

#### Demo 2: "From Weak Form to K*a = F" (Week 3)
- **Workbench functions**: `split_linear_weak_form()`, `LinearElement1D`, `local_trial_expansion()`, `substitute_fe()`, `extract_coefficient_matrix()`
- **What students see**: The bilinear form with symbolic shape functions substituted, then coefficient extraction producing the element stiffness matrix.
- **Interactive element**: Compare P1 shape functions on different element lengths.

#### Demo 3: "Exact vs. Gauss Integration" (Week 5)
- **Workbench functions**: `integrate_reference_triangle_exact()`, `quadrilateral_gauss_rule()`, `_gauss_legendre_1d()`
- **What students see**: Side-by-side comparison of exact symbolic integral vs. quadrature approximation for polynomial and non-polynomial integrands.
- **Interactive element**: Vary quadrature order and watch the error shrink.

#### Demo 4: "Reference Element Gallery" (Week 8)
- **Workbench functions**: All `Reference*` classes from `reference.py`
- **What students see**: Node locations, shape functions, partition of unity, and Kronecker property for each element type.
- **Interactive element**: Students pick an element type and inspect its properties. Matplotlib visualizations of 2D shape functions.

#### Demo 5: "2D Triangle Assembly from Scratch" (Week 9-10)
- **Workbench functions**: `build_poisson_triangle_p1_local_problem()`, `assembly.py` functions
- **What students see**: The `manual_assembly_square_4tri.py` example, but as an annotated notebook with each step explained.
- **Interactive element**: Students modify the mesh (move the center node, change triangle connectivity) and observe how the global system changes.

#### Demo 6: "Gradient Pullback Visualized" (Week 9)
- **Workbench functions**: `AffineTriangleMap2D`, `pullback_gradient_2d()`, `grad_2d()`
- **What students see**: How J^{-T} transforms reference gradients to physical gradients on distorted triangles.
- **Interactive element**: Drag triangle vertices and see how the Jacobian and gradient fields change.

### 3.2 Embedding Workbench Outputs in Lecture Slides

The Quarto slides already use Python code cells with `{python}` blocks. Workbench outputs can be embedded directly:

```qmd
## Local Stiffness Matrix Derivation

```{python}
#| code-fold: true
from symbolic_fem_workbench.workflow import build_bar_1d_local_problem
result = build_bar_1d_local_problem()
display(Math(sp.latex(result["Ke"])))
```

This works because Quarto renders Python cells during build. The symbolic SymPy expressions produce LaTeX that KaTeX renders in the slides.

**Specific embedding opportunities:**

| Slide Topic | What to Embed | Output Format |
|---|---|---|
| Weak form derivation | `weighted_residual()` output | LaTeX equation |
| Integration by parts | `integrate_divergence_1d()` result showing domain + boundary terms | LaTeX with annotations |
| Element stiffness matrix (1D) | `result["Ke"]` from `build_bar_1d_local_problem()` | LaTeX 2x2 matrix |
| Element stiffness matrix (2D) | `result["Ke"]` from `build_poisson_triangle_p1_local_problem()` | LaTeX 3x3 matrix |
| Reference element mapping | `AffineTriangleMap2D` properties | LaTeX for J, det(J), J^{-T} |
| Shape function plots | Matplotlib output from `reference.py` shape functions | PNG figures |
| Assembly visualization | Print intermediate global K during assembly loop | Formatted matrix output |


---

## 4. Required Modifications to the Workbench

### 4.1 API Cleanup (Priority: High)

1. **Unify 1D element representation.** Currently 1D elements live in `fe_spaces.py` while 2D/3D elements are in `reference.py`. Create a `ReferenceIntervalP1` class in `reference.py` that mirrors the 2D/3D pattern, keeping `LinearElement1D` as a convenience alias.

2. **Add docstrings to all public functions.** Every function in `transforms.py`, `extract.py`, and `workflow.py` needs a one-paragraph docstring explaining: what it does, what it returns, and what FEM step it corresponds to. Students will read docstrings.

3. **Improve error messages.** `substitute_field()` is the most complex function and produces opaque SymPy errors on misuse. Add input validation with clear messages (e.g., "Expected a Function application like u(x), got Symbol").

4. **`__init__.py` exports are already clean.** The package already has comprehensive `__all__` with 30+ symbols. Consider adding the remaining reference elements (`ReferenceQuadrilateralQ1`, `ReferenceTetrahedronP1`, etc.) to the top-level exports as the course expands to cover those element types.

### 4.2 New Workflows (Priority: Medium)

5. **Add `build_bar_1d_mass_matrix()` to workflow.py.** Derive the consistent mass matrix M_e for a 1D bar element. This bridges to the time-dependent content in Lectures.qmd.

6. **Add `build_poisson_quad_q1_local_problem()` to workflow.py.** The Q1 reference element exists but has no workflow. This would demonstrate bilinear shape functions and tensor-product integration.

7. **Add triangle Gauss quadrature rules.** Currently only 1-point rule on triangle. Add 3-point and 6-point rules to `quadrature.py` for the "exact vs. approximate" comparison demo.

### 4.3 Documentation (Priority: Medium)

8. **Create a `docs/` directory** with a module-by-module API reference. Can be generated from docstrings with sphinx or pdoc.

9. **Add a GETTING_STARTED.md** that walks through the `bar_1d_workflow.py` example line-by-line with explanations.

### 4.4 Notebook Wrappers (Priority: High)

10. **Create `notebooks/` directory** with Jupyter notebooks wrapping each example into a teaching-ready format:
    - `01_strong_to_weak_form.ipynb`
    - `02_galerkin_discretization.ipynb`
    - `03_1d_element_stiffness.ipynb`
    - `04_exact_vs_gauss_integration.ipynb`
    - `05_reference_element_gallery.ipynb`
    - `06_2d_triangle_poisson.ipynb`
    - `07_manual_assembly_2d.ipynb`

    Each notebook should have:
    - Markdown cells explaining the FEM theory
    - Code cells using the workbench with step-by-step output
    - Exercises (empty cells for students to fill in)
    - Verification cells that check student work against known answers

### 4.5 Visualization Helpers (Priority: Low-Medium)

11. **Add a `viz.py` module** with matplotlib-based helpers for:
    - Plotting 2D shape functions on reference and physical triangles
    - Showing the affine mapping visually (reference -> physical)
    - Visualizing sparsity patterns of assembled matrices
    - Plotting FEM solutions on triangular meshes

12. **Integrate with existing L3.py visualization code.** The `FEMShapeFunctions` class in `src/Lectures/L3.py` already has good shape function plotting. Factor this into the workbench package or create a bridge module.


---

## 5. VS Code Extensions and Tooling

### 5.1 For Working with Lecture HTML/QMD Files

| Extension | Purpose | Required? |
|---|---|---|
| **Quarto** (`quarto-dev.quarto`) | Render .qmd files, preview reveal.js slides, run embedded Python cells | Yes |
| **Jupyter** (`ms-toolsai.jupyter`) | Run notebook cells within .qmd files | Yes |
| **Live Preview** (`ms-vscode.live-server`) | View rendered HTML slides locally | Recommended |
| **Marp** (alternative) | Only if migrating away from Quarto | No |

### 5.2 For Working with symbolic_fem_workbench

| Extension | Purpose | Required? |
|---|---|---|
| **Python** (`ms-python.python`) | Python language support, IntelliSense | Yes |
| **Pylance** (`ms-python.vscode-pylance`) | Type checking, go-to-definition | Recommended |
| **Jupyter** (`ms-toolsai.jupyter`) | Run .ipynb notebooks for demos | Yes |
| **Ruff** (`charliermarsh.ruff`) | Linting (matches pyproject.toml config) | Recommended |
| **mypy** (via Pylance or `matangover.mypy`) | Type checking (configured in pyproject.toml) | Optional |

### 5.3 For Math Rendering

| Extension | Purpose | Required? |
|---|---|---|
| **Markdown+Math** (`goessner.mdmath`) | Preview KaTeX in markdown files | Recommended |
| **LaTeX Workshop** (`james-yu.latex-workshop`) | Only if exporting to LaTeX/PDF | Optional |

### 5.4 Build Tooling

| Tool | Purpose | Installation |
|---|---|---|
| **Quarto CLI** | Render .qmd to HTML slides | `brew install quarto` or quarto.org |
| **uv** | Python package manager (project uses uv.lock) | `curl -LsSf https://astral.sh/uv/install.sh \| sh` |
| **SymPy >= 1.13** | Core dependency | `uv sync` |
| **NumPy >= 2.0** | Numerical assembly | `uv sync` |
| **pytest >= 8.0** | Test runner | `uv sync --group dev` |
| **ruff >= 0.6** | Linter | `uv sync --group dev` |


---

## 6. Integration Sequencing

### Phase 1: Foundation (Weeks 1-2 of development)

**Goal:** Make the workbench importable from lecture code cells and create the first two teaching notebooks.

1. Review `__init__.py` exports (already clean with `__all__`; add missing reference elements as needed).
2. Add docstrings to all public functions in `transforms.py`, `extract.py`, `workflow.py`.
3. Create `notebooks/01_strong_to_weak_form.ipynb` wrapping the 1D weak form derivation.
4. Create `notebooks/02_galerkin_discretization.ipynb` wrapping the shape function substitution and K extraction.
5. Verify that `uv sync` installs cleanly and all tests pass.
6. Add workbench imports to `Lecture1.qmd` and `Lecture2.qmd` code cells (replace or supplement the ad-hoc L1.py/L2.py symbolic code).

### Phase 2: 1D Completion (Weeks 3-4)

**Goal:** Full 1D coverage with integration into Lecture 1, 2, L3.

7. Create `notebooks/03_1d_element_stiffness.ipynb`.
8. Create `notebooks/04_exact_vs_gauss_integration.ipynb`.
9. Add triangle Gauss quadrature rules (3-point, 6-point) to `quadrature.py`.
10. Add `build_bar_1d_mass_matrix()` to `workflow.py`.
11. Embed workbench outputs in `L3.qmd` code cells for element derivation slides.
12. Update `L3.py` visualization code to use workbench shape functions where applicable.

### Phase 3: 2D Integration (Weeks 5-7)

**Goal:** Full 2D P1 triangle coverage with assembly demo.

13. Create `notebooks/05_reference_element_gallery.ipynb`.
14. Create `notebooks/06_2d_triangle_poisson.ipynb` (gradient pullback, local tensor extraction).
15. Create `notebooks/07_manual_assembly_2d.ipynb` (the crown jewel: full pipeline on 4-triangle mesh).
16. Add visualization helpers for 2D shape functions and affine mapping.
17. Embed 2D workbench outputs in `Lectures.qmd` code cells for the 2D formulation section.
18. Add Q1 quad workflow to `workflow.py` (optional, for comparison with triangles).

### Phase 4: Polish and Documentation (Weeks 8-9)

**Goal:** Production-ready teaching materials.

19. Write `GETTING_STARTED.md` for students.
20. Generate API reference docs (pdoc or sphinx).
21. Add input validation and clear error messages to `substitute_field()` and other complex functions.
22. Create a `ReferenceIntervalP1` class in `reference.py` for 1D consistency.
23. Run all tests, fix any regressions, verify notebook outputs.
24. Review all lecture QMD files for consistency with workbench notation.

### Phase 5: Advanced Features (Optional, Weeks 10+)

25. Add 2D boundary integral helpers for edge Neumann conditions.
26. Add vector Poisson / elasticity support (plane strain/stress).
27. Create I-heart-LA export demo notebook.
28. Add interactive widgets (ipywidgets) to notebooks for parametric exploration.
29. Consider Quarto dashboard format for the reference element gallery.


---

## 7. File Structure After Integration

```
IntroFem/
  src/
    Lectures/
      Lecture1.qmd          # Updated: imports workbench for live derivation
      Lecture2.qmd          # Updated: imports workbench for Galerkin demo
      L3.qmd                # Updated: imports workbench for element stiffness
      Lectures.qmd          # Updated: imports workbench for 2D content
      Function_Spaces.qmd   # Unchanged
      L3.py                 # Updated: uses workbench shape functions
    symbolic_fem_workbench/
      __init__.py            # Updated: clean __all__ exports
      symbols.py             # Docstrings added
      forms.py               # Docstrings added
      transforms.py          # Docstrings + input validation added
      fe_spaces.py           # Docstrings added
      reference.py           # ReferenceIntervalP1 added, docstrings
      quadrature.py          # Triangle 3/6-point rules added
      extract.py             # Docstrings added
      validate.py            # Docstrings added
      workflow.py            # Mass matrix + Q1 quad workflows added
      assembly.py            # Docstrings added
      viz.py                 # NEW: visualization helpers
      printers/
      codegen/
  notebooks/                 # NEW: teaching notebooks
    01_strong_to_weak_form.ipynb
    02_galerkin_discretization.ipynb
    03_1d_element_stiffness.ipynb
    04_exact_vs_gauss_integration.ipynb
    05_reference_element_gallery.ipynb
    06_2d_triangle_poisson.ipynb
    07_manual_assembly_2d.ipynb
  examples/                  # Unchanged (kept as minimal scripts)
  tests/                     # Unchanged + new tests for added features
  docs/
    GETTING_STARTED.md       # NEW: student-facing quick start
    API_REFERENCE.md         # NEW: generated module docs
  COURSE_INTEGRATION_NOTES.md
  SYLLABUS_REVISION_DRAFT.md
  TEACHING_NOTE_MANUAL_ASSEMBLY_2D.md
  README.md
  pyproject.toml
```


---

## 8. Risk Assessment

| Risk | Impact | Mitigation |
|---|---|---|
| SymPy rendering breaks in Quarto slides | High | Test early with a minimal QMD cell that imports workbench; KaTeX may need `sp.latex()` wrapping |
| Workbench API changes break existing examples | Medium | Run `pytest` in CI after every change; examples double as integration tests |
| Notebook cells are too slow for live demos | Medium | Pre-compute heavy symbolic results; use `sp.simplify()` sparingly (it's slow) |
| Students can't install the workbench | Low | Provide a `uv sync` one-liner; consider a requirements.txt fallback |
| 2D visualization is too complex for matplotlib | Low | Keep visualizations simple; use plotly only if matplotlib is insufficient |
| Scope creep into vector problems delays core work | Medium | Defer elasticity (Phase 5) until scalar pipeline is fully integrated |
