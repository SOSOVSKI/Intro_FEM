# Build Verification Report

**Date:** 2026-04-12
**Environment:** Linux sandbox (Python 3.10.12, SymPy 1.13+, NumPy 2.0+)
**Quarto CLI:** Not available in sandbox (GitHub/quarto.org downloads blocked by network egress). All verification below uses static QMD analysis plus full Python cell execution.

## 1. Verification Methodology

Since Quarto CLI cannot be installed in the sandbox environment, this report uses three complementary approaches:

1. **Static analysis** of all 18 QMD files (YAML headers, LaTeX math, code cells, image references, cross-file dependencies)
2. **Python cell execution** — every `{python}` code block extracted from QMDs and executed with the correct PYTHONPATH
3. **Notebook execution** — all 12 Jupyter notebooks executed end-to-end
4. **pytest test suite** — all 34 unit tests run

## 2. Issues Found and Fixed

### 2.1 ERRORS (would prevent compilation)

| # | File | Issue | Fix Applied |
|---|------|-------|-------------|
| 1 | `Lecture1.qmd` line 33 | `sys.path.insert(0, "../symbolic_fem_workbench/../../")` resolves to project root (wrong), cannot import workbench | Changed to `sys.path.insert(0, os.path.abspath(".."))` which correctly resolves to `src/` |
| 2 | `Lecture2.qmd` line 377 | Same broken sys.path | Same fix |
| 3 | `L3.qmd` line 375 | Same broken sys.path | Same fix |
| 4 | `Lectures.qmd` line 442 | Same broken sys.path | Same fix |
| 5 | `Function_Spaces.qmd` line 224 | LaTeX typo: `$L^2$$` (double closing dollar) parses as `$L^2$` followed by an empty `$$` block | Changed to `$L^2$` |
| 6 | Notebook `09_boundary_conditions_2d.ipynb` cell 3 | `ReferenceTriangleP1()` called without required `xi, eta` arguments; symbols defined after constructor | Reordered: define `xi, eta` first, then `ReferenceTriangleP1(xi, eta)` |
| 7 | Notebook `08_boundary_conditions_1d.ipynb` cell 16 | Non-homogeneous lifting case: `expand_reduced_solution(..., prescribed_values=[0.0, 1.0])` double-counts prescribed values when added to `u_lift` | Changed to `prescribed_values=[0.0, 0.0]` with comment explaining the lift already carries prescribed values |

### 2.2 WARNINGS (compile but display incorrectly)

| # | File | Issue | Status |
|---|------|-------|--------|
| 1 | `ch10-elasticity-3d.qmd` line 15 | `![](images/3dElast.png)` — initially flagged as missing | Verified: file exists in `New/images/` (false alarm) |
| 2 | Multiple QMD files | Long equations (6x6 D-matrix, Jacobian chains) may overflow PDF margins | **Not fixed** — recommend adding `code-overflow: wrap` and `\small` in PDF format when Quarto is available |
| 3 | `ch04-1d-fem.qmd`, `ch07-b-matrix.qmd` | Code cells with lines >80 chars may overflow in PDF | **Not fixed** — same recommendation |

### 2.3 INFO (cosmetic)

| # | File | Issue |
|---|------|-------|
| 1 | `Lecture2.qmd` | Unused import `Q` from sympy |
| 2 | Multiple legacy QMDs | Missing `#| fig-cap` labels on code cell outputs |
| 3 | `ch01`–`ch04` | `%run _common.py` used repeatedly (only needed once) |

## 3. Verification Results by Component

### 3.1 Unit Test Suite

```
34 passed in 12.41s
```

All test files:
- `test_3d_elasticity.py` (7 tests) — PASS
- `test_elasticity.py` (7 tests) — PASS
- `test_manual_assembly_square_4tri.py` (1 test) — PASS
- `test_new_features.py` (10 tests) — PASS
- `test_reference_elements_extended.py` (8 tests) — PASS
- `test_triangle_p1.py` (1 test) — PASS

### 3.2 Jupyter Notebooks (12 total)

| Notebook | Cells | Result | Notes |
|----------|-------|--------|-------|
| 01_strong_to_weak_form.ipynb | 7 code | PASS | All symbolic derivations produce correct output |
| 02_galerkin_discretization.ipynb | 8 code | PASS | Shape function substitution and extraction correct |
| 03_1d_element_stiffness.ipynb | 6 code | PASS | Lambdify, numerical verification, parameter study all work |
| 04_exact_vs_gauss_integration.ipynb | 6 code | PASS | 3-point triangle rule exact for quadratic confirmed |
| 05_reference_element_gallery.ipynb | 8 code | PASS | All 6 element types display correctly |
| 06_2d_triangle_poisson.ipynb | 8 code | PASS | Gradient pullback and Ke extraction verified |
| 07_manual_assembly_2d.ipynb | 7 code | PASS | Center node u=1/12 confirmed |
| 08_boundary_conditions_1d.ipynb | 10 code | PASS | Row substitution and lifting match (both homogeneous and non-homogeneous) |
| 09_boundary_conditions_2d.ipynb | 12 code | PASS | 2D mixed BCs, both methods match |
| 10_2d_elasticity.ipynb | 8 code | PASS | Plane stress/strain D-matrices, B-matrix, Ke symmetry confirmed |
| 11_3d_elasticity.ipynb | 7 code | SKIP | Exceeds sandbox compute limits (symbolic 12x12 matrix); library code verified via unit tests |
| 12_iheartla_and_widgets.ipynb | 8 code | PASS | I-heart-LA export, interactive widgets (ipywidgets detected) |

**Overall: 11/12 PASS, 1 SKIP (resource limit only)**

### 3.3 Legacy Lecture QMDs (Workbench Code Cells)

After fixing the sys.path issues, all workbench imports resolve correctly:

| File | Workbench cells | Result |
|------|----------------|--------|
| Lecture1.qmd | 2 cells (setup + 2 new slides) | PASS — `build_bar_1d_local_problem()` executes, Ke and fe display correctly |
| Lecture2.qmd | 1 cell (Galerkin slide) | PASS — shape functions, Ke, fe all display |
| L3.qmd | 1 cell (element derivation + mass matrix) | PASS — `build_bar_1d_mass_matrix()` produces correct Me |
| Lectures.qmd | 2 cells (reference elements + 2D stiffness) | PASS — Q1 quad and P1 triangle shape functions, gradient pullback, Ke_unit all display |

### 3.4 New Book Chapters (New/ directory)

| File | Code cells | LaTeX | Images | Status |
|------|-----------|-------|--------|--------|
| index.qmd | 0 | N/A | N/A | PASS (pure markdown) |
| ch01-introduction.qmd | 7 | OK | OK | PASS (depends on python/L1.py via _common.py) |
| ch02-algebraic-system.qmd | 9 | OK | N/A | PASS |
| ch03-function-spaces.qmd | 2 | OK | N/A | PASS |
| ch04-1d-fem.qmd | 13 | OK | OK (shape function PNGs) | PASS |
| ch05-time-dependent.qmd | 0 | OK | N/A | PASS (pure math markdown) |
| ch06-math-prelim-2d3d.qmd | 0 | OK | 1 (GLmap.png) | PASS |
| ch07-b-matrix.qmd | 0 | OK | 5 (hex20, elementscases, etc.) | PASS — all images verified present |
| ch08-heat-transfer-2d.qmd | 0 | OK | N/A | PASS |
| ch09-continuum-mechanics.qmd | 0 | OK | N/A | PASS |
| ch10-elasticity-3d.qmd | 0 | OK | 1 (3dElast.png) | PASS — image verified in New/images/ |
| ch11-fem-implementation.qmd | 0 | OK | N/A | PASS |
| ch12-error-estimates.qmd | 0 | OK | 1 (mref.png) | PASS |

### 3.5 `_quarto.yml` Validation

- All 12 chapter files referenced in `_quarto.yml` exist
- Format configs (html + pdf) are syntactically valid
- `execute: freeze: auto` correctly configured for incremental builds
- LaTeX preamble includes amsmath, amssymb, booktabs
- Custom CSS file (`custom.css`) exists

## 4. PDF-Specific Concerns (Predicted)

Since Quarto PDF rendering is not possible in this environment, these are predicted issues based on static analysis of the content:

| Concern | Files Affected | Severity | Recommended Fix |
|---------|---------------|----------|-----------------|
| Wide 6x6 constitutive matrices | ch10-elasticity-3d.qmd, 10_2d_elasticity.ipynb | MEDIUM | Add `\small` or `\footnotesize` environment in LaTeX header, or use `\resizebox` |
| Long derivative chains | ch07-b-matrix.qmd lines 344-369 | LOW | Already uses aligned environment; should fit A4 with 25mm margins |
| Code blocks with long lines | ch04-1d-fem.qmd, notebooks with matrix prints | MEDIUM | Add `code-overflow: wrap` to `_quarto.yml` PDF format section |
| Algorithm pseudo-code boxes | ch11-fem-implementation.qmd | LOW | Current callout box formatting should render; verify font size |

**Recommended addition to `_quarto.yml` for PDF:**
```yaml
format:
  pdf:
    code-overflow: wrap
    code-block-font-size: \small
```

## 5. Summary

| Component | Items | Pass | Fail | Skip |
|-----------|-------|------|------|------|
| Unit tests | 34 | 34 | 0 | 0 |
| Notebooks | 12 | 11 | 0 | 1 |
| Legacy QMD workbench cells | 6 | 6 | 0 | 0 |
| New book chapters | 12 | 12 | 0 | 0 |
| _quarto.yml config | 1 | 1 | 0 | 0 |

**Critical bugs fixed in this verification pass:** 7 (4 sys.path, 1 LaTeX typo, 1 notebook constructor, 1 notebook lifting logic)

**Remaining action items for when Quarto CLI is available:**
1. Run `quarto render src/Lectures/New/ --to html` and verify no runtime errors
2. Run `quarto render src/Lectures/New/ --to pdf` and inspect for overflow
3. Add `code-overflow: wrap` to `_quarto.yml` PDF config if overflow is found
4. Verify KaTeX rendering of complex equations in HTML output
5. Run `quarto render` on legacy lecture QMDs individually (they are standalone, not part of the book project)
