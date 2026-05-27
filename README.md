# Introduction to the Finite Element Method

Course materials for the *Introduction to the Finite Element Method* graduate course at the Technion — Israel Institute of Technology.

The repository contains three interconnected components:

| Component | What it is | Where |
|---|---|---|
| **Lecture notes** | Quarto book (HTML + PDF) covering theory through implementation | `book/` |
| **Companion notebooks** | 14 Jupyter notebooks for hands-on practice | `notebooks/` |
| **`symbolic_fem_workbench`** | Teaching Python package used by the notebooks | `src/symbolic_fem_workbench/` |

---

## Quick start

### 1. Set up the environment

FEniCSx (chapters 14 and 16) is only available through conda-forge. Everything else is pure Python.

```bash
conda env create -f environment.yml
conda activate fem-env
uv pip install -e ".[lectures,dev]"
```

### 2. Build everything

```bash
make all
```

This renders the lecture-note book (HTML + PDF), the companion notebooks (HTML + PDF), and the revealjs lecture slides.

### 3. Open the materials

| Output | Path |
|---|---|
| Book (HTML) | `book/_output/index.html` |
| Book (PDF) | `book/_output/Introduction-to-the-Finite-Element-Method.pdf` |
| Notebooks (HTML) | `notebooks/_output/index.html` |
| Notebooks (PDF) | `notebooks/_output/FEM-Companion-Notebooks.pdf` |
| Lecture slides | `book/_slides/ch01-introduction.html` … |
| Package docs | `docs/_output/index.html` |

---

## Repository layout

```
IntroFem/
├── book/                  # Lecture note source (Quarto book)
│   ├── _quarto.yml        #   Book configuration
│   ├── index.qmd          #   Preface
│   ├── ch01–ch16 *.qmd    #   Chapters
│   ├── _slides/           #   Generated revealjs slides (gitignored)
│   └── _output/           #   Generated HTML + PDF (gitignored)
│
├── notebooks/             # Companion Jupyter notebooks
│   ├── _quarto.yml        #   Book configuration
│   ├── 01–14 *.ipynb      #   Exercises (14 notebooks)
│   └── _output/           #   Generated HTML + PDF (gitignored)
│
├── src/
│   └── symbolic_fem_workbench/   # Teaching Python package
│
├── docs/                  # Package documentation (Quarto website)
├── examples/              # Standalone usage examples
├── tests/                 # Test suite (pytest)
├── environment.yml        # conda environment (Python 3.11 + FEniCSx)
├── pyproject.toml         # Python package metadata
└── Makefile               # All build targets
```

---

## Build targets

```bash
make all               # build everything
make book              # lecture notes HTML + PDF
make book-html         # HTML only
make book-pdf          # PDF only
make notebooks         # companion notebooks HTML + PDF
make slides            # revealjs lecture slides
make docs              # package documentation website
make purge             # remove all build outputs and freeze caches
```

---

## Course structure

The lecture notes cover FEM from first principles through 3D implementation in FEniCSx.

**Part I — Foundations:** Strong and weak forms, algebraic system, function spaces

**Part II — 1D Finite Elements:** Shape functions, assembly, time-dependent problems

**Part III — Multi-Dimensional FEM:** Reference elements, isoparametric mapping, 2D heat transfer

**Part IV — Elasticity:** Continuum mechanics, 3D elasticity

**Part V — Implementation & Error Analysis:** FEM software, error estimates, FEniCSx tutorial

**Part VI — Structural Elements:** Beam theory (Euler-Bernoulli and Timoshenko), FEniCSx implementation

---

## The `symbolic_fem_workbench` package

A small, teaching-first symbolic FEM library built on SymPy. Every step from strong form to element matrix is explicit and inspectable. It is **not** a production FEM solver — see `docs/` for the full rationale and API reference.

```python
from symbolic_fem_workbench.workflow import build_bar_1d_local_problem
result = build_bar_1d_local_problem()
print(result["Ke"])   # 2×2 symbolic stiffness matrix
```

---

## Requirements

- Python 3.11
- conda / mamba (for FEniCSx)
- [Quarto](https://quarto.org) ≥ 1.4
- pdflatex (TeX Live 2023 or later)


---

## Interactive Workbench (Streamlit)

This repository can be paired with a small Streamlit UI so students can run symbolic FEM workflows interactively (without changing the core symbolic package).

### Startup command

```bash
# from repository root
uv pip install streamlit
streamlit run apps/streamlit/app.py
```

If you do not use `uv`, a plain `pip install streamlit` is also fine.

### Expected environment

- Python `3.11` (same as the course environment)
- The package installed in editable mode:
  - `uv pip install -e ".[lectures,dev]"`
- Streamlit installed in the same environment
- Launch from repo root so imports like `symbolic_fem_workbench.workflow` resolve consistently

### Common troubleshooting

- **`ModuleNotFoundError: symbolic_fem_workbench`**
  - Reinstall editable package: `uv pip install -e ".[lectures,dev]"`
  - Confirm you are running Streamlit from the project root.
- **`streamlit: command not found`**
  - Install Streamlit in the active environment: `uv pip install streamlit`.
- **Port already in use (default 8501)**
  - Run: `streamlit run apps/streamlit/app.py --server.port 8502`.
- **Slow first run / stale cache**
  - Clear Streamlit cache: `streamlit cache clear`, then restart.

### Architecture notes

- **Core symbolic logic stays in `src/symbolic_fem_workbench/`**
  - The UI should only call public helpers such as:
    - `build_bar_1d_local_problem()`
    - `build_poisson_triangle_p1_local_problem()`
    - `build_elasticity_triangle_p1_2d()`
    - `build_elasticity_tetra_p1_3d()`
- **UI adapters live in Streamlit app modules**
  - Keep adapter code in `apps/streamlit/` (for example, `app.py`, `pages/*.py`, `adapters/*.py`).
  - Adapters transform UI inputs (dropdowns, sliders, checkboxes) into arguments for `workflow.py` functions, then format outputs for display.
- **Presets map to course examples**
  - A preset should correspond to one teaching workflow/example file:
    - `bar_1d` → `examples/bar_1d_workflow.py`
    - `triangle_poisson_p1` → `examples/triangle_p1_poisson_workflow.py`
    - `manual_assembly_square_4tri` → `examples/manual_assembly_square_4tri.py`
  - Keep preset definitions declarative (e.g., one dict/list) so UI pages remain thin.
- **Adding new PDE workflows**
  1. Add/extend a pure symbolic builder in `src/symbolic_fem_workbench/workflow.py`.
  2. Add a scriptable reference usage in `examples/`.
  3. Add a Streamlit adapter/preset entry that calls the new builder.
  4. Display returned symbolic objects using existing formatting utilities (no symbolic math in UI layer).

### Contributor checklist (UI modules)

When adding a new Streamlit UI module:

- [ ] Create/update UI files under `apps/streamlit/` only.
- [ ] Use existing public API from `symbolic_fem_workbench` (avoid reaching into private internals).
- [ ] Do **not** move symbolic derivations, transformations, or extraction logic into UI code.
- [ ] Add/update a preset-to-example mapping entry.
- [ ] Add/extend tests for any new adapter logic (pure functions preferred).
- [ ] Keep outputs reproducible with a matching `examples/*.py` script.

