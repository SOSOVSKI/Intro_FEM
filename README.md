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
