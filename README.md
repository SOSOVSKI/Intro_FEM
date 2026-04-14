# symbolic-fem-workbench

A pedagogical symbolic FEM workbench built on top of SymPy.

This package is designed for teaching finite elements without turning the derivation into a black box. It supports explicit symbolic steps from strong form to weak form, finite-dimensional substitution at the element level, extraction of local tensors, and optional export of finalized local algebra to I❤️LA-ready text.

## Design goals

- keep every analytical step visible,
- reduce sign, derivative, and indexing mistakes,
- separate calculus from topology and global assembly,
- support teaching-first workflows,
- stay narrow in scope for the first phase.

## Current scope

### 1D

- 1D linear scalar problems in divergence form,
- weighted residual construction,
- integration by parts for divergence terms,
- explicit Dirichlet and Neumann boundary handling,
- P1 interval finite elements,
- element matrix/vector extraction.

### 2D phase 1

- local P1 triangle geometry helpers,
- affine reference-to-physical mapping,
- gradient pullback using `J^{-T}`,
- exact integration on the reference triangle,
- local stiffness/load extraction for scalar Poisson-type forms.

## Package layout

```text
src/symbolic_fem_workbench/
    symbols.py
    forms.py
    transforms.py
    fe_spaces.py
    reference.py
    quadrature.py
    extract.py
    validate.py
    workflow.py
    assembly.py
    printers/
        iheartla_printer.py
    codegen/
        iheartla_backend.py
```

## Quick start

See:

- `examples/bar_1d_workflow.py`
- `examples/triangle_p1_poisson_workflow.py`
- `examples/manual_assembly_square_4tri.py`
- `examples/triangle_p1_poisson_workflow.py`

## What this package does not do

- global assembly,
- mesh generation,
- full boundary-condition management in 2D,
- nonlinear automation,
- replacement of a real FEM framework.

That is intentional. Students should still write their own assembly and topology code, because that is where most of the method actually becomes real.
