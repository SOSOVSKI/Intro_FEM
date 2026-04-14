# Teaching Note: Manual 2D Assembly on a Tiny Triangle Mesh

This note adds the missing bridge between local element algebra and global systems.
The package now includes a worked example:

- `examples/manual_assembly_square_4tri.py`

## Why this example matters

Students often understand the weak form and can follow the derivation of a local
triangle stiffness matrix, but they still treat global assembly as mysterious bookkeeping.
This example isolates exactly that step.

The example uses:

- a unit square,
- split into four triangles around a center node,
- scalar Poisson problem,
- constant source term,
- homogeneous Dirichlet data on the outer boundary.

This setup is small enough to inspect by hand, but large enough to include:

- repeated local-to-global contributions,
- a nontrivial interior degree of freedom,
- a real reduced solve after Dirichlet elimination.

## What students should do with it

1. Inspect each local `K_local` and `F_local`.
2. Track how each element connectivity tuple maps local rows/columns into the global system.
3. Verify that the center node receives contributions from all four elements.
4. Print the assembled `K_global` and `F_global`.
5. Apply Dirichlet reduction explicitly.
6. Solve the one-degree-of-freedom reduced system.
7. Compare the center-node value with the expected result for this tiny mesh.

## What this example should teach

- Local matrices do not become global matrices automatically.
- Assembly is topological, not analytical.
- Boundary conditions modify the algebra after assembly, not before.
- A symbolic local derivation can feed directly into a numerical element loop.

## Recommended teaching use

Use this example immediately after the first 2D P1 triangle derivation.
A good sequence is:

1. Derive the weak form in 2D.
2. Derive the local P1 triangle matrix symbolically.
3. Evaluate the local tensor numerically on a specific element.
4. Assemble the four-triangle square mesh by hand in Python.
5. Only then show students how larger FEM frameworks automate the same process.
