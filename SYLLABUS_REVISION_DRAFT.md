# Draft syllabus amendment: symbolic FEM workbench thread

This draft keeps the original course spine intact and inserts a lightweight symbolic workbench thread into the existing lecture-discussion-homework structure.

## Guiding principle

The package is used to shorten routine analytical work and reduce common symbolic mistakes. It is not presented as an automatic FEM solver and it does not replace student implementation of assembly, connectivity, boundary-condition handling, or solver setup.

## Proposed week-by-week amendment

### Week 1

**Keep as is**
- introduction,
- essence of FEM,
- standard solution procedure.

**Small addition in discussion**
- show the full stack once: strong form -> weak form -> local tensor -> assembly.

### Week 2

**Current topic**
- derivation of the weak form,
- the test function.

**Add symbolic workbench task**
- define trial and test fields symbolically,
- build the weighted residual,
- stress the use of divergence form rather than raw second derivatives.

### Week 3

**Current topic**
- Galerkin method,
- shape functions.

**Add symbolic workbench task**
- split weak expressions into bilinear and linear forms,
- introduce `u_h = N d` and `v_h = N w` in 1D.

### Week 4

**Current topic**
- development of 1D finite elements,
- numerical solution using 1D finite elements.

**Add symbolic workbench task**
- derive the local 2-node bar/Poisson element,
- extract `K_e` and `f_e` symbolically,
- compare with hand derivation.

### Week 5

**Current topic**
- numerical integration using Gauss weights,
- generic 1D FE solution.

**Add symbolic workbench task**
- compare exact symbolic integration with quadrature,
- identify what changes when exact integration is replaced by quadrature.

### Week 6

**Current topic**
- natural and essential boundary conditions.

**Add symbolic workbench task**
- inspect boundary terms explicitly,
- separate Dirichlet test-space vanishing from Neumann flux substitution.

### Week 7

**Current topic**
- global and local error,
- characteristics of the numerical solution.

**Add symbolic workbench task**
- verify consistency of local tensors against patch tests or simple manufactured solutions.

### Week 8

**Current topic**
- time dependent problems in 1D,
- beginning of 2D weak form.

**Add symbolic workbench task**
- introduce the reference triangle,
- define P1 triangle shape functions,
- derive the affine map and Jacobian.

### Week 9

**Current topic**
- 2D formulation,
- divergence theorem,
- development of 2D finite elements.

**Add symbolic workbench task**
- compute physical gradients using `J^{-T}`,
- derive the local scalar Poisson triangle stiffness matrix symbolically.

### Week 10

**Current topic**
- shape functions in 2D,
- meshing and node bookkeeping,
- enforcing BCs in 2D.

**Add symbolic workbench task**
- derive the local triangle load vector,
- hand off local kernels to plain Python assembly,
- keep node bookkeeping and assembly fully manual.

### Week 11 onward

**Keep course emphasis**
- implementation details,
- characteristics of finite elements in 2D,
- broader comparisons to higher-level FEM tools.

**Optional bridge**
- map the package workflow to UFL/FEniCS notation,
- map finalized local formulas to I❤️LA and generated NumPy.

## Suggested homework amendments

### Homework 1
- weighted residual and integration by parts for a 1D bar/Poisson problem.

### Homework 2
- symbolic extraction of the 1D element stiffness and load vector.

### Homework 3
- plain Python assembly using symbolically derived local tensors.

### Homework 4
- local P1 triangle stiffness/load derivation and assembly on a very small mesh.

## Teaching note

The symbolic thread should stay in the discussion and homework ecosystem more than in the lecture core. That preserves the theoretical emphasis of the course while giving students a rigorous and visible path from calculus to code.
