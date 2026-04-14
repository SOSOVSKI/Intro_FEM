# Course integration notes for the symbolic FEM workbench

## Why this fits the existing course

The current syllabus already emphasizes two pillars:

1. theory and mathematics in the lectures,
2. student-written FEM code in the discussions and homework.

That is a strong fit for the symbolic workbench. The tool should be introduced as a bridge between those two pillars, not as a replacement for either one.

## Recommended amendments to the course flow

### Existing topics that map directly to the workbench

- weak form derivation,
- test functions and Galerkin method,
- shape functions,
- 1D finite elements,
- numerical integration,
- natural and essential boundary conditions,
- 2D weak form and divergence theorem,
- node bookkeeping and implementation.

### Suggested additions

#### Add a symbolic workflow thread in weeks 2 to 10

- Weeks 2 to 3: weighted residual, divergence form, and integration by parts in symbolic form.
- Weeks 3 to 4: bilinear and linear form splitting.
- Weeks 4 to 5: local trial and test substitutions with shape functions in 1D.
- Weeks 5 to 6: element stiffness and load extraction in 1D.
- Weeks 8 to 9: 2D reference triangle, affine mapping, and gradient pullback.
- Weeks 9 to 10: local P1 triangle stiffness/load extraction and comparison with hand derivation.
- Optional later step: export finalized local formulas to I❤️LA text and NumPy.

#### Keep global assembly fully manual

Students should still write:

- connectivity maps,
- assembly loops,
- boundary condition enforcement,
- solve and post-processing.

That preserves the educational value of the implementation sessions.

## A practical teaching sequence

### Module A: strong to weak

Students start from divergence form, not raw second derivatives.

### Module B: local discretization in 1D

Students substitute `u_h = N d` and `v_h = N w` and inspect the resulting algebra.

### Module C: local tensors in 1D

Students extract `K_e` and `f_e` symbolically and verify them against hand derivations.

### Module D: local tensors in 2D

Students build the P1 triangle on the reference element, apply the affine map, and derive the physical gradients using `J^{-T}`.

### Module E: code generation

Students manually transcribe one small result into I❤️LA once, then later use a printer helper.

### Module F: topology and assembly

Students assemble global matrices in plain Python.

## Assessment ideas

- homework on symbolic derivation steps,
- homework on 1D element extraction and verification,
- homework on 2D triangle local tensor extraction,
- homework on assembly using the generated local tensors,
- one comparison exercise that maps the same pipeline to FEniCS/UFL notation.
