# Streamlit App Plan for `src/symbolic_fem_workbench`

## Goal

Create a Streamlit app around the `src/symbolic_fem_workbench` package that exposes the symbolic FEM workflow as an interactive, guided workbench without duplicating core FEM logic in the UI layer.

## 1. Create App Skeleton and Dependency Wiring

1. Create `apps/streamlit_app.py` or `streamlit_app.py` at the repository root as the runnable entrypoint.
2. Add Streamlit to project dependencies in `pyproject.toml`, or document an optional install path if the app should remain an extra.
3. In app bootstrap, import from `src.symbolic_fem_workbench` modules used in demos: `workflow`, `forms`, `fe_spaces`, `reference`, `assembly`, `elasticity`, and `viz`.
4. Add page config with `st.set_page_config`.
5. Add a top-level navigation model using `st.sidebar.radio`, tabs, or Streamlit multipage routing.
6. Keep app code separate from core library logic. The app should wrap and orchestrate existing APIs, not reimplement FEM behavior.

## 2. Build Guided Workflow Tabs

Implement workflow-first UI sections that match the symbolic FEM pipeline:

1. Problem Setup
2. FE Space
3. Form Definition
4. Assembly
5. Results and Export

For each step, map UI inputs to existing module APIs:

1. Setup: element type, dimension, and PDE flavor.
2. FE Space: basis/order using `fe_spaces` and `reference`.
3. Form: bilinear and linear form definitions via `forms` and `elasticity`.
4. Assembly: symbolic assembly via `assembly` and `workflow`.

Persist intermediate objects in `st.session_state` so users can move between sections without losing work or recomputing unnecessarily.

Each step should show both a human-readable summary and raw symbolic expressions for inspection.

## 3. Add Runnable Presets from Existing Examples

Create a "Load Preset" section that mirrors repository examples:

1. `examples/bar_1d_workflow.py`
2. `examples/triangle_p1_poisson_workflow.py`
3. `examples/manual_assembly_square_4tri.py`

Implementation notes:

1. Add preset loader functions in a small adapter module such as `apps/presets.py`.
2. Each preset should return a normalized parameter dictionary for the UI.
3. On preset load, populate `st.session_state` values and render outputs immediately.
4. Add "Reset to default" and "Clone preset to editable config" actions.

## 4. Add Visualization and Symbolic Inspection Panels

Add expandable panels for:

1. Symbolic form text using `sympy.pretty` or plain `repr`.
2. Assembled matrix and vector previews.
3. Reference or mesh visual hints using `viz` where available.

Use main-panel tabs such as:

1. Math View
2. Table View
3. Codegen/Export

The app should gracefully fall back to textual representations if a visualization backend is unavailable.

Include numeric truncation toggles and full-precision downloads for reproducibility.

## 5. Implement Export and Reproducibility Features

Add download buttons for:

1. JSON or YAML config of chosen parameters.
2. Symbolic outputs as text or Markdown.
3. Optional generated backend snippets via the `codegen` module if present.

Define a deterministic config schema in `apps/schema.py` or a similar adapter module so saved sessions can be reloaded.

Add a "Share run metadata" section that includes:

1. Package version.
2. Timestamp.
3. Selected preset.
4. Relevant app/runtime details.

Ensure all exports are generated from in-memory state rather than ad-hoc recomputation.

## 6. Add Validation, Error Handling, and UX Safety Rails

1. Validate user inputs before calling symbolic routines, including dimension/order compatibility, required coefficients, and boundary selections.
2. Wrap each compute step in exception handling with user-readable errors and expandable traceback details.
3. Add a lightweight "sanity checks" panel using the existing validation-related module, such as `validate.py`, where applicable.
4. Disable downstream buttons until prerequisite artifacts exist in `st.session_state`.

## 7. Document Running and Extension Workflow

Add a `docs/` section or README subsection titled "Interactive Workbench (Streamlit)".

Document:

1. Startup command.
2. Expected environment.
3. Common troubleshooting.
4. Architecture notes explaining where UI adapters live.
5. How presets map to examples.
6. How to add new PDE workflows.
7. A short contributor checklist for adding new UI modules without modifying core symbolic logic.

## Suggested MVP Order

1. App skeleton and dependency wiring.
2. One end-to-end preset, preferably `triangle_p1_poisson_workflow`.
3. Session-state pipeline tabs.
4. Export and documentation.
5. Remaining presets and advanced visualization.

## Initial Checks From Planning Pass

1. `pwd; rg --files | head -n 80`
2. `rg --files -g 'AGENTS.md'`
3. `rg --files | rg -n "AGENTS.md|symbolic_work_bench|Src"`

Notes:

1. No `AGENTS.md` was found in repository scope during the planning pass.
2. The referenced path/name appeared to differ from the repository package name. The repository uses `src/symbolic_fem_workbench`.
