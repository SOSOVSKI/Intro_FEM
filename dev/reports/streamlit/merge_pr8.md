# Merge Report: PR #8

Generated: 2026-05-27

## Confirmation

Merged PR #8 (`pr-8`) into local `main` with merge commit message `Merge PR #8 streamlit export schema`.

## What landed

- Added deterministic app schema/export helpers under `src/symbolic_fem_workbench/apps/`.
- Added config JSON/YAML exports, symbolic text/Markdown exports, codegen snippet export handling, and run metadata export helpers.
- Added `pyyaml>=6.0` to runtime dependencies.
- Added `tests/test_app_schema_exports.py`.

## Fixes Required

- Resolved `pyproject.toml` dependency conflict by retaining both `streamlit>=1.36` from PR #3 and `pyyaml>=6.0` from PR #8.
- Added a local fallback for `__version__` so source-tree imports do not fail before the package is installed.
- Wired the Streamlit `Codegen/Export` tab to construct a `RunState` from current in-memory `st.session_state` and assembled symbolic outputs.
- Added Streamlit download buttons for deterministic config JSON/YAML, symbolic Markdown/text, codegen snippets, and run metadata.

## Verification

- `python3 -m py_compile streamlit_app.py src/symbolic_fem_workbench/__init__.py src/symbolic_fem_workbench/apps/schema.py src/symbolic_fem_workbench/apps/__init__.py tests/test_app_schema_exports.py` passed.
- `python3 -m pytest -q tests/test_app_schema_exports.py` did not run to completion because package import reaches `assembly.py`, which requires missing dependency `numpy`.

## Going Forward

- Full export tests remain blocked until dependencies are installed.
- PR #7 should be used to gate downstream workflow/export actions and present user-readable compute errors.
