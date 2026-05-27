from __future__ import annotations

import sympy as sp

from symbolic_fem_workbench.apps.schema import (
    SCHEMA_VERSION,
    RunConfig,
    RunState,
    export_codegen_snippets,
    export_config_json,
    export_config_yaml,
    export_run_metadata_markdown,
    export_symbolic_markdown,
    export_symbolic_text,
    session_to_dict,
)


def _state() -> RunState:
    return RunState(
        config=RunConfig(
            preset="triangle_p1",
            parameters={"E": 210e9, "nu": 0.3},
            symbolic_targets=["Ke", "fe"],
            include_codegen=True,
        ),
        symbolic_outputs={"fe": sp.Matrix([1, 2]), "Ke": sp.Matrix([[1, 0], [0, 1]])},
        codegen_outputs={"python": "def assemble():\n    pass"},
    )


def test_session_schema_deterministic_version_and_preset():
    state = _state()
    data = session_to_dict(state)
    assert data["schema"] == SCHEMA_VERSION
    assert data["config"]["preset"] == "triangle_p1"
    assert data["metadata"]["selected_preset"] == "triangle_p1"


def test_json_and_yaml_exports_include_schema_and_parameters():
    state = _state()
    json_text = export_config_json(state)
    yaml_text = export_config_yaml(state)

    assert '"schema": "1.0.0"' in json_text
    assert '"E": 210000000000.0' in json_text
    assert "schema: 1.0.0" in yaml_text
    assert "nu: 0.3" in yaml_text


def test_symbolic_and_codegen_exports_from_in_memory_state():
    state = _state()
    md = export_symbolic_markdown(state)
    txt = export_symbolic_text(state)
    codegen = export_codegen_snippets(state)

    assert "## Ke" in md and "## fe" in md
    assert "[Ke]" in txt and "[fe]" in txt
    assert "def assemble()" in codegen


def test_metadata_section_and_no_codegen_case():
    state = _state()
    state.codegen_outputs = None

    metadata_md = export_run_metadata_markdown(state)
    codegen = export_codegen_snippets(state)

    assert "## Share run metadata" in metadata_md
    assert "package_version" in metadata_md
    assert "No codegen outputs available" in codegen
