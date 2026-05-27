from __future__ import annotations

import copy

import streamlit as st

from presets import load_all_presets
from symbolic_fem_workbench.workflow import (
    build_bar_1d_local_problem,
    build_poisson_triangle_p1_local_problem,
)

DEFAULT_STATE = {
    "selected_preset": "bar_1d",
    "active_config": {},
    "editable_config": {},
}


def _ensure_state() -> None:
    for key, value in DEFAULT_STATE.items():
        st.session_state.setdefault(key, copy.deepcopy(value))


def _reset_to_default() -> None:
    for key, value in DEFAULT_STATE.items():
        st.session_state[key] = copy.deepcopy(value)


def _apply_preset(preset_key: str, presets: dict[str, object]) -> None:
    preset = presets[preset_key]
    st.session_state["selected_preset"] = preset_key
    st.session_state["active_config"] = copy.deepcopy(preset.parameters)


def _clone_preset_to_editable() -> None:
    st.session_state["editable_config"] = copy.deepcopy(st.session_state["active_config"])


def _render_output(config: dict[str, object]) -> None:
    problem_type = config.get("problem_type")

    if problem_type == "bar_1d":
        problem = build_bar_1d_local_problem()
        st.subheader("Bar 1D local tensors")
        st.write("Ke")
        st.code(str(problem["Ke"]))
        st.write("fe")
        st.code(str(problem["fe"]))
    elif problem_type == "triangle_p1_poisson":
        data = build_poisson_triangle_p1_local_problem()
        st.subheader("Triangle P1 Poisson local tensors")
        st.write("Ke")
        st.code(str(data["Ke"]))
        st.write("fe")
        st.code(str(data["fe"]))
    elif problem_type == "manual_assembly_square_4tri":
        st.subheader("Manual Assembly / Square 4 Triangles")
        st.json(config)
        st.info("Configuration mirrors examples/manual_assembly_square_4tri.py")
    else:
        st.warning("Unknown preset configuration.")


def main() -> None:
    st.title("Symbolic FEM Workbench")
    _ensure_state()
    presets = load_all_presets()

    st.header("Load Preset")
    with st.container(border=True):
        preset_key = st.selectbox(
            "Repository example preset",
            options=list(presets.keys()),
            index=list(presets.keys()).index(st.session_state["selected_preset"]),
            format_func=lambda key: presets[key].title,
        )
        st.caption(presets[preset_key].description)

        col1, col2, col3 = st.columns(3)
        with col1:
            if st.button("Load preset", use_container_width=True):
                _apply_preset(preset_key, presets)
        with col2:
            if st.button("Reset to default", use_container_width=True):
                _reset_to_default()
        with col3:
            if st.button("Clone preset to editable config", use_container_width=True):
                if not st.session_state["active_config"]:
                    _apply_preset(preset_key, presets)
                _clone_preset_to_editable()

    if not st.session_state["active_config"]:
        _apply_preset(st.session_state["selected_preset"], presets)

    st.header("Active Preset Config")
    st.json(st.session_state["active_config"])

    st.header("Rendered Output")
    _render_output(st.session_state["active_config"])

    st.header("Editable Cloned Config")
    st.json(st.session_state["editable_config"])


if __name__ == "__main__":
    main()
