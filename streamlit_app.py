"""Streamlit UI wrapper for symbolic_fem_workbench demos.

This module intentionally keeps FEM derivations in the core library.
"""

from __future__ import annotations

import streamlit as st

from src.symbolic_fem_workbench import assembly, elasticity, fe_spaces, forms, reference, viz, workflow


st.set_page_config(
    page_title="Symbolic FEM Workbench",
    page_icon="🧮",
    layout="wide",
)


PAGES = [
    "Overview",
    "1D Bar Workflow",
    "2D Triangle Workflow",
    "Elasticity",
    "Module Surface",
]


def _render_expr(title: str, expr: object) -> None:
    st.subheader(title)
    st.code(str(expr), language="text")


def page_overview() -> None:
    st.title("Symbolic FEM Workbench")
    st.write(
        "Interactive UI around the existing symbolic FEM library modules. "
        "All FEM logic stays in `src/symbolic_fem_workbench` and is imported here."
    )
    st.markdown(
        "**Imported demo modules:** `workflow`, `forms`, `fe_spaces`, `reference`, "
        "`assembly`, `elasticity`, `viz`."
    )


def page_bar_1d() -> None:
    st.title("1D Bar Local Element")
    data = workflow.build_bar_1d_local_problem()
    col1, col2 = st.columns(2)
    with col1:
        _render_expr("Weak bilinear form", data["weak_bilinear"])
        _render_expr("Weak linear form", data["weak_linear"])
    with col2:
        _render_expr("Element stiffness Ke", data["Ke"])
        _render_expr("Element load fe", data["fe"])


def page_triangle() -> None:
    st.title("2D Poisson on P1 Triangle")
    data = workflow.build_poisson_triangle_p1_local_problem()
    col1, col2 = st.columns(2)
    with col1:
        _render_expr("Reference bilinear integrand", data["bilinear_integrand_reference"])
        _render_expr("Reference linear integrand", data["linear_integrand_reference"])
    with col2:
        _render_expr("Ke on unit right triangle", data["Ke_unit_right_triangle"])
        _render_expr("fe on unit right triangle", data["fe_unit_right_triangle"])


def page_elasticity() -> None:
    st.title("2D Elasticity (P1 Triangle)")
    formulation = st.selectbox("Constitutive model", ["plane_stress", "plane_strain"], index=0)
    data = workflow.build_elasticity_triangle_p1_2d(formulation=formulation)
    col1, col2 = st.columns(2)
    with col1:
        _render_expr("Constitutive matrix D", data["D"])
        _render_expr("B-matrix", data["B"])
    with col2:
        _render_expr("Ke (symbolic)", data["Ke"])
        _render_expr("Ke on unit right triangle", data["Ke_unit_right_triangle"])


def page_module_surface() -> None:
    st.title("Module Surface / Quick Introspection")
    modules = {
        "forms": forms,
        "fe_spaces": fe_spaces,
        "reference": reference,
        "assembly": assembly,
        "elasticity": elasticity,
        "viz": viz,
    }
    selected = st.selectbox("Module", list(modules.keys()))
    names = [n for n in dir(modules[selected]) if not n.startswith("_")]
    st.write(f"Public names in `{selected}`: {len(names)}")
    st.code("\n".join(names[:200]), language="text")


def main() -> None:
    with st.sidebar:
        st.header("Navigation")
        page = st.radio("Go to", PAGES)

    if page == "Overview":
        page_overview()
    elif page == "1D Bar Workflow":
        page_bar_1d()
    elif page == "2D Triangle Workflow":
        page_triangle()
    elif page == "Elasticity":
        page_elasticity()
    elif page == "Module Surface":
        page_module_surface()


if __name__ == "__main__":
    main()
