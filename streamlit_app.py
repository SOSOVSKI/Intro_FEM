"""Streamlit UI wrapper for symbolic_fem_workbench workflows.

This module keeps FEM derivations in the core library and only orchestrates
interactive state, rendering, and inspection.
"""

from __future__ import annotations

import sympy as sp
import streamlit as st

from symbolic_fem_workbench import assembly, elasticity, fe_spaces, forms, reference, viz, workflow
from symbolic_fem_workbench.elasticity import plane_stress_D, plane_strain_D
from symbolic_fem_workbench.reference import (
    ReferenceIntervalP1,
    ReferenceTetrahedronP1,
    ReferenceTriangleP1,
)


WORKFLOW_STEPS = [
    "Problem Setup",
    "FE Space",
    "Form Definition",
    "Assembly",
    "Results/Export",
]

PAGES = [
    "Guided Workflow",
    "Example Demos",
    "Module Surface",
]


def _init_state() -> None:
    st.session_state.setdefault(
        "setup",
        {
            "element_type": "Triangle P1",
            "dimension": 2,
            "pde_flavor": "Poisson",
        },
    )
    st.session_state.setdefault("fe_space", {})
    st.session_state.setdefault("form", {})
    st.session_state.setdefault("assembly", {})


def _render_expr(title: str, expr: object) -> None:
    st.subheader(title)
    st.caption("Readable")
    if isinstance(expr, (sp.Expr, sp.MatrixBase)):
        st.latex(sp.latex(expr))
    else:
        st.write(expr)
    st.caption("Raw symbolic")
    st.code(repr(expr), language="python")


def _reference_from_setup() -> object:
    xi, eta, zeta = sp.symbols("xi eta zeta", real=True)
    element_type = st.session_state.setup.get("element_type", "Triangle P1")
    if element_type == "Interval P1":
        return ReferenceIntervalP1(xi=xi)
    if element_type == "Tetrahedron P1":
        return ReferenceTetrahedronP1(xi=xi, eta=eta, zeta=zeta)
    return ReferenceTriangleP1(xi=xi, eta=eta)


def _build_current_problem() -> dict[str, object]:
    setup = st.session_state.setup
    element_type = setup.get("element_type", "Triangle P1")
    pde_flavor = setup.get("pde_flavor", "Poisson")

    if element_type == "Interval P1":
        return workflow.build_bar_1d_local_problem()
    if pde_flavor == "Elasticity":
        formulation = setup.get("elasticity_formulation", "plane_stress")
        return workflow.build_elasticity_triangle_p1_2d(formulation=formulation)
    return workflow.build_poisson_triangle_p1_local_problem()


def _setup_step() -> None:
    st.header("1. Problem Setup")
    current = st.session_state.setup
    element_type = st.selectbox(
        "Element type",
        ["Interval P1", "Triangle P1", "Tetrahedron P1"],
        index=["Interval P1", "Triangle P1", "Tetrahedron P1"].index(
            current.get("element_type", "Triangle P1")
        ),
    )
    dimension = st.selectbox(
        "Dimension",
        [1, 2, 3],
        index=[1, 2, 3].index(int(current.get("dimension", 2))),
    )
    pde_flavor = st.selectbox(
        "PDE flavor",
        ["Poisson", "Elasticity"],
        index=["Poisson", "Elasticity"].index(current.get("pde_flavor", "Poisson")),
    )
    setup = {
        "element_type": element_type,
        "dimension": dimension,
        "pde_flavor": pde_flavor,
    }
    if pde_flavor == "Elasticity":
        setup["elasticity_formulation"] = st.selectbox(
            "Elasticity formulation",
            ["plane_stress", "plane_strain"],
            index=0,
        )

    st.session_state.setup = setup
    _render_expr("Setup summary", setup)


def _fe_space_step() -> None:
    st.header("2. FE Space")
    basis = st.selectbox("Basis", ["Lagrange"])
    order = st.selectbox("Order", [1])
    ref = _reference_from_setup()
    dofs = sp.symbols(f"d0:{len(ref.shape_functions)}", real=True)
    weights = sp.symbols(f"w0:{len(ref.shape_functions)}", real=True)

    st.session_state.fe_space = {
        "basis": basis,
        "order": order,
        "reference": ref,
        "nodes": ref.nodes,
        "shape_functions": ref.shape_functions,
        "shape_gradients_reference": ref.shape_gradients_reference,
        "trial_expansion": fe_spaces.local_trial_expansion(ref.shape_functions, dofs),
        "test_expansion": fe_spaces.local_test_expansion(ref.shape_functions, weights),
    }

    _render_expr(
        "FE space summary",
        {
            "basis": basis,
            "order": order,
            "nodes": ref.nodes,
            "shape_functions": ref.shape_functions,
        },
    )
    _render_expr("Trial expansion", st.session_state.fe_space["trial_expansion"])
    _render_expr("Test expansion", st.session_state.fe_space["test_expansion"])


def _form_step() -> None:
    st.header("3. Form Definition")
    problem = _build_current_problem()
    pde_flavor = st.session_state.setup.get("pde_flavor", "Poisson")

    if st.session_state.setup.get("element_type") == "Interval P1":
        form_obj = forms.WeakForm(
            bilinear=problem["weak_bilinear"],
            linear=problem["weak_linear"],
            trial=problem["u"],
            test=problem["v"],
        )
    elif pde_flavor == "Elasticity":
        E, nu = sp.symbols("E nu", positive=True)
        formulation = st.session_state.setup.get("elasticity_formulation", "plane_stress")
        form_obj = {
            "constitutive_D": plane_stress_D(E, nu)
            if formulation == "plane_stress"
            else plane_strain_D(E, nu),
            "formulation": formulation,
            "B": problem["B"],
        }
    else:
        form_obj = forms.WeakForm(
            bilinear=problem["bilinear_integrand_reference"],
            linear=problem["linear_integrand_reference"],
            trial=problem["u"],
            test=problem["v"],
        )

    st.session_state.form = {
        "pde_flavor": pde_flavor,
        "dimension": st.session_state.setup.get("dimension", 2),
        "form": form_obj,
        "forms_module": forms.__name__,
        "elasticity_module": elasticity.__name__,
    }
    _render_expr("Form summary", form_obj)


def _assembly_step() -> None:
    st.header("4. Assembly")
    assembled = _build_current_problem()
    Ke = assembled["Ke"]
    fe = assembled.get("fe", sp.zeros(Ke.shape[0], 1))

    st.session_state.assembly = {
        "assembled": assembled,
        "Ke": Ke,
        "fe": fe,
        "assembly_module": assembly.__name__,
        "workflow_module": workflow.__name__,
    }

    _render_expr("Element stiffness Ke", Ke)
    _render_expr("Element load fe", fe)


def _results_step() -> None:
    st.header("5. Results / Export")
    for key in ["setup", "fe_space", "form", "assembly"]:
        with st.expander(key.replace("_", " ").title(), expanded=key in {"setup", "assembly"}):
            data = st.session_state.get(key, {})
            st.json({str(k): str(v) for k, v in data.items()})
            st.code(repr(data), language="python")


def page_guided_workflow() -> None:
    st.title("Symbolic FEM Workbench")
    _init_state()
    selected_step = st.sidebar.radio("Workflow", WORKFLOW_STEPS)
    if selected_step == "Problem Setup":
        _setup_step()
    elif selected_step == "FE Space":
        _fe_space_step()
    elif selected_step == "Form Definition":
        _form_step()
    elif selected_step == "Assembly":
        _assembly_step()
    else:
        _results_step()


def page_example_demos() -> None:
    st.title("Example Demos")
    demo = st.radio("Demo", ["1D Bar Workflow", "2D Triangle Workflow", "Elasticity"], horizontal=True)
    if demo == "1D Bar Workflow":
        data = workflow.build_bar_1d_local_problem()
        col1, col2 = st.columns(2)
        with col1:
            _render_expr("Weak bilinear form", data["weak_bilinear"])
            _render_expr("Weak linear form", data["weak_linear"])
        with col2:
            _render_expr("Element stiffness Ke", data["Ke"])
            _render_expr("Element load fe", data["fe"])
    elif demo == "2D Triangle Workflow":
        data = workflow.build_poisson_triangle_p1_local_problem()
        col1, col2 = st.columns(2)
        with col1:
            _render_expr("Reference bilinear integrand", data["bilinear_integrand_reference"])
            _render_expr("Reference linear integrand", data["linear_integrand_reference"])
        with col2:
            _render_expr("Ke on unit right triangle", data["Ke_unit_right_triangle"])
            _render_expr("fe on unit right triangle", data["fe_unit_right_triangle"])
    else:
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
    st.title("Module Surface")
    modules = {
        "forms": forms,
        "fe_spaces": fe_spaces,
        "reference": reference,
        "assembly": assembly,
        "elasticity": elasticity,
        "viz": viz,
        "workflow": workflow,
    }
    selected = st.selectbox("Module", list(modules.keys()))
    names = [name for name in dir(modules[selected]) if not name.startswith("_")]
    st.write(f"Public names in `{selected}`: {len(names)}")
    st.code("\n".join(names[:200]), language="text")


def main() -> None:
    st.set_page_config(
        page_title="Symbolic FEM Workbench",
        page_icon=":material/functions:",
        layout="wide",
    )
    with st.sidebar:
        st.header("Navigation")
        page = st.radio("Go to", PAGES)

    if page == "Guided Workflow":
        page_guided_workflow()
    elif page == "Example Demos":
        page_example_demos()
    else:
        page_module_surface()


if __name__ == "__main__":
    main()
