from __future__ import annotations

import streamlit as st
import sympy as sp

from symbolic_fem_workbench import assembly, fe_spaces, forms, workflow
from symbolic_fem_workbench.elasticity import plane_stress_D, plane_strain_D
from symbolic_fem_workbench.reference import ReferenceIntervalP1, ReferenceTriangleP1, ReferenceTetrahedronP1


STEPS = ["Problem Setup", "FE Space", "Form Definition", "Assembly", "Results/Export"]


def _init_state() -> None:
    st.session_state.setdefault("setup", {})
    st.session_state.setdefault("fe_space", {})
    st.session_state.setdefault("form", {})
    st.session_state.setdefault("assembly", {})


def _expr_block(title: str, expr: object) -> None:
    st.subheader(title)
    st.write("**Readable**")
    st.latex(sp.latex(expr) if isinstance(expr, (sp.Expr, sp.MatrixBase)) else str(expr))
    st.write("**Raw symbolic**")
    st.code(str(expr), language="python")


def _setup_step() -> None:
    st.header("1) Problem Setup")
    element_type = st.selectbox("Element type", ["Interval P1", "Triangle P1", "Tetrahedron P1"])
    dimension = st.selectbox("Dimension", [1, 2, 3], index=1)
    pde_flavor = st.selectbox("PDE flavor", ["Poisson", "Elasticity"])

    st.session_state.setup = {
        "element_type": element_type,
        "dimension": dimension,
        "pde_flavor": pde_flavor,
    }
    _expr_block("Setup summary", st.session_state.setup)


def _fe_space_step() -> None:
    st.header("2) FE Space")
    basis = st.selectbox("Basis", ["Lagrange"])
    order = st.selectbox("Order", [1])

    xi, eta, zeta = sp.symbols("xi eta zeta", real=True)
    element_type = st.session_state.setup.get("element_type", "Triangle P1")

    if element_type == "Interval P1":
        ref = ReferenceIntervalP1(xi=xi)
    elif element_type == "Triangle P1":
        ref = ReferenceTriangleP1(xi=xi, eta=eta)
    else:
        ref = ReferenceTetrahedronP1(xi=xi, eta=eta, zeta=zeta)

    st.session_state.fe_space = {
        "basis": basis,
        "order": order,
        "reference": ref,
        "shape_functions": ref.shape_functions,
        "shape_gradients_reference": ref.shape_gradients_reference,
        "local_trial_api": fe_spaces.local_trial_expansion,
        "local_test_api": fe_spaces.local_test_expansion,
    }

    _expr_block("FE space summary", {
        "basis": basis,
        "order": order,
        "shape_functions": ref.shape_functions,
    })


def _form_step() -> None:
    st.header("3) Form Definition")
    pde_flavor = st.session_state.setup.get("pde_flavor", "Poisson")
    dimension = st.session_state.setup.get("dimension", 2)

    if pde_flavor == "Poisson":
        result = workflow.build_poisson_triangle_p1_local_problem()
        bilinear = result["bilinear_integrand_reference"]
        linear = result["linear_integrand_reference"]
        form_obj = forms.WeakForm(bilinear=bilinear, linear=linear, trial=result["u"], test=result["v"])
    else:
        E, nu = sp.symbols("E nu", positive=True)
        formulation = st.selectbox("Elasticity formulation", ["plane_stress", "plane_strain"])
        D = plane_stress_D(E, nu) if formulation == "plane_stress" else plane_strain_D(E, nu)
        bilinear = sp.Symbol("B")
        linear = sp.Symbol("f")
        form_obj = {
            "constitutive_D": D,
            "bilinear_placeholder": bilinear,
            "linear_placeholder": linear,
            "module": "elasticity",
        }

    st.session_state.form = {
        "pde_flavor": pde_flavor,
        "dimension": dimension,
        "form": form_obj,
        "forms_module": forms,
    }

    _expr_block("Form summary", form_obj)


def _assembly_step() -> None:
    st.header("4) Assembly")
    pde_flavor = st.session_state.setup.get("pde_flavor", "Poisson")

    if pde_flavor == "Poisson":
        assembled = workflow.build_poisson_triangle_p1_local_problem()
        Ke = assembled["Ke"]
        fe = assembled["fe"]
    else:
        assembled = workflow.build_elasticity_triangle_p1_2d()
        Ke = assembled["Ke"]
        fe = sp.zeros(Ke.shape[0], 1)

    st.session_state.assembly = {
        "assembled": assembled,
        "Ke": Ke,
        "fe": fe,
        "assembly_module": assembly,
        "workflow_module": workflow,
    }

    _expr_block("Assembly summary", {"Ke": Ke, "fe": fe})


def _results_step() -> None:
    st.header("5) Results / Export")
    st.write("Session data is preserved in `st.session_state` across workflow steps.")

    for key in ["setup", "fe_space", "form", "assembly"]:
        st.markdown(f"### {key}")
        data = st.session_state.get(key, {})
        st.write("**Readable**")
        st.json({k: str(v) for k, v in data.items()})
        st.write("**Raw symbolic**")
        st.code(repr(data), language="python")


st.set_page_config(page_title="Symbolic FEM Workflow", layout="wide")
st.title("Symbolic FEM Workflow Explorer")
_init_state()
selected_step = st.sidebar.radio("Workflow", STEPS)

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
