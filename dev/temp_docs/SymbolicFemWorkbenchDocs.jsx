import React, { useMemo, useState } from "react";

const tabs = [
  "Overview",
  "Workflow",
  "Modules",
  "Examples",
  "Teaching",
  "Roadmap",
];

const workflowSteps = [
  {
    title: "1. Define the problem",
    body:
      "Create domains and symbolic fields with make_field_1d or make_field_2d. Keep the unknown, the test field, and the coefficients separate so the mathematics stays legible.",
  },
  {
    title: "2. Build the weak form explicitly",
    body:
      "Use divergence-form helpers, boundary-term handling, and split_linear_weak_form to move from the strong statement to a bilinear form and a linear form without pretending the package should make formulation choices for the student.",
  },
  {
    title: "3. Insert finite-dimensional approximations",
    body:
      "Construct local trial and test expansions with explicit degree-of-freedom vectors. Substitute u_h and v_h at the element level rather than hiding the approximation step.",
  },
  {
    title: "4. Extract local tensors",
    body:
      "Use extract_coefficient_matrix and extract_coefficient_vector to recover K_e and f_e from the symbolic weak form after substitution and integration.",
  },
  {
    title: "5. Assemble globally in plain Python",
    body:
      "Keep assembly, connectivity, and Dirichlet elimination visible. The package helps with calculus and local algebra, not with replacing the parts students need to understand.",
  },
  {
    title: "6. Export final local algebra if desired",
    body:
      "Use the narrow I❤️LA printer bridge only after the local algebra is settled. That way code generation is the end of the reasoning chain rather than a substitute for it.",
  },
];

const modules = [
  {
    name: "symbols.py",
    purpose: "Defines teaching-first symbolic building blocks for 1D and 2D domains and fields.",
    exports: ["Domain1D", "Domain2D", "make_field_1d", "make_field_2d"],
    notes:
      "Use this module first. It keeps the distinction between the physical setting and the later algebraic manipulations from dissolving into generic SymPy clutter.",
  },
  {
    name: "forms.py",
    purpose: "Structured containers for domain integrals, boundary contributions, weighted residuals, and split weak forms.",
    exports: ["DomainIntegral", "BoundaryContribution", "WeightedResidual", "WeakForm"],
    notes:
      "This is where the package refuses to flatten everything into an anonymous expression blob. A rare act of discipline in scientific code.",
  },
  {
    name: "transforms.py",
    purpose: "Symbolic transformations from strong form to weak form, FE substitution, 2D gradients, and linearization.",
    exports: [
      "weighted_residual",
      "integrate_divergence_1d",
      "drop_dirichlet_boundary",
      "apply_neumann_flux",
      "split_linear_weak_form",
      "grad_2d",
      "pullback_gradient_2d",
      "substitute_field",
      "substitute_fe",
      "gateaux_derivative",
    ],
    notes:
      "This is the analytical engine. It should remain explicit and composable, not drift into an automation ritual that picks formulations behind the user's back.",
  },
  {
    name: "fe_spaces.py",
    purpose: "Finite element spaces and local trial/test expansions for the current phase-1 elements.",
    exports: ["LinearElement1D", "local_trial_expansion", "local_test_expansion"],
    notes:
      "This module should grow carefully. It is where convenience can easily mutate into hidden assumptions.",
  },
  {
    name: "reference.py",
    purpose: "Reference-element definitions and affine mapping support for the current teaching elements in 2D and 3D.",
    exports: ["ReferenceTriangleP1", "ReferenceQuadrilateralQ1", "ReferenceQuadrilateralQ2", "ReferenceHexahedronQ1", "ReferenceTetrahedronP1", "ReferenceTetrahedronP2", "AffineTriangleMap2D"],
    notes:
      "Use the course-facing names in class: triangle P1, quadrilateral degree 1 (Q1), quadrilateral degree 2 (Q2), Hex8, Tet4, and Tet10. The docs keep the naming simple on purpose.",
  },
  {
    name: "quadrature.py",
    purpose: "Reference-triangle integration and quadrature helpers for local 2D forms.",
    exports: ["integrate_reference_triangle_exact", "triangle_one_point_rule"],
    notes:
      "Phase 1 keeps this modest. The goal is transparency, not a heroic quadrature library nobody asked for.",
  },
  {
    name: "extract.py",
    purpose: "Extraction of coefficient matrices and vectors from symbolic bilinear and linear expressions.",
    exports: ["extract_coefficient_matrix", "extract_coefficient_vector"],
    notes:
      "This is the bridge from symbolic weak forms to the local algebra students actually recognize as K_e and f_e.",
  },
  {
    name: "validate.py",
    purpose: "Reserved for correctness checks and workflow guards.",
    exports: ["Validation helpers to expand over time"],
    notes:
      "Even when small, validation matters. Scientific code loves letting bad assumptions survive until the final minus sign explodes during grading week.",
  },
  {
    name: "workflow.py",
    purpose: "Course-ready workflow constructors that assemble the symbolic steps into reusable teaching examples.",
    exports: ["build_bar_1d_local_problem", "build_poisson_triangle_p1_local_problem"],
    notes:
      "Use these as reference implementations, not as a replacement for students doing the mathematics themselves.",
  },
  {
    name: "assembly.py",
    purpose: "Explicit dense assembly helpers and Dirichlet reduction utilities for manual teaching examples.",
    exports: [
      "assemble_dense_matrix",
      "assemble_dense_vector",
      "apply_dirichlet_by_reduction",
      "expand_reduced_solution",
    ],
    notes:
      "This keeps the topology layer visible. Local tensors do not become a global system by metaphysical uplift.",
  },
  {
    name: "printers/iheartla_printer.py",
    purpose: "Narrow printer layer for final local expressions that are ready to be exported toward I❤️LA syntax.",
    exports: ["I❤️LA-oriented formatting helpers"],
    notes:
      "Keep this narrow. Code generation belongs at the end of the chain, after the mathematics is settled.",
  },
];

const examples = [
  {
    title: "bar_1d_workflow.py",
    focus: "Strong form to weak form to local bar-element tensors.",
    teaches: [
      "weighted residuals",
      "integration by parts in divergence form",
      "Dirichlet and Neumann treatment",
      "P1 interval discretization",
      "symbolic extraction of K_e and f_e",
    ],
  },
  {
    title: "triangle_p1_poisson_workflow.py",
    focus: "Reference-triangle formulation for a scalar 2D Poisson-type element.",
    teaches: [
      "reference coordinates",
      "affine map and Jacobian",
      "gradient pullback via J^{-T}",
      "local 2D bilinear and linear forms",
    ],
  },
  {
    title: "manual_assembly_square_4tri.py",
    focus: "Local-to-global assembly on a tiny 2D mesh with boundary elimination and solve.",
    teaches: [
      "connectivity",
      "global accumulation",
      "Dirichlet reduction",
      "recovering the full solution vector",
    ],
  },
];

const teachingTopics = [
  {
    week: "Weeks 1-2",
    topic: "1D weak form as a symbolic workflow",
    packageUse:
      "Students define the strong form in divergence form, build the weighted residual, integrate by parts, and handle boundary contributions explicitly.",
  },
  {
    week: "Weeks 3-4",
    topic: "Element-level discretization in 1D",
    packageUse:
      "Students substitute u_h and v_h with independent coefficient vectors, integrate symbolically, and extract K_e and f_e before writing any assembly code.",
  },
  {
    week: "Weeks 5-6",
    topic: "Assembly as a separate intellectual task",
    packageUse:
      "Students use manual assembly helpers sparingly or replicate them by hand. The point is to isolate topology from calculus, not to hide topology.",
  },
  {
    week: "Weeks 7-9",
    topic: "2D formulation on triangles",
    packageUse:
      "The package supports reference-triangle shape functions, affine geometry, gradient pullback, and local 2D tensor extraction for scalar Poisson-type examples.",
  },
  {
    week: "Later bridge",
    topic: "Connecting to FEniCS/UFL or I❤️LA",
    packageUse:
      "Students compare their own symbolic and local-tensor pipeline against higher-level tools, which stops those tools from seeming like black magic dropped from the sky.",
  },
];

const roadmap = [
  "Add additional teaching cases carefully, using the course labels Quadrilateral degree 1 (Q1) and Quadrilateral degree 2 (Q2) without turning the docs into an element-taxonomy detour.",
  "Expand validate.py with checks for linearity, test/trial dependence, and common gradient/Jacobian mistakes.",
  "Add a restricted SymPy to I❤️LA printer for small matrices, vectors, and scalar expressions to reduce manual transcription errors.",
  "Introduce a nonlinear layer through residual forms and Gateaux derivatives rather than by inventing a separate workflow universe.",
  "Add notebook-based teaching assets that mirror the lecture order without turning the package into a notebook graveyard.",
];

function SectionCard({ title, children, subtitle }) {
  return (
    <div className="rounded-2xl border border-slate-200 bg-white p-5 shadow-sm">
      <div className="mb-3">
        <h3 className="text-lg font-semibold text-slate-900">{title}</h3>
        {subtitle ? <p className="mt-1 text-sm text-slate-500">{subtitle}</p> : null}
      </div>
      <div className="text-sm leading-6 text-slate-700">{children}</div>
    </div>
  );
}

function Pill({ children }) {
  return (
    <span className="inline-flex rounded-full border border-slate-300 bg-slate-50 px-3 py-1 text-xs font-medium text-slate-700">
      {children}
    </span>
  );
}

export default function SymbolicFemWorkbenchDocs() {
  const [activeTab, setActiveTab] = useState("Overview");

  const activeContent = useMemo(() => {
    switch (activeTab) {
      case "Overview":
        return (
          <div className="grid gap-5 lg:grid-cols-2">
            <SectionCard title="What the package is">
              <p>
                <strong>symbolic-fem-workbench</strong> is a teaching-first symbolic FEM package built on top of SymPy.
                It is meant to make the derivation pipeline explicit: strong form, weighted residual, weak form,
                finite-dimensional substitution, local tensor extraction, and finally manual assembly.
              </p>
              <p className="mt-3">
                It is intentionally narrow. It does not try to be a full FEM framework, because students need to see
                where the method lives instead of outsourcing the whole intellectual problem to a convenience layer.
              </p>
            </SectionCard>

            <SectionCard title="Design rules">
              <ul className="list-disc space-y-2 pl-5">
                <li>Keep analytical steps visible.</li>
                <li>Keep local algebra separate from global topology.</li>
                <li>Use structured weak-form objects instead of anonymous expression soup.</li>
                <li>Use simple course-facing element names like Quadrilateral degree 1 (Q1) and degree 2 (Q2) instead of sending students into classification side quests.</li>
                <li>Prefer explicit transformations over hidden automation.</li>
                <li>Use code generation only after the mathematics is settled.</li>
              </ul>
            </SectionCard>

            <SectionCard title="Current scope" subtitle="Phase 1 package coverage">
              <div className="space-y-3">
                <div>
                  <div className="font-medium text-slate-900">1D</div>
                  <p>Linear scalar problems in divergence form, boundary handling, P1 interval elements, and local tensor extraction.</p>
                </div>
                <div>
                  <div className="font-medium text-slate-900">2D</div>
                  <p>P1 triangles, affine reference-to-physical mapping, gradient pullback, reference integration, and local Poisson-type tensors.</p>
                </div>
                <div>
                  <div className="font-medium text-slate-900">Assembly</div>
                  <p>Minimal dense assembly helpers and explicit Dirichlet reduction for teaching examples.</p>
                </div>
              </div>
            </SectionCard>

            <SectionCard title="Public imports" subtitle="Main entry points exposed from __init__.py">
              <div className="flex flex-wrap gap-2">
                {[
                  "Domain1D",
                  "Domain2D",
                  "make_field_1d",
                  "make_field_2d",
                  "WeakForm",
                  "integrate_divergence_1d",
                  "split_linear_weak_form",
                  "LinearElement1D",
                  "ReferenceTriangleP1",
                  "AffineTriangleMap2D",
                  "extract_coefficient_matrix",
                  "extract_coefficient_vector",
                  "build_bar_1d_local_problem",
                  "build_poisson_triangle_p1_local_problem",
                  "assemble_dense_matrix",
                  "apply_dirichlet_by_reduction",
                ].map((item) => (
                  <Pill key={item}>{item}</Pill>
                ))}
              </div>
            </SectionCard>
          </div>
        );

      case "Workflow":
        return (
          <div className="space-y-4">
            {workflowSteps.map((step) => (
              <SectionCard key={step.title} title={step.title}>
                <p>{step.body}</p>
              </SectionCard>
            ))}
          </div>
        );

      case "Modules":
        return (
          <div className="grid gap-5 lg:grid-cols-2">
            {modules.map((module) => (
              <SectionCard key={module.name} title={module.name} subtitle={module.purpose}>
                <div>
                  <div className="mb-2 text-xs font-semibold uppercase tracking-wide text-slate-500">Key exports</div>
                  <div className="mb-3 flex flex-wrap gap-2">
                    {module.exports.map((item) => (
                      <Pill key={item}>{item}</Pill>
                    ))}
                  </div>
                </div>
                <p>{module.notes}</p>
              </SectionCard>
            ))}
          </div>
        );

      case "Examples":
        return (
          <div className="grid gap-5 lg:grid-cols-3">
            {examples.map((example) => (
              <SectionCard key={example.title} title={example.title} subtitle={example.focus}>
                <div className="mb-2 text-xs font-semibold uppercase tracking-wide text-slate-500">What students see</div>
                <ul className="list-disc space-y-2 pl-5">
                  {example.teaches.map((item) => (
                    <li key={item}>{item}</li>
                  ))}
                </ul>
              </SectionCard>
            ))}
          </div>
        );

      case "Teaching":
        return (
          <div className="space-y-4">
            {teachingTopics.map((topic) => (
              <SectionCard key={topic.week} title={`${topic.week}: ${topic.topic}`}>
                <p>{topic.packageUse}</p>
              </SectionCard>
            ))}
          </div>
        );

      case "Roadmap":
        return (
          <SectionCard title="Suggested next steps" subtitle="Build carefully, not theatrically.">
            <ol className="list-decimal space-y-3 pl-5">
              {roadmap.map((item) => (
                <li key={item}>{item}</li>
              ))}
            </ol>
          </SectionCard>
        );

      default:
        return null;
    }
  }, [activeTab]);

  return (
    <div className="min-h-screen bg-slate-50 text-slate-900">
      <div className="mx-auto max-w-7xl px-4 py-8 sm:px-6 lg:px-8">
        <header className="mb-8 rounded-3xl border border-slate-200 bg-white p-6 shadow-sm">
          <div className="flex flex-col gap-4 lg:flex-row lg:items-end lg:justify-between">
            <div>
              <p className="text-sm font-semibold uppercase tracking-[0.2em] text-slate-500">
                Package documentation
              </p>
              <h1 className="mt-2 text-3xl font-bold tracking-tight text-slate-950 sm:text-4xl">
                symbolic-fem-workbench
              </h1>
              <p className="mt-3 max-w-3xl text-sm leading-6 text-slate-600 sm:text-base">
                A tabbed documentation page for the teaching-first symbolic FEM package. It is built to explain the
                package structure, the intended workflow, the module boundaries, and the teaching strategy without
                burying everything under generic API sludge.
              </p>
            </div>
            <div className="flex flex-wrap gap-2">
              <Pill>SymPy-first</Pill>
              <Pill>Teaching-first</Pill>
              <Pill>Element-level</Pill>
              <Pill>Manual assembly visible</Pill>
            </div>
          </div>
        </header>

        <nav className="mb-6 overflow-x-auto rounded-2xl border border-slate-200 bg-white p-2 shadow-sm">
          <div className="flex min-w-max gap-2">
            {tabs.map((tab) => {
              const selected = activeTab === tab;
              return (
                <button
                  key={tab}
                  type="button"
                  onClick={() => setActiveTab(tab)}
                  className={[
                    "rounded-xl px-4 py-2 text-sm font-medium transition",
                    selected
                      ? "bg-slate-900 text-white shadow"
                      : "bg-slate-50 text-slate-700 hover:bg-slate-100",
                  ].join(" ")}
                >
                  {tab}
                </button>
              );
            })}
          </div>
        </nav>

        <main>{activeContent}</main>
      </div>
    </div>
  );
}
