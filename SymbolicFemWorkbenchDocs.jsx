import React, { useMemo, useState } from "react";

const tabs = [
  "Overview",
  "Structure",
  "Workflow",
  "Modules",
  "Examples",
  "Teaching",
  "Roadmap",
];

const packageTree = `src/symbolic_fem_workbench/
  symbols.py
  forms.py
  transforms.py
  fe_spaces.py
  reference.py
  quadrature.py
  extract.py
  validate.py
  workflow.py
  assembly.py
  printers/
    iheartla_printer.py
  codegen/
    iheartla_backend.py`;

const workflowSteps = [
  {
    title: "1. Define symbolic objects",
    body:
      "Start with domains, coordinates, fields, test functions, and coefficients. The package should make these explicit instead of flattening everything into anonymous SymPy symbols.",
  },
  {
    title: "2. Build the weighted residual and weak form",
    body:
      "Use the transform helpers to move from the strong statement to a weighted residual and then to a weak form through explicit divergence-form integration by parts and boundary handling.",
  },
  {
    title: "3. Insert local FE approximations",
    body:
      "Use the FE-space helpers and reference-element definitions to substitute local trial and test expansions. This is where continuous fields become element-level algebra.",
  },
  {
    title: "4. Integrate and extract local tensors",
    body:
      "Use quadrature or exact symbolic integration where appropriate, then extract coefficient matrices and vectors to recover K_e and f_e.",
  },
  {
    title: "5. Assemble globally in plain Python",
    body:
      "Assembly remains visible. Connectivity, local-to-global indexing, Dirichlet elimination, and solving should stay understandable to students rather than disappearing behind a framework wall.",
  },
  {
    title: "6. Export local algebra only at the end",
    body:
      "Use the printer/codegen layer only after the mathematics is settled. I❤️LA should be a narrow backend bridge, not a substitute for formulation.",
  },
];

const modules = [
  {
    name: "symbols.py",
    purpose: "Teaching-first symbolic building blocks for domains, coordinates, and fields.",
    exports: ["Domain objects", "field constructors", "coordinate helpers"],
  },
  {
    name: "forms.py",
    purpose: "Structured containers for domain integrals, boundary terms, weighted residuals, and split weak forms.",
    exports: ["DomainIntegral", "BoundaryContribution", "WeightedResidual", "WeakForm"],
  },
  {
    name: "transforms.py",
    purpose: "Explicit symbolic transformations: weighting, integration by parts, BC handling, FE substitution, pullback, and linearization.",
    exports: [
      "weighted_residual",
      "integrate_divergence_*",
      "boundary helpers",
      "substitute_fe",
      "gateaux_derivative",
    ],
  },
  {
    name: "fe_spaces.py",
    purpose: "Finite-element-space helpers and local trial/test expansions.",
    exports: ["1D element helpers", "local_trial_expansion", "local_test_expansion"],
  },
  {
    name: "reference.py",
    purpose: "Reference-element definitions, shape functions, gradients, and affine maps.",
    exports: ["Q1", "Q2", "Hex8", "Tet4", "Tet10", "triangle P1 helpers"],
  },
  {
    name: "quadrature.py",
    purpose: "Reference-domain quadrature and exact symbolic integration helpers.",
    exports: ["triangle rules", "quad rules", "hex rules", "tet rules"],
  },
  {
    name: "extract.py",
    purpose: "Extraction of coefficient matrices and vectors from symbolic FE expressions.",
    exports: ["extract_coefficient_matrix", "extract_coefficient_vector"],
  },
  {
    name: "validate.py",
    purpose: "Validation helpers for workflow sanity checks and future guards.",
    exports: ["linearity checks", "dependency checks", "consistency checks"],
  },
  {
    name: "workflow.py",
    purpose: "Course-ready worked pipelines that connect the lower-level pieces.",
    exports: ["1D bar workflow", "triangle Poisson workflow"],
  },
  {
    name: "assembly.py",
    purpose: "Manual global assembly helpers kept explicit for teaching purposes.",
    exports: [
      "assemble_dense_matrix",
      "assemble_dense_vector",
      "apply_dirichlet_by_reduction",
      "expand_reduced_solution",
    ],
  },
  {
    name: "printers/iheartla_printer.py",
    purpose: "Narrow formatting layer for local expressions that are ready to become I❤️LA text.",
    exports: ["small matrix/vector/scalar printer helpers"],
  },
  {
    name: "codegen/iheartla_backend.py",
    purpose: "Backend glue for exporting or compiling finalized local algebra through I❤️LA.",
    exports: ["backend bridge helpers"],
  },
];

const examples = [
  {
    title: "bar_1d_workflow.py",
    focus: "Strong form to weak form to element tensors in 1D.",
  },
  {
    title: "triangle_p1_poisson_workflow.py",
    focus: "Reference-triangle Poisson element pipeline in 2D.",
  },
  {
    title: "manual_assembly_square_4tri.py",
    focus: "Manual local-to-global assembly on a tiny 2D mesh.",
  },
  {
    title: "reference_element_library_demo.py",
    focus: "Inspection of reference elements, shape functions, and gradients.",
  },
];

const teachingTopics = [
  {
    topic: "Weak form before coding",
    body:
      "Use the package to show that weighting, integration by parts, and boundary-term reasoning are mathematical steps, not software magic.",
  },
  {
    topic: "Discretization as a visible substitution",
    body:
      "Students should see u_h and v_h inserted explicitly and watch the local matrix emerge from symbolic algebra.",
  },
  {
    topic: "Assembly as its own topic",
    body:
      "Keep local tensor derivation separate from topology and bookkeeping so students can debug one layer at a time.",
  },
  {
    topic: "2D and 3D as geometry extensions",
    body:
      "Reference elements, mappings, Jacobians, and quadrature become the new ideas, not a total rewrite of the method.",
  },
];

const roadmap = [
  "Expand validate.py so it catches linearity mistakes, missing Jacobians, and trial/test misuse early.",
  "Broaden workflow examples carefully instead of turning the package into a giant abstraction pile.",
  "Keep the I❤️LA bridge narrow and local-expression focused.",
  "Add notebook teaching assets that mirror lecture progression without becoming the main codebase.",
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
            <SectionCard title="What this package is">
              <p>
                <strong>symbolic-fem-workbench</strong> is a teaching-first symbolic FEM package built around a very
                unfashionable idea: students should see the mathematics before the software hides it.
              </p>
              <p className="mt-3">
                The package covers symbolic weak-form construction, FE substitution, local tensor extraction,
                reference-element definitions, quadrature, and manual assembly helpers. It is not trying to be a
                full FEM framework because that would defeat the point.
              </p>
            </SectionCard>
            <SectionCard title="Design rules">
              <ul className="list-disc space-y-2 pl-5">
                <li>Keep analytical steps explicit.</li>
                <li>Keep local algebra separate from global topology.</li>
                <li>Use structured weak-form objects instead of expression soup.</li>
                <li>Let code generation come last.</li>
                <li>Prefer clarity over framework theatrics.</li>
              </ul>
            </SectionCard>
          </div>
        );

      case "Structure":
        return (
          <SectionCard title="Current package structure" subtitle="This matches the aligned package layout.">
            <pre className="overflow-x-auto rounded-xl bg-slate-50 p-4 text-xs text-slate-800">{packageTree}</pre>
          </SectionCard>
        );

      case "Workflow":
        return (
          <div className="grid gap-5 lg:grid-cols-2">
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
                <div className="mb-3 flex flex-wrap gap-2">
                  {module.exports.map((item) => (
                    <Pill key={item}>{item}</Pill>
                  ))}
                </div>
              </SectionCard>
            ))}
          </div>
        );

      case "Examples":
        return (
          <div className="grid gap-5 lg:grid-cols-2">
            {examples.map((example) => (
              <SectionCard key={example.title} title={example.title}>
                <p>{example.focus}</p>
              </SectionCard>
            ))}
          </div>
        );

      case "Teaching":
        return (
          <div className="grid gap-5 lg:grid-cols-2">
            {teachingTopics.map((item) => (
              <SectionCard key={item.topic} title={item.topic}>
                <p>{item.body}</p>
              </SectionCard>
            ))}
          </div>
        );

      case "Roadmap":
        return (
          <SectionCard title="Near-term roadmap">
            <ul className="list-disc space-y-2 pl-5">
              {roadmap.map((item) => (
                <li key={item}>{item}</li>
              ))}
            </ul>
          </SectionCard>
        );

      default:
        return null;
    }
  }, [activeTab]);

  return (
    <div className="min-h-screen bg-slate-50 p-4 md:p-8">
      <div className="mx-auto max-w-7xl">
        <header className="mb-6 rounded-3xl border border-slate-200 bg-white p-6 shadow-sm md:p-8">
          <div className="flex flex-col gap-4 lg:flex-row lg:items-end lg:justify-between">
            <div>
              <p className="text-xs font-bold uppercase tracking-[0.2em] text-slate-500">Package documentation</p>
              <h1 className="mt-2 text-3xl font-bold tracking-tight text-slate-950 md:text-5xl">
                symbolic-fem-workbench
              </h1>
              <p className="mt-4 max-w-3xl text-sm leading-6 text-slate-600 md:text-base">
                Documentation for the teaching-first symbolic FEM package, with tabs for structure, workflow,
                modules, examples, teaching use, and roadmap. Because code and docs occasionally deserve to agree.
              </p>
            </div>
            <div className="flex flex-wrap gap-2">
              <Pill>SymPy-first</Pill>
              <Pill>Teaching-first</Pill>
              <Pill>Local tensors visible</Pill>
              <Pill>Assembly still manual</Pill>
              <Pill>uv-managed project</Pill>
            </div>
          </div>
        </header>

        <div className="mb-6 overflow-x-auto rounded-2xl border border-slate-200 bg-white p-2 shadow-sm">
          <div className="flex min-w-max gap-2">
            {tabs.map((tab) => {
              const active = activeTab === tab;
              return (
                <button
                  key={tab}
                  type="button"
                  onClick={() => setActiveTab(tab)}
                  className={`rounded-xl px-4 py-2 text-sm font-semibold transition ${
                    active ? "bg-slate-900 text-white" : "bg-slate-100 text-slate-700 hover:bg-slate-200"
                  }`}
                >
                  {tab}
                </button>
              );
            })}
          </div>
        </div>

        {activeContent}
      </div>
    </div>
  );
}
