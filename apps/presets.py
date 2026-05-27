from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Preset:
    """Normalized preset payload for the UI layer."""

    key: str
    title: str
    description: str
    parameters: dict[str, object]


def load_bar_1d_preset() -> Preset:
    return Preset(
        key="bar_1d",
        title="1D Bar Local Workflow",
        description="Mirrors examples/bar_1d_workflow.py.",
        parameters={
            "problem_type": "bar_1d",
            "source_q": "q",
            "neumann_load": "P",
            "length": "L",
            "youngs_modulus": "E",
            "area": "A",
        },
    )


def load_triangle_p1_poisson_preset() -> Preset:
    return Preset(
        key="triangle_p1_poisson",
        title="Triangle P1 Poisson Local Workflow",
        description="Mirrors examples/triangle_p1_poisson_workflow.py.",
        parameters={
            "problem_type": "triangle_p1_poisson",
            "source_f": "f",
            "x1": "x1",
            "y1": "y1",
            "x2": "x2",
            "y2": "y2",
            "x3": "x3",
            "y3": "y3",
        },
    )


def load_manual_assembly_square_4tri_preset() -> Preset:
    return Preset(
        key="manual_assembly_square_4tri",
        title="Manual Assembly: Square / 4 Triangles",
        description="Mirrors examples/manual_assembly_square_4tri.py.",
        parameters={
            "problem_type": "manual_assembly_square_4tri",
            "source_value": 1.0,
            "boundary_nodes": [0, 1, 2, 3],
            "nodes": [
                [0.0, 0.0],
                [1.0, 0.0],
                [1.0, 1.0],
                [0.0, 1.0],
                [0.5, 0.5],
            ],
            "elements": [
                [0, 1, 4],
                [1, 2, 4],
                [2, 3, 4],
                [3, 0, 4],
            ],
        },
    )


def load_all_presets() -> dict[str, Preset]:
    presets = (
        load_bar_1d_preset(),
        load_triangle_p1_poisson_preset(),
        load_manual_assembly_square_4tri_preset(),
    )
    return {preset.key: preset for preset in presets}
