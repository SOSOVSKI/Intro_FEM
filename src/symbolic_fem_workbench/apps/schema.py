"""Deterministic schema and export helpers for interactive app sessions."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import UTC, datetime
import json
from typing import Any

import sympy as sp
import yaml

from symbolic_fem_workbench import __version__


SCHEMA_VERSION = "1.0.0"


@dataclass(slots=True)
class RunConfig:
    """Serializable run configuration for reproducible app sessions."""

    preset: str
    parameters: dict[str, Any]
    symbolic_targets: list[str] = field(default_factory=list)
    include_codegen: bool = False


@dataclass(slots=True)
class RunMetadata:
    """Metadata shared with exported session artifacts."""

    package_version: str
    timestamp_utc: str
    selected_preset: str


@dataclass(slots=True)
class RunState:
    """In-memory canonical state from which all exports are derived."""

    config: RunConfig
    symbolic_outputs: dict[str, sp.Expr | sp.Matrix]
    codegen_outputs: dict[str, str] | None = None
    metadata: RunMetadata | None = None

    def ensure_metadata(self) -> RunMetadata:
        if self.metadata is None:
            self.metadata = RunMetadata(
                package_version=__version__,
                timestamp_utc=datetime.now(UTC).isoformat(),
                selected_preset=self.config.preset,
            )
        return self.metadata


@dataclass(slots=True)
class SessionEnvelope:
    """Top-level deterministic schema for persisted app sessions."""

    schema: str
    config: RunConfig
    metadata: RunMetadata


def _stable_dumps(data: dict[str, Any]) -> str:
    return json.dumps(data, sort_keys=True, indent=2)


def session_to_dict(state: RunState) -> dict[str, Any]:
    metadata = state.ensure_metadata()
    envelope = SessionEnvelope(schema=SCHEMA_VERSION, config=state.config, metadata=metadata)
    return asdict(envelope)


def export_config_json(state: RunState) -> str:
    """Export deterministic JSON config payload from in-memory state."""
    return _stable_dumps(session_to_dict(state))


def export_config_yaml(state: RunState) -> str:
    """Export deterministic YAML config payload from in-memory state."""
    data = session_to_dict(state)
    return yaml.safe_dump(data, sort_keys=True)


def export_symbolic_markdown(state: RunState) -> str:
    """Export symbolic outputs as markdown generated directly from current state."""
    lines = ["# Symbolic Outputs", ""]
    for key in sorted(state.symbolic_outputs):
        lines.append(f"## {key}")
        lines.append("```")
        lines.append(str(state.symbolic_outputs[key]))
        lines.append("```")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def export_symbolic_text(state: RunState) -> str:
    """Export symbolic outputs as plain text generated from current state."""
    chunks: list[str] = []
    for key in sorted(state.symbolic_outputs):
        chunks.append(f"[{key}]\n{state.symbolic_outputs[key]}")
    return "\n\n".join(chunks) + ("\n" if chunks else "")


def export_codegen_snippets(state: RunState) -> str:
    """Export generated backend snippets when codegen outputs are available."""
    if not state.codegen_outputs:
        return "# No codegen outputs available.\n"
    lines = ["# Backend snippets", ""]
    for key in sorted(state.codegen_outputs):
        lines.append(f"## {key}")
        lines.append("```")
        lines.append(state.codegen_outputs[key])
        lines.append("```")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def export_run_metadata_markdown(state: RunState) -> str:
    """Shareable markdown section for run metadata."""
    metadata = state.ensure_metadata()
    return (
        "## Share run metadata\n"
        f"- package_version: `{metadata.package_version}`\n"
        f"- timestamp_utc: `{metadata.timestamp_utc}`\n"
        f"- selected_preset: `{metadata.selected_preset}`\n"
    )
