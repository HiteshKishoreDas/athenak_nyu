"""Compatibility wrapper for the standalone chunked center reducer."""

from __future__ import annotations

from pathlib import Path
import json
from copy import deepcopy

from .provenance import atomic_write_json
from .local_pass import build_reporter

from .chunking.centers import (
    build_base_summary_from_output_store,
    finalize_chunked_shock_outputs,
)
from .chunking.layout import load_chunk_layout
from .chunking.writer import NpzChunkedOutput
from .config import ShockFinderConfig


def reduce_chunked_shock_outputs(
    output_root: str | Path,
    config: ShockFinderConfig,
    *,
    progress: bool = True,
) -> dict[str, object]:
    """Recompute public masks from immutable full fields in a completed result.

    Interrupted retries use incomplete.json; ordinary opening must reject that
    marker even if summary.json also exists. Simple-score centers need at least
    28 bytes/global cell before temporaries. Exact medians retain detections.
    """
    config.validate()
    config = deepcopy(config)
    build_reporter(progress)
    output_path = Path(output_root)
    marker = output_path / "incomplete.json"
    summary_path = output_path / "summary.json"
    if marker.exists():
        summary = json.loads(marker.read_text()).get("previous_summary")
        if not isinstance(summary, dict) or not summary:
            raise ValueError("Incomplete reduction marker has no valid previous summary.")
    elif summary_path.exists():
        summary = json.loads(summary_path.read_text())
        if not isinstance(summary, dict) or not summary:
            raise ValueError("Standalone reduction requires a complete output summary.")
    else:
        raise ValueError("Standalone reduction requires a complete output summary.")
    layout = load_chunk_layout(output_path / "full_shock_mask")
    output_store = NpzChunkedOutput(output_path, layout)
    atomic_write_json(marker, {"previous_summary": summary})
    summary_path.unlink(missing_ok=True)
    result = finalize_chunked_shock_outputs(output_store, config, summary, progress=progress)
    marker.unlink()
    return result
