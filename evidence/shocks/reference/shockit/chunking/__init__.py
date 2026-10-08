"""Explicit chunked shock-finder workflow helpers."""

from .layout import save_chunked_field
from .reader import NpzChunkedInput
from .writer import NpzChunkedOutput, join_all_fields, join_chunked_field

__all__ = [
    "NpzChunkedInput",
    "NpzChunkedOutput",
    "join_all_fields",
    "join_chunked_field",
    "run_chunked_shock_finder",
    "save_chunked_field",
]


def __getattr__(name: str):
    """Resolve compatibility exports lazily to avoid import cycles."""
    if name == "run_chunked_shock_finder":
        from .runner import run_chunked_shock_finder

        return run_chunked_shock_finder
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
