"""Compatibility module for core concepts."""

from .config import ShockFinderConfig
from .fields import FluidCube
from .finder import ShockFinder
from .result import ShockFinderResult

__all__ = [
    "FluidCube",
    "ShockFinder",
    "ShockFinderConfig",
    "ShockFinderResult",
]
