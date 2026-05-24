"""BAMCP evaluation harness.

A MARRVEL-MCP-compatible runner that exercises BAMCP tools through an LLM
provider, captures the resulting tool-call telemetry, and grades each case
deterministically (when possible) plus via an LLM judge (fallback).

See :mod:`bamcp.eval.runner` for the main entry point.
"""

from .schema import EvalCase, EvalResult, GraderVerdict, RunConfig, load_cases

__all__ = ["EvalCase", "EvalResult", "GraderVerdict", "RunConfig", "load_cases"]
