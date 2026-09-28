"""Shared filename helpers for state/seed outputs and TC-current runs."""

from __future__ import annotations

import math
from pathlib import Path


def filename_float(value: float) -> str:
    """Match the deterministic float encoding used by ``config.py``."""
    if value == 0.0:
        value = 0.0
    return (
        format(value, ".12g")
        .replace("-", "m")
        .replace("+", "p")
        .replace(".", "p")
    )


def tc_current_tag(current_uA_cm2: float | None, ramp_ms: float) -> str:
    """Return the filename tag for an optional TC-current intervention."""
    if current_uA_cm2 is None:
        return ""
    if not math.isfinite(current_uA_cm2) or not -1.0 <= current_uA_cm2 <= 0.0:
        raise ValueError("TC current must be finite and between -1.0 and 0.0")
    if not math.isfinite(ramp_ms) or ramp_ms < 0.0:
        raise ValueError("TC-current ramp must be finite and nonnegative")
    if current_uA_cm2 == 0.0:
        return ""
    return (
        f"tcI{filename_float(current_uA_cm2)}_"
        f"r{filename_float(ramp_ms)}ms"
    )


def state_output_name(
    prefix: str,
    state: int,
    seed: int,
    nhost: int,
    extension: str,
    run_tag: str = "",
) -> str:
    """Build one exact seed-tagged output filename."""
    tag = f"_{run_tag}" if run_tag else ""
    return (
        f"{prefix}_state{state}_seed{seed}{tag}_"
        f"nhost={nhost}{extension}"
    )


def resolve_seeded_output(
    root: Path,
    prefix: str,
    state: int,
    seed: int,
    nhost: int,
    extension: str,
    run_tag: str = "",
    allow_seed1_legacy: bool = False,
) -> Path:
    """Resolve an exact tagged output, optionally allowing old seed-1 files."""
    seeded = root / state_output_name(
        prefix, state, seed, nhost, extension, run_tag
    )
    if seeded.exists():
        return seeded
    if run_tag:
        raise FileNotFoundError(f"No exact current-tagged output: {seeded}")
    legacy = root / f"{prefix}_state{state}_nhost={nhost}{extension}"
    if allow_seed1_legacy and seed == 1 and legacy.exists():
        return legacy
    raise FileNotFoundError(
        f"No output found for state {state}, seed {seed}: {seeded}"
    )
