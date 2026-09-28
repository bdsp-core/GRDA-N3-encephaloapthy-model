"""Bout-level screen for RDA-like dynamics in saved model outputs.

This is a transparent hypothesis screen, not an automated clinical EEG
classifier.  It operationalizes the explicit ACNS RDA requirements (0.5--4 Hz,
at least six cycles, relative cycle regularity, similar morphology, and no
interwave interval) and keeps declared engineering thresholds separate from
those requirements.  Every detected bout is also aligned to the saved cortical
and thalamic raster so the population recruitment order can be inspected.

Default outputs are written below ``output/rda_bout_screen``.  The continuous
canonical N3 interval (180--210 s) is analyzed as the primary internal control.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy import signal


RAW_FS_HZ = 40000.0
VIEW_FS_HZ = 200.0
BROAD_BAND_HZ = (0.5, 30.0)
DEFAULT_CYCLE_BAND_HZ = (0.5, 4.0)
MIN_CYCLES = 6
CYCLE_POINTS = 128
FIXED_Q_CYCLES = 8
MIN_PHASE_PLV = 0.20
MIN_PHASE_CYCLE_FRACTION = 0.50

# The first four entries are fixed by the model's gid scheme.
POPULATIONS = {
    "pyr": (0, 500),
    "inh": (500, 600),
    "re": (600, 700),
    "tc": (700, 800),
}

LFP_RE = re.compile(
    r"^lfp_state(?P<state>\d+)"
    r"(?:_seed(?P<seed>\d+))?"
    r"(?:_(?P<tag>tcI[^_]+_r[^_]+ms))?"
    r"_nhost=(?P<nhost>\d+)\.txt(?:\.gz)?$"
)


@dataclass(frozen=True)
class RunSpec:
    """One unique saved simulation readout."""

    state: int | None
    seed: int
    run_tag: str
    nhost: int
    lfp_path: Path
    raster_path: Path
    duration_s: float
    analysis_start_s: float
    analysis_end_s: float
    condition_id: str
    label: str
    is_control: bool = False


@dataclass
class Bout:
    """One merged sequence of overlapping six-cycle passing windows."""

    start_s: float
    end_s: float
    boundary_s: np.ndarray
    metrics: dict[str, object]
    normalized_cycles: np.ndarray
    template: np.ndarray
    acf_lag_s: np.ndarray
    acf: np.ndarray
    psd_hz: np.ndarray
    psd: np.ndarray


@dataclass
class RunResult:
    """Processed trace, bouts, and one run-level summary."""

    spec: RunSpec
    fs: float
    raw_fs: float
    broad: np.ndarray
    delta: np.ndarray
    analysis_slice: slice
    bouts: list[Bout]
    row: dict[str, object]
    median_template: np.ndarray | None = None


def finite_or_blank(value: object) -> object:
    """Replace non-finite scalar values with an empty CSV field."""
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return ""
    return value


def write_csv(rows: Sequence[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"No rows available for {output}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: finite_or_blank(row.get(key, "")) for key in fields})


def load_text_vector(path: Path) -> np.ndarray:
    """Load one whitespace-delimited vector, including optional gzip files."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            values = np.loadtxt(stream, dtype=np.float64)
    else:
        values = np.fromfile(path, sep=" ", dtype=np.float64)
    values = np.ravel(values)
    if values.size == 0:
        raise ValueError(f"Empty signal: {path}")
    if not np.isfinite(values).all():
        raise ValueError(f"Non-finite signal: {path}")
    return values


def infer_duration_s(path: Path) -> tuple[int, float, float]:
    """Load enough metadata to infer a nominal 30/120/360-second duration."""
    opener = gzip.open if path.suffix == ".gz" else open
    samples = 0
    final_byte = b""
    with opener(path, "rb") as stream:
        while True:
            block = stream.read(8 * 1024 * 1024)
            if not block:
                break
            samples += block.count(b"\n")
            final_byte = block[-1:]
    if final_byte and final_byte != b"\n":
        samples += 1
    estimated = samples / RAW_FS_HZ
    nominal = float(round(estimated)) if estimated >= 1.0 else float(estimated)
    if nominal <= 0:
        raise ValueError(f"Cannot infer duration for {path}")
    raw_fs = samples / nominal
    if not np.isclose(raw_fs, RAW_FS_HZ, atol=1.0, rtol=0.0):
        raise ValueError(
            f"Unexpected sampling for {path}: {samples} samples over {nominal:g} s"
        )
    return samples, nominal, raw_fs


def resolve_raster_for_lfp(lfp_path: Path) -> Path:
    """Resolve the raster sharing an LFP filename suffix."""
    name = lfp_path.name
    if not name.startswith("lfp_"):
        raise ValueError(f"Unexpected LFP filename: {name}")
    stem = "raster_" + name[len("lfp_") :]
    candidates = [lfp_path.with_name(stem)]
    if stem.endswith(".txt"):
        candidates.append(lfp_path.with_name(stem + ".gz"))
    elif stem.endswith(".txt.gz"):
        candidates.append(lfp_path.with_name(stem[:-3]))
    existing = [path for path in candidates if path.exists()]
    if not existing:
        raise FileNotFoundError(f"No raster paired with {lfp_path}")
    # Prefer plain text when both an identical plain and compressed copy exist.
    return sorted(existing, key=lambda path: path.suffix == ".gz")[0]


def discover_runs(root: Path, nhost: int = 8) -> list[RunSpec]:
    """Discover and deduplicate complete static 120-second simulations."""
    chosen: dict[tuple[int, int, str, int], tuple[Path, bool]] = {}
    for path in sorted(root.glob("lfp_state*_nhost=*.txt*")):
        match = LFP_RE.match(path.name)
        if not match or int(match.group("nhost")) != nhost:
            continue
        state = int(match.group("state"))
        explicit_seed = match.group("seed") is not None
        seed = int(match.group("seed") or 1)
        tag = match.group("tag") or ""
        key = (state, seed, tag, nhost)
        prior = chosen.get(key)
        if prior is None or (explicit_seed and not prior[1]):
            chosen[key] = (path, explicit_seed)

    runs: list[RunSpec] = []
    for (state, seed, tag, host_count), (path, _) in sorted(chosen.items()):
        samples, duration_s, _ = infer_duration_s(path)
        if samples < int(119.5 * RAW_FS_HZ) or not np.isclose(
            duration_s, 120.0, atol=0.5
        ):
            continue
        raster = resolve_raster_for_lfp(path)
        condition = f"state_{state}" + (f"_{tag}" if tag else "")
        label = f"State {state}, seed {seed}" + (f", {tag}" if tag else "")
        runs.append(
            RunSpec(
                state=state,
                seed=seed,
                run_tag=tag,
                nhost=host_count,
                lfp_path=path,
                raster_path=raster,
                duration_s=duration_s,
                analysis_start_s=10.0,
                analysis_end_s=119.0,
                condition_id=condition,
                label=label,
            )
        )
    return runs


def canonical_n3_spec(root: Path, nhost: int = 8) -> RunSpec:
    canonical = root / "tmp" / "canonical_fink"
    lfp_path = canonical / f"lfp_nhost={nhost}.txt"
    raster_path = canonical / f"raster_nhost={nhost}.txt"
    if not lfp_path.exists() or not raster_path.exists():
        raise FileNotFoundError("Canonical continuous Fink outputs are missing")
    _, duration_s, _ = infer_duration_s(lfp_path)
    return RunSpec(
        state=None,
        seed=1,
        run_tag="",
        nhost=nhost,
        lfp_path=lfp_path,
        raster_path=raster_path,
        duration_s=duration_s,
        analysis_start_s=180.0,
        analysis_end_s=210.0,
        condition_id="canonical_continuous_n3",
        label="Canonical continuous N3, 180–210 s",
        is_control=True,
    )


def preprocess_lfp(
    path: Path,
    duration_s: float,
    view_fs: float = VIEW_FS_HZ,
    broad_band: tuple[float, float] = BROAD_BAND_HZ,
    cycle_band: tuple[float, float] = DEFAULT_CYCLE_BAND_HZ,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Median-center, anti-alias, resample, and zero-phase filter an LFP."""
    raw = load_text_vector(path)
    raw_fs = raw.size / duration_s
    if not np.isclose(raw_fs, RAW_FS_HZ, atol=1.0, rtol=0.0):
        raise ValueError(f"Unexpected raw rate {raw_fs:g} Hz for {path}")
    raw -= np.median(raw)
    ratio = Fraction(view_fs / raw_fs).limit_denominator(10000)
    down = signal.resample_poly(
        raw, ratio.numerator, ratio.denominator, padtype="line"
    )

    def filt(band: tuple[float, float]) -> np.ndarray:
        sos = signal.butter(4, band, btype="bandpass", fs=view_fs, output="sos")
        padlen = min(down.size - 1, int(round(3.0 * view_fs / band[0])))
        y = signal.sosfiltfilt(sos, down, padtype="even", padlen=padlen)
        return y - np.median(y)

    return filt(broad_band), filt(cycle_band), view_fs, raw_fs


def analytic_cycle_boundaries(x: np.ndarray, fs: float) -> np.ndarray:
    """Return monotonically ordered phase-to-phase boundaries in samples."""
    if x.size < int(3 * fs):
        return np.empty(0, dtype=np.float64)
    phase = np.unwrap(np.angle(signal.hilbert(x)))
    # A noisy analytic phase can move briefly backwards.  Its running maximum
    # provides deterministic phase crossings; implausibly short cycles are
    # still rejected explicitly later.
    monotonic = np.maximum.accumulate(phase)
    first_turn = math.ceil(monotonic[0] / (2 * np.pi))
    last_turn = math.floor(monotonic[-1] / (2 * np.pi))
    if last_turn - first_turn < MIN_CYCLES:
        return np.empty(0, dtype=np.float64)
    levels = 2 * np.pi * np.arange(first_turn, last_turn + 1)
    indices = np.searchsorted(monotonic, levels, side="left")
    valid = (indices > 0) & (indices < monotonic.size)
    indices = indices[valid]
    levels = levels[valid]
    lo = indices - 1
    denom = monotonic[indices] - monotonic[lo]
    fraction = np.divide(
        levels - monotonic[lo], denom, out=np.zeros_like(levels), where=denom > 0
    )
    boundaries = lo.astype(np.float64) + fraction
    return boundaries[np.concatenate(([True], np.diff(boundaries) > 1e-6))]


def resample_cycle(x: np.ndarray, start: float, end: float) -> np.ndarray:
    """Interpolate one cycle to a fixed phase grid and robustly normalize it."""
    if end <= start + 2:
        return np.full(CYCLE_POINTS, np.nan)
    sample_positions = np.linspace(start, end, CYCLE_POINTS, endpoint=False)
    indices = np.arange(x.size, dtype=np.float64)
    cycle = np.interp(sample_positions, indices, x)
    cycle -= np.median(cycle)
    amplitude = np.percentile(cycle, 95) - np.percentile(cycle, 5)
    if amplitude <= np.finfo(float).eps:
        return np.full(CYCLE_POINTS, np.nan)
    return cycle / amplitude


def safe_corr(a: np.ndarray, b: np.ndarray) -> float:
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        return np.nan
    if np.std(a) <= 1e-12 or np.std(b) <= 1e-12:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def robust_cv(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=np.float64)
    median = np.median(values)
    if values.size == 0 or abs(median) <= np.finfo(float).eps:
        return np.nan
    mad = np.median(np.abs(values - median))
    return float(1.4826 * mad / abs(median))


def bandpowers(x: np.ndarray, fs: float) -> dict[str, object]:
    """Welch band powers used for screening and sharp/fast flags."""
    if x.size < int(2 * fs):
        return {
            "frequency_hz": np.array([]),
            "power": np.array([]),
            "delta_fraction": np.nan,
            "fast_delta_ratio": np.nan,
        }
    nperseg = min(x.size, int(round(4.0 * fs)))
    frequency, power = signal.welch(
        x,
        fs=fs,
        nperseg=nperseg,
        noverlap=nperseg // 2,
        detrend="constant",
    )

    def area(low: float, high: float) -> float:
        use = (frequency >= low) & (frequency <= high)
        if np.count_nonzero(use) < 2:
            return 0.0
        return float(np.trapezoid(power[use], frequency[use]))

    delta = area(0.5, 4.0)
    total = area(0.5, 30.0)
    fast = area(8.0, 30.0)
    return {
        "frequency_hz": frequency,
        "power": power,
        "delta_fraction": delta / total if total > 0 else np.nan,
        "fast_delta_ratio": fast / delta if delta > 0 else np.nan,
    }


def morphology_metrics(cycles: np.ndarray) -> tuple[dict[str, float], np.ndarray]:
    """Leave-one-out template similarity and adjacent-cycle similarity."""
    if cycles.ndim != 2 or cycles.shape[0] < MIN_CYCLES:
        raise ValueError("Morphology metrics require at least six cycles")
    leave_one_out = []
    for index, cycle in enumerate(cycles):
        template = np.nanmedian(np.delete(cycles, index, axis=0), axis=0)
        leave_one_out.append(safe_corr(cycle, template))
    adjacent = [safe_corr(a, b) for a, b in zip(cycles[:-1], cycles[1:])]
    template = np.nanmedian(cycles, axis=0)
    loo = np.asarray(leave_one_out)
    adj = np.asarray(adjacent)
    return (
        {
            "template_corr_median": float(np.nanmedian(loo)),
            "template_corr_pass_fraction": float(np.mean(loo >= 0.50)),
            "adjacent_corr_median": float(np.nanmedian(adj)),
        },
        template,
    )


def longest_true_run(mask: np.ndarray) -> int:
    padded = np.concatenate(([False], mask, [False])).astype(np.int8)
    edges = np.flatnonzero(np.diff(padded))
    return int(np.max(edges[1::2] - edges[::2])) if edges.size else 0


def quiet_interval_fraction(
    x: np.ndarray,
    boundaries: np.ndarray,
    fs: float,
    median_period_s: float,
) -> float:
    """Operational flag for a discernible resting interval between waves."""
    start = max(0, int(math.floor(boundaries[0])))
    end = min(x.size, int(math.ceil(boundaries[-1])))
    bout = x[start:end]
    if bout.size < 3:
        return np.nan
    derivative = np.abs(np.gradient(bout) * fs)
    low_derivative = derivative <= np.percentile(derivative, 20)
    baseline = float(np.median(bout[low_derivative]))
    flags = []
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        lo = max(0, int(math.floor(left)))
        hi = min(x.size, int(math.ceil(right)))
        cycle = x[lo:hi]
        period_s = (right - left) / fs
        if cycle.size < 3:
            continue
        amplitude = np.percentile(cycle, 95) - np.percentile(cycle, 5)
        if amplitude <= np.finfo(float).eps:
            flags.append(True)
            continue
        slope = np.abs(np.gradient(cycle) * fs)
        quiet = (np.abs(cycle - baseline) < 0.10 * amplitude) & (
            slope < 0.15 * amplitude / max(period_s, 1 / fs)
        )
        required = int(round(max(0.05, 0.20 * period_s) * fs))
        flags.append(longest_true_run(quiet) >= max(required, 1))
    return float(np.mean(flags)) if flags else np.nan


def acf_metrics(
    x: np.ndarray, fs: float, period_s: float
) -> tuple[dict[str, float], np.ndarray, np.ndarray]:
    """Overlap-normalized ACF metrics through three cycles.

    Normalizing each lag by the samples that actually overlap prevents an
    otherwise identical short bout from receiving a lower ACF merely because
    fewer samples remain at later lags.
    """
    centered = x - np.mean(x)
    total_energy = float(np.dot(centered, centered))
    if total_energy <= np.finfo(float).eps:
        return (
            {"acf_first_peak": np.nan, "acf_persistence": np.nan},
            np.array([]),
            np.array([]),
        )
    max_lag = min(x.size - 2, int(round(3.25 * period_s * fs)))
    lag_samples = np.arange(max_lag + 1, dtype=np.int64)
    overlap = x.size - lag_samples
    raw = signal.correlate(centered, centered, mode="full", method="fft")
    raw = raw[x.size - 1 : x.size + max_lag]
    prefix = np.concatenate(([0.0], np.cumsum(centered)))
    prefix_sq = np.concatenate(([0.0], np.cumsum(centered**2)))
    sum_left = prefix[overlap]
    sum_right = prefix[-1] - prefix[lag_samples]
    sum_sq_left = prefix_sq[overlap]
    sum_sq_right = prefix_sq[-1] - prefix_sq[lag_samples]
    covariance = raw - sum_left * sum_right / overlap
    variance_left = sum_sq_left - sum_left**2 / overlap
    variance_right = sum_sq_right - sum_right**2 / overlap
    denominator = np.sqrt(np.maximum(variance_left * variance_right, 0.0))
    corr = np.divide(
        covariance,
        denominator,
        out=np.full_like(covariance, np.nan, dtype=np.float64),
        where=denominator > np.finfo(float).eps,
    )
    corr[0] = 1.0
    lags = lag_samples / fs

    def peak_near(multiplier: int) -> float:
        low = (multiplier - 0.25) * period_s
        high = (multiplier + 0.25) * period_s
        use = (lags >= low) & (lags <= high)
        return float(np.max(corr[use])) if np.any(use) else np.nan

    peaks = np.asarray([peak_near(k) for k in range(1, 4)])
    return (
        {
            "acf_first_peak": float(peaks[0]),
            "acf_persistence": float(np.nanmean(peaks[1:])),
        },
        lags,
        corr,
    )


def multitaper_q(
    x: np.ndarray, fs: float, duration_s: float, cycles: int
) -> tuple[dict[str, object], np.ndarray, np.ndarray]:
    """Return a three-taper PSD and fixed-eight-cycle Q estimate."""
    centered = x - np.mean(x)
    tapers = signal.windows.dpss(centered.size, NW=2.5, Kmax=3, sym=False)
    spectra = []
    for taper in np.atleast_2d(tapers):
        fft = np.fft.rfft(centered * taper)
        spectra.append(np.abs(fft) ** 2)
    power = np.mean(spectra, axis=0)
    frequency = np.fft.rfftfreq(centered.size, 1 / fs)
    use = (frequency >= 0.5) & (frequency <= 4.0)
    if not np.any(use):
        return {"dominant_hz": np.nan, "spectral_q": np.nan, "q_reliable": False}, frequency, power
    eligible_indices = np.flatnonzero(use)
    peak_index = int(eligible_indices[np.argmax(power[use])])
    dominant = float(frequency[peak_index])
    reliable = bool(cycles == FIXED_Q_CYCLES)
    q = np.nan
    if reliable and 0 < peak_index < power.size - 1:
        prominence = signal.peak_prominences(power, [peak_index])[0][0]
        if prominence > 0:
            with np.errstate(invalid="ignore"):
                width_samples = float(
                    signal.peak_widths(power, [peak_index], rel_height=0.5)[0][0]
                )
            df = fs / centered.size
            bandwidth = max(width_samples * df, 2.0 / max(duration_s, 1 / fs))
            q = dominant / bandwidth if bandwidth > 0 else np.nan
    return (
        {"dominant_hz": dominant, "spectral_q": float(q), "q_reliable": reliable},
        frequency,
        power,
    )


def sharp_fast_metrics(
    x: np.ndarray,
    boundaries: np.ndarray,
    fs: float,
    period_s: float,
    fast_delta_ratio: float,
) -> dict[str, object]:
    sharpness = []
    event_phase = []
    event_present = []
    fast_sos = signal.butter(4, (8.0, 30.0), btype="bandpass", fs=fs, output="sos")
    fast = signal.sosfiltfilt(fast_sos, x, padtype="even")
    envelope = np.abs(signal.hilbert(fast))
    threshold = np.median(envelope) + 2.0 * 1.4826 * np.median(
        np.abs(envelope - np.median(envelope))
    )
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        lo = max(0, int(math.floor(left)))
        hi = min(x.size, int(math.ceil(right)))
        cycle = x[lo:hi]
        if cycle.size < 3:
            continue
        amplitude = np.percentile(cycle, 95) - np.percentile(cycle, 5)
        current_period = (right - left) / fs
        if amplitude <= np.finfo(float).eps:
            continue
        sharpness.append(float(np.max(np.abs(np.gradient(cycle) * fs)) * current_period / amplitude))
        local_envelope = envelope[lo:hi]
        peak = int(np.argmax(local_envelope))
        event_phase.append(2 * np.pi * peak / max(local_envelope.size, 1))
        event_present.append(bool(local_envelope[peak] >= threshold))
    phase = np.asarray(event_phase)
    phase_locking = float(abs(np.mean(np.exp(1j * phase)))) if phase.size else np.nan
    median_sharpness = float(np.median(sharpness)) if sharpness else np.nan
    event_fraction = float(np.mean(event_present)) if event_present else np.nan
    criteria = [
        median_sharpness > 8.0,
        fast_delta_ratio > 0.15,
        event_fraction >= 0.80,
        phase_locking >= 0.70,
    ]
    return {
        "sharpness_median": median_sharpness,
        "fast_event_cycle_fraction": event_fraction,
        "fast_event_phase_plv": phase_locking,
        "sharp_fast_phase_locked": bool(sum(bool(item) for item in criteria) >= 3),
    }


def window_passes(
    broad: np.ndarray,
    fs: float,
    boundaries: np.ndarray,
    normalized_cycles: np.ndarray,
    run_scale: float,
) -> tuple[bool, dict[str, float]]:
    periods = np.diff(boundaries) / fs
    in_band_fraction = float(np.mean((periods >= 0.25) & (periods <= 2.0)))
    ratios = np.maximum(periods[:-1], periods[1:]) / np.minimum(
        periods[:-1], periods[1:]
    )
    regular_fraction = float(np.mean(ratios < 1.5)) if ratios.size else 0.0
    morphology, _ = morphology_metrics(normalized_cycles)
    lo = max(0, int(math.floor(boundaries[0])))
    hi = min(broad.size, int(math.ceil(boundaries[-1])))
    segment = broad[lo:hi]
    power = bandpowers(segment, fs)
    amplitudes = []
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        cycle = broad[int(math.floor(left)) : int(math.ceil(right))]
        amplitudes.append(np.percentile(cycle, 95) - np.percentile(cycle, 5))
    amplitude_ratio = float(np.median(amplitudes) / max(run_scale, 1e-12))
    metrics = {
        "median_frequency_hz": float(1.0 / np.median(periods)),
        "period_regular_pair_fraction": regular_fraction,
        "period_cv_robust": robust_cv(periods),
        "in_frequency_band_fraction": in_band_fraction,
        "delta_fraction": float(power["delta_fraction"]),
        "amplitude_vs_run_scale": amplitude_ratio,
        **morphology,
    }
    passes = bool(
        0.5 <= metrics["median_frequency_hz"] <= 4.0
        and in_band_fraction >= 0.80
        and regular_fraction > 0.50
        and metrics["delta_fraction"] >= 0.35
        and metrics["template_corr_median"] >= 0.65
        and metrics["template_corr_pass_fraction"] >= (2 / 3)
        and metrics["adjacent_corr_median"] >= 0.50
    )
    return passes, metrics


def globally_valid_bout_ranges(
    broad: np.ndarray,
    fs: float,
    boundaries: np.ndarray,
    normalized_cycles: np.ndarray,
    run_scale: float,
    windows: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """Greedily extend overlapping windows only while their union still passes.

    Touching but non-overlapping six-cycle windows are deliberately separate.
    When adding the next overlapping window makes the full union fail, the
    current valid bout is closed and that failed extension becomes the start of
    the next candidate. This may retain overlapping bouts, whose burden is
    therefore calculated from the union of their time intervals.
    """
    if not windows:
        return []
    starts = [start for start, _ in windows]
    ranges: list[tuple[int, int]] = []
    cursor = 0
    while cursor < len(starts):
        start = starts[cursor]
        end = start + MIN_CYCLES
        cursor += 1
        while cursor < len(starts) and starts[cursor] < end:
            candidate_end = starts[cursor] + MIN_CYCLES
            passed, _ = window_passes(
                broad,
                fs,
                boundaries[start : candidate_end + 1],
                normalized_cycles[start:candidate_end],
                run_scale,
            )
            if not passed:
                break
            end = candidate_end
            cursor += 1
        ranges.append((start, end))
    return ranges


def harmonic_ratios(
    frequency: np.ndarray, power: np.ndarray, fundamental: float
) -> dict[str, float]:
    def local_power(center: float) -> float:
        if center <= 0 or center > 30:
            return np.nan
        half_width = max(0.15, 0.20 * fundamental)
        use = (frequency >= center - half_width) & (frequency <= center + half_width)
        return float(np.trapezoid(power[use], frequency[use])) if np.count_nonzero(use) >= 2 else np.nan

    first = local_power(fundamental)
    second = local_power(2 * fundamental)
    third = local_power(3 * fundamental)
    return {
        "second_harmonic_ratio": second / first if first and np.isfinite(first) else np.nan,
        "third_harmonic_ratio": third / first if first and np.isfinite(first) else np.nan,
    }


def build_bout(
    broad: np.ndarray,
    fs: float,
    all_boundaries: np.ndarray,
    cycle_start: int,
    cycle_end: int,
    offset_s: float,
) -> Bout:
    boundaries = all_boundaries[cycle_start : cycle_end + 1]
    lo = max(0, int(math.floor(boundaries[0])))
    hi = min(broad.size, int(math.ceil(boundaries[-1])))
    segment = broad[lo:hi]
    local_boundaries = boundaries - lo
    periods = np.diff(boundaries) / fs
    cycles = np.vstack(
        [resample_cycle(broad, left, right) for left, right in zip(boundaries[:-1], boundaries[1:])]
    )
    morphology, template = morphology_metrics(cycles)
    power = bandpowers(segment, fs)
    median_period = float(np.median(periods))
    acf_values, acf_lag, acf = acf_metrics(segment, fs, median_period)
    if periods.size >= FIXED_Q_CYCLES:
        q_start = (periods.size - FIXED_Q_CYCLES) // 2
        q_lo = max(0, int(math.floor(boundaries[q_start])))
        q_hi = min(
            broad.size,
            int(math.ceil(boundaries[q_start + FIXED_Q_CYCLES])),
        )
        q_segment = broad[q_lo:q_hi]
        q_cycles = FIXED_Q_CYCLES
    else:
        q_segment = segment
        q_cycles = int(periods.size)
    q_values, psd_hz, psd = multitaper_q(
        q_segment, fs, q_segment.size / fs, q_cycles
    )
    quiet_fraction = quiet_interval_fraction(
        broad, boundaries, fs, median_period
    )
    sharp = sharp_fast_metrics(
        broad, boundaries, fs, median_period, float(power["fast_delta_ratio"])
    )
    ratios = np.maximum(periods[:-1], periods[1:]) / np.minimum(
        periods[:-1], periods[1:]
    )
    amplitudes = []
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        cycle = broad[int(math.floor(left)) : int(math.ceil(right))]
        amplitudes.append(np.percentile(cycle, 95) - np.percentile(cycle, 5))
    start_s = offset_s + boundaries[0] / fs
    end_s = offset_s + boundaries[-1] / fs
    metrics: dict[str, object] = {
        "start_s": float(start_s),
        "end_s": float(end_s),
        "duration_s": float(end_s - start_s),
        "n_cycles": int(periods.size),
        "median_frequency_hz": float(1.0 / median_period),
        "period_regular_pair_fraction": float(np.mean(ratios < 1.5)),
        "period_cv_robust": robust_cv(periods),
        "amplitude_cv_robust": robust_cv(np.asarray(amplitudes)),
        "delta_fraction": float(power["delta_fraction"]),
        "fast_delta_ratio": float(power["fast_delta_ratio"]),
        "quiet_cycle_fraction": quiet_fraction,
        "continuity_class": (
            "RDA-compatible"
            if quiet_fraction <= 0.20
            else "indeterminate"
            if quiet_fraction < 0.50
            else "interval/PD-like"
        ),
        **morphology,
        **acf_values,
        **q_values,
        **sharp,
    }
    metrics.update(harmonic_ratios(psd_hz, psd, float(q_values["dominant_hz"])))
    return Bout(
        start_s=float(start_s),
        end_s=float(end_s),
        boundary_s=offset_s + boundaries / fs,
        metrics=metrics,
        normalized_cycles=cycles,
        template=template,
        acf_lag_s=acf_lag,
        acf=acf,
        psd_hz=psd_hz,
        psd=psd,
    )


def detect_bouts(
    broad: np.ndarray,
    delta: np.ndarray,
    fs: float,
    offset_s: float,
) -> list[Bout]:
    """Detect and merge transparent six-cycle RDA-screen windows."""
    boundaries = analytic_cycle_boundaries(delta, fs)
    if boundaries.size < MIN_CYCLES + 1:
        return []
    normalized = np.vstack(
        [resample_cycle(broad, left, right) for left, right in zip(boundaries[:-1], boundaries[1:])]
    )
    run_scale = float(np.percentile(broad, 95) - np.percentile(broad, 5))
    passing: list[tuple[int, int]] = []
    for start in range(0, normalized.shape[0] - MIN_CYCLES + 1):
        end = start + MIN_CYCLES
        cycle_block = normalized[start:end]
        if not np.isfinite(cycle_block).all():
            continue
        passed, _ = window_passes(
            broad,
            fs,
            boundaries[start : end + 1],
            cycle_block,
            run_scale,
        )
        if passed:
            passing.append((start, end))
    merged = globally_valid_bout_ranges(
        broad, fs, boundaries, normalized, run_scale, passing
    )
    return [
        build_bout(broad, fs, boundaries, start, end, offset_s)
        for start, end in merged
        if end - start >= MIN_CYCLES
    ]


def load_raster(path: Path) -> np.ndarray:
    """Load and time-sort a two-column ``time_ms gid`` raster."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            raster = np.loadtxt(stream, dtype=np.float64)
    else:
        raster = np.loadtxt(path, dtype=np.float64)
    if raster.size == 0:
        return np.empty((0, 2), dtype=np.float64)
    raster = np.atleast_2d(raster)
    if raster.shape[1] != 2 or not np.isfinite(raster).all():
        raise ValueError(f"Invalid raster: {path}")
    return raster[np.argsort(raster[:, 0], kind="stable")]


def circular_summary(phases: np.ndarray) -> tuple[float, float]:
    if phases.size == 0:
        return np.nan, np.nan
    vector = np.mean(np.exp(1j * phases))
    return float(np.mod(np.angle(vector), 2 * np.pi)), float(abs(vector))


def phases_within_cycles(
    spikes_s: np.ndarray, boundaries_s: np.ndarray
) -> np.ndarray:
    cycle_index = np.searchsorted(boundaries_s, spikes_s, side="right") - 1
    valid = (cycle_index >= 0) & (cycle_index < boundaries_s.size - 1)
    if not np.any(valid):
        return np.empty(0, dtype=np.float64)
    index = cycle_index[valid]
    period = boundaries_s[index + 1] - boundaries_s[index]
    return 2 * np.pi * (spikes_s[valid] - boundaries_s[index]) / period


def sector_entropy(gids: np.ndarray, low: int, high: int, sectors: int = 10) -> float:
    """Normalized entropy across the model's periodic one-dimensional topology."""
    if gids.size == 0:
        return np.nan
    coordinate = (gids - low) / (high - low)
    count, _ = np.histogram(coordinate, bins=np.linspace(0, 1, sectors + 1))
    probability = count[count > 0] / np.sum(count)
    return float(-np.sum(probability * np.log(probability)) / np.log(sectors))


def population_bout_metrics(raster: np.ndarray, bout: Bout) -> dict[str, object]:
    """Align one detected bout to cortical and thalamic population spikes."""
    start_ms = bout.start_s * 1000.0
    end_ms = bout.end_s * 1000.0
    in_bout = (raster[:, 0] >= start_ms) & (raster[:, 0] < end_ms)
    data = raster[in_bout]
    duration_s = bout.end_s - bout.start_s
    phases: dict[str, np.ndarray] = {}
    metrics: dict[str, object] = {}
    n_cycles = max(1, bout.boundary_s.size - 1)
    for name, (low, high) in POPULATIONS.items():
        use = (data[:, 1] >= low) & (data[:, 1] < high)
        population = data[use]
        spike_s = population[:, 0] / 1000.0
        phases[name] = phases_within_cycles(spike_s, bout.boundary_s)
        cycle_index = np.searchsorted(
            bout.boundary_s, spike_s, side="right"
        ) - 1
        valid_cycle = (cycle_index >= 0) & (cycle_index < n_cycles)
        participating_cycle_fraction = float(
            np.unique(cycle_index[valid_cycle]).size / n_cycles
        )
        mean_phase, plv = circular_summary(phases[name])
        metrics[f"{name}_rate_hz_per_cell"] = float(
            population.shape[0] / ((high - low) * duration_s)
        )
        metrics[f"{name}_active_cell_fraction"] = float(
            np.unique(population[:, 1]).size / (high - low)
        )
        metrics[f"{name}_lfp_phase_rad"] = mean_phase
        metrics[f"{name}_lfp_phase_plv"] = plv
        metrics[f"{name}_phase_spike_count"] = int(phases[name].size)
        metrics[f"{name}_participating_cycle_fraction"] = (
            participating_cycle_fraction
        )
        metrics[f"{name}_sector_entropy"] = sector_entropy(
            population[:, 1], low, high
        )

    pyr_phase = metrics["pyr_lfp_phase_rad"]
    tc_phase = metrics["tc_lfp_phase_rad"]
    phase_is_supported = bool(
        np.isfinite(pyr_phase)
        and np.isfinite(tc_phase)
        and float(metrics["pyr_lfp_phase_plv"]) >= MIN_PHASE_PLV
        and float(metrics["tc_lfp_phase_plv"]) >= MIN_PHASE_PLV
        and int(metrics["pyr_phase_spike_count"]) >= n_cycles
        and int(metrics["tc_phase_spike_count"]) >= n_cycles
        and float(metrics["pyr_participating_cycle_fraction"])
        >= MIN_PHASE_CYCLE_FRACTION
        and float(metrics["tc_participating_cycle_fraction"])
        >= MIN_PHASE_CYCLE_FRACTION
    )
    if phase_is_supported:
        difference = float(np.angle(np.exp(1j * (pyr_phase - tc_phase))))
        # Positive values mean that the TC mean phase precedes the PYR phase.
        metrics["tc_to_pyr_phase_lag_ms"] = (
            difference
            * float(bout.metrics["median_frequency_hz"]) ** -1
            * 1000.0
            / (2 * np.pi)
        )
    else:
        metrics["tc_to_pyr_phase_lag_ms"] = np.nan

    recruitment_lags = []
    for left, right in zip(bout.boundary_s[:-1], bout.boundary_s[1:]):
        pyr_spikes = data[
            (data[:, 1] >= POPULATIONS["pyr"][0])
            & (data[:, 1] < POPULATIONS["pyr"][1])
            & (data[:, 0] >= left * 1000)
            & (data[:, 0] < right * 1000),
            0,
        ]
        tc_spikes = data[
            (data[:, 1] >= POPULATIONS["tc"][0])
            & (data[:, 1] < POPULATIONS["tc"][1])
            & (data[:, 0] >= left * 1000)
            & (data[:, 0] < right * 1000),
            0,
        ]
        # A 10th-percentile recruitment time is more stable than the first
        # isolated spike.  Require enough spikes to make it interpretable.
        if pyr_spikes.size >= 5 and tc_spikes.size >= 2:
            recruitment_lags.append(
                float(np.quantile(pyr_spikes, 0.10) - np.quantile(tc_spikes, 0.10))
            )
    metrics["tc_to_pyr_recruitment_lag_ms"] = (
        float(np.median(recruitment_lags)) if recruitment_lags else np.nan
    )
    metrics["tc_first_among_corecruited_cycles_fraction"] = (
        float(np.mean(np.asarray(recruitment_lags) > 0))
        if recruitment_lags
        else np.nan
    )
    metrics["cycles_with_tc_pyr_recruitment"] = int(len(recruitment_lags))
    metrics["tc_pyr_corecruited_cycle_fraction"] = float(
        len(recruitment_lags) / n_cycles
    )
    return metrics


def median_numeric(rows: Sequence[dict[str, object]], key: str) -> float:
    values = []
    for row in rows:
        value = row.get(key, np.nan)
        if isinstance(value, (int, float, np.integer, np.floating)) and np.isfinite(value):
            values.append(float(value))
    return float(np.median(values)) if values else np.nan


def choose_representative_bout(bouts: Sequence[Bout]) -> Bout | None:
    if not bouts:
        return None

    compatible = [
        bout
        for bout in bouts
        if bout.metrics["continuity_class"] == "RDA-compatible"
    ]
    pool = compatible if compatible else list(bouts)

    def rank(bout: Bout) -> tuple[float, ...]:
        return (
            float(not bool(bout.metrics["sharp_fast_phase_locked"])),
            float(bout.metrics["template_corr_median"]),
            float(bout.metrics["acf_first_peak"]),
            -float(bout.metrics["period_cv_robust"]),
            float(bout.metrics["duration_s"]),
        )

    return max(pool, key=rank)


def interval_union_duration(bouts: Sequence[Bout]) -> float:
    """Return covered seconds without double-counting overlapping bouts."""
    if not bouts:
        return 0.0
    intervals = sorted((bout.start_s, bout.end_s) for bout in bouts)
    total = 0.0
    current_start, current_end = intervals[0]
    for start, end in intervals[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            total += current_end - current_start
            current_start, current_end = start, end
    return float(total + current_end - current_start)


def summarize_run(
    spec: RunSpec,
    bouts: list[Bout],
    raw_fs: float,
    view_fs: float,
    analysis_duration_s: float,
) -> tuple[dict[str, object], np.ndarray | None]:
    compatible = [
        bout for bout in bouts if bout.metrics["continuity_class"] == "RDA-compatible"
    ]
    clean_compatible = [
        bout
        for bout in compatible
        if not bool(bout.metrics["sharp_fast_phase_locked"])
    ]
    indeterminate = [
        bout for bout in bouts if bout.metrics["continuity_class"] == "indeterminate"
    ]
    total_duration = interval_union_duration(bouts)
    compatible_duration = interval_union_duration(compatible)
    clean_compatible_duration = interval_union_duration(clean_compatible)
    representative = choose_representative_bout(bouts)
    median_template = (
        np.median(np.vstack([bout.template for bout in compatible]), axis=0)
        if compatible
        else None
    )
    row: dict[str, object] = {
        "condition_id": spec.condition_id,
        "state": "" if spec.state is None else spec.state,
        "seed": spec.seed,
        "run_tag": spec.run_tag,
        "is_control": spec.is_control,
        "source_file": str(spec.lfp_path),
        "raster_file": str(spec.raster_path),
        "duration_s": spec.duration_s,
        "analysis_start_s": spec.analysis_start_s,
        "analysis_end_s": spec.analysis_end_s,
        "analysis_duration_s": analysis_duration_s,
        "raw_fs_hz": raw_fs,
        "view_fs_hz": view_fs,
        "n_screen_bouts": len(bouts),
        "n_rda_compatible_bouts": len(compatible),
        "n_clean_rda_compatible_bouts": len(clean_compatible),
        "n_indeterminate_bouts": len(indeterminate),
        "screen_burden_pct": 100.0 * total_duration / analysis_duration_s,
        "rda_compatible_burden_pct": 100.0
        * compatible_duration
        / analysis_duration_s,
        "clean_rda_compatible_burden_pct": 100.0
        * clean_compatible_duration
        / analysis_duration_s,
        "longest_bout_s": max(
            (float(bout.metrics["duration_s"]) for bout in bouts), default=0.0
        ),
        "longest_rda_compatible_bout_s": max(
            (float(bout.metrics["duration_s"]) for bout in compatible),
            default=0.0,
        ),
        "sharp_fast_bout_fraction": (
            float(
                np.mean(
                    [bool(bout.metrics["sharp_fast_phase_locked"]) for bout in bouts]
                )
            )
            if bouts
            else np.nan
        ),
        "sharp_fast_compatible_bout_fraction": (
            float(
                np.mean(
                    [
                        bool(bout.metrics["sharp_fast_phase_locked"])
                        for bout in compatible
                    ]
                )
            )
            if compatible
            else np.nan
        ),
    }
    summary_keys = (
        "duration_s",
        "n_cycles",
        "median_frequency_hz",
        "period_regular_pair_fraction",
        "period_cv_robust",
        "amplitude_cv_robust",
        "delta_fraction",
        "fast_delta_ratio",
        "quiet_cycle_fraction",
        "template_corr_median",
        "adjacent_corr_median",
        "acf_first_peak",
        "acf_persistence",
        "dominant_hz",
        "spectral_q",
        "sharpness_median",
        "fast_event_phase_plv",
        "second_harmonic_ratio",
        "third_harmonic_ratio",
        "tc_to_pyr_phase_lag_ms",
        "tc_to_pyr_recruitment_lag_ms",
        "tc_first_among_corecruited_cycles_fraction",
        "tc_pyr_corecruited_cycle_fraction",
        "pyr_active_cell_fraction",
        "tc_active_cell_fraction",
        "pyr_lfp_phase_plv",
        "tc_lfp_phase_plv",
        "pyr_sector_entropy",
        "tc_sector_entropy",
    )
    for prefix, selected in (
        ("median_screen_bout", bouts),
        ("median_compatible_bout", compatible),
    ):
        selected_metrics = [bout.metrics for bout in selected]
        for key in summary_keys:
            row[f"{prefix}_{key}"] = median_numeric(selected_metrics, key)
    if representative is not None:
        row["representative_bout_start_s"] = representative.start_s
        row["representative_bout_end_s"] = representative.end_s
        row["representative_continuity_class"] = representative.metrics[
            "continuity_class"
        ]
    else:
        row["representative_bout_start_s"] = ""
        row["representative_bout_end_s"] = ""
        row["representative_continuity_class"] = ""
    return row, median_template


def analyze_spec(
    spec: RunSpec,
    view_fs: float,
    broad_band: tuple[float, float],
    cycle_band: tuple[float, float],
) -> RunResult:
    broad_full, delta_full, fs, raw_fs = preprocess_lfp(
        spec.lfp_path, spec.duration_s, view_fs, broad_band, cycle_band
    )
    start = int(round(spec.analysis_start_s * fs))
    end = min(broad_full.size, int(round(spec.analysis_end_s * fs)))
    if end - start < int(8 * fs):
        raise ValueError(f"Analysis epoch is too short for {spec.label}")
    broad = broad_full[start:end]
    delta = delta_full[start:end]
    bouts = detect_bouts(broad, delta, fs, spec.analysis_start_s)
    if bouts:
        raster = load_raster(spec.raster_path)
        for bout in bouts:
            bout.metrics.update(population_bout_metrics(raster, bout))
    row, median_template = summarize_run(
        spec, bouts, raw_fs, fs, (end - start) / fs
    )
    return RunResult(
        spec=spec,
        fs=fs,
        raw_fs=raw_fs,
        broad=broad,
        delta=delta,
        analysis_slice=slice(start, end),
        bouts=bouts,
        row=row,
        median_template=median_template,
    )


def pairwise_template_corr(templates: Sequence[np.ndarray]) -> float:
    correlations = []
    for index, first in enumerate(templates):
        for second in templates[index + 1 :]:
            correlations.append(safe_corr(first, second))
    return float(np.nanmedian(correlations)) if correlations else np.nan


def summarize_conditions(
    results: Sequence[RunResult], control: RunResult
) -> list[dict[str, object]]:
    groups: dict[str, list[RunResult]] = {}
    for result in results:
        groups.setdefault(result.spec.condition_id, []).append(result)

    n3 = control.row
    condition_rows: list[dict[str, object]] = []
    for condition_id, group in groups.items():
        group = sorted(group, key=lambda result: result.spec.seed)
        run_rows = [result.row for result in group]
        with_bouts = [result for result in group if result.bouts]
        compatible_seeds = [
            result
            for result in group
            if int(result.row["n_rda_compatible_bouts"]) > 0
        ]
        clean_compatible_seeds = [
            result
            for result in group
            if int(result.row["n_clean_rda_compatible_bouts"]) > 0
        ]
        seed_frequency = np.asarray(
            [
                float(result.row["median_compatible_bout_median_frequency_hz"])
                for result in compatible_seeds
                if np.isfinite(
                    result.row["median_compatible_bout_median_frequency_hz"]
                )
            ]
        )
        templates = [
            result.median_template
            for result in compatible_seeds
            if result.median_template is not None
        ]
        frequency_cv = (
            float(np.std(seed_frequency, ddof=1) / np.mean(seed_frequency))
            if seed_frequency.size >= 2 and np.mean(seed_frequency) > 0
            else np.nan
        )
        template_corr = pairwise_template_corr(templates)
        n_seeds = len(group)
        n_with = len(with_bouts)
        n_compatible = len(compatible_seeds)
        replicated = bool(
            n_seeds >= 2
            and n_compatible >= 2
            and np.isfinite(frequency_cv)
            and frequency_cv <= 0.20
            and np.isfinite(template_corr)
            and template_corr >= 0.60
        )
        if replicated:
            evidence = "replicated RDA-compatible bouts"
            evidence_tier = 3
        elif n_seeds >= 2 and n_compatible >= 1:
            evidence = "multi-seed, not replicated"
            evidence_tier = 2
        elif n_compatible >= 1:
            evidence = "one-seed lead"
            evidence_tier = 1
        elif n_with >= 1:
            evidence = "screen bouts only; none RDA-compatible"
            evidence_tier = 0
        else:
            evidence = "no eligible bout"
            evidence_tier = 0

        state = group[0].spec.state
        row: dict[str, object] = {
            "condition_id": condition_id,
            "state": "" if state is None else state,
            "run_tag": group[0].spec.run_tag,
            "seeds_available": ",".join(str(result.spec.seed) for result in group),
            "n_seeds": n_seeds,
            "seeds_with_screen_bout": n_with,
            "seeds_with_rda_compatible_bout": n_compatible,
            "seeds_with_clean_rda_compatible_bout": len(
                clean_compatible_seeds
            ),
            "compatible_seed_fraction": n_compatible / n_seeds,
            "compatible_seed_result": f"{n_compatible}/{n_seeds}",
            "replicated": replicated,
            "evidence": evidence,
            "evidence_tier": evidence_tier,
            "across_seed_frequency_cv": frequency_cv,
            "across_seed_template_corr": template_corr,
        }
        for key in (
            "screen_burden_pct",
            "rda_compatible_burden_pct",
            "clean_rda_compatible_burden_pct",
            "longest_bout_s",
            "longest_rda_compatible_bout_s",
            "sharp_fast_bout_fraction",
            "sharp_fast_compatible_bout_fraction",
            "median_compatible_bout_duration_s",
            "median_compatible_bout_n_cycles",
            "median_compatible_bout_median_frequency_hz",
            "median_compatible_bout_period_cv_robust",
            "median_compatible_bout_delta_fraction",
            "median_compatible_bout_quiet_cycle_fraction",
            "median_compatible_bout_template_corr_median",
            "median_compatible_bout_adjacent_corr_median",
            "median_compatible_bout_acf_first_peak",
            "median_compatible_bout_acf_persistence",
            "median_compatible_bout_spectral_q",
            "median_compatible_bout_sharpness_median",
            "median_compatible_bout_fast_delta_ratio",
            "median_compatible_bout_tc_to_pyr_phase_lag_ms",
            "median_compatible_bout_tc_to_pyr_recruitment_lag_ms",
            "median_compatible_bout_tc_first_among_corecruited_cycles_fraction",
            "median_compatible_bout_tc_pyr_corecruited_cycle_fraction",
            "median_compatible_bout_pyr_active_cell_fraction",
            "median_compatible_bout_tc_active_cell_fraction",
            "median_compatible_bout_pyr_lfp_phase_plv",
            "median_compatible_bout_tc_lfp_phase_plv",
            "median_compatible_bout_pyr_sector_entropy",
            "median_compatible_bout_tc_sector_entropy",
        ):
            row[key] = median_numeric(run_rows, key)

        row["sharp_fast_warning"] = any(
            isinstance(
                run_row.get("sharp_fast_compatible_bout_fraction"),
                (int, float, np.integer, np.floating),
            )
            and np.isfinite(
                run_row.get("sharp_fast_compatible_bout_fraction", np.nan)
            )
            and float(run_row["sharp_fast_compatible_bout_fraction"]) > 0
            for run_row in run_rows
        )

        # N3 is a dependent internal control, not a biological population.
        # Report transparent differences rather than a composite score.
        for key in (
            "median_compatible_bout_template_corr_median",
            "median_compatible_bout_period_cv_robust",
            "median_compatible_bout_acf_first_peak",
            "median_compatible_bout_acf_persistence",
            "median_compatible_bout_quiet_cycle_fraction",
            "median_compatible_bout_spectral_q",
            "rda_compatible_burden_pct",
        ):
            reference = n3.get(key, np.nan)
            current = row.get(key, np.nan)
            row[f"{key}_minus_n3"] = (
                float(current) - float(reference)
                if isinstance(current, (int, float, np.integer, np.floating))
                and isinstance(reference, (int, float, np.integer, np.floating))
                and np.isfinite(current)
                and np.isfinite(reference)
                else np.nan
            )
        if templates and control.median_template is not None:
            row["template_corr_to_continuous_n3"] = safe_corr(
                np.median(np.vstack(templates), axis=0), control.median_template
            )
        else:
            row["template_corr_to_continuous_n3"] = np.nan
        condition_rows.append(row)

    def order(row: dict[str, object]) -> tuple[float, ...]:
        def number(key: str, default: float) -> float:
            value = row.get(key, default)
            return float(value) if isinstance(value, (int, float)) and np.isfinite(value) else default

        return (
            -number("evidence_tier", 0),
            -number("seeds_with_clean_rda_compatible_bout", 0),
            number("sharp_fast_compatible_bout_fraction", 1),
            -number("median_compatible_bout_template_corr_median", -1),
            -number("median_compatible_bout_acf_first_peak", -1),
            number("median_compatible_bout_period_cv_robust", 99),
            -number("rda_compatible_burden_pct", 0),
            int(row["state"]) if row["state"] != "" else 9999,
        )

    condition_rows.sort(key=order)
    for rank, row in enumerate(condition_rows, start=1):
        row["review_order"] = rank
    return condition_rows


def flatten_bouts(results: Sequence[RunResult]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        for index, bout in enumerate(result.bouts, start=1):
            row: dict[str, object] = {
                "condition_id": result.spec.condition_id,
                "state": "" if result.spec.state is None else result.spec.state,
                "seed": result.spec.seed,
                "run_tag": result.spec.run_tag,
                "is_control": result.spec.is_control,
                "bout_index": index,
                "source_file": str(result.spec.lfp_path),
            }
            row.update(bout.metrics)
            rows.append(row)
    return rows


def representative_result(group: Sequence[RunResult]) -> RunResult | None:
    eligible = [result for result in group if result.bouts]
    if not eligible:
        return None

    def rank(result: RunResult) -> tuple[float, ...]:
        row = result.row
        clean = float(row["n_clean_rda_compatible_bouts"])
        compatible = float(row["n_rda_compatible_bouts"])
        template = float(row["median_compatible_bout_template_corr_median"])
        acf = float(row["median_compatible_bout_acf_first_peak"])
        cv = float(row["median_compatible_bout_period_cv_robust"])
        return (
            clean,
            compatible,
            template if np.isfinite(template) else -1,
            acf if np.isfinite(acf) else -1,
            -cv if np.isfinite(cv) else -99,
            float(row["longest_bout_s"]),
        )

    return max(eligible, key=rank)


def add_panel_letter(ax: plt.Axes, letter: str) -> None:
    ax.text(
        0.008,
        0.96,
        letter,
        transform=ax.transAxes,
        fontweight="bold",
        fontsize=11,
        va="top",
        ha="left",
        zorder=10,
    )


def make_review_page(pdf: PdfPages, result: RunResult) -> None:
    bout = choose_representative_bout(result.bouts)
    if bout is None:
        return
    fs = result.fs
    absolute_time = result.spec.analysis_start_s + np.arange(result.broad.size) / fs
    context_start = max(result.spec.analysis_start_s, bout.start_s - 5.0)
    context_end = min(result.spec.analysis_end_s, bout.end_s + 5.0)
    context = (absolute_time >= context_start) & (absolute_time <= context_end)
    within = (absolute_time >= bout.start_s) & (absolute_time <= bout.end_s)

    raster = load_raster(result.spec.raster_path)
    raster_s = raster[:, 0] / 1000.0
    raster_context = (raster_s >= context_start) & (raster_s <= context_end)

    fig = plt.figure(figsize=(11, 8.5), facecolor="white")
    grid = fig.add_gridspec(
        4,
        2,
        height_ratios=(0.90, 1.0, 0.92, 0.82),
        left=0.075,
        right=0.975,
        top=0.855,
        bottom=0.075,
        hspace=0.62,
        wspace=0.28,
    )
    m = bout.metrics
    pdf_label = result.spec.label.replace("–", "-").replace("·", "/")
    fig.suptitle(pdf_label, x=0.075, y=0.965, ha="left", fontsize=16)
    compatible_status = (
        "RDA-compatible screen bout"
        if m["continuity_class"] == "RDA-compatible"
        else "No RDA-compatible bout; highest-ranked screen bout shown"
    )
    if bool(m["sharp_fast_phase_locked"]):
        compatible_status += " / sharp-fast phase-locked warning"
    run_warning = result.row.get("sharp_fast_compatible_bout_fraction", np.nan)
    if (
        isinstance(run_warning, (int, float, np.integer, np.floating))
        and np.isfinite(run_warning)
        and float(run_warning) > 0
        and not bool(m["sharp_fast_phase_locked"])
    ):
        compatible_status += f" / {float(run_warning):.0%} of compatible bouts flagged sharp-fast"
    fig.text(
        0.075,
        0.925,
        (
            f"Representative bout {bout.start_s:.1f}-{bout.end_s:.1f} s / "
            f"{m['median_frequency_hz']:.2f} Hz / {m['n_cycles']} cycles / "
            f"template r={m['template_corr_median']:.2f} / "
            f"period CV={m['period_cv_robust']:.2f}"
        ),
        ha="left",
        fontsize=9.5,
        color="#424b55",
    )
    fig.text(
        0.075,
        0.897,
        compatible_status,
        ha="left",
        fontsize=9.5,
        color="#424b55",
    )

    ax = fig.add_subplot(grid[0, :])
    ax.plot(absolute_time[context], result.broad[context], color="#1c4e80", lw=0.65)
    ax.plot(absolute_time[context], result.delta[context], color="#d17c32", lw=1.0, alpha=0.8)
    ax.axvspan(bout.start_s, bout.end_s, color="#5b8f5a", alpha=0.12)
    ax.set(xlim=(context_start, context_end), ylabel="LFP (a.u.)", xlabel="Time (s)")
    ax.grid(axis="x", color="#d8dde2", lw=0.5)
    ax.legend(["0.5-30 Hz", "0.5-4 Hz"], frameon=False, ncol=2, loc="upper right")
    add_panel_letter(ax, "A")

    ax = fig.add_subplot(grid[1, 0])
    ax.plot(absolute_time[within], result.broad[within], color="#1c4e80", lw=0.75)
    for boundary in bout.boundary_s:
        ax.axvline(boundary, color="#9aa3aa", lw=0.45)
    ax.set(xlim=(bout.start_s, bout.end_s), ylabel="LFP (a.u.)", xlabel="Time (s)")
    ax.grid(axis="x", color="#e1e5e8", lw=0.45)
    add_panel_letter(ax, "B")

    ax = fig.add_subplot(grid[1, 1])
    phase = np.linspace(0, 1, CYCLE_POINTS, endpoint=False)
    for cycle in bout.normalized_cycles:
        ax.plot(phase, cycle, color="#9aa3aa", lw=0.6, alpha=0.45)
    ax.plot(phase, bout.template, color="#1c4e80", lw=2.0, label="Median template")
    ax.set(xlim=(0, 1), ylabel="Normalized LFP", xlabel="Cycle phase")
    ax.grid(color="#e1e5e8", lw=0.45)
    ax.legend(frameon=False, loc="upper right")
    add_panel_letter(ax, "C")

    ax = fig.add_subplot(grid[2, 0])
    ax.plot(bout.acf_lag_s, bout.acf, color="#1c4e80", lw=1.0)
    for multiple in range(1, 6):
        ax.axvline(
            multiple / float(m["median_frequency_hz"]),
            color="#9aa3aa",
            lw=0.45,
        )
    ax.set(
        xlim=(0, bout.acf_lag_s[-1] if bout.acf_lag_s.size else 1),
        ylim=(-1, 1.05),
        ylabel="Autocorrelation",
        xlabel="Lag (s)",
    )
    ax.grid(color="#e1e5e8", lw=0.45)
    add_panel_letter(ax, "D")

    ax = fig.add_subplot(grid[2, 1])
    use_psd = (bout.psd_hz >= 0.5) & (bout.psd_hz <= 15.0)
    ax.semilogy(
        bout.psd_hz[use_psd],
        np.maximum(bout.psd[use_psd], np.finfo(float).tiny),
        color="#1c4e80",
        lw=1.0,
    )
    ax.axvspan(0.5, 4.0, color="#d17c32", alpha=0.10)
    ax.axvline(float(m["dominant_hz"]), color="#d17c32", lw=0.8)
    ax.set(xlim=(0.5, 15.0), ylabel="Power (a.u.)", xlabel="Frequency (Hz)")
    ax.grid(color="#e1e5e8", lw=0.45)
    add_panel_letter(ax, "E")

    ax = fig.add_subplot(grid[3, :])
    raster_subset = raster[raster_context]
    if raster_subset.size:
        colors = np.where(
            raster_subset[:, 1] < 500,
            "#1c4e80",
            np.where(
                raster_subset[:, 1] < 600,
                "#aa4c3c",
                np.where(raster_subset[:, 1] < 700, "#7b5d96", "#3b7d62"),
            ),
        )
        ax.scatter(
            raster_subset[:, 0] / 1000.0,
            raster_subset[:, 1],
            s=0.6,
            c=colors,
            linewidths=0,
            rasterized=True,
        )
    ax.axvspan(bout.start_s, bout.end_s, color="#5b8f5a", alpha=0.10)
    for boundary in (500, 600, 700):
        ax.axhline(boundary, color="#a8afb5", lw=0.4)
    ax.set(
        xlim=(context_start, context_end),
        ylim=(800, 0),
        ylabel="Cell gid (PYR / INH / RE / TC)",
        xlabel="Time (s)",
    )
    ax.text(
        0.006,
        0.92,
        "F",
        transform=ax.transAxes,
        fontweight="bold",
        fontsize=11,
        va="top",
    )
    pdf.savefig(fig, dpi=180)
    plt.close(fig)


def make_no_bout_page(pdf: PdfPages, result: RunResult) -> None:
    """Show a negative seed rather than silently omitting it from review."""
    fs = result.fs
    absolute_time = result.spec.analysis_start_s + np.arange(result.broad.size) / fs
    fig, ax = plt.subplots(figsize=(11, 8.5), facecolor="white")
    fig.suptitle(
        result.spec.label.replace("–", "-").replace("·", "/"),
        x=0.075,
        y=0.965,
        ha="left",
        fontsize=16,
    )
    fig.text(
        0.075,
        0.91,
        "No qualifying six-cycle screen bout was detected in this seed.",
        ha="left",
        fontsize=11,
        color="#424b55",
    )
    ax.plot(absolute_time, result.broad, color="#1c4e80", lw=0.45)
    ax.set(
        xlim=(result.spec.analysis_start_s, result.spec.analysis_end_s),
        xlabel="Time (s)",
        ylabel="LFP (a.u.)",
    )
    ax.grid(axis="x", color="#d8dde2", lw=0.5)
    fig.subplots_adjust(left=0.09, right=0.975, top=0.84, bottom=0.12)
    pdf.savefig(fig, dpi=180)
    plt.close(fig)


def make_review_pdf(
    results: Sequence[RunResult],
    control: RunResult,
    condition_rows: Sequence[dict[str, object]],
    output: Path,
    top_n: int,
) -> list[str]:
    groups: dict[str, list[RunResult]] = {}
    for result in results:
        groups.setdefault(result.spec.condition_id, []).append(result)
    selected_conditions: list[str] = []
    selected_ids = {control.spec.condition_id}

    def add_rows(rows: Sequence[dict[str, object]], limit: int) -> None:
        added = 0
        for row in rows:
            if len(selected_conditions) >= top_n or added >= limit:
                return
            condition_id = str(row["condition_id"])
            candidate = representative_result(groups[condition_id])
            if candidate is not None and condition_id not in selected_ids:
                selected_conditions.append(condition_id)
                selected_ids.add(condition_id)
                added += 1

    repeated = [row for row in condition_rows if bool(row["replicated"])]
    unstable_multiseed = [
        row
        for row in condition_rows
        if int(row["n_seeds"]) >= 2
        and not bool(row["replicated"])
        and int(row["seeds_with_rda_compatible_bout"]) > 0
    ]
    clean_one_seed = [
        row
        for row in condition_rows
        if int(row["n_seeds"]) == 1
        and int(row["seeds_with_clean_rda_compatible_bout"]) > 0
        and not bool(row["sharp_fast_warning"])
    ]
    repeated_limit = min(len(repeated), max(1, int(math.ceil(0.55 * top_n))))
    add_rows(repeated, repeated_limit)
    add_rows(unstable_multiseed, 1)
    add_rows(clean_one_seed, top_n - len(selected_conditions))
    add_rows(condition_rows, top_n - len(selected_conditions))
    output.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(output) as pdf:
        make_review_page(pdf, control)
        for condition_id in selected_conditions:
            for result in sorted(groups[condition_id], key=lambda item: item.spec.seed):
                if result.bouts:
                    make_review_page(pdf, result)
                else:
                    make_no_bout_page(pdf, result)
    return [control.spec.condition_id, *selected_conditions]


def make_overview_plot(
    condition_rows: Sequence[dict[str, object]],
    control: RunResult,
    output: Path,
    label_n: int = 24,
) -> None:
    eligible = [
        row
        for row in condition_rows
        if isinstance(
            row.get("median_compatible_bout_period_cv_robust"), (int, float)
        )
        and np.isfinite(row["median_compatible_bout_period_cv_robust"])
        and isinstance(
            row.get("median_compatible_bout_template_corr_median"),
            (int, float),
        )
        and np.isfinite(row["median_compatible_bout_template_corr_median"])
    ]
    fig, ax = plt.subplots(figsize=(10.5, 7.0), facecolor="white")
    palette = {
        "replicated RDA-compatible bouts": "#2d6a4f",
        "multi-seed, not replicated": "#cc7a00",
        "one-seed lead": "#406c9b",
        "screen bouts only; none RDA-compatible": "#999999",
        "no eligible bout": "#999999",
    }
    for evidence, color in palette.items():
        subset = [row for row in eligible if row["evidence"] == evidence]
        if not subset:
            continue
        for warning, marker in ((False, "o"), (True, "X")):
            marked = [
                row for row in subset if bool(row["sharp_fast_warning"]) == warning
            ]
            if not marked:
                continue
            x = [
                float(row["median_compatible_bout_period_cv_robust"])
                for row in marked
            ]
            y = [
                float(row["median_compatible_bout_template_corr_median"])
                for row in marked
            ]
            size = [
                24
                + 2.2
                * min(float(row["rda_compatible_burden_pct"]), 50.0)
                for row in marked
            ]
            label = evidence + (" / sharp-fast warning" if warning else "")
            ax.scatter(
                x,
                y,
                s=size,
                color=color,
                marker=marker,
                alpha=0.70,
                linewidths=0,
                label=label,
            )

    n3_x = float(control.row["median_compatible_bout_period_cv_robust"])
    n3_y = float(control.row["median_compatible_bout_template_corr_median"])
    if np.isfinite(n3_x) and np.isfinite(n3_y):
        ax.axvline(n3_x, color="#7b5d96", lw=1.0, alpha=0.8)
        ax.axhline(n3_y, color="#7b5d96", lw=1.0, alpha=0.8)
        ax.scatter([n3_x], [n3_y], marker="*", s=220, color="#7b5d96", zorder=5)
        ax.annotate("continuous N3", (n3_x, n3_y), xytext=(6, 6), textcoords="offset points")

    for row in eligible[:label_n]:
        short_tag = str(row["run_tag"]).replace("tcI", "I").replace(
            "_r5000ms", ""
        )
        label = f"{row['state']}" + (f"/{short_tag}" if short_tag else "")
        ax.annotate(
            label,
            (
                float(row["median_compatible_bout_period_cv_robust"]),
                float(row["median_compatible_bout_template_corr_median"]),
            ),
            xytext=(4, 3),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set(
        xlabel="Robust cycle-period CV (lower is more regular)",
        ylabel="Median cycle-template correlation (higher is more monomorphic)",
        title="Bout-level RDA screen across completed simulations",
    )
    ax.grid(color="#e1e5e8", lw=0.55)
    ax.legend(frameon=False, loc="lower left")
    fig.text(
        0.99,
        0.01,
        "Point size = RDA-compatible burden; X = sharp-fast warning; * = TC-current intervention",
        ha="right",
        fontsize=8.5,
        color="#505861",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=190, bbox_inches="tight")
    plt.close(fig)


def synthetic_filter(x: np.ndarray, fs: float) -> tuple[np.ndarray, np.ndarray]:
    def filt(band: tuple[float, float]) -> np.ndarray:
        sos = signal.butter(4, band, btype="bandpass", fs=fs, output="sos")
        return signal.sosfiltfilt(sos, x, padtype="even")

    return filt(BROAD_BAND_HZ), filt(DEFAULT_CYCLE_BAND_HZ)


def run_self_tests() -> None:
    """Deterministic synthetic checks for frequency and scale behavior."""
    fs = VIEW_FS_HZ
    time = np.arange(int(24 * fs)) / fs
    rng = np.random.default_rng(42)

    sine = np.sin(2 * np.pi * 1.5 * time) + 0.02 * rng.standard_normal(time.size)
    broad, delta = synthetic_filter(sine, fs)
    bouts = detect_bouts(broad[400:-400], delta[400:-400], fs, 2.0)
    if not bouts or not any(1.35 <= float(b.metrics["median_frequency_hz"]) <= 1.65 for b in bouts):
        raise AssertionError("1.5-Hz sinusoidal control did not pass")

    scaled_broad, scaled_delta = synthetic_filter(7.5 * sine, fs)
    scaled = detect_bouts(
        scaled_broad[400:-400], scaled_delta[400:-400], fs, 2.0
    )
    if not scaled:
        raise AssertionError("Amplitude-scaled sinusoid did not pass")
    if abs(
        float(bouts[0].metrics["median_frequency_hz"])
        - float(scaled[0].metrics["median_frequency_hz"])
    ) > 0.02:
        raise AssertionError("Detector is not amplitude invariant")

    theta = np.sin(2 * np.pi * 7.0 * time)
    theta_broad, theta_delta = synthetic_filter(theta, fs)
    if detect_bouts(theta_broad[400:-400], theta_delta[400:-400], fs, 2.0):
        raise AssertionError("7-Hz theta control falsely passed")

    noise = rng.standard_normal(time.size)
    noise_broad, noise_delta = synthetic_filter(noise, fs)
    if detect_bouts(noise_broad[400:-400], noise_delta[400:-400], fs, 2.0):
        raise AssertionError("White-noise control falsely passed")

    acf_first = []
    for cycle_count in (6, 8, 12, 20):
        local_time = np.arange(int(cycle_count * fs)) / fs
        acf_value, _, _ = acf_metrics(
            np.sin(2 * np.pi * local_time), fs, 1.0
        )
        acf_first.append(float(acf_value["acf_first_peak"]))
    if np.ptp(acf_first) > 0.01:
        raise AssertionError("ACF remains confounded by bout length")
    print("synthetic detector checks passed")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--nhost", type=int, default=8)
    parser.add_argument("--view-fs", type=float, default=VIEW_FS_HZ)
    parser.add_argument("--broad-band", nargs=2, type=float, default=BROAD_BAND_HZ)
    parser.add_argument(
        "--cycle-band", nargs=2, type=float, default=DEFAULT_CYCLE_BAND_HZ
    )
    parser.add_argument("--states", nargs="*", type=int)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("output/rda_bout_screen")
    )
    parser.add_argument(
        "--review-pdf",
        type=Path,
        default=Path("output/pdf/rda_bout_screen_review.pdf"),
    )
    parser.add_argument(
        "--top-review-pages",
        type=int,
        default=14,
        help="Number of shortlisted conditions; every available seed is shown",
    )
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.self_test:
        run_self_tests()
        return
    root = args.root.resolve()
    output_dir = (
        args.output_dir
        if args.output_dir.is_absolute()
        else (root / args.output_dir)
    )
    review_pdf = (
        args.review_pdf
        if args.review_pdf.is_absolute()
        else (root / args.review_pdf)
    )
    runs = discover_runs(root, args.nhost)
    if args.states:
        requested = set(args.states)
        runs = [run for run in runs if run.state in requested]
    if not runs:
        raise SystemExit("No complete 120-second runs found")
    control_spec = canonical_n3_spec(root, args.nhost)
    print(f"found {len(runs)} complete 120-second runs")

    results: list[RunResult] = []
    for index, spec in enumerate(runs, start=1):
        print(f"[{index:03d}/{len(runs):03d}] {spec.label}", flush=True)
        results.append(
            analyze_spec(
                spec,
                args.view_fs,
                tuple(args.broad_band),
                tuple(args.cycle_band),
            )
        )
    print("analyzing canonical continuous N3 control", flush=True)
    control = analyze_spec(
        control_spec,
        args.view_fs,
        tuple(args.broad_band),
        tuple(args.cycle_band),
    )

    run_rows = [result.row for result in results] + [control.row]
    bout_rows = flatten_bouts([*results, control])
    condition_rows = summarize_conditions(results, control)
    write_csv(run_rows, output_dir / "rda_run_metrics.csv")
    if bout_rows:
        write_csv(bout_rows, output_dir / "rda_bout_metrics.csv")
    write_csv(condition_rows, output_dir / "rda_condition_summary.csv")
    make_overview_plot(
        condition_rows, control, output_dir / "rda_screen_overview.png"
    )
    reviewed = make_review_pdf(
        results,
        control,
        condition_rows,
        review_pdf,
        args.top_review_pages,
    )
    provenance = {
        "method": "transparent bout-level RDA screen",
        "acns_fixed_requirements": {
            "frequency_hz": [0.5, 4.0],
            "minimum_cycles": MIN_CYCLES,
            "adjacent_period_ratio": "<1.5 in >50% of pairs",
            "morphology": "relatively uniform",
            "continuity": "no interwave interval",
        },
        "declared_engineering_thresholds": {
            "delta_fraction_min": 0.35,
            "median_leave_one_out_template_corr_min": 0.65,
            "cycle_template_pass_fraction_min": 2 / 3,
            "median_adjacent_cycle_corr_min": 0.50,
            "amplitude_vs_run_scale": "reported for QC; not used as a gate",
            "quiet_cycle_fraction_rda_compatible_max": 0.20,
            "quiet_cycle_fraction_pd_like_min": 0.50,
        },
        "preprocessing": {
            "anti_alias_resampling": "scipy.signal.resample_poly",
            "view_fs_hz": args.view_fs,
            "broad_band_hz": list(args.broad_band),
            "cycle_locator_band_hz": list(args.cycle_band),
            "filter": "4th-order Butterworth, zero phase",
        },
        "runs_analyzed": len(results),
        "conditions_analyzed": len(condition_rows),
        "continuous_n3_window_s": [
            control_spec.analysis_start_s,
            control_spec.analysis_end_s,
        ],
        "reviewed_conditions": reviewed,
        "important_limit": (
            "This is a model-field-potential screen, not an automated clinical "
            "EEG or GRDA diagnosis. Generalized spatial distribution is not tested."
        ),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "rda_screen_provenance.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(provenance, stream, indent=2)
        stream.write("\n")
    print(f"wrote bout screen to {output_dir}")


if __name__ == "__main__":
    main()
