"""Reproducible summary metrics for static-state LFP simulations.

The analysis used for the project handoff is intentionally simple and fully
declared here: median centering, polyphase anti-alias resampling to 200 Hz,
zero-phase 0.5-30 Hz filtering, and Welch spectral estimates on an edge-free
physiological epoch.

Beyond the amplitude and band-power summary, three metric families separate
GRDA-like activity from the N3 slow oscillation along the axes that define the
clinical pattern:

1. Rhythmicity (``delta_peak_q``, ``acf_first_peak``, ``acf_sidelobes``,
   ``cycle_cv``, ``inst_freq_sd_hz``). A thalamic delta clock is narrowband and
   repeats; the cortical slow oscillation is broadband and irregular.
2. Absence of faster activity (``fast_slow_ratio``, ``sef95_hz``, plus the
   existing band percentages).
3. Absence of spindles (``spindle_per_min``, ``sigma_env_cv``, ``sigma_pct``).

Every metric here is computable on scalp EEG as well as on the model LFP, and
all of them are amplitude-scale-free, so model and patient recordings can be
compared without calibrating arbitrary units. Metrics that require the network
state rather than the field potential live in ``analyze_network_state.py``.
"""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy import signal, stats


def resolve_lfp_file(root: Path, state: int, seed: int, nhost: int) -> Path:
    """Resolve new seed-tagged outputs, retaining seed-1 legacy support."""
    seeded = root / f"lfp_state{state}_seed{seed}_nhost={nhost}.txt"
    if seeded.exists():
        return seeded
    legacy = root / f"lfp_state{state}_nhost={nhost}.txt"
    if seed == 1 and legacy.exists():
        return legacy
    raise FileNotFoundError(f"No LFP output found for state {state}, seed {seed}")


def load_eegband(
    path: Path,
    duration_s: float = 120.0,
    view_fs: float = 200.0,
    band: tuple[float, float] = (0.5, 30.0),
) -> tuple[np.ndarray, float, float]:
    """Load, anti-alias, resample, and band-pass one LFP trace."""
    raw = np.loadtxt(path, dtype=np.float64)
    if not np.isfinite(raw).all():
        raise ValueError(f"Non-finite samples in {path}")
    raw_fs = len(raw) / duration_s
    raw = raw - np.median(raw)
    ratio = Fraction(view_fs / raw_fs).limit_denominator(10000)
    down = signal.resample_poly(
        raw, ratio.numerator, ratio.denominator, padtype="line"
    )
    sos = signal.butter(4, band, btype="bandpass", fs=view_fs, output="sos")
    padlen = min(len(down) - 1, int(3 * view_fs / band[0]))
    eeg = signal.sosfiltfilt(
        sos, down, padtype="even", padlen=padlen
    )
    return eeg - np.median(eeg), view_fs, raw_fs


def bandpass(x: np.ndarray, fs: float, low: float, high: float) -> np.ndarray:
    """Zero-phase Butterworth band-pass, matching the loader's conventions."""
    sos = signal.butter(4, (low, high), btype="bandpass", fs=fs, output="sos")
    padlen = min(len(x) - 1, int(3 * fs / low))
    return signal.sosfiltfilt(sos, x, padtype="even", padlen=padlen)


def rhythmicity_metrics(
    x: np.ndarray,
    fs: float,
    band: tuple[float, float] = (0.5, 4.0),
    segment_s: float = 20.0,
    max_lag_s: float = 5.0,
) -> dict[str, float]:
    """Quantify how clock-like the slow activity is.

    GRDA is described clinically as more rhythmic than the slow activity of N3.
    The distinction is generator-level: a thalamic Ih/I_T delta oscillation is
    narrowband and repeats at a near-constant interval, whereas the cortical
    slow oscillation is a bistable Up/Down process whose cycle length varies
    widely. Each metric below is scale-free, so it transfers between the model
    LFP and a scalp recording unchanged.

    The delta peak uses a longer Welch segment than the band-power summary
    because peak width is the measurement: at the 5.12 s segment used for band
    power a sharp 2 Hz peak occupies about two frequency bins, which makes any
    width estimate resolution-limited rather than signal-limited.
    """
    xd = bandpass(x, fs, band[0], band[1])

    nperseg = min(len(x), int(round(segment_s * fs)))
    freq, power = signal.welch(
        x, fs=fs, nperseg=nperseg, noverlap=nperseg // 2, detrend="constant"
    )
    # Lightly smooth before measuring peak width. On an unsmoothed estimate the
    # half-power walk below stops at the first noise dip, which reports a
    # narrow peak -- and so a high Q -- for broadband signals that have no
    # rhythm at all.
    kernel = np.ones(3) / 3.0
    power = np.convolve(power, kernel, mode="same")
    in_band = (freq >= band[0]) & (freq <= band[1])
    wide = (freq >= 0.5) & (freq <= 30.0)
    peak_index = int(np.argmax(np.where(in_band, power, -np.inf)))
    peak_hz = float(freq[peak_index])
    peak_power = float(power[peak_index])
    background = float(np.median(power[wide]))

    # Half-power width of the peak, walked outward until the spectrum first
    # drops below half the peak. Q = centre frequency / width, as for a filter:
    # a high Q means the energy sits in a narrow band, which is what makes a
    # trace look monomorphic.
    half = peak_power / 2.0
    lo = peak_index
    while lo > 0 and power[lo] > half:
        lo -= 1
    hi = peak_index
    while hi < len(power) - 1 and power[hi] > half:
        hi += 1
    width_hz = float(freq[hi] - freq[lo])
    peak_q = float(peak_hz / width_hz) if width_hz > 0 else float("nan")

    # Autocorrelation of the band-limited signal. A clock keeps returning to
    # itself, so side lobes persist over many cycles; an irregular oscillation
    # decorrelates after roughly one.
    centered = xd - np.mean(xd)
    acf_full = signal.correlate(centered, centered, mode="full", method="fft")
    acf = acf_full[len(centered) - 1 :]
    acf = acf / acf[0]
    max_lag = min(len(acf) - 1, int(round(max_lag_s * fs)))
    acf = acf[: max_lag + 1]
    # Ignore the central lobe: no peak can precede half a cycle of the band's
    # upper edge.
    min_lag = int(round(fs / (2.0 * band[1])))
    peaks, _ = signal.find_peaks(acf[min_lag:])
    peaks = peaks + min_lag
    if len(peaks):
        first_peak = float(acf[peaks[0]])
        first_peak_lag_s = float(peaks[0] / fs)
        sidelobes = int(np.sum(acf[peaks] > 0.30))
    else:
        first_peak = 0.0
        first_peak_lag_s = float("nan")
        sidelobes = 0

    # Cycle-to-cycle interval variability from positive-going zero crossings.
    crossings = np.flatnonzero((xd[:-1] < 0) & (xd[1:] >= 0))
    if len(crossings) > 2:
        intervals = np.diff(crossings) / fs
        cycle_cv = float(np.std(intervals) / np.mean(intervals))
        cycle_hz = float(1.0 / np.mean(intervals))
    else:
        cycle_cv = float("nan")
        cycle_hz = float("nan")

    # Instantaneous frequency, restricted to samples where the band-limited
    # envelope is above its own median. Phase is meaningless where there is no
    # oscillation to track, and including those samples would report the noise
    # floor rather than the rhythm.
    analytic = signal.hilbert(xd)
    envelope = np.abs(analytic)
    inst_freq = np.diff(np.unwrap(np.angle(analytic))) * fs / (2.0 * np.pi)
    strong = envelope[:-1] > np.median(envelope)
    valid = strong & (inst_freq > band[0]) & (inst_freq < band[1])
    inst_freq_sd = float(np.std(inst_freq[valid])) if np.any(valid) else float("nan")

    return {
        "delta_peak_hz": peak_hz,
        "delta_peak_q": peak_q,
        "delta_peak_prominence": float(peak_power / background)
        if background > 0
        else float("nan"),
        "acf_first_peak": first_peak,
        "acf_first_peak_lag_s": first_peak_lag_s,
        "acf_sidelobes": sidelobes,
        "cycle_hz": cycle_hz,
        "cycle_cv": cycle_cv,
        "inst_freq_sd_hz": inst_freq_sd,
    }


def spindle_metrics(
    x: np.ndarray,
    fs: float,
    band: tuple[float, float] = (7.0, 15.0),
    smooth_s: float = 0.2,
    high_factor: float = 2.5,
    low_factor: float = 1.5,
    min_dur_s: float = 0.4,
    max_dur_s: float = 3.0,
) -> dict[str, float]:
    """Detect spindles and describe the sigma envelope.

    Model LFP is in arbitrary units, so the detector thresholds each run
    against its own median sigma envelope rather than against a fixed
    microvolt criterion. That makes the count comparable across states with
    very different delta amplitudes, but it also means a spindle-free run can
    still register occasional events from background fluctuations.

    ``sigma_env_cv`` is the threshold-free companion measure and should carry
    the weight when the two disagree: spindles make the sigma envelope bursty
    (high CV) whereas a spindle-free background of any amplitude keeps it
    steady (low CV).
    """
    band = (band[0], min(band[1], 0.45 * fs))
    sigma = bandpass(x, fs, band[0], band[1])
    envelope = np.abs(signal.hilbert(sigma))
    window = max(1, int(round(smooth_s * fs)))
    envelope = np.convolve(envelope, np.ones(window) / window, mode="same")

    median_env = float(np.median(envelope))
    duration_s = len(x) / fs
    if median_env <= 0:
        return {
            "sigma_pct_of_broadband": 0.0,
            "sigma_env_cv": float("nan"),
            "spindle_count": 0,
            "spindle_per_min": 0.0,
            "spindle_mean_dur_s": float("nan"),
            "spindle_peak_ratio": float("nan"),
        }

    high = high_factor * median_env
    low = low_factor * median_env
    above_low = envelope > low
    # Candidate events are runs above the lower bound; each must also reach the
    # upper bound somewhere inside, which is the usual two-threshold scheme for
    # keeping event boundaries at the envelope's foot rather than its shoulder.
    edges = np.diff(above_low.astype(np.int8))
    starts = list(np.flatnonzero(edges == 1) + 1)
    ends = list(np.flatnonzero(edges == -1) + 1)
    if above_low[0]:
        starts.insert(0, 0)
    if above_low[-1]:
        ends.append(len(envelope))

    durations = []
    peaks = []
    for start, end in zip(starts, ends):
        dur = (end - start) / fs
        if not (min_dur_s <= dur <= max_dur_s):
            continue
        segment_peak = float(np.max(envelope[start:end]))
        if segment_peak < high:
            continue
        durations.append(dur)
        peaks.append(segment_peak)

    sigma_power = float(np.mean(sigma * sigma))
    broadband_power = float(np.mean(x * x))
    return {
        "sigma_pct_of_broadband": 100.0 * sigma_power / broadband_power
        if broadband_power > 0
        else float("nan"),
        "sigma_env_cv": float(np.std(envelope) / np.mean(envelope)),
        "spindle_count": len(durations),
        "spindle_per_min": 60.0 * len(durations) / duration_s,
        "spindle_mean_dur_s": float(np.mean(durations)) if durations else float("nan"),
        "spindle_peak_ratio": float(np.mean(peaks) / median_env)
        if peaks
        else float("nan"),
    }


def compute_metrics(
    eeg: np.ndarray,
    fs: float,
    start_s: float = 10.0,
    end_s: float = 119.0,
) -> dict[str, float]:
    """Compute amplitude and spectral metrics on the declared epoch."""
    start = int(round(start_s * fs))
    end = min(len(eeg), int(round(end_s * fs)))
    x = eeg[start:end]
    if len(x) < int(5.12 * fs):
        raise ValueError("Analysis epoch is too short")

    nperseg = min(len(x), int(round(5.12 * fs)))
    freq, power = signal.welch(
        x,
        fs=fs,
        nperseg=nperseg,
        noverlap=nperseg // 2,
        detrend="constant",
    )
    use = (freq >= 0.5) & (freq <= 30.0)
    f = freq[use]
    p = power[use]
    total = np.trapezoid(p, f)

    def band_percent(low: float, high: float) -> float:
        mask = (f >= low) & (f <= high)
        return 100.0 * np.trapezoid(p[mask], f[mask]) / total

    delta_pct = band_percent(0.5, 4.0)
    theta_pct = band_percent(4.0, 8.0)
    alpha_pct = band_percent(8.0, 13.0)
    beta_pct = band_percent(13.0, 30.0)

    # Spectral edge: the frequency below which 95% of the power falls.
    cumulative = np.concatenate([[0.0], np.cumsum(np.diff(f) * (p[:-1] + p[1:]) / 2.0)])
    sef95 = float(np.interp(0.95 * cumulative[-1], cumulative, f))

    metrics = {
        "epoch_start_s": start_s,
        "epoch_end_s": end_s,
        "rms_au": float(np.sqrt(np.mean(x * x))),
        "p99_9_abs_au": float(np.percentile(np.abs(x), 99.9)),
        "peak_abs_au": float(np.max(np.abs(x))),
        "delta_pct": delta_pct,
        "theta_pct": theta_pct,
        "alpha_pct": alpha_pct,
        "beta_pct": beta_pct,
        "dominant_hz": float(f[np.argmax(p)]),
        "spectral_centroid_hz": float(np.trapezoid(f * p, f) / total),
        "excess_kurtosis": float(stats.kurtosis(x)),
        "theta_delta_ratio": theta_pct / delta_pct if delta_pct > 0 else float("nan"),
        "alpha_delta_ratio": alpha_pct / delta_pct if delta_pct > 0 else float("nan"),
        "beta_delta_ratio": beta_pct / delta_pct if delta_pct > 0 else float("nan"),
        "fast_slow_ratio": (theta_pct + alpha_pct + beta_pct) / delta_pct
        if delta_pct > 0
        else float("nan"),
        "sef95_hz": sef95,
    }
    metrics.update(rhythmicity_metrics(x, fs))
    metrics.update(spindle_metrics(x, fs))
    return metrics


def analyze_run(
    root: Path,
    state: int,
    seed: int,
    nhost: int,
    duration_s: float,
    view_fs: float,
    band: tuple[float, float],
    start_s: float,
    end_s: float,
) -> tuple[dict[str, object], np.ndarray, float]:
    source = resolve_lfp_file(root, state, seed, nhost)
    eeg, fs, raw_fs = load_eegband(source, duration_s, view_fs, band)
    row: dict[str, object] = {
        "state": state,
        "seed": seed,
        "source_file": source.name,
        "samples_raw": int(round(raw_fs * duration_s)),
        "duration_s": duration_s,
        "raw_fs_hz": raw_fs,
        "view_fs_hz": fs,
        "band_low_hz": band[0],
        "band_high_hz": band[1],
    }
    row.update(compute_metrics(eeg, fs, start_s, end_s))
    return row, eeg, fs


def write_csv(rows: list[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("states", nargs="+", type=int)
    parser.add_argument("--seeds", nargs="+", type=int, default=[1])
    parser.add_argument("--nhost", type=int, default=8)
    parser.add_argument("--duration-s", type=float, default=120.0)
    parser.add_argument("--view-fs", type=float, default=200.0)
    parser.add_argument("--band", nargs=2, type=float, default=(0.5, 30.0))
    parser.add_argument("--start-s", type=float, default=10.0)
    parser.add_argument("--end-s", type=float, default=119.0)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/metrics/lfp_state_metrics.csv"),
    )
    args = parser.parse_args()

    rows = []
    for state in args.states:
        for seed in args.seeds:
            row, _, _ = analyze_run(
                args.root,
                state,
                seed,
                args.nhost,
                args.duration_s,
                args.view_fs,
                tuple(args.band),
                args.start_s,
                args.end_s,
            )
            rows.append(row)
    write_csv(rows, args.output)
    print(f"wrote {args.output} ({len(rows)} runs)")


if __name__ == "__main__":
    main()
