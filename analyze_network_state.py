"""Network-state metrics computed from the spike raster.

These are the measurements the field potential cannot supply. An EEG shows
that activity is slow; it cannot show *which* generator is producing the
slowing. The raster can:

- ``down_state_pct`` and ``py_rate_bimodality`` measure cortical bistability.
  The N3 slow oscillation is an alternation between an active Up state and a
  silent Down state, so cortical firing is bimodal and stops outright for tens
  to hundreds of milliseconds at a time. A thalamically driven delta rhythm
  modulates cortical firing without silencing it, so the same LFP slowing
  should appear with a near-zero Down-state fraction.
- ``tc_burst_pct`` measures whether relay cells are in T-current burst mode,
  which is the signature of the hyperpolarized Ih/I_T delta clock.

Together these separate two states that can look similar in the LFP but arise
from different mechanisms, which is the point of running the model at all.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy import signal, stats


def resolve_raster_file(root: Path, state: int, seed: int, nhost: int) -> Path:
    """Resolve seed-tagged raster outputs, retaining seed-1 legacy support."""
    seeded = root / f"raster_state{state}_seed{seed}_nhost={nhost}.txt"
    if seeded.exists():
        return seeded
    legacy = root / f"raster_state{state}_nhost={nhost}.txt"
    if seed == 1 and legacy.exists():
        return legacy
    raise FileNotFoundError(f"No raster output found for state {state}, seed {seed}")


def load_raster(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load a raster file as (spike times in seconds, gids).

    Parallel runs write spikes grouped by rank, so the file is not sorted in
    time. Nothing here depends on the order.
    """
    data = np.loadtxt(path, dtype=np.float64)
    if data.size == 0:
        return np.array([]), np.array([], dtype=int)
    data = np.atleast_2d(data)
    return data[:, 0] / 1000.0, data[:, 1].astype(int)


def burst_fraction(times_s: np.ndarray, gids: np.ndarray, max_isi_s: float) -> float:
    """Fraction of spikes belonging to a high-frequency burst.

    A spike counts as part of a burst when it sits within ``max_isi_s`` of the
    previous or next spike from the same cell. Low-threshold Ca spikes carry
    2-6 action potentials at 200-400 Hz, well inside the default 10 ms bound,
    while tonic firing at the rates seen here is far outside it.
    """
    if len(times_s) == 0:
        return float("nan")
    in_burst = 0
    total = 0
    for gid in np.unique(gids):
        cell_times = np.sort(times_s[gids == gid])
        total += len(cell_times)
        if len(cell_times) < 2:
            continue
        isi = np.diff(cell_times)
        close = isi <= max_isi_s
        flagged = np.zeros(len(cell_times), dtype=bool)
        flagged[:-1] |= close
        flagged[1:] |= close
        in_burst += int(np.sum(flagged))
    return 100.0 * in_burst / total if total else float("nan")


def bimodality_coefficient(x: np.ndarray) -> float:
    """Sarle's bimodality coefficient; above about 0.555 suggests two modes."""
    n = len(x)
    if n < 4 or np.std(x) == 0:
        return float("nan")
    skew = stats.skew(x)
    kurt = stats.kurtosis(x)  # excess kurtosis
    correction = 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    return float((skew**2 + 1.0) / (kurt + correction))


def raster_metrics(
    path: Path,
    npyr: int = 500,
    ninh: int = 100,
    nre: int = 100,
    ntc: int = 100,
    start_s: float = 10.0,
    end_s: float = 119.0,
    bin_ms: float = 10.0,
    min_down_ms: float = 50.0,
    burst_isi_ms: float = 10.0,
    down_frac_of_median: float = 0.25,
) -> dict[str, float]:
    """Population rates, cortical bistability, and thalamic burst mode."""
    times_s, gids = load_raster(path)
    keep = (times_s >= start_s) & (times_s < end_s)
    times_s, gids = times_s[keep], gids[keep]
    span_s = end_s - start_s

    bounds = {
        "py": (0, npyr),
        "in": (npyr, npyr + ninh),
        "re": (npyr + ninh, npyr + ninh + nre),
        "tc": (npyr + ninh + nre, npyr + ninh + nre + ntc),
    }
    metrics: dict[str, float] = {}
    for name, (low, high) in bounds.items():
        count = int(np.sum((gids >= low) & (gids < high)))
        metrics[f"{name}_rate_hz"] = count / ((high - low) * span_s)

    for name in ("tc", "re"):
        low, high = bounds[name]
        mask = (gids >= low) & (gids < high)
        metrics[f"{name}_burst_pct"] = burst_fraction(
            times_s[mask], gids[mask], burst_isi_ms / 1000.0
        )

    # Cortical population rate on a fine grid, used for both the Down-state
    # detection and the bimodality test.
    py_low, py_high = bounds["py"]
    py_times = times_s[(gids >= py_low) & (gids < py_high)]
    n_bins = max(1, int(round(span_s * 1000.0 / bin_ms)))
    counts, _ = np.histogram(py_times, bins=n_bins, range=(start_s, end_s))
    rate_hz = counts / (npyr * bin_ms / 1000.0)

    # A Down state is a period of strongly suppressed cortical firing, not one
    # of absolute silence. Requiring every PY cell to be quiet at once is a far
    # stricter test than the physiology: the model's background miniature EPSPs
    # keep a handful of cells firing through a Down state, so an
    # all-cells-silent criterion reports 0% even for a network whose population
    # rate swings sixty-fold. The threshold is therefore relative to each run's
    # own median rate, which also keeps the measure comparable across states
    # with very different mean rates. The strict version is kept separately as
    # ``silent_pct``.
    min_bins = max(1, int(round(min_down_ms / bin_ms)))
    threshold = down_frac_of_median * float(np.median(rate_hz))
    quiet = rate_hz <= threshold

    def runs_above_min(flags: np.ndarray) -> list[float]:
        edges = np.diff(flags.astype(np.int8))
        starts = list(np.flatnonzero(edges == 1) + 1)
        ends = list(np.flatnonzero(edges == -1) + 1)
        if flags[0]:
            starts.insert(0, 0)
        if flags[-1]:
            ends.append(len(flags))
        return [
            (end - start) * bin_ms / 1000.0
            for start, end in zip(starts, ends)
            if (end - start) >= min_bins
        ]

    down_durations = runs_above_min(quiet)
    metrics["silent_pct"] = 100.0 * float(np.mean(counts == 0))
    metrics["down_threshold_hz"] = threshold
    metrics["down_state_pct"] = 100.0 * float(np.sum(down_durations)) / span_s
    metrics["down_state_per_min"] = 60.0 * len(down_durations) / span_s
    metrics["down_mean_dur_s"] = (
        float(np.mean(down_durations)) if down_durations else float("nan")
    )
    metrics["py_rate_bimodality"] = bimodality_coefficient(rate_hz)
    metrics["py_rate_cv"] = (
        float(np.std(rate_hz) / np.mean(rate_hz)) if np.mean(rate_hz) > 0 else float("nan")
    )

    # Frequency at which the cortical population rate itself is modulated. For
    # a thalamic delta clock this should land in the delta band even when the
    # cortex never falls silent.
    fs_rate = 1000.0 / bin_ms
    if len(rate_hz) > 16 and np.std(rate_hz) > 0:
        nperseg = min(len(rate_hz), int(round(20.0 * fs_rate)))
        freq, power = signal.welch(
            rate_hz - np.mean(rate_hz), fs=fs_rate, nperseg=nperseg, detrend="constant"
        )
        band = (freq >= 0.2) & (freq <= 30.0)
        metrics["py_rate_peak_hz"] = float(freq[band][np.argmax(power[band])])
    else:
        metrics["py_rate_peak_hz"] = float("nan")

    return metrics


def analyze_run(
    root: Path, state: int, seed: int, nhost: int, **kwargs
) -> dict[str, object]:
    source = resolve_raster_file(root, state, seed, nhost)
    row: dict[str, object] = {
        "state": state,
        "seed": seed,
        "raster_file": source.name,
    }
    row.update(raster_metrics(source, **kwargs))
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("states", nargs="+", type=int)
    parser.add_argument("--seeds", nargs="+", type=int, default=[1])
    parser.add_argument("--nhost", type=int, default=8)
    parser.add_argument("--start-s", type=float, default=10.0)
    parser.add_argument("--end-s", type=float, default=119.0)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--output", type=Path, default=Path("output/metrics/network_state_metrics.csv")
    )
    args = parser.parse_args()

    rows = [
        analyze_run(
            args.root,
            state,
            seed,
            args.nhost,
            start_s=args.start_s,
            end_s=args.end_s,
        )
        for state in args.states
        for seed in args.seeds
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {args.output} ({len(rows)} runs)")


if __name__ == "__main__":
    main()
