"""Compare states across the GRDA-vs-N3 discriminator metrics.

One row per state: the band-limited LFP, its power spectrum, its
autocorrelation, and the discriminator panel. The autocorrelation column is
the one that matters most for this project, because "more rhythmic than N3"
is the defining feature of GRDA and side-lobes that persist over many cycles
are what a clock looks like.

The bottom strip carries the network-state metrics from the raster, which is
where cortical bistability shows up. Two states can produce similar-looking
slow LFP and still differ completely there: an N3 slow oscillation silences
cortex in every cycle, a thalamically driven delta rhythm does not.

Usage:
    .venv/bin/python plot_state_comparison.py --states 0 1 2 42 --nhost 8
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

from analyze_lfp_states import analyze_run, bandpass
import analyze_network_state

STATE_LABELS = {
    0: "Wake",
    1: "N2",
    2: "N3",
    3: "REM",
    35: "state 35 (irregular slowing)",
    40: "state 40",
    42: "state 42 (diffuse-slowing candidate)",
    44: "state 44",
}

ACCENT = "#123b69"


def state_label(state: int) -> str:
    return STATE_LABELS.get(state, f"state {state}")


def build(states: list[int], seed: int, nhost: int, root: Path, output: Path) -> None:
    available = []
    for state in states:
        try:
            row, eeg, fs = analyze_run(
                root, state, seed, nhost, 120.0, 200.0, (0.5, 30.0), 10.0, 119.0
            )
        except FileNotFoundError:
            print(f"skipping state {state}: no LFP output yet")
            continue
        try:
            net = analyze_network_state.raster_metrics(
                analyze_network_state.resolve_raster_file(root, state, seed, nhost)
            )
        except FileNotFoundError:
            net = {}
        available.append((state, row, eeg, fs, net))

    if not available:
        raise SystemExit("no state outputs found")

    nrows = len(available)
    fig, axes = plt.subplots(
        nrows, 4, figsize=(19.5, 3.05 * nrows), facecolor="white", squeeze=False
    )

    for r, (state, row, eeg, fs, net) in enumerate(available):
        time = np.arange(len(eeg)) / fs
        window = (time >= 100) & (time < 115)

        ax = axes[r][0]
        ax.plot(time[window] - 100, eeg[window], color=ACCENT, linewidth=0.7)
        ax.set_ylabel(f"{state_label(state)}\nLFP (a.u.)", fontsize=10, fontweight="bold")
        ax.set_xlim(0, 15)
        ax.grid(color="#eef1f4", linewidth=0.6)
        if r == 0:
            ax.set_title("Band-limited LFP, 15 s", fontsize=11, fontweight="bold")
        if r == nrows - 1:
            ax.set_xlabel("Time (s)")

        ax = axes[r][1]
        x = eeg[int(10 * fs) : int(119 * fs)]
        nperseg = min(len(x), int(round(20.0 * fs)))
        freq, power = signal.welch(x, fs=fs, nperseg=nperseg, noverlap=nperseg // 2)
        show = (freq >= 0.3) & (freq <= 25)
        ax.semilogy(freq[show], power[show], color=ACCENT, linewidth=1.1)
        ax.axvspan(0.5, 4, color="#c0392b", alpha=0.10)
        ax.axvspan(7, 15, color="#1f6f4a", alpha=0.10)
        ax.set_xlim(0, 25)
        ax.grid(color="#eef1f4", linewidth=0.6)
        if r == 0:
            ax.set_title(
                "Power spectrum (delta red, spindle green)", fontsize=11, fontweight="bold"
            )
        if r == nrows - 1:
            ax.set_xlabel("Frequency (Hz)")

        ax = axes[r][2]
        xd = bandpass(x, fs, 0.5, 4.0)
        centered = xd - np.mean(xd)
        acf = signal.correlate(centered, centered, mode="full", method="fft")
        acf = acf[len(centered) - 1 :]
        acf = acf / acf[0]
        lags = np.arange(len(acf)) / fs
        keep = lags <= 5.0
        ax.plot(lags[keep], acf[keep], color=ACCENT, linewidth=1.0)
        ax.axhline(0, color="#98a4ae", linewidth=0.7)
        ax.axhline(0.30, color="#c0392b", linewidth=0.8, linestyle="--")
        ax.set_ylim(-1.05, 1.05)
        ax.grid(color="#eef1f4", linewidth=0.6)
        ax.text(
            0.97,
            0.93,
            f"side-lobes>0.3: {int(row['acf_sidelobes'])}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9.5,
            fontweight="bold",
            color="#c0392b",
        )
        if r == 0:
            ax.set_title(
                "Autocorrelation of 0.5-4 Hz (rhythmicity)", fontsize=11, fontweight="bold"
            )
        if r == nrows - 1:
            ax.set_xlabel("Lag (s)")

        ax = axes[r][3]
        ax.axis("off")
        lines = [
            ("delta power", f"{row['delta_pct']:.1f} %"),
            ("fast/slow ratio", f"{row['fast_slow_ratio']:.2f}"),
            ("spectral centroid", f"{row['spectral_centroid_hz']:.2f} Hz"),
            ("delta peak", f"{row['delta_peak_hz']:.2f} Hz  (Q {row['delta_peak_q']:.1f})"),
            ("cycle CV", f"{row['cycle_cv']:.2f}"),
            ("spindles/min", f"{row['spindle_per_min']:.1f}"),
            ("sigma envelope CV", f"{row['sigma_env_cv']:.2f}"),
        ]
        if net:
            lines += [
                ("PY rate", f"{net['py_rate_hz']:.2f} Hz"),
                ("Down-state time", f"{net['down_state_pct']:.1f} %"),
                ("PY rate bimodality", f"{net['py_rate_bimodality']:.2f}"),
                ("TC spikes in bursts", f"{net['tc_burst_pct']:.0f} %"),
            ]
        for i, (name, value) in enumerate(lines):
            y = 0.96 - i * 0.088
            ax.text(0.02, y, name, fontsize=9.6, color="#4a545c", transform=ax.transAxes)
            ax.text(
                0.98,
                y,
                value,
                fontsize=9.6,
                fontweight="bold",
                ha="right",
                color="#16324a",
                transform=ax.transAxes,
            )
        if r == 0:
            ax.set_title("Discriminator metrics", fontsize=11, fontweight="bold")

    fig.suptitle(
        "Discriminating GRDA-like activity from the N3 slow oscillation",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.005,
        0.004,
        "Analysis epoch 10-119 s. GRDA should show high delta power AND many "
        "autocorrelation side-lobes AND near-zero spindles AND a low Down-state "
        "fraction. N3 shows high delta power with few side-lobes and a large "
        "Down-state fraction, because its slow oscillation is a bistable "
        "cortical Up/Down process rather than a clock.",
        fontsize=8.4,
        color="#4a545c",
    )
    fig.tight_layout(rect=(0, 0.022, 1, 0.975))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170, facecolor="white")
    plt.close(fig)
    print(f"wrote {output} ({len(available)} states)")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--states", nargs="+", type=int, default=[0, 1, 2, 42])
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--nhost", type=int, default=8)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument(
        "--output", type=Path, default=Path("docs/figures/state_discriminators.png")
    )
    args = parser.parse_args()
    build(args.states, args.seed, args.nhost, args.root, args.output)


if __name__ == "__main__":
    main()
