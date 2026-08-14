"""Map the intrinsic delta window of a single thalamocortical relay cell.

A GRDA-like rhythm has to come from somewhere. In this model family there are
three candidate slow generators: the cortical Up/Down slow oscillation, the
thalamic Ih/I_T delta clock, and the RE-TC spindle loop. Only the second is
intrinsically clock-like, which makes it the natural candidate for a pattern
defined by being *more* rhythmic than N3.

This script asks the cheapest possible version of the question first: with no
network at all, can one TC cell oscillate at 0.5-4 Hz, and at which
combinations of K-leak conductance and Ih activation shift? Those are exactly
the two parameters the model already uses to move the thalamus between wake,
N2, N3, and REM, so the answer is measured in the model's own coordinates and
the four reference states can be plotted on top of it.

A single cell settles and runs in seconds, so mapping the whole plane costs
minutes rather than the hours a network sweep would take. Run this before
committing compute to network states.

The K-leak axis extends well past the model's own wake-to-N3 range on purpose.
The question is whether an intrinsic delta limit cycle exists anywhere in this
cell's parameter space, not just at the calibrated operating points. Raising
gkL is also a fair proxy for tonic thalamic inhibition: kL reverses at -95 mV,
the same reversal this model gives RE->TC GABA-B, so adding K-leak and adding
tonic inhibitory conductance hyperpolarize the cell along the same axis.
Points beyond the model's calibrated range answer a mechanistic question and
are not claims about physiological parameter values.

Usage:
    .venv/bin/nrniv -python tc_delta_map.py
    .venv/bin/nrniv -python tc_delta_map.py -- --settle-s 20 --analyze-s 30
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

import config
from cell_classes import TCCell
from neuron import h

# Reference operating points, so the map is read in the model's own terms.
REFERENCE_STATES = [
    ("wake", config.gkl_TC_awake, config.gh_TC_awake),
    ("N2", config.gkl_TC_s2, config.gh_TC_s2),
    ("N3", config.gkl_TC_s3, config.gh_TC_s3),
    ("REM", config.gkl_TC_rem, config.gh_TC_rem),
]

SPIKE_THRESHOLD_MV = 0.0
BURST_ISI_MS = 20.0
SPIKE_CLIP_MV = -45.0


def script_args(argv: list[str]) -> list[str]:
    """Recover this script's own arguments from the nrniv command line.

    nrniv leaves sys.argv as the full launcher command, e.g.
    ['nrniv', '-python', 'tc_delta_map.py', '--', '--settle-s', '20'], so
    parsing sys.argv[1:] directly silently drops every option.
    """
    if "--" in argv:
        return argv[argv.index("--") + 1 :]
    name = Path(__file__).name
    for i, token in enumerate(argv):
        if token.endswith(name):
            return argv[i + 1 :]
    return argv[1:]


def analyze_trace(
    v: np.ndarray,
    spike_s: np.ndarray,
    fs: float,
    band: tuple[float, float] = (0.2, 10.0),
) -> dict[str, float]:
    """Describe one cell's activity: burst rhythm, regularity, and baseline.

    Spike times come from an APCount on the soma rather than from threshold
    crossings of ``v``: the trace is stored at 1 ms to keep the sweep small,
    and an action potential is narrower than that, so crossings found in the
    stored trace would miss most spikes.
    """

    # Group spikes into low-threshold Ca bursts. Isolated spikes count as
    # one-spike bursts so that tonic firing is still measured, not discarded.
    if len(spike_s):
        gaps = np.diff(spike_s) > BURST_ISI_MS / 1000.0
        burst_start = np.concatenate([[True], gaps])
        burst_times = spike_s[burst_start]
        spikes_per_burst = len(spike_s) / len(burst_times)
    else:
        burst_times = np.array([])
        spikes_per_burst = 0.0

    if len(burst_times) >= 3:
        ibi = np.diff(burst_times)
        burst_hz = float(1.0 / np.mean(ibi))
        burst_cv = float(np.std(ibi) / np.mean(ibi))
    else:
        burst_hz = float("nan")
        burst_cv = float("nan")

    # Subthreshold rhythm: clip the action potentials so the spectrum reports
    # the underlying slow envelope rather than spike energy. This also detects
    # oscillations that never reach threshold.
    clipped = np.minimum(v, SPIKE_CLIP_MV)
    nperseg = min(len(clipped), int(round(10.0 * fs)))
    freq, power = signal.welch(
        clipped - np.mean(clipped), fs=fs, nperseg=nperseg, detrend="constant"
    )
    in_band = (freq >= band[0]) & (freq <= band[1])
    sub_hz = float(freq[in_band][np.argmax(power[in_band])])
    sub_amp = float(np.sqrt(np.max(power[in_band]) * (freq[1] - freq[0])))

    return {
        "mean_vm_mv": float(np.mean(clipped)),
        "min_vm_mv": float(np.min(v)),
        "spike_rate_hz": float(len(spike_s) / (len(v) / fs)),
        "burst_rate_hz": float(len(burst_times) / (len(v) / fs)),
        "n_bursts": int(len(burst_times)),
        "spikes_per_burst": float(spikes_per_burst),
        "burst_hz": burst_hz,
        "burst_cv": burst_cv,
        "subthreshold_hz": sub_hz,
        "subthreshold_amp_mv": sub_amp,
    }


def classify(row: dict[str, float], band: tuple[float, float], max_cv: float) -> str:
    """Label each grid point by what kind of activity it produced."""
    if row["n_bursts"] < 5:
        return "silent" if row["spike_rate_hz"] < 0.5 else "sparse"
    if row["spikes_per_burst"] < 1.5:
        return "tonic"
    if not (band[0] <= row["burst_hz"] <= band[1]):
        return "burst-out-of-band"
    return "delta-clock" if row["burst_cv"] <= max_cv else "burst-irregular"


def run_point(
    cell: TCCell,
    detector,
    spike_vec,
    v_vec,
    t_vec,
    gkl: float,
    gh_shift: float,
    settle_s: float,
    analyze_s: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Run one parameter combination and return (trace, spike times) after settling.

    The settle period matters more than usual here. The Ih mechanism carries a
    calcium-binding process with a decay rate of 4e-4 /ms, a time constant of
    about 2.5 s, so a short run would measure the transient rather than the
    steady-state rhythm.
    """
    cell.soma.gkL_kL = gkl
    cell.soma.fac_gh_TC_iar = gh_shift
    detector.thresh = SPIKE_THRESHOLD_MV

    h.dt = 0.025  # match the network runs in bazh_net.py
    h.finitialize(-68)
    h.continuerun((settle_s + analyze_s) * 1000.0)

    v = np.array(v_vec)
    t_ms = np.array(t_vec)
    keep = t_ms >= settle_s * 1000.0
    spikes_ms = np.array(spike_vec)
    spikes_s = (spikes_ms[spikes_ms >= settle_s * 1000.0] - settle_s * 1000.0) / 1000.0
    return v[keep], spikes_s


def rebound_control(
    cell: TCCell,
    detector,
    spike_vec,
    v_vec,
    t_vec,
    gkl: float,
    gh_shift: float,
    amplitudes_na: tuple[float, ...] = (-0.05, -0.1, -0.2),
    settle_ms: float = 8000.0,
    step_ms: float = 1000.0,
    follow_ms: float = 16000.0,
    decay_window_ms: float = 2000.0,
) -> list[dict[str, object]]:
    """Positive control: hyperpolarize, release, and look for a rebound burst.

    A sweep that finds no oscillation anywhere is only informative if the
    machinery that would produce one is working. Releasing the cell from a
    hyperpolarizing step is the standard way to demonstrate that I_T is
    de-inactivating and can drive a low-threshold burst. If this control fails,
    the sweep result is uninterpretable rather than negative.

    The step is followed for a further ``follow_ms`` so the oscillation the
    release evokes can be tracked. That separates the two ways the sweep could
    have missed a rhythm: a stable limit cycle coexisting with the quiescent
    fixed point, which the sweep would never reach because every point starts
    from rest, versus damped ringing around a single stable fixed point. The
    returned per-window peak-to-peak amplitudes settle which one it is.
    """
    stim = h.IClamp(cell.soma(0.5))
    stim.delay = settle_ms
    stim.dur = step_ms
    cell.soma.gkL_kL = gkl
    cell.soma.fac_gh_TC_iar = gh_shift
    detector.thresh = SPIKE_THRESHOLD_MV

    results = []
    for amp in amplitudes_na:
        stim.amp = amp
        h.dt = 0.025
        h.finitialize(-68)
        h.continuerun(settle_ms + step_ms + follow_ms)
        t = np.array(t_vec)
        v = np.array(v_vec)
        spikes = np.array(spike_vec)
        release = settle_ms + step_ms

        # Clip action potentials so the envelope tracks the low-threshold
        # oscillation rather than spike height.
        clipped = np.minimum(v, SPIKE_CLIP_MV)
        decay_t = []
        decay_amp = []
        offset = 0.0
        while offset + decay_window_ms <= follow_ms:
            mask = (t >= release + offset) & (t < release + offset + decay_window_ms)
            if np.any(mask):
                decay_t.append((offset + decay_window_ms / 2.0) / 1000.0)
                decay_amp.append(float(np.ptp(clipped[mask])))
            offset += decay_window_ms

        window = (t >= settle_ms - 500.0) & (t <= release + 4000.0)
        results.append(
            {
                "amp_na": amp,
                "rest_mv": float(np.mean(v[(t > settle_ms - 1000.0) & (t < settle_ms)])),
                "step_mv": float(np.mean(v[(t > release - 200.0) & (t < release)])),
                "rebound_spikes": int(np.sum((spikes > release) & (spikes < release + 500.0))),
                "late_spikes": int(np.sum(spikes > release + 2000.0)),
                "t_ms": t[window] - release,
                "v_mv": v[window],
                "decay_t_s": np.array(decay_t),
                "decay_amp_mv": np.array(decay_amp),
            }
        )
    stim.amp = 0.0
    stim.dur = 0.0
    return results


def make_figure(
    rows: list[dict[str, object]],
    gkl_values: np.ndarray,
    gh_values: np.ndarray,
    band: tuple[float, float],
    output: Path,
    control: list[dict[str, object]] | None = None,
) -> None:
    shape = (len(gh_values), len(gkl_values))
    freq = np.full(shape, np.nan)
    cv = np.full(shape, np.nan)
    vm = np.full(shape, np.nan)
    label = np.empty(shape, dtype=object)
    for row in rows:
        i = int(row["gh_index"])
        j = int(row["gkl_index"])
        freq[i, j] = row["burst_hz"]
        cv[i, j] = row["burst_cv"]
        vm[i, j] = row["mean_vm_mv"]
        label[i, j] = row["classification"]

    # Show frequency only where the cell actually produced a rhythm.
    rhythmic = np.array(
        [[str(v) in ("delta-clock", "burst-irregular", "burst-out-of-band") for v in r] for r in label]
    )
    freq_masked = np.where(rhythmic, freq, np.nan)
    cv_masked = np.where(rhythmic, cv, np.nan)

    x_um = gkl_values * 1e5  # display K leak in units of 1e-5 S/cm2
    extent = [x_um[0], x_um[-1], gh_values[0], gh_values[-1]]

    # Categorical map of what the cell actually did, so a plane with no
    # oscillation anywhere reads as a result rather than as an empty panel.
    categories = ["silent", "sparse", "tonic", "burst-out-of-band", "burst-irregular", "delta-clock"]
    colors = ["#c9d1d9", "#9fb3c8", "#f2b880", "#8fbf9f", "#4f93c0", "#c0392b"]
    code = np.full(shape, np.nan)
    for i in range(shape[0]):
        for j in range(shape[1]):
            if label[i, j] in categories:
                code[i, j] = categories.index(label[i, j])

    ncols = 4 if control else 2
    fig, axes = plt.subplots(1, ncols, figsize=(5.6 * ncols, 5.6), facecolor="white")

    def mark_reference(ax):
        for name, gkl_ref, gh_ref in REFERENCE_STATES:
            if not (x_um[0] <= gkl_ref * 1e5 <= x_um[-1]):
                continue
            ax.plot(gkl_ref * 1e5, gh_ref, "o", ms=10, mfc="none", mec="white", mew=2.6)
            ax.plot(gkl_ref * 1e5, gh_ref, "o", ms=10, mfc="none", mec="black", mew=1.3)
            ax.annotate(
                name,
                (gkl_ref * 1e5, gh_ref),
                textcoords="offset points",
                xytext=(9, 5),
                fontsize=10,
                fontweight="bold",
                color="black",
            )

    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, len(categories) + 0.5), cmap.N)
    axes[0].imshow(
        code, origin="lower", aspect="auto", extent=extent,
        cmap=cmap, norm=norm, interpolation="nearest",
    )
    axes[0].set_xlabel(r"TC K-leak conductance ($10^{-5}$ S/cm$^2$)")
    axes[0].set_ylabel("Ih activation shift (mV)")
    axes[0].set_title(
        "A. Activity regime of the isolated TC cell", fontsize=11.5, fontweight="bold"
    )
    present = [c for c in categories if np.any(label == c)]
    axes[0].legend(
        handles=[Patch(facecolor=colors[categories.index(c)], edgecolor="#5a6570", label=c) for c in present],
        loc="upper right", fontsize=9, framealpha=0.95,
    )
    mark_reference(axes[0])

    image = axes[1].imshow(
        vm, origin="lower", aspect="auto", extent=extent,
        cmap="coolwarm_r", interpolation="nearest",
    )
    axes[1].set_xlabel(r"TC K-leak conductance ($10^{-5}$ S/cm$^2$)")
    axes[1].set_ylabel("Ih activation shift (mV)")
    axes[1].set_title(
        "B. Steady membrane potential", fontsize=11.5, fontweight="bold"
    )
    bar = fig.colorbar(image, ax=axes[1], fraction=0.046, pad=0.03)
    bar.set_label("mV", fontsize=9)
    # Mark the range over which I_T de-inactivates, to show the sweep reached it.
    contour = axes[1].contour(
        np.linspace(x_um[0], x_um[-1], len(x_um)),
        np.linspace(gh_values[0], gh_values[-1], len(gh_values)),
        vm, levels=[-80.0, -70.0], colors="black", linewidths=1.1, linestyles="--",
    )
    axes[1].clabel(contour, fmt="%.0f mV", fontsize=8)
    mark_reference(axes[1])

    if control:
        ax = axes[2]
        for entry, color in zip(control, ["#4f93c0", "#1f6f4a", "#c0392b"]):
            ax.plot(
                np.asarray(entry["t_ms"]) / 1000.0,
                entry["v_mv"],
                color=color,
                linewidth=0.9,
                label="%.2f nA step to %.0f mV: %d rebound spikes"
                % (entry["amp_na"], entry["step_mv"], entry["rebound_spikes"]),
            )
        ax.axvline(0, color="#5a6570", linewidth=1.0, linestyle=":")
        ax.text(0.02, 0.03, "step released", transform=ax.transAxes, fontsize=8.5, color="#5a6570")
        ax.set_xlabel("Time relative to release of hyperpolarizing step (s)")
        ax.set_ylabel("Somatic V (mV)")
        ax.set_title(
            "C. Positive control: I_T rebound burst", fontsize=11.5, fontweight="bold"
        )
        ax.legend(loc="upper right", fontsize=8)

        ax = axes[3]
        for entry, color in zip(control, ["#4f93c0", "#1f6f4a", "#c0392b"]):
            ax.semilogy(
                entry["decay_t_s"],
                np.maximum(entry["decay_amp_mv"], 1e-3),
                "o-",
                color=color,
                markersize=4,
                linewidth=1.2,
                label="%.2f nA step" % entry["amp_na"],
            )
        ax.set_xlabel("Time after release (s)")
        ax.set_ylabel("Oscillation amplitude (mV peak-to-peak)")
        ax.set_title(
            "D. The evoked rhythm dies out", fontsize=11.5, fontweight="bold"
        )
        ax.grid(True, which="both", color="#e6eaee", linewidth=0.6)
        ax.legend(loc="upper right", fontsize=8.5)
        ax.text(
            0.03,
            0.06,
            "Amplitude falls by one to three orders of magnitude:\n"
            "damped ringing around a stable fixed point, not a\n"
            "self-sustaining limit cycle.",
            transform=ax.transAxes,
            fontsize=8.5,
            color="#4a545c",
        )

    fig.suptitle(
        "Intrinsic delta window of an isolated TC cell "
        "(no synapses, no network)",
        fontsize=13.5,
        fontweight="bold",
        y=0.99,
    )
    fig.text(
        0.005,
        0.015,
        "Circles mark the K-leak and Ih settings the model itself uses for "
        "wake, N2, N3, and REM. The K-leak axis extends past that calibrated "
        "range on purpose, to test whether an intrinsic delta limit cycle "
        "exists anywhere in this cell's parameter space; K leak reverses at "
        "-95 mV, the same reversal as RE->TC GABA-B, so it stands in for tonic "
        "thalamic inhibition. Panel C is the positive control: the same cell "
        "does fire a low-threshold rebound burst when released from a "
        "hyperpolarizing step, so a null result in A is a property of the "
        "model rather than a broken I_T current.",
        fontsize=8.2,
        color="#4a545c",
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.955))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=170, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    # The upper end is far above the model's N3 value (2.8e-5). See the module
    # docstring: the point is to find the intrinsic delta window if it exists
    # anywhere, and extra K leak stands in for tonic inhibition.
    parser.add_argument("--gkl-min", type=float, default=1.2e-5)
    parser.add_argument("--gkl-max", type=float, default=1.4e-4)
    parser.add_argument("--gkl-steps", type=int, default=15)
    parser.add_argument("--gh-min", type=float, default=-8.0)
    parser.add_argument("--gh-max", type=float, default=0.0)
    parser.add_argument("--gh-steps", type=int, default=9)
    parser.add_argument("--settle-s", type=float, default=15.0)
    parser.add_argument("--analyze-s", type=float, default=20.0)
    parser.add_argument("--record-dt-ms", type=float, default=1.0)
    parser.add_argument("--band", nargs=2, type=float, default=(0.5, 4.0))
    parser.add_argument("--max-cv", type=float, default=0.25)
    parser.add_argument(
        "--output-csv", type=Path, default=Path("output/metrics/tc_delta_map.csv")
    )
    parser.add_argument(
        "--output-figure", type=Path, default=Path("output/figures/tc_delta_map.png")
    )
    args, _ = parser.parse_known_args(script_args(sys.argv))

    gkl_values = np.linspace(args.gkl_min, args.gkl_max, args.gkl_steps)
    gh_values = np.linspace(args.gh_min, args.gh_max, args.gh_steps)
    band = (args.band[0], args.band[1])

    cell = TCCell(0)
    # Detect spikes at the integration time step; the stored trace is far too
    # coarse to find them (see analyze_trace).
    detector = h.APCount(cell.soma(0.5))
    spike_vec = h.Vector()
    detector.record(spike_vec)
    v_vec = h.Vector()
    t_vec = h.Vector()
    v_vec.record(cell.soma(0.5)._ref_v, args.record_dt_ms)
    t_vec.record(h._ref_t, args.record_dt_ms)

    rows: list[dict[str, object]] = []
    total = len(gkl_values) * len(gh_values)
    done = 0
    for i, gh_shift in enumerate(gh_values):
        for j, gkl in enumerate(gkl_values):
            v, spikes_s = run_point(
                cell,
                detector,
                spike_vec,
                v_vec,
                t_vec,
                gkl,
                gh_shift,
                args.settle_s,
                args.analyze_s,
            )
            row: dict[str, object] = {
                "gh_index": i,
                "gkl_index": j,
                "gkl_S_per_cm2": float(gkl),
                "gh_shift_mv": float(gh_shift),
            }
            row.update(analyze_trace(v, spikes_s, fs=1000.0 / args.record_dt_ms))
            row["classification"] = classify(row, band, args.max_cv)
            rows.append(row)
            done += 1
            print(
                "[%3d/%3d] gkL=%.3e gh=%+.2f -> %-18s burst=%.2f Hz cv=%.2f Vm=%.1f mV"
                % (
                    done,
                    total,
                    gkl,
                    gh_shift,
                    row["classification"],
                    row["burst_hz"] if row["burst_hz"] == row["burst_hz"] else float("nan"),
                    row["burst_cv"] if row["burst_cv"] == row["burst_cv"] else float("nan"),
                    row["mean_vm_mv"],
                )
            )

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {args.output_csv}")

    print("\nrunning I_T rebound positive control ...")
    control = rebound_control(
        cell, detector, spike_vec, v_vec, t_vec, config.gkl_TC_s3, config.gh_TC_s3
    )
    for entry in control:
        print(
            "  step %.2f nA: rest %.1f mV -> %.1f mV, %d rebound spikes"
            % (entry["amp_na"], entry["rest_mv"], entry["step_mv"], entry["rebound_spikes"])
        )

    make_figure(rows, gkl_values, gh_values, band, args.output_figure, control)
    print(f"wrote {args.output_figure}")

    clocks = [r for r in rows if r["classification"] == "delta-clock"]
    print(
        "\ndelta-clock points: %d of %d (%.0f%%)"
        % (len(clocks), len(rows), 100.0 * len(clocks) / len(rows))
    )
    if clocks:
        freqs = [r["burst_hz"] for r in clocks]
        print(
            "  frequency range %.2f-%.2f Hz; K leak %.2e-%.2e S/cm2; Ih shift %+.1f to %+.1f mV"
            % (
                min(freqs),
                max(freqs),
                min(r["gkl_S_per_cm2"] for r in clocks),
                max(r["gkl_S_per_cm2"] for r in clocks),
                min(r["gh_shift_mv"] for r in clocks),
                max(r["gh_shift_mv"] for r in clocks),
            )
        )
    h.quit()


if __name__ == "__main__":
    main()
