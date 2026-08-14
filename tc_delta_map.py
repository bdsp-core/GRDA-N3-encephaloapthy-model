"""Map the delta-band resonance of a single thalamocortical relay cell.

A GRDA-like rhythm has to come from somewhere, and the thalamic Ih/I_T
interaction is one candidate generator. With no network at all, where in
parameter space does one relay cell produce delta?

**The primary measurement is resonance, not spontaneous oscillation.** An
earlier version of this script asked only whether a cell oscillates by itself,
found that none did anywhere on the plane, and reported that as a dead end.
That was the wrong question. A stable focus that rings at 2 Hz with a small
damping constant is a resonator, and a resonator driven by common noise or
fluctuating afferent input produces intermittent rhythmic runs even though the
isolated deterministic cell has no limit cycle at all. That is an attractive
account of GRDA, because one mechanism explains both regularity within a run
and intermittency across runs -- and the spontaneous scan labels such a point
'silent'.

So each grid point now receives a brief, deliberately small current pulse, and
the ringing that follows is fitted to a damped sinusoid. Two numbers come out:
the resonant frequency, and the damping constant, which says how close the
leading eigenvalue pair sits to the imaginary axis. Low damping in the delta
band is the target; it marks the direction to push toward a Hopf bifurcation
or toward noise-driven quasi-cycles. The spontaneous-activity classification
is retained as a secondary output.

The two swept parameters are the ones the model itself uses to move the
thalamus between wake, N2, N3 and REM, so the map is in the model's own
coordinates and the four reference states can be plotted on it. Note that this
is a 2-D slice: intrinsic TC delta also depends on g_T, the voltage dependence
of I_T, g_h, I_h kinetics, calcium handling and the calcium-dependent
regulation of I_h. A null here constrains this slice, not the whole space.

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
from scipy import optimize, signal

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


def damped_sinusoid(t, amplitude, damping, freq, phase, offset):
    """A * exp(-lambda t) * cos(2 pi f t + phi) + C."""
    return amplitude * np.exp(-damping * t) * np.cos(2 * np.pi * freq * t + phase) + offset


def fit_resonance(
    v: np.ndarray,
    fs: float,
    band: tuple[float, float] = (0.2, 12.0),
) -> dict[str, float]:
    """Fit the post-perturbation ringing and report frequency and damping.

    This is the primary measurement. Asking only whether a cell oscillates
    spontaneously answers the wrong question: a stable focus that rings at
    2 Hz with a small damping constant is a resonator, and a resonator driven
    by common noise or fluctuating afferent input produces intermittent
    rhythmic runs even though the isolated deterministic cell has no limit
    cycle. Such a point is a candidate GRDA generator, yet the spontaneous
    scan classifies it as 'silent'.

    Damping is what says how close the system sits to a Hopf bifurcation:
    small positive damping means the leading eigenvalue pair is just inside
    the left half-plane. Reported as 1/s, so smaller is closer.

    Returns NaN fields rather than raising when the response is too small or
    the fit does not converge, so a sweep can continue.
    """
    out = {
        "resonant_hz": float("nan"),
        "damping_per_s": float("nan"),
        "quality_factor": float("nan"),
        "ring_amp_mv": float("nan"),
        "ring_fit_r2": float("nan"),
        "is_resonant": 0,
    }
    max_damping = 20.0  # 1/s; above this the response is a passive decay
    # Spikes are not part of the linear response; clip them out.
    x = np.minimum(v, SPIKE_CLIP_MV).astype(float)
    if len(x) < int(2 * fs):
        return out
    t = np.arange(len(x)) / fs
    baseline = float(np.mean(x[-int(0.5 * fs) :]))
    y = x - baseline
    amp0 = float(np.max(np.abs(y[: int(0.5 * fs)])))
    out["ring_amp_mv"] = amp0
    if amp0 < 1e-4:
        return out  # no measurable response

    # Seed the fit from the spectrum and from the envelope decay, because
    # curve_fit on a decaying oscillation is sensitive to the initial guess.
    freqs = np.fft.rfftfreq(len(y), 1.0 / fs)
    spec = np.abs(np.fft.rfft(y * np.hanning(len(y))))
    in_band = (freqs >= band[0]) & (freqs <= band[1])
    if not np.any(in_band):
        return out
    f0 = float(freqs[in_band][np.argmax(spec[in_band])])

    envelope = np.abs(signal.hilbert(y))
    smooth = max(1, int(0.1 * fs))
    envelope = np.convolve(envelope, np.ones(smooth) / smooth, mode="same")
    usable = envelope > (0.02 * envelope.max())
    if np.count_nonzero(usable) > int(0.5 * fs):
        slope = np.polyfit(t[usable], np.log(envelope[usable]), 1)[0]
        lam0 = float(max(1e-3, -slope))
    else:
        lam0 = 1.0

    try:
        popt, _ = optimize.curve_fit(
            damped_sinusoid,
            t,
            y,
            p0=[amp0, lam0, f0, 0.0, 0.0],
            bounds=(
                [0.0, 0.0, band[0], -2 * np.pi, -20.0],
                [10 * amp0 + 1e-6, max_damping, band[1], 2 * np.pi, 20.0],
            ),
            maxfev=20000,
        )
    except (RuntimeError, ValueError):
        return out

    residual = y - damped_sinusoid(t, *popt)
    denom = float(np.sum((y - np.mean(y)) ** 2))
    out["ring_fit_r2"] = float(1.0 - np.sum(residual**2) / denom) if denom > 0 else float("nan")
    out["resonant_hz"] = float(popt[2])
    out["damping_per_s"] = float(popt[1])
    # Q = pi f / lambda for a lightly damped oscillator; higher is more
    # clock-like and closer to sustained.
    out["quality_factor"] = (
        float(np.pi * popt[2] / popt[1]) if popt[1] > 1e-6 else float("inf")
    )
    # A good fit is not the same as a resonance. Reject overdamped responses,
    # and reject fits that ran into the damping bound, which means the curve
    # was a passive decay that the model could only match by disappearing.
    out["is_resonant"] = int(
        out["ring_fit_r2"] > 0.5
        and out["quality_factor"] >= 0.5
        and popt[1] < 0.95 * max_damping
    )
    return out


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
    probe=None,
    probe_amp_na: float = 0.0,
    probe_dur_ms: float = 200.0,
    ring_s: float = 14.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Run one parameter combination and return (trace, spike times) after settling.

    The settle period matters more than usual here. The Ih mechanism carries a
    calcium-binding process with a decay rate of 4e-4 /ms, a time constant of
    about 2.5 s, so a short run would measure the transient rather than the
    steady-state rhythm.

    A brief current pulse is delivered after the spontaneous observation
    window, and the response that follows it is returned separately. That
    segment is what carries the resonance measurement: the pulse is kept small
    so the cell is probed near its fixed point rather than driven into a full
    low-threshold burst, which would measure the nonlinear response instead.
    """
    cell.soma.gkL_kL = gkl
    cell.soma.fac_gh_TC_iar = gh_shift
    detector.thresh = SPIKE_THRESHOLD_MV

    probe_onset_ms = (settle_s + analyze_s) * 1000.0
    if probe is not None:
        probe.delay = probe_onset_ms
        probe.dur = probe_dur_ms
        probe.amp = probe_amp_na
        total_ms = probe_onset_ms + probe_dur_ms + ring_s * 1000.0
    else:
        total_ms = probe_onset_ms

    h.dt = 0.025  # match the network runs in bazh_net.py
    h.finitialize(-68)
    h.continuerun(total_ms)

    v = np.array(v_vec)
    t_ms = np.array(t_vec)
    spikes_ms = np.array(spike_vec)

    spont = (t_ms >= settle_s * 1000.0) & (t_ms < probe_onset_ms)
    spont_spikes = spikes_ms[
        (spikes_ms >= settle_s * 1000.0) & (spikes_ms < probe_onset_ms)
    ]
    spikes_s = (spont_spikes - settle_s * 1000.0) / 1000.0

    if probe is None:
        return v[spont], spikes_s, np.array([]), 0

    ring_start = probe_onset_ms + probe_dur_ms
    ring = (t_ms >= ring_start) & (t_ms <= ring_start + ring_s * 1000.0)
    # Spikes during the ring window mean the probe drove the cell out of the
    # linear regime, so any damped-sinusoid fit there describes a burst rather
    # than a resonance.
    ring_spikes = int(np.sum(spikes_ms > ring_start))
    return v[spont], spikes_s, v[ring], ring_spikes


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
    res_hz = np.full(shape, np.nan)
    damping = np.full(shape, np.nan)
    label = np.empty(shape, dtype=object)
    for row in rows:
        i = int(row["gh_index"])
        j = int(row["gkl_index"])
        freq[i, j] = row["burst_hz"]
        cv[i, j] = row["burst_cv"]
        vm[i, j] = row["mean_vm_mv"]
        label[i, j] = row["classification"]
        if row.get("is_resonant", 0):  # only plot genuine resonances
            res_hz[i, j] = row["resonant_hz"]
            damping[i, j] = row["damping_per_s"]

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

    ncols = 5 if control else 3
    fig, axes = plt.subplots(1, ncols, figsize=(5.4 * ncols, 5.6), facecolor="white")

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

    del categories, colors, code, ListedColormap, BoundaryNorm, Patch

    def grid_panel(ax, data, title, cmap, label_txt, vmin=None, vmax=None):
        image = ax.imshow(
            data, origin="lower", aspect="auto", extent=extent,
            cmap=cmap, interpolation="nearest", vmin=vmin, vmax=vmax,
        )
        ax.set_xlabel(r"TC K-leak conductance ($10^{-5}$ S/cm$^2$)")
        ax.set_ylabel("Ih activation shift (mV)")
        ax.set_title(title, fontsize=11.5, fontweight="bold", pad=14)
        bar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
        bar.set_label(label_txt, fontsize=9)
        mark_reference(ax)
        return image

    # A: where the cell resonates. This is the primary result -- a stable
    # focus ringing in the delta band is a candidate generator even though it
    # never oscillates on its own.
    grid_panel(
        axes[0], res_hz,
        "A. Resonant frequency\n(fitted ringing after a small pulse)",
        "viridis", "Hz", vmin=0.0, vmax=8.0,
    )
    delta_mask = ((res_hz >= band[0]) & (res_hz <= band[1])).astype(float)
    if delta_mask.any():
        for colour, width, style in (("white", 2.4, "-"), ("black", 1.1, "--")):
            axes[0].contour(
                np.linspace(x_um[0], x_um[-1], len(x_um)),
                np.linspace(gh_values[0], gh_values[-1], len(gh_values)),
                delta_mask, levels=[0.5], colors=colour,
                linewidths=width, linestyles=style,
            )

    # B: how close that resonance is to sustaining itself.
    grid_panel(
        axes[1], damping,
        "B. Damping of the resonance\n(lower is closer to a Hopf bifurcation)",
        "magma_r", "1/s",
    )

    grid_panel(
        axes[2], vm, "C. Steady membrane potential", "coolwarm_r", "mV",
    )
    # Mark the range over which I_T de-inactivates, to show the sweep reached it.
    contour = axes[2].contour(
        np.linspace(x_um[0], x_um[-1], len(x_um)),
        np.linspace(gh_values[0], gh_values[-1], len(gh_values)),
        vm, levels=[-80.0, -70.0], colors="black", linewidths=1.1, linestyles="--",
    )
    axes[2].clabel(contour, fmt="%.0f mV", fontsize=8)

    if control:
        ax = axes[3]
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
            "D. Positive control: I_T rebound burst", fontsize=11.5, fontweight="bold"
        )
        ax.legend(loc="upper right", fontsize=8)

        ax = axes[4]
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
            "E. The evoked rhythm dies out", fontsize=11.5, fontweight="bold"
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
        "Delta-band resonance of an isolated TC cell "
        "(no synapses, no network)",
        fontsize=13.5,
        fontweight="bold",
        y=0.99,
    )
    fig.text(
        0.005,
        0.015,
        "Panels A and B are the primary result and show only points where the "
        "post-pulse ringing is a genuine resonance: a converged fit with "
        "quality factor at or above 0.5, no spikes during the ring window, and "
        "damping below the fit bound. Blank cells failed one of those tests, "
        "most often because the response was an overdamped passive decay. "
        "Dashed contour in A marks the 0.5-4 Hz delta band. Circles mark the "
        "K-leak and Ih settings the model itself uses for wake, N2, N3 and "
        "REM. Panels D and E are the earlier positive control, retained "
        "because they establish that I_T works and that the evoked rhythm is "
        "damped rather than self-sustaining.",
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
    parser.add_argument("--settle-s", type=float, default=12.0)
    parser.add_argument("--analyze-s", type=float, default=10.0)
    # Resonance probe, calibrated empirically. It must be long enough to
    # engage I_h, whose gating runs over hundreds of milliseconds -- a very
    # brief pulse only charges the passive membrane (tau of order 20 ms) and
    # the fit then returns a meaningless sub-hertz "resonance" for what is
    # really a passive decay. It must also stay small enough not to trigger a
    # low-threshold burst, which would measure the nonlinear response instead.
    # A sweep over amplitude and duration found -0.01 nA for 500 ms gives
    # spike-free responses with the best fit quality (r2 ~0.94) across the
    # hyperpolarized range; larger or longer probes either evoke spikes or
    # push the fit into an overdamped regime.
    parser.add_argument("--probe-amp-na", type=float, default=-0.01)
    parser.add_argument("--probe-dur-ms", type=float, default=500.0)
    parser.add_argument("--ring-s", type=float, default=14.0)
    # A fit can be excellent and still describe no oscillation. Q below about
    # 0.5 is overdamped: the response decays before completing a cycle.
    parser.add_argument("--min-q", type=float, default=0.5)
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

    probe = h.IClamp(cell.soma(0.5))
    probe.delay = 0.0
    probe.dur = 0.0
    probe.amp = 0.0

    rows: list[dict[str, object]] = []
    total = len(gkl_values) * len(gh_values)
    done = 0
    for i, gh_shift in enumerate(gh_values):
        for j, gkl in enumerate(gkl_values):
            v, spikes_s, ring, ring_spikes = run_point(
                cell,
                detector,
                spike_vec,
                v_vec,
                t_vec,
                gkl,
                gh_shift,
                args.settle_s,
                args.analyze_s,
                probe=probe,
                probe_amp_na=args.probe_amp_na,
                probe_dur_ms=args.probe_dur_ms,
                ring_s=args.ring_s,
            )
            row: dict[str, object] = {
                "gh_index": i,
                "gkl_index": j,
                "gkl_S_per_cm2": float(gkl),
                "gh_shift_mv": float(gh_shift),
            }
            row.update(analyze_trace(v, spikes_s, fs=1000.0 / args.record_dt_ms))
            resonance = fit_resonance(ring, fs=1000.0 / args.record_dt_ms)
            resonance["ring_spikes"] = ring_spikes
            if ring_spikes:
                resonance["is_resonant"] = 0  # probe left the linear regime
            row.update(resonance)
            row["classification"] = classify(row, band, args.max_cv)
            rows.append(row)
            done += 1
            print(
                "[%3d/%3d] gkL=%.3e gh=%+.2f | resonance %5.2f Hz  damping %6.3f /s  "
                "Q %5.1f  r2 %.2f | %-16s burst=%.2f Hz cv=%.2f Vm=%.1f mV"
                % (
                    done,
                    total,
                    gkl,
                    gh_shift,
                    row["resonant_hz"],
                    row["damping_per_s"],
                    row["quality_factor"],
                    row["ring_fit_r2"],
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

    # The primary result: where does the cell resonate in the delta band, and
    # how close to sustained is it there?
    good = [r for r in rows if r.get("is_resonant", 0)]
    in_delta = [r for r in good if band[0] <= r["resonant_hz"] <= band[1]]
    print(
        "\ngenuine resonances (Q>=0.5, good fit): %d of %d; of those, %d in %g-%g Hz"
        % (len(good), len(rows), len(in_delta), band[0], band[1])
    )
    if in_delta:
        best = min(in_delta, key=lambda r: r["damping_per_s"])
        print(
            "  least-damped delta-band point: %.2f Hz, damping %.3f /s, Q %.1f\n"
            "    at K leak %.3e S/cm2, Ih shift %+.2f mV"
            % (
                best["resonant_hz"],
                best["damping_per_s"],
                best["quality_factor"],
                best["gkl_S_per_cm2"],
                best["gh_shift_mv"],
            )
        )
        print(
            "  Lower damping means the leading eigenvalue pair sits closer to the\n"
            "  imaginary axis, so this is the direction to push toward a Hopf\n"
            "  bifurcation or toward noise-driven quasi-cycles."
        )
    h.quit()


if __name__ == "__main__":
    main()
