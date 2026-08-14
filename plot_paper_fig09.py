"""Reproduce Figure 9 of Krishnan et al. 2016 eLife (propofol anesthesia).

Panels follow the paper:

  A  propofol: ACh and HA reduced to the cat natural-sleep (SWS) level and the
     GABA-A decay time constant lengthened          (model state 49)
  B  the GABA-A decay change alone, from awake      (model state 50)
  C  control, cat natural sleep (SWS)               (model state 51)
  D  power in the 0.2-4 Hz and 8-15 Hz bands for all three conditions

Each of panels A-C shows, from top to bottom, the PY spike raster, the
simulated cortical LFP, and that LFP filtered to the spindle band.

Required input files, produced by `bazh_net.py` in the working directory:

    raster_state{S}_seed{SEED}_nhost={N}.txt
    vcort_state{S}_seed{SEED}_nhost={N}.txt   (or lfp_... with --signal lfp)

for S in the three (or five) requested states.  Example production runs:

    for S in 49 50 51; do
        FINK_STATE=$S FINK_SEED=1 FINK_DURATION_MS=120000 \
            mpiexec -n 8 nrniv -mpi -python bazh_net.py \
            2>&1 | tee run_state${S}_seed1.log
    done
    python plot_paper_fig09.py --nhost 8 --duration-ms 120000

Two deliberate deviations from the paper are noted on the figure itself.
First, the paper's panels A-C show colour-coded PY membrane voltages, which
this port does not record; the spike raster is shown instead.  Second, the
paper's panel D averages membrane voltage within 10 groups of 50 PY neurons
and shows the standard deviation across those groups, while this port records
only the summed cortical voltage, so the error bars here are the standard
deviation across 10 consecutive epochs of the same trace.
"""

import argparse
from fractions import Fraction
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy import signal

# Cortical cells contributing to the summed voltage trace: 500 PY and 100 IN
# cells, each with a somatic and a dendritic segment (see cell_classes.py).
N_CORT_SEGMENTS = (500 + 100) * 2
N_PYR = 500

PANEL_TITLES = {
    49: "A   Propofol - ACh, HA reduction and\n"
        "     GABA decay time constant increase",
    50: "B   Only GABA decay time constant increase",
    51: "C   Control - Natural Sleep (SWS)",
    52: "A'   Propofol, GABA decay x2.5 reading",
    53: "B'   Only GABA decay x2.5 reading",
}
BAR_LABELS = {
    49: "Propofol",
    50: "Only\nGABA",
    51: "SWS",
    52: "Propofol\n(x2.5)",
    53: "Only GABA\n(x2.5)",
}

parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument(
    "--states",
    nargs="+",
    type=int,
    default=[49, 50, 51],
    help="state IDs for panels A, B, C (default: 49 50 51)",
)
parser.add_argument("--seed", type=int, default=1, help="FINK_SEED of the runs")
parser.add_argument("--nhost", type=int, default=8, help="MPI rank count of the runs")
parser.add_argument(
    "--duration-ms",
    type=float,
    default=120000.0,
    help="FINK_DURATION_MS of the runs (default: 120000)",
)
parser.add_argument(
    "--skip-s",
    type=float,
    default=10.0,
    help="seconds of initialization transient to discard (default: 10)",
)
parser.add_argument(
    "--window-s",
    type=float,
    default=30.0,
    help="length of the trace window shown in panels A-C (default: 30)",
)
parser.add_argument(
    "--signal",
    choices=("vcort", "lfp"),
    default="vcort",
    help="vcort (mean cortical membrane voltage, as in the paper) or the "
         "biophysical extracellular LFP (default: vcort)",
)
parser.add_argument(
    "--view-fs",
    type=float,
    default=200.0,
    help="analysis/display sampling rate in Hz (default: 200)",
)
parser.add_argument(
    "--slow-band",
    type=float,
    nargs=2,
    metavar=("LOW_HZ", "HIGH_HZ"),
    default=(0.2, 4.0),
    help="slow band for panel D (default: 0.2 4)",
)
parser.add_argument(
    "--spindle-band",
    type=float,
    nargs=2,
    metavar=("LOW_HZ", "HIGH_HZ"),
    default=(8.0, 15.0),
    help="spindle band for the filtered trace and panel D (default: 8 15)",
)
parser.add_argument(
    "--n-epochs",
    type=int,
    default=10,
    help="epochs used for the panel D error bars (default: 10)",
)
parser.add_argument(
    "--output",
    default="docs/figures/paper_fig09_propofol.png",
    help="output PNG path (default: docs/figures/paper_fig09_propofol.png)",
)
args = parser.parse_args()

if len(args.states) < 2:
    parser.error("provide at least two states, for example: --states 49 50 51")
if args.skip_s < 0 or args.skip_s >= args.duration_ms / 1000.0:
    parser.error("--skip-s must be within the simulated duration")
if args.window_s <= 0:
    parser.error("--window-s must be positive")
if args.n_epochs < 2:
    parser.error("--n-epochs must be at least 2")
if not 0 < args.spindle_band[0] < args.spindle_band[1] < args.view_fs / 2:
    parser.error("--spindle-band must satisfy 0 < LOW < HIGH < view-fs/2")
if not 0 < args.slow_band[0] < args.slow_band[1] < args.view_fs / 2:
    parser.error("--slow-band must satisfy 0 < LOW < HIGH < view-fs/2")


def state_file(prefix, state):
    return Path(f"{prefix}_state{state}_seed{args.seed}_nhost={args.nhost}.txt")


def load_column(path):
    """Load a one-column text file, preferring pandas for large files."""
    try:
        import pandas as pd

        return pd.read_csv(path, header=None, sep=r"\s+").iloc[:, 0].to_numpy(float)
    except ImportError:
        return np.loadtxt(path)


def load_raster(path):
    """Return spike times (s) and gids of PY cells."""
    if path.stat().st_size == 0:
        return np.empty(0), np.empty(0)
    try:
        import pandas as pd

        data = pd.read_csv(path, header=None, sep=r"\s+").to_numpy(float)
    except ImportError:
        data = np.loadtxt(path)
    if data.size == 0:
        return np.empty(0), np.empty(0)
    data = np.atleast_2d(data)
    pyr = data[data[:, 1] < N_PYR]
    return pyr[:, 0] / 1000.0, pyr[:, 1]


def resample_to_view(y, raw_fs):
    """Anti-aliased resampling of a raw trace to the analysis rate."""
    y = y - np.median(y)
    ratio = Fraction(args.view_fs / raw_fs).limit_denominator(10000)
    return signal.resample_poly(
        y, ratio.numerator, ratio.denominator, padtype="line"
    )


def bandpass(y, band):
    sos = signal.butter(4, band, btype="bandpass", fs=args.view_fs, output="sos")
    padlen = min(len(y) - 1, int(3 * args.view_fs / band[0]))
    return signal.sosfiltfilt(sos, y, padtype="even", padlen=padlen)


def band_power(y, band):
    """Mean-square power in a band, from Welch's periodogram."""
    nperseg = int(min(len(y), max(args.view_fs * 10, 4 * args.view_fs / band[0])))
    freqs, psd = signal.welch(y, fs=args.view_fs, nperseg=nperseg)
    mask = (freqs >= band[0]) & (freqs <= band[1])
    if not mask.any():
        return float("nan")
    return float(np.trapezoid(psd[mask], freqs[mask])
                 if hasattr(np, "trapezoid") else np.trapz(psd[mask], freqs[mask]))


missing = [
    str(state_file(prefix, state))
    for state in args.states
    for prefix in (args.signal, "raster")
    if not state_file(prefix, state).exists()
]
if missing:
    raise SystemExit(
        "Missing simulation output:\n  "
        + "\n  ".join(missing)
        + "\nRun bazh_net.py for each state with the matching FINK_SEED, "
          "FINK_DURATION_MS and rank count, or pass --seed/--nhost/"
          "--duration-ms to match existing files."
    )

traces = {}
for state in args.states:
    raw = load_column(state_file(args.signal, state))
    if args.signal == "vcort":
        # Summed cortical voltage -> mean cortical membrane voltage, the
        # closest available analogue of the paper's "average activity in the
        # PY network (simulated LFP)".
        raw = raw / N_CORT_SEGMENTS
    raw_fs = len(raw) / (args.duration_ms / 1000.0)
    y = resample_to_view(raw, raw_fs)
    t = np.arange(len(y)) / args.view_fs
    keep = t >= args.skip_s
    y, t = y[keep], t[keep]
    y = y - np.median(y)
    spike_t, spike_gid = load_raster(state_file("raster", state))
    traces[state] = {"t": t, "y": y, "spike_t": spike_t, "spike_gid": spike_gid}

# Panel D: band power over the whole analysis epoch, with the spread across
# consecutive sub-epochs as an error estimate.
powers = {}
for state in args.states:
    y = traces[state]["y"]
    chunks = np.array_split(y, args.n_epochs)
    entry = {}
    for name, band in (("slow", args.slow_band), ("spindle", args.spindle_band)):
        entry[name] = band_power(y, band)
        entry[name + "_sd"] = float(
            np.std([band_power(c, band) for c in chunks if len(c) > 8])
        )
    powers[state] = entry

# Panels A-C show a common window near the end of the analysis epoch.
window_start = max(
    args.skip_s, args.duration_ms / 1000.0 - args.window_s
)
window_end = min(args.duration_ms / 1000.0, window_start + args.window_s)

n_trace_panels = len(args.states)
n_rows = int(np.ceil((n_trace_panels + 1) / 2))
fig = plt.figure(figsize=(13.5, 3.9 * n_rows))
outer = gridspec.GridSpec(n_rows, 2, figure=fig, hspace=0.42, wspace=0.16)

trace_axes = []
for idx, state in enumerate(args.states):
    inner = gridspec.GridSpecFromSubplotSpec(
        3, 1, subplot_spec=outer[idx // 2, idx % 2], hspace=0.12,
        height_ratios=[1.5, 1.0, 0.8],
    )
    ax_raster = fig.add_subplot(inner[0])
    ax_lfp = fig.add_subplot(inner[1], sharex=ax_raster)
    ax_band = fig.add_subplot(inner[2], sharex=ax_raster)
    trace_axes.append((state, ax_raster, ax_lfp, ax_band))

    d = traces[state]
    sel = (d["spike_t"] >= window_start) & (d["spike_t"] <= window_end)
    ax_raster.plot(
        d["spike_t"][sel], d["spike_gid"][sel], "k.", ms=0.9, alpha=0.55
    )
    ax_raster.set_ylim(0, N_PYR)
    ax_raster.set_ylabel("PY cell")
    ax_raster.set_title(
        PANEL_TITLES.get(state, f"state {state}"), loc="left", fontsize=11
    )
    plt.setp(ax_raster.get_xticklabels(), visible=False)

    sel = (d["t"] >= window_start) & (d["t"] <= window_end)
    ax_lfp.plot(d["t"][sel], d["y"][sel], lw=0.6, color="navy")
    ax_lfp.set_ylabel("LFP (mV)" if args.signal == "vcort" else "LFP (a.u.)")
    plt.setp(ax_lfp.get_xticklabels(), visible=False)

    filtered = bandpass(d["y"], args.spindle_band)
    ax_band.plot(d["t"][sel], filtered[sel], lw=0.6, color="darkred")
    ax_band.set_ylabel(
        "%g-%g Hz" % (args.spindle_band[0], args.spindle_band[1])
    )
    ax_band.set_xlabel("time (s)")
    ax_band.set_xlim(window_start, window_end)

# Share the voltage scale across conditions, as the paper does with one scale
# bar per row.
for row_axes in (2, 3):
    lims = [
        max(abs(np.asarray(ax[row_axes].get_ylim()))) for ax in trace_axes
    ]
    for ax in trace_axes:
        ax[row_axes].set_ylim(-max(lims), max(lims))

# Panel D of the paper is two bar charts with independent scales, ordered
# propofol, SWS, only GABA.
bar_order = [s for s in (49, 51, 50, 52, 53) if s in args.states]
bar_states = bar_order + [s for s in args.states if s not in bar_order]
bar_grid = gridspec.GridSpecFromSubplotSpec(
    1, 2, subplot_spec=outer[n_trace_panels // 2, n_trace_panels % 2], wspace=0.35
)
unit = "mV$^2$" if args.signal == "vcort" else "a.u."
x = np.arange(len(bar_states))
for col, (key, band, color, title) in enumerate(
    (
        ("slow", args.slow_band, "#4c72b0", "Power at %g-%g Hz"),
        ("spindle", args.spindle_band, "#c44e52", "Power in %g-%g Hz"),
    )
):
    ax_bar = fig.add_subplot(bar_grid[col])
    ax_bar.bar(
        x,
        [powers[s][key] for s in bar_states],
        0.62,
        yerr=[powers[s][key + "_sd"] for s in bar_states],
        capsize=3,
        color=color,
    )
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(
        [BAR_LABELS.get(s, f"state {s}") for s in bar_states], fontsize=9
    )
    ax_bar.set_title(title % (band[0], band[1]), fontsize=10)
    if col == 0:
        ax_bar.set_ylabel("power (%s)" % unit)
        ax_bar.text(
            -0.22, 1.16, "D", transform=ax_bar.transAxes, fontsize=11,
            fontweight="bold",
        )

fig.suptitle(
    "Krishnan et al. 2016, Figure 9: mechanism of action of propofol "
    "induced anesthesia (NEURON port)",
    fontsize=12,
)
fig.text(
    0.01,
    0.005,
    "Rasters replace the paper's colour-coded PY membrane voltages, which "
    "this port does not record. Panel D error bars are the standard "
    "deviation across %d consecutive epochs of the %.0f s analysis window; "
    "the paper instead used 10 spatial groups of 50 PY neurons. "
    "seed %d, %g s simulated, %g s discarded as initialization."
    % (
        args.n_epochs,
        args.duration_ms / 1000.0 - args.skip_s,
        args.seed,
        args.duration_ms / 1000.0,
        args.skip_s,
    ),
    fontsize=7.5,
    color="0.35",
)

output = Path(args.output)
output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output, dpi=200, bbox_inches="tight")
print("wrote %s" % output)
for state in args.states:
    print(
        "state %d: %g-%g Hz power = %.4g, %g-%g Hz power = %.4g"
        % (
            state,
            args.slow_band[0],
            args.slow_band[1],
            powers[state]["slow"],
            args.spindle_band[0],
            args.spindle_band[1],
            powers[state]["spindle"],
        )
    )
