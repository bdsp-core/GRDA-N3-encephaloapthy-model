#!/usr/bin/env python
"""Reproduce Figure 10 of Krishnan et al. 2016 eLife (e18607) from this NEURON port.

Figure 10, "Effect of GABA on the network state transitions", is three rows,
each one complete Wake -> N2 -> N3 -> REM -> N2 cycle run under a different
GABA model:

  A  "Fixed GABA conductance" -- both phasic and tonic GABA conductances fixed
     for the entire simulation.
  B  "Tonic inhibition based on microdialysis measurements" -- tonic GABA
     conductance (the conductance of the miniature IPSPs) at 115%, 130% and 75%
     of awake for N2, N3 and REM; phasic GABA fixed.
  C  "Tonic inhibition based on co-release of GABA with ACh and HA" -- tonic
     GABA conductance at 85%, 70% and 125% of awake; phasic GABA fixed.

Each row of the paper has three blocks, reproduced here in the same order:
a full-cycle population display, the full-cycle LFP with a 10 s scale bar, and
three zoomed windows (N2, N3, REM) with a 2 s scale bar.

IMPORTANT DEVIATION FROM THE PAPER
----------------------------------
The paper's population display is a *voltage image*: one row per PY cell,
colour-coded V_m between -75 and -45 mV.  This port does not save the membrane
voltage of all 500 PY cells for a 360 s run, so the honest equivalent that the
saved data supports is a spike raster, exactly as in plot_paper_fig03.py.  A
raster shows the spatiotemporal structure of firing (spindle waxing/waning,
Up/Down states, spatially localised alpha-bursts) but not sub-threshold
polarisation, so there is no -75/-45 mV colour bar here.

The trace under each raster is this port's summed cortical membrane voltage
(`vcort_*`), divided by the cell count and mean-subtracted; the paper's LFP is
the mean voltage of the 500 PY cells, so this is the closest available
quantity.  `--signal lfp` switches to the biophysical extracellular LFP
(`lfp_*`) instead.

Inputs (produced by bazh_net.py, written to the repo root):
    raster_fig10{A,B,C}_sweep_seed{SEED}_nhost={N}.txt   "t_ms  gid"
    vcort_fig10{A,B,C}_sweep_seed{SEED}_nhost={N}.txt    one sample per h.dt
    lfp_fig10{A,B,C}_sweep_seed{SEED}_nhost={N}.txt      (only for --signal lfp)

Produce them with (360 s each, so run them on MPI, not serially):

    FINK_SLEEPSTATES=1 FINK_FIG10=A FINK_SEED=1 FINK_DURATION_MS=360000 \
        mpiexec -n 8 ./.venv/bin/nrniv -mpi -python bazh_net.py
    ... and the same with FINK_FIG10=B and FINK_FIG10=C.

Usage:
    .venv/bin/python plot_paper_fig10.py
    .venv/bin/python plot_paper_fig10.py --seed 1 --nhost 8
    .venv/bin/python plot_paper_fig10.py --zoom N3:190:10
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

REPO = os.path.dirname(os.path.abspath(__file__))
OUT_PATH = os.path.join(REPO, "docs", "figures", "paper_fig10_gaba.png")

# ---------------------------------------------------------------------------
# The three GABA models.  Percentages are of the awake mIPSP conductance, in
# stage order N2, N3, REM, and match config.py's fig10_tonic_* tuples.
# ---------------------------------------------------------------------------
PANELS = [
    ("A", "Fixed GABA conductance", (1.00, 1.00, 1.00)),
    ("B", "Tonic inhibition based on microdialysis measurements", (1.15, 1.30, 0.75)),
    ("C", "Tonic inhibition based on co-release of GABA with ACh and HA", (0.85, 0.70, 1.25)),
]

# Sleep-cycle schedule, from config.py's do_sleepstates block: each transition
# is a 17.5 s ramp, and the run ends at 360 s.
TRANSITIONS = [(80.0, 97.5), (150.0, 167.5), (220.0, 237.5), (290.0, 307.5)]
CYCLE_END_S = 360.0
STAGE_LABELS = ["Awake", "N2", "N3", "REM", "N2"]

# Zoom windows (start_s, duration_s), one per inset column, chosen to sit well
# inside each stage plateau and away from the ramps.
DEFAULT_ZOOMS = {
    "N2": (120.0, 10.0),
    "N3": (190.0, 10.0),
    "REM": (260.0, 10.0),
}
ZOOM_ORDER = ["N2", "N3", "REM"]

NPYR = 500        # config.Npyr; PY cells are gids 0..499
NCORT = 600       # PY + IN, the cells summed into vcort
DT_MS = 0.025     # h.dt in bazh_net.py; one vcort/lfp sample per step
VIEW_FS = 200.0   # display sampling rate after anti-aliased decimation, Hz

INK = "#0b0b0b"
MUTED = "#6b6a66"
FAINT = "#d8d7d2"
RASTER_C = "#2a78d6"
TRACE_C = "#101010"


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def run_tag(panel, seed, nhost):
    return "fig10%s_sweep_seed%d_nhost=%d" % (panel, seed, nhost)


def data_path(root, kind, panel, seed, nhost):
    return os.path.join(root, "%s_%s.txt" % (kind, run_tag(panel, seed, nhost)))


def load_raster(path, gid_hi=NPYR):
    """Return (t_s, gid) for spikes from cells below gid_hi.

    The raster file is written by every MPI rank, so it is not time sorted; we
    only ever mask, never assume order.
    """
    dat = np.loadtxt(path)
    if dat.size == 0:
        return np.empty(0), np.empty(0)
    dat = np.atleast_2d(dat)
    m = dat[:, 1] < gid_hi
    return dat[m, 0] / 1000.0, dat[m, 1]


def load_signal(path, sample_dt_ms=DT_MS, view_fs=VIEW_FS):
    """Load a one-sample-per-timestep trace and decimate it to view_fs.

    The files hold one ASCII value per line at 1/h.dt = 40 kHz, i.e. 14.4e6
    values for a 360 s cycle, so this reads with np.fromfile (C parser) and
    box-car averages before subsampling to avoid aliasing the display.
    """
    x = np.fromfile(path, dtype=float, sep=" ")
    if x.size == 0:
        return np.empty(0), np.empty(0)
    step = max(1, int(round((1000.0 / sample_dt_ms) / view_fs)))
    n = (x.size // step) * step
    if n == 0:
        return np.arange(x.size) * sample_dt_ms / 1000.0, x
    y = x[:n].reshape(-1, step).mean(axis=1)
    t = (np.arange(y.size) + 0.5) * step * sample_dt_ms / 1000.0
    return t, y


def panel_available(root, panel, seed, nhost, signal_kind):
    return os.path.exists(data_path(root, "raster", panel, seed, nhost)) and os.path.exists(
        data_path(root, signal_kind, panel, seed, nhost)
    )


# ---------------------------------------------------------------------------
# Quantitative summary (printed, not plotted)
# ---------------------------------------------------------------------------
def stage_windows():
    """Plateau windows (start_s, end_s) for the five cycle blocks."""
    edges = [0.0]
    for a, b in TRANSITIONS:
        edges += [a, b]
    edges.append(CYCLE_END_S)
    # plateaus are [0,80], [97.5,150], [167.5,220], [237.5,290], [307.5,360]
    return [(edges[i], edges[i + 1]) for i in range(0, len(edges) - 1, 2)]


def band_power(t, y, lo, hi):
    if t.size < 16:
        return 0.0
    fs = 1.0 / np.median(np.diff(t))
    x = y - y.mean()
    freqs = np.fft.rfftfreq(x.size, d=1.0 / fs)
    psd = np.abs(np.fft.rfft(x * np.hanning(x.size))) ** 2
    b = (freqs >= lo) & (freqs < hi)
    tot = (freqs >= 0.3) & (freqs < 30.0)
    denom = psd[tot].sum()
    return float(psd[b].sum() / denom) if denom > 0 else 0.0


def describe_panel(letter, t_sig, y_sig, t_spk, tonic):
    print(
        "\n--- Panel %s  (tonic mIPSP N2/N3/REM = %d%%/%d%%/%d%% of awake) ---"
        % (letter, round(100 * tonic[0]), round(100 * tonic[1]), round(100 * tonic[2]))
    )
    for (a, b), name in zip(stage_windows(), STAGE_LABELS):
        m = (t_sig >= a) & (t_sig < b)
        ms = (t_spk >= a) & (t_spk < b)
        rate = ms.sum() / float(NPYR) / max(b - a, 1e-9)
        if m.sum() < 16:
            print("    %-5s  no signal samples in %.0f-%.0f s" % (name, a, b))
            continue
        print(
            "    %-5s %6.1f-%6.1f s | PY rate %5.2f Hz | slow 0.5-1.5 Hz %4.1f%% | "
            "delta 1.5-4 Hz %4.1f%% | spindle/alpha 8-15 Hz %4.1f%%"
            % (
                name, a, b, rate,
                100 * band_power(t_sig[m], y_sig[m], 0.5, 1.5),
                100 * band_power(t_sig[m], y_sig[m], 1.5, 4.0),
                100 * band_power(t_sig[m], y_sig[m], 8.0, 15.0),
            )
        )


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------
def bare(ax, panel=False):
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    if panel:
        ax.set_facecolor("#f6f6f4")
    else:
        ax.patch.set_alpha(0.0)


def draw_raster(ax, t_spk, gids, t0, t1, show_stage_lines):
    bare(ax, panel=True)
    m = (t_spk >= t0) & (t_spk < t1)
    ax.plot(
        t_spk[m], gids[m], linestyle="none", marker="o", markersize=0.35,
        markeredgewidth=0, color=RASTER_C, alpha=0.8, rasterized=True,
    )
    if m.sum() == 0:
        ax.text(0.5, 0.5, "no PY spikes", transform=ax.transAxes, ha="center",
                va="center", fontsize=7.5, color=MUTED, style="italic")
    ax.set_xlim(t0, t1)
    ax.set_ylim(0, NPYR)
    if show_stage_lines:
        for a, b in TRANSITIONS:
            ax.axvline(0.5 * (a + b), color=INK, lw=1.1)


def draw_trace(ax, t_sig, y_sig, t0, t1, ylim=None):
    bare(ax)
    m = (t_sig >= t0) & (t_sig < t1)
    if m.sum():
        y = y_sig[m] - y_sig[m].mean()
        ax.plot(t_sig[m], y, lw=0.45, color=TRACE_C, rasterized=True)
    ax.set_xlim(t0, t1)
    if ylim is not None:
        ax.set_ylim(*ylim)
    return ax


def scale_bar(ax, seconds, label, y=0.02):
    """Horizontal scale bar in axes-fraction y, data-coordinate x."""
    x0, x1 = ax.get_xlim()
    ax.plot([x0, x0 + seconds], [y, y], transform=ax.get_xaxis_transform(),
            color=INK, lw=1.5, solid_capstyle="butt", clip_on=False)
    ax.text(x0 + seconds / 2.0, y + 0.04, label, transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=7.5, color=INK)


def draw_missing(fig, gs, letter, title, panel_tag):
    ax = fig.add_subplot(gs)
    bare(ax)
    ax.text(
        0.5, 0.5,
        "%s  %s\n\nno data: expected raster_%s.txt and the matching signal file\n"
        "in the repo root (see this script's docstring for the run command)"
        % (letter, title, panel_tag),
        ha="center", va="center", fontsize=9, color=MUTED, linespacing=1.6,
        transform=ax.transAxes,
    )
    ax.add_patch(plt.Rectangle((0.01, 0.05), 0.98, 0.9, transform=ax.transAxes,
                               fill=False, ec=FAINT, lw=0.8, ls=(0, (4, 4))))


def draw_panel(fig, gs_row, letter, title, tonic, t_spk, gids, t_sig, y_sig, zooms, ylim):
    """One paper row: full-cycle raster + LFP, then three zoom columns."""
    inner = GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs_row, height_ratios=[1.55, 1.0], hspace=0.42
    )
    top = GridSpecFromSubplotSpec(2, 1, subplot_spec=inner[0], height_ratios=[2.5, 1.0], hspace=0.06)
    bot = GridSpecFromSubplotSpec(1, 3, subplot_spec=inner[1], wspace=0.075)

    # ---- full cycle -------------------------------------------------------
    ax_r = fig.add_subplot(top[0])
    draw_raster(ax_r, t_spk, gids, 0.0, CYCLE_END_S, show_stage_lines=True)
    ax_r.text(-0.022, 0.5, "PY", transform=ax_r.transAxes, ha="right", va="center",
              fontsize=9.5, color=INK, fontweight="bold")
    ax_r.set_title("%s   %s" % (letter, title), fontsize=10.5, color=INK, pad=14, loc="left")
    for (a, b), name in zip(stage_windows(), STAGE_LABELS):
        ax_r.text(0.5 * (a + b), 1.012, name, transform=ax_r.get_xaxis_transform(),
                  ha="center", va="bottom", fontsize=8.5, color=INK)

    ax_t = fig.add_subplot(top[1])
    draw_trace(ax_t, t_sig, y_sig, 0.0, CYCLE_END_S, ylim)
    scale_bar(ax_t, 10.0, "10 s", y=-0.08)

    # ---- zoom columns -----------------------------------------------------
    for col, name in enumerate(ZOOM_ORDER):
        z0, zdur = zooms[name]
        cell = GridSpecFromSubplotSpec(2, 1, subplot_spec=bot[col],
                                       height_ratios=[1.7, 1.0], hspace=0.08)
        axz = fig.add_subplot(cell[0])
        draw_raster(axz, t_spk, gids, z0, z0 + zdur, show_stage_lines=False)
        axz.set_title("%s  %.0f-%.0f s" % (name, z0, z0 + zdur), fontsize=8,
                      color=MUTED, pad=3, loc="left")
        axzt = fig.add_subplot(cell[1])
        draw_trace(axzt, t_sig, y_sig, z0, z0 + zdur, ylim)
        if col == 0:
            scale_bar(axzt, 2.0, "2 s", y=-0.12)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--nhost", type=int, default=8, help="MPI rank count in the filenames")
    ap.add_argument("--root", default=REPO, help="directory holding the run output files")
    ap.add_argument("--signal", choices=("vcort", "lfp"), default="vcort",
                    help="vcort = summed cortical Vm (paper-like LFP, default); "
                         "lfp = biophysical extracellular LFP")
    ap.add_argument("--out", default=OUT_PATH)
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument("--sample-dt-ms", type=float, default=DT_MS,
                    help="sample spacing of the signal files (default h.dt = 0.025 ms)")
    ap.add_argument("--zoom", action="append", default=[], metavar="STAGE:START_S:DUR_S",
                    help="override a zoom window, e.g. --zoom N3:190:10")
    args = ap.parse_args()

    zooms = dict(DEFAULT_ZOOMS)
    for spec in args.zoom:
        name, a, b = spec.split(":")
        if name not in zooms:
            ap.error("--zoom stage must be one of %s" % ", ".join(ZOOM_ORDER))
        zooms[name] = (float(a), float(b))

    fig = plt.figure(figsize=(13.0, 15.0), facecolor="white", dpi=args.dpi)
    gs = GridSpec(3, 1, figure=fig, left=0.06, right=0.985, top=0.900, bottom=0.045,
                  hspace=0.30)

    fig.text(0.06, 0.973, "Effect of GABA on the network state transitions",
             fontsize=13.5, color=INK, fontweight="semibold", ha="left", va="center")
    fig.text(
        0.06, 0.952,
        "NEURON port of Krishnan et al. 2016 eLife, Figure 10  -  seed %d, %d MPI ranks  -  "
        "one Wake-N2-N3-REM-N2 cycle per row, phasic GABA fixed at its awake value in all three rows"
        % (args.seed, args.nhost),
        fontsize=8.2, color=MUTED, ha="left", va="center",
    )
    fig.text(
        0.06, 0.936,
        "Each row: PY spike raster over the cycle (the paper shows a membrane-voltage image; this "
        "port does not store Vm for 500 cells), the %s below it, then N2/N3/REM zooms"
        % ("mean cortical membrane voltage" if args.signal == "vcort"
           else "biophysical extracellular LFP"),
        fontsize=8.2, color=MUTED, ha="left", va="center",
    )

    print("Figure 10 reproduction -- seed %d, nhost %d, signal %s" % (args.seed, args.nhost, args.signal))

    # Load everything first so all three rows share one signal scale.
    loaded, missing = {}, []
    for letter, _title, _tonic in PANELS:
        if not panel_available(args.root, letter, args.seed, args.nhost, args.signal):
            missing.append(letter)
            continue
        t_spk, gids = load_raster(data_path(args.root, "raster", letter, args.seed, args.nhost))
        t_sig, y_sig = load_signal(
            data_path(args.root, args.signal, letter, args.seed, args.nhost),
            sample_dt_ms=args.sample_dt_ms,
        )
        if args.signal == "vcort":
            y_sig = y_sig / float(NCORT)  # summed -> mean cortical membrane voltage
        loaded[letter] = (t_spk, gids, t_sig, y_sig)

    ylim = None
    if loaded:
        spread = max(np.percentile(np.abs(v[3] - v[3].mean()), 99.5) for v in loaded.values()
                     if v[3].size)
        ylim = (-1.35 * spread, 1.35 * spread)

    for row, (letter, title, tonic) in enumerate(PANELS):
        if letter not in loaded:
            draw_missing(fig, gs[row], letter, title, run_tag(letter, args.seed, args.nhost))
            print("\n--- Panel %s: no data on disk ---" % letter)
            continue
        t_spk, gids, t_sig, y_sig = loaded[letter]
        draw_panel(fig, gs[row], letter, title, tonic, t_spk, gids, t_sig, y_sig, zooms, ylim)
        describe_panel(letter, t_sig, y_sig, t_spk, tonic)

    fig.text(
        0.06, 0.012,
        "Rasters: each dot is one PY spike, y is PY index 0-499. Vertical lines mark the midpoint of "
        "each 17.5 s neuromodulator transition. Traces are mean-subtracted and share one scale.",
        fontsize=7.8, color=MUTED, ha="left", va="center",
    )

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    fig.savefig(args.out, dpi=args.dpi, facecolor="white")
    print("\nWrote %s" % args.out)
    if missing:
        print("Panels missing: %s" % ", ".join(missing))
        print("Run each missing panel with FINK_SLEEPSTATES=1 FINK_FIG10=<panel> "
              "FINK_DURATION_MS=360000 (see the docstring).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
