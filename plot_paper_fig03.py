#!/usr/bin/env python
"""Reproduce Figure 3 of Krishnan et al. 2016 eLife (e18607) from this NEURON port.

Figure 3 of the paper shows, in four columns (Awake / N2 / N3 / REM):

  * top block  -- network activity of PY, IN, TC and RE cells over a short window
                  (5 s for awake, 10 s for the sleep stages), y = neuron index,
                  colour = membrane voltage;
  * bottom block -- membrane voltage of four representative cells,
                  PY #250, IN #50, TC #50, RE #50.

IMPORTANT DEVIATION FROM THE PAPER
----------------------------------
The paper's top block is a *voltage image*: one row per cell, colour-coded V_m
between -80 and -40 mV.  This port does not save the membrane voltage of all 800
cells (it would be ~800 x 4.8e6 samples per run), so the honest equivalent that
the saved data supports is a **spike raster** for each population over the same
window.  A raster shows the spatiotemporal structure of firing (Up/Down states,
spindle waxing/waning, spatial localisation of alpha-bursts) but it cannot show
sub-threshold membrane polarisation, so the depolarised-but-silent TC cells of
N3 and the tonic hyperpolarisation of N2 are visible only in the bottom block
of single-cell traces.  There is consequently no -80/-40 mV colour bar; colour
in the top block encodes population identity, not voltage.

Everything else follows the paper: same populations, same row order
(PY, IN, TC, RE), same representative cells, same window lengths, and a single
shared voltage scale across all four columns so the columns are comparable.

Inputs (produced by bazh_net.py, all in the repo root):
    raster_state{S}_seed{SEED}_nhost={N}.txt        "t_ms  gid", not time sorted
    trace_PY_gid250_state{S}_seed{SEED}_nhost={N}.txt   "# header" then "t_ms  V_mV"
    trace_IN_gid550_...  trace_RE_gid650_...  trace_TC_gid750_...

Usage:
    .venv/bin/python plot_paper_fig03.py
    .venv/bin/python plot_paper_fig03.py --win 2:60:10      # N3 window 60-70 s
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
OUT_PATH = os.path.join(REPO, "docs", "figures", "paper_fig03_stage_activity.png")

# ----------------------------------------------------------------------------
# Populations.  gid layout is fixed by config.py: Npyr=500, Ninh=100, Nre=100,
# Ntc=100 assigned in that order, so PY = 0-499, IN = 500-599, RE = 600-699,
# TC = 700-799.  The paper's "neuron #k" is the k-th cell within its population,
# hence IN #50 = gid 550, RE #50 = gid 650, TC #50 = gid 750.
# Row order below is the paper's: PY, IN, TC, RE.
# ----------------------------------------------------------------------------
POPS = [
    # label, gid_lo, gid_hi (exclusive), n_cells, trace gid, trace file tag, colour
    ("PY", 0, 500, 500, 250, "PY_gid250", "#2a78d6"),
    ("IN", 500, 600, 100, 550, "IN_gid550", "#eb6834"),
    ("TC", 700, 800, 100, 750, "TC_gid750", "#1baf7a"),
    ("RE", 600, 700, 100, 650, "RE_gid650", "#4a3aa7"),
]

STATES = [
    (0, "Awake", "A"),
    (1, "N2", "B"),
    (2, "N3", "C"),
    (3, "REM", "D"),
]

# Display windows, (start_s, duration_s).  The first 10 s of every run is
# initialisation and must be excluded (README_FINK_ET_AL.md); these windows sit
# deep in the steady state.  Window lengths match the paper: 5 s awake, 10 s
# for the sleep stages.
DEFAULT_WINDOWS = {
    0: (100.0, 5.0),
    1: (100.0, 10.0),
    2: (100.0, 10.0),
    3: (100.0, 10.0),
}

# One shared voltage scale for every single-cell trace panel in the figure, so
# the four columns (and the four cell types) are directly comparable.
VLIM = (-100.0, 50.0)
SCALE_MV = 40.0  # length of the voltage scale bar, as in the paper

INK = "#0b0b0b"
MUTED = "#6b6a66"
FAINT = "#d8d7d2"


# ----------------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------------
def raster_path(state: int, seed: int, nhost: int) -> str:
    return os.path.join(REPO, "raster_state%d_seed%d_nhost=%d.txt" % (state, seed, nhost))


def trace_path(tag: str, state: int, seed: int, nhost: int) -> str:
    return os.path.join(
        REPO, "trace_%s_state%d_seed%d_nhost=%d.txt" % (tag, state, seed, nhost)
    )


def load_raster(state, seed, nhost, t0_ms, t1_ms):
    """Return (times_ms, gids) for spikes inside [t0_ms, t1_ms).

    The raster file is written by all MPI ranks and is therefore not time
    sorted; we only ever mask on time, so sorting is unnecessary.
    """
    path = raster_path(state, seed, nhost)
    if not os.path.exists(path):
        return None
    dat = np.loadtxt(path)
    if dat.size == 0:
        return np.empty(0), np.empty(0)
    dat = np.atleast_2d(dat)
    t, g = dat[:, 0], dat[:, 1]
    m = (t >= t0_ms) & (t < t1_ms)
    return t[m], g[m]


def load_trace(tag, state, seed, nhost, t0_ms, t1_ms):
    path = trace_path(tag, state, seed, nhost)
    if not os.path.exists(path):
        return None
    dat = np.loadtxt(path, skiprows=1)
    if dat.size == 0:
        return None
    dat = np.atleast_2d(dat)
    t, v = dat[:, 0], dat[:, 1]
    m = (t >= t0_ms) & (t < t1_ms)
    return t[m], v[m]


def state_available(state, seed, nhost):
    """A state is plottable if its raster and all four traces exist."""
    if not os.path.exists(raster_path(state, seed, nhost)):
        return False
    return all(
        os.path.exists(trace_path(p[5], state, seed, nhost)) for p in POPS
    )


# ----------------------------------------------------------------------------
# Quantitative summary (printed, not plotted) -- backs the qualitative claims
# ----------------------------------------------------------------------------
def describe_state(name, t_ms, gids, t0_ms, t1_ms, state=None, seed=1, nhost=8):
    dur_s = (t1_ms - t0_ms) / 1000.0
    print("\n--- %s  (%.1f-%.1f s, %.0f s window) ---" % (name, t0_ms / 1e3, t1_ms / 1e3, dur_s))
    rates = {}
    for label, lo, hi, n, _tg, _tag, _c in POPS:
        m = (gids >= lo) & (gids < hi)
        rate = m.sum() / float(n) / dur_s
        rates[label] = rate
        print("    %-3s mean firing rate %6.2f Hz  (%d spikes)" % (label, rate, m.sum()))

    # Cortical Up/Down structure: bin PY spikes at 10 ms and measure the
    # fraction of bins that are completely silent (Down states) plus the
    # dominant frequency of the resulting population-rate signal.
    m_py = (gids >= 0) & (gids < 500)
    bins = np.arange(t0_ms, t1_ms + 1e-9, 10.0)
    cnt, _ = np.histogram(t_ms[m_py], bins=bins)
    silent = float((cnt == 0).mean())
    print("    PY population-rate silence (10 ms bins, 0 spikes/500 cells): %.1f%%" % (100 * silent))
    if cnt.size > 8 and cnt.sum() > 0:
        x = cnt - cnt.mean()
        fs = 100.0  # 10 ms bins
        freqs = np.fft.rfftfreq(x.size, d=1.0 / fs)
        psd = np.abs(np.fft.rfft(x * np.hanning(x.size))) ** 2

        def band(lo, hi):
            b = (freqs >= lo) & (freqs < hi)
            return psd[b].sum()

        tot = band(0.2, 30.0)
        if tot > 0:
            print(
                "    PY rate power  slow 0.5-1.5 Hz %5.1f%% | delta 1.5-4 Hz %5.1f%% | "
                "spindle 7-15 Hz %5.1f%%"
                % (100 * band(0.5, 1.5) / tot, 100 * band(1.5, 4.0) / tot, 100 * band(7.0, 15.0) / tot)
            )
            pk = freqs[1:][np.argmax(psd[1:])]
            print("    PY rate peak frequency %.2f Hz" % pk)

    # Sub-threshold membrane polarisation of the four traced cells -- this is
    # what the paper's voltage image would have shown and the raster cannot.
    if state is not None:
        for label, _lo, _hi, _n, _tg, tag, _c in POPS:
            tr = load_trace(tag, state, seed, nhost, t0_ms, t1_ms)
            if tr is None or tr[1].size == 0:
                continue
            v = tr[1]
            sub = v[v < -40.0]  # exclude spikes to get the sub-threshold level
            print(
                "    %-3s traced cell  median Vm %6.1f mV  (sub-threshold %6.1f, "
                "min %6.1f, max %6.1f)"
                % (label, np.median(v), np.median(sub) if sub.size else np.nan, v.min(), v.max())
            )
    return rates


# ----------------------------------------------------------------------------
# Plotting helpers
# ----------------------------------------------------------------------------
def bare(ax, panel=False):
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    if panel:
        # A barely-there ground so an empty panel reads as "this population was
        # silent" rather than "this panel failed to draw".
        ax.set_facecolor("#f6f6f4")
    else:
        ax.patch.set_alpha(0.0)


def draw_column(fig, gs_top, gs_bot, state, name, letter, seed, nhost, win, first_col):
    """Draw one column (raster block + trace block).  Returns True if data drawn."""
    t0_ms, t1_ms = win[0] * 1000.0, (win[0] + win[1]) * 1000.0

    top = GridSpecFromSubplotSpec(
        4, 1, subplot_spec=gs_top, height_ratios=[3.4, 1.0, 1.0, 1.0], hspace=0.10
    )
    bot = GridSpecFromSubplotSpec(4, 1, subplot_spec=gs_bot, hspace=0.10)

    have = state_available(state, seed, nhost)

    if not have:
        ax = fig.add_subplot(gs_top)
        bare(ax)
        ax.text(
            0.5, 0.5,
            "%s\nnot available\n\nno completed state-%d run was on disk\nwhen this figure was built.\n"
            "Re-run this script once\ntrace_*_state%d_* exist and this\ncolumn fills in automatically."
            % (name, state, state),
            ha="center", va="center", fontsize=8.5, color=MUTED, linespacing=1.5,
            transform=ax.transAxes,
        )
        ax.set_title("%s   %s" % (letter, name), fontsize=11, color=MUTED, pad=6, loc="left")
        ax2 = fig.add_subplot(gs_bot)
        bare(ax2)
        for a in (ax, ax2):
            a.add_patch(
                plt.Rectangle(
                    (0.02, 0.02), 0.96, 0.96, transform=a.transAxes, fill=False,
                    ec=FAINT, lw=0.8, ls=(0, (4, 4)),
                )
            )
        return False, None

    t_ms, gids = load_raster(state, seed, nhost, t0_ms, t1_ms)

    # ---- top block: one raster per population -----------------------------
    for i, (label, lo, hi, n, _tg, _tag, colour) in enumerate(POPS):
        ax = fig.add_subplot(top[i])
        bare(ax, panel=True)
        m = (gids >= lo) & (gids < hi)
        ax.plot(
            (t_ms[m] - t0_ms) / 1000.0, gids[m] - lo,
            linestyle="none", marker="o", markersize=0.55, markeredgewidth=0,
            color=colour, alpha=0.85, rasterized=True,
        )
        if m.sum() == 0:
            ax.text(
                0.5, 0.5, "no spikes", transform=ax.transAxes, ha="center", va="center",
                fontsize=7.5, color=MUTED, style="italic",
            )
        ax.set_xlim(0, win[1])
        ax.set_ylim(-0.02 * n, 1.02 * n)
        if first_col:
            ax.text(
                -0.045, 0.5, label, transform=ax.transAxes, ha="right", va="center",
                fontsize=9.5, color=INK, fontweight="bold",
            )
        if i == 0:
            ax.set_title(
                "%s   %s" % (letter, name), fontsize=11, color=INK, pad=17, loc="left"
            )
            ax.text(
                0.0, 1.035,
                "%.0f–%.0f s of the run  ·  %.0f s window"
                % (win[0], win[0] + win[1], win[1]),
                transform=ax.transAxes, ha="left", va="bottom",
                fontsize=7.8, color=MUTED,
            )

    # ---- bottom block: representative single-cell voltages -----------------
    for i, (label, _lo, _hi, _n, tgid, tag, colour) in enumerate(POPS):
        ax = fig.add_subplot(bot[i])
        bare(ax, panel=True)
        tr = load_trace(tag, state, seed, nhost, t0_ms, t1_ms)
        if tr is not None:
            tt, vv = tr
            ax.plot((tt - t0_ms) / 1000.0, vv, lw=0.5, color=colour, rasterized=True)
        ax.set_xlim(0, win[1])
        ax.set_ylim(*VLIM)
        # -60 mV fiducial, identical in every panel, so the reader can see the
        # tonic hyperpolarisation across stages without axis ticks.
        ax.axhline(-60.0, color="#e2e1dc", lw=0.6, zorder=0)
        if first_col:
            ax.text(
                -0.045, 0.5, "%s #%d" % (label, tgid - _lo), transform=ax.transAxes,
                ha="right", va="center", fontsize=8.5, color=INK,
            )
    return True, (t_ms, gids)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--nhost", type=int, default=8)
    ap.add_argument("--out", default=OUT_PATH)
    ap.add_argument("--dpi", type=int, default=200)
    ap.add_argument(
        "--win", action="append", default=[],
        metavar="STATE:START_S:DUR_S",
        help="override a display window, e.g. --win 2:60:10",
    )
    args = ap.parse_args()

    windows = dict(DEFAULT_WINDOWS)
    for spec in args.win:
        s, a, b = spec.split(":")
        windows[int(s)] = (float(a), float(b))

    # width of each column is proportional to its window length, so 1 s of
    # simulated time is the same number of pixels in every column (the paper
    # uses a single 1 s scale bar for all four panels).
    wr = [windows[s][1] for s, _, _ in STATES]

    fig = plt.figure(figsize=(15.0, 8.6), facecolor="white", dpi=args.dpi)
    gs = GridSpec(
        2, 4, figure=fig, height_ratios=[1.10, 1.0], width_ratios=wr,
        left=0.055, right=0.985, top=0.857, bottom=0.085, wspace=0.075, hspace=0.115,
    )

    fig.text(
        0.055, 0.965,
        "Characteristic patterns of electrical activity during different sleep stages",
        fontsize=13.5, color=INK, fontweight="bold", ha="left", va="center",
    )
    fig.text(
        0.055, 0.938,
        "NEURON port of Krishnan et al. 2016 eLife, Figure 3  ·  seed %d, %d MPI ranks  ·  "
        "120 s static-state runs, windows taken from the steady state (first 10 s discarded)"
        % (args.seed, args.nhost),
        fontsize=8.6, color=MUTED, ha="left", va="center",
    )

    print("Figure 3 reproduction -- seed %d, nhost %d" % (args.seed, args.nhost))
    print("Window lengths follow the paper: 5 s awake, 10 s for N2/N3/REM.")

    drawn = []
    missing = []
    for col, (state, name, letter) in enumerate(STATES):
        ok, dat = draw_column(
            fig, gs[0, col], gs[1, col], state, name, letter,
            args.seed, args.nhost, windows[state], first_col=(col == 0),
        )
        if ok:
            drawn.append(name)
            t0 = windows[state][0] * 1000.0
            t1 = t0 + windows[state][1] * 1000.0
            describe_state(name, dat[0], dat[1], t0, t1, state, args.seed, args.nhost)
        else:
            missing.append(name)
            print("\n--- %s: no data (state %d not on disk) ---" % (name, state))

    # ---- scale bars, bottom-right, as in the paper ------------------------
    # 1 s horizontally, 40 mV vertically, drawn in figure coordinates against
    # the axes of the last column.
    last_state = STATES[-1][0]
    dur = windows[last_state][1]
    BAR_W, BAR_H = 0.085, 0.058
    bar_ax = fig.add_axes([0.895, 0.014, BAR_W, BAR_H])
    bare(bar_ax)
    bar_ax.set_xlim(0, 1)
    bar_ax.set_ylim(0, 1)
    # horizontal 1 s bar, scaled to the width of one column's time axis
    col_w = gs[1, 3].get_position(fig).width
    frac_1s = (1.0 / dur) * col_w / BAR_W
    bar_ax.plot([0.0, min(frac_1s, 1.0)], [0.72, 0.72], color=INK, lw=1.6, solid_capstyle="butt")
    bar_ax.text(min(frac_1s, 1.0) / 2, 0.86, "1 s", ha="center", va="bottom", fontsize=8.5, color=INK)
    # vertical 40 mV bar, scaled to the height of one trace panel
    trace_h = gs[1, 3].get_position(fig).height / 4.0
    frac_mv = (SCALE_MV / (VLIM[1] - VLIM[0])) * trace_h / BAR_H
    bar_ax.plot([0.0, 0.0], [0.05, 0.05 + min(frac_mv, 0.6)], color=INK, lw=1.6, solid_capstyle="butt")
    bar_ax.text(0.04, 0.05 + min(frac_mv, 0.6) / 2, "%.0f mV" % SCALE_MV, ha="left", va="center",
                fontsize=8.5, color=INK)

    for y, line in [
        (0.058, "Top block — spike raster per population, one dot per spike; y = neuron index "
                "within its population (PY 0–499, IN / TC / RE 0–99)."),
        (0.036, "The paper shows a membrane-voltage image here; this port does not store Vₘ for "
                "all 800 cells, so sub-threshold polarisation appears only in the traces below."),
        (0.014, "Bottom block — Vₘ of the paper's cells (PY #250, IN #50, TC #50, RE #50), one "
                "shared −%.0f to +%.0f mV scale in every panel; grey line = −60 mV."
                % (abs(VLIM[0]), VLIM[1])),
    ]:
        fig.text(0.055, y, line, fontsize=8.2, color=MUTED, ha="left", va="center")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=args.dpi, facecolor="white")
    print("\nWrote %s" % args.out)
    print("Columns with data : %s" % (", ".join(drawn) if drawn else "none"))
    if missing:
        print("Columns missing   : %s" % ", ".join(missing))
    return 0


if __name__ == "__main__":
    sys.exit(main())
