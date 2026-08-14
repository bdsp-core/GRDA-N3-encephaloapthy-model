"""Reproduce Figure 11 of Krishnan et al. 2016 eLife (DOI 10.7554/eLife.18607),
"Effect of additional parameter changes on the stage transitions."

Panel A  network activity of the extended model, in which ACh scales the
         pyramidal dendritic KCa, Km, and NaP conductances in addition to the
         baseline K-leak / Ih / GABA / AMPA controls
         (FINK_SLEEPSTATES=1 FINK_FIG11=1)
Panel B  the reduced model, in which only those three dendritic conductances
         follow the sleep schedule while K leak, Ih, GABA, and AMPA stay at
         their awake values (FINK_SLEEPSTATES=1 FINK_FIG11=2).  In the paper
         this model fails to produce stage transitions.
Panel C  zoom of panel A during a spindle (default: inside N2)
Panel D  zoom of panel A during a slow oscillation (default: inside N3)

Each network panel stacks, from top to bottom, PY / IN / TC / RE spike rasters
and the mean cortical membrane potential, matching the paper's ordering.

Required input files, produced by bazh_net.py in sleep-cycle mode:

    raster_sweep_seed<SEED>_nhost=<N>.txt
    vcort_sweep_seed<SEED>_nhost=<N>.txt

The cycle needs the full 360000 ms run; the last scheduled transition ends at
307500 ms.  Because both variants write the same `sweep_seed<SEED>` file names,
run them in separate directories (or rename the outputs) and pass the two
directories with --extended-dir and --reduced-dir.

Example:

    mkdir -p out_fig11a out_fig11b
    (cd out_fig11a && FINK_SLEEPSTATES=1 FINK_FIG11=1 FINK_SEED=1 \
        FINK_DURATION_MS=360000 mpiexec -n 8 ../.venv/bin/nrniv -mpi \
        -python ../bazh_net.py)
    (cd out_fig11b && FINK_SLEEPSTATES=1 FINK_FIG11=2 FINK_SEED=1 \
        FINK_DURATION_MS=360000 mpiexec -n 8 ../.venv/bin/nrniv -mpi \
        -python ../bazh_net.py)
    ./.venv/bin/python plot_paper_fig11.py \
        --extended-dir out_fig11a --reduced-dir out_fig11b --nhost 8
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

# Cell-type gid ranges, from Net.setGids in network_class.py: pyramidal cells
# come first, then inhibitory interneurons, then RE, then TC.
NPYR, NINH, NRE, NTC = 500, 100, 100, 100
GID_RANGES = {
    "PY": (0, NPYR),
    "IN": (NPYR, NPYR + NINH),
    "RE": (NPYR + NINH, NPYR + NINH + NRE),
    "TC": (NPYR + NINH + NRE, NPYR + NINH + NRE + NTC),
}
# The paper plots PY, IN, TC, RE from top to bottom.
ROW_ORDER = ("PY", "IN", "TC", "RE")
ROW_COLORS = {"PY": "#1f4e79", "IN": "#a03030", "TC": "#256d3a", "RE": "#6a3d9a"}

# The callback in config.py sums seg.v over every segment of every cortical
# cell, so dividing by the segment count recovers a mean membrane potential.
# PyrCell and InhCell each have a 1-segment soma and a 1-segment dendrite.
N_CORT_SEGMENTS = 2 * (NPYR + NINH)

DT_MS = 0.025  # h.dt in bazh_net.py

# Sleep-cycle schedule; these mirror the do_sleepstates block of config.py.
SCHEDULE = [
    ("Awake", 0.0, 80000.0),
    ("N2", 97500.0, 150000.0),
    ("N3", 167500.0, 220000.0),
    ("REM", 237500.0, 290000.0),
    ("N2", 307500.0, 360000.0),
]


def load_column(path):
    """Read a one-number-per-line output file quickly."""
    with open(path) as handle:
        return np.array(handle.read().split(), dtype=float)


def load_raster(path):
    """Read a 'time  gid' raster file into two arrays."""
    with open(path) as handle:
        raw = np.array(handle.read().split(), dtype=float)
    if raw.size == 0:
        return np.zeros(0), np.zeros(0)
    raw = raw.reshape(-1, 2)
    return raw[:, 0], raw[:, 1]


class RunData(object):
    """Spike times and mean cortical Vm for one simulation."""

    def __init__(self, directory, seed, nhost, label):
        directory = Path(directory)
        self.label = label
        self.raster_path = directory / ("raster_sweep_seed%d_nhost=%g.txt" % (seed, nhost))
        self.vcort_path = directory / ("vcort_sweep_seed%d_nhost=%g.txt" % (seed, nhost))
        missing = [p for p in (self.raster_path, self.vcort_path) if not p.exists()]
        if missing:
            raise SystemExit(
                "Missing input file(s):\n  %s\nRun the sleep cycle first; see the "
                "module docstring of this script." % "\n  ".join(str(p) for p in missing)
            )
        self.spike_t, self.spike_gid = load_raster(self.raster_path)
        vcort = load_column(self.vcort_path)
        self.vm = vcort / N_CORT_SEGMENTS
        self.vm_t = np.arange(self.vm.size) * DT_MS

    def spikes_in(self, kind, t0, t1):
        lo, hi = GID_RANGES[kind]
        sel = (
            (self.spike_gid >= lo)
            & (self.spike_gid < hi)
            & (self.spike_t >= t0)
            & (self.spike_t <= t1)
        )
        return self.spike_t[sel], self.spike_gid[sel] - lo

    def vm_in(self, t0, t1):
        i0 = max(0, int(t0 / DT_MS))
        i1 = min(self.vm.size, int(t1 / DT_MS) + 1)
        return self.vm_t[i0:i1], self.vm[i0:i1]


def decimate_for_plot(t, y, max_points):
    """Envelope-preserving reduction so long traces stay legible and light."""
    if t.size <= max_points or max_points <= 0:
        return t, y
    step = int(np.ceil(t.size / float(max_points)))
    usable = (t.size // step) * step
    tt = t[:usable].reshape(-1, step)
    yy = y[:usable].reshape(-1, step)
    # interleave the per-bin min and max so oscillation envelopes survive
    out_t = np.repeat(tt[:, 0], 2)
    out_y = np.empty(out_t.size)
    out_y[0::2] = yy.min(axis=1)
    out_y[1::2] = yy.max(axis=1)
    return out_t, out_y


def draw_network_panel(fig, spec, run, t0, t1, title, show_stage_labels,
                       max_spikes, max_points, vm_lim):
    """Draw the PY/IN/TC/RE rasters plus mean Vm inside one gridspec cell."""
    inner = GridSpecFromSubplotSpec(
        5, 1, subplot_spec=spec, hspace=0.12, height_ratios=[1.6, 0.8, 0.8, 0.8, 1.2]
    )
    axes = []
    rng = np.random.default_rng(0)  # only used to thin dense rasters, deterministically

    for row, kind in enumerate(ROW_ORDER):
        ax = fig.add_subplot(inner[row])
        st, sg = run.spikes_in(kind, t0, t1)
        if max_spikes and st.size > max_spikes:
            keep = rng.choice(st.size, max_spikes, replace=False)
            keep.sort()
            st, sg = st[keep], sg[keep]
        # dots have to shrink as the displayed window grows or the raster
        # saturates into a solid block
        markersize = 0.35 if (t1 - t0) > 60000.0 else 1.4
        ax.plot(
            st / 1000.0,
            sg,
            linestyle="none",
            marker=".",
            markersize=markersize,
            alpha=0.7,
            color=ROW_COLORS[kind],
            rasterized=True,
        )
        span = GID_RANGES[kind][1] - GID_RANGES[kind][0]
        ax.set_ylim(-1, span)
        ax.set_xlim(t0 / 1000.0, t1 / 1000.0)
        ax.set_yticks([])
        ax.set_xticks([])
        for side in ("top", "right", "bottom", "left"):
            ax.spines[side].set_visible(False)
        ax.set_ylabel(kind, rotation=0, ha="right", va="center", fontsize=8)
        axes.append(ax)

    ax = fig.add_subplot(inner[4])
    vt, vy = run.vm_in(t0, t1)
    vt, vy = decimate_for_plot(vt, vy, max_points)
    ax.plot(vt / 1000.0, vy, color="0.15", linewidth=0.4, rasterized=True)
    ax.set_xlim(t0 / 1000.0, t1 / 1000.0)
    if vm_lim is not None:
        ax.set_ylim(*vm_lim)
    ax.set_ylabel("mean\nV$_m$ (mV)", rotation=0, ha="right", va="center", fontsize=8)
    ax.tick_params(labelsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_xlabel("time (s)", fontsize=8)
    axes.append(ax)

    if show_stage_labels:
        for name, s0, s1 in SCHEDULE:
            mid = 0.5 * (max(s0, t0) + min(s1, t1))
            if not (t0 <= mid <= t1):
                continue
            axes[0].text(
                mid / 1000.0,
                1.06,
                name,
                transform=axes[0].get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=8,
            )
        for _, s0, s1 in SCHEDULE:
            for boundary in (s0, s1):
                if t0 < boundary < t1:
                    for a in axes:
                        a.axvline(boundary / 1000.0, color="0.75", linewidth=0.5, zorder=0)

    axes[0].set_title(title, fontsize=9, pad=16 if show_stage_labels else 4)
    return axes


def default_window(stage, width_ms):
    """Centre a window of the requested width inside the named stage."""
    for name, s0, s1 in SCHEDULE:
        if name == stage:
            mid = 0.5 * (s0 + s1)
            return mid - width_ms / 2.0, mid + width_ms / 2.0
    raise ValueError("unknown stage %r" % stage)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--extended-dir", default="out_fig11a",
                        help="directory holding the FINK_FIG11=1 cycle output")
    parser.add_argument("--reduced-dir", default="out_fig11b",
                        help="directory holding the FINK_FIG11=2 cycle output")
    parser.add_argument("--seed", type=int, default=1, help="FINK_SEED used for the runs")
    parser.add_argument("--nhost", type=int, default=8, help="MPI rank count of the runs")
    parser.add_argument("--duration-ms", type=float, default=360000.0,
                        help="simulated duration to display (default 360000)")
    parser.add_argument("--spindle-window-ms", type=float, nargs=2, default=None,
                        help="start and end of the panel C zoom (default: 5 s inside N2)")
    parser.add_argument("--slow-window-ms", type=float, nargs=2, default=None,
                        help="start and end of the panel D zoom (default: 5 s inside N3)")
    parser.add_argument("--zoom-width-ms", type=float, default=5000.0,
                        help="width of the default zoom windows (default 5000)")
    parser.add_argument("--vm-lim", type=float, nargs=2, default=(-75.0, -50.0),
                        help="y limits of the mean-Vm traces (default -75 -50)")
    parser.add_argument("--max-spikes", type=int, default=200000,
                        help="thin each raster to at most this many dots (0 disables)")
    parser.add_argument("--max-points", type=int, default=6000,
                        help="envelope-decimate each Vm trace to this many bins")
    parser.add_argument("--out", default="docs/figures/paper_fig11_intrinsic.png",
                        help="output image path")
    parser.add_argument("--dpi", type=int, default=200)
    args = parser.parse_args()

    extended = RunData(args.extended_dir, args.seed, args.nhost, "extended")
    reduced = RunData(args.reduced_dir, args.seed, args.nhost, "reduced")

    spindle = tuple(args.spindle_window_ms) if args.spindle_window_ms else default_window(
        "N2", args.zoom_width_ms
    )
    slow = tuple(args.slow_window_ms) if args.slow_window_ms else default_window(
        "N3", args.zoom_width_ms
    )

    fig = plt.figure(figsize=(13.0, 10.0))
    outer = GridSpec(2, 2, figure=fig, hspace=0.34, wspace=0.20,
                     left=0.08, right=0.98, top=0.855, bottom=0.06)

    draw_network_panel(fig, outer[0, 0], extended, 0.0, args.duration_ms,
                       "A  Extended model\n(ACh also scales dendritic K$_{Ca}$, K$_m$, Na$_p$)",
                       True, args.max_spikes, args.max_points, tuple(args.vm_lim))
    draw_network_panel(fig, outer[0, 1], reduced, 0.0, args.duration_ms,
                       "B  Reduced model\n(only K$_{Ca}$, K$_m$, Na$_p$ change; core fixed at awake)",
                       True, args.max_spikes, args.max_points, tuple(args.vm_lim))
    draw_network_panel(fig, outer[1, 0], extended, spindle[0], spindle[1],
                       "C  Extended model, spindle epoch (%.1f-%.1f s)"
                       % (spindle[0] / 1000.0, spindle[1] / 1000.0),
                       False, args.max_spikes, args.max_points, tuple(args.vm_lim))
    draw_network_panel(fig, outer[1, 1], extended, slow[0], slow[1],
                       "D  Extended model, slow-oscillation epoch (%.1f-%.1f s)"
                       % (slow[0] / 1000.0, slow[1] / 1000.0),
                       False, args.max_spikes, args.max_points, tuple(args.vm_lim))

    fig.suptitle(
        "Krishnan et al. 2016 Figure 11: effect of ACh acting on dendritic "
        "intrinsic currents\n(seed %d, %d ranks, NEURON port)" % (args.seed, args.nhost),
        fontsize=11,
        y=0.985,
        va="top",
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=args.dpi)
    print("wrote %s" % out_path)


if __name__ == "__main__":
    main()
