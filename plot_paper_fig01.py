"""Reproduce Figure 1 of Krishnan et al. 2016 eLife 5:e18607 from this port.

Figure 1 ("Model design") is a picture of the model definition, not of a
simulation, so nothing here is simulated.  Every number is read out of
``config.py`` (the per-state parameter values and the transition times) and
every schedule shape is taken from the ``h.Vector.play`` calls in
``bazh_net.py``, which is what actually drives the Wake -> N2 -> N3 -> REM -> N2
cycle during a run.

Panels
------
A  Network schematic: 500 PY / 100 IN in cortex, 100 TC / 100 RE in thalamus,
   with the connection classes drawn the way the paper draws them.  The paper's
   panel is hand-drawn; this is a programmatic equivalent with the same
   connectivity, and the pathway list is taken from ``Net.connectCells`` in
   ``network_class.py``.
B  Ih activation curve per vigilance state.  ``mod/Ih.mod`` defines
       h_inf = 1 / (1 + exp((v + 75 + fac_gh_TC) / 5.5))
   and ``fac_gh_TC`` is the Shift_HA knob, played from ``gh_TC_*`` in config.
C  Cortical conductance schedule (K-leak, AMPA, GABA), normalized to awake.
D  Thalamic conductance schedule (K-leak-TC, K-leak-RE, GABA), same.

Run with::

    .venv/bin/python plot_paper_fig01.py

Importing ``config`` instantiates a NEURON ParallelContext; that is cheap and
runs no simulation.  ``FINK_SLEEPSTATES`` is forced on before the import
because the transition times only exist inside config's ``do_sleepstates``
branch.
"""

import os
import re

# config.py only defines awake_to_s2_start / ... inside `if do_sleepstates:`,
# so ask for the full sleep-stage cycle before importing it.
os.environ["FINK_SLEEPSTATES"] = "1"

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch

import config

HERE = os.path.dirname(os.path.abspath(__file__))
OUTFILE = os.path.join(HERE, "docs", "figures", "paper_fig01_model_design.png")

# ---------------------------------------------------------------------------
# Palette chosen to match the printed figure.
# ---------------------------------------------------------------------------
C_KLEAK = "#1f6fb4"  # blue
C_AMPA = "#bf3f22"  # brick red
C_GABA = "#f2a900"  # amber
C_PORT = "#7a7a7a"  # gray, for traces that are port-specific

STATE_COLORS = {
    "Awake": "#000000",
    "N2": "#e8192c",
    "N3": "#c2186b",
    "REM": "#2fb14a",
}

CELL_COLORS = {
    "PY": "#3aa63a",
    "IN": "#e03030",
    "TC": "#2b3a9c",
    "RE": "#ec86b4",
}
BOX_FACE = "#d7e7ef"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.linewidth": 1.3,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 8.5,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.facecolor": "white",
    }
)


# ---------------------------------------------------------------------------
# Panel B: Ih activation
# ---------------------------------------------------------------------------
def ih_activation(v, fac_gh_tc):
    """Steady-state Ih activation, mod/Ih.mod PROCEDURE evaluate_fct."""
    return 1.0 / (1.0 + np.exp((v + 75.0 + fac_gh_tc) / 5.5))


def verify_ih_formula():
    """Fail loudly if mod/Ih.mod no longer matches the formula hard-coded above."""
    path = os.path.join(HERE, "mod", "Ih.mod")
    with open(path) as fh:
        text = fh.read()
    match = re.search(
        r"h_inf\s*=\s*1\s*/\s*\(\s*1\s*\+\s*exp\(\(v\+(\d+(?:\.\d+)?)"
        r"\+fac_gh_TC\)/(\d+(?:\.\d+)?)\)\s*\)",
        text.replace(" ", "").replace("\t", ""),
    )
    if match is None:
        raise RuntimeError("could not find the h_inf expression in mod/Ih.mod")
    v_off, slope = float(match.group(1)), float(match.group(2))
    if (v_off, slope) != (75.0, 5.5):
        raise RuntimeError(
            "mod/Ih.mod h_inf constants changed to (%g, %g); update this script"
            % (v_off, slope)
        )


# ---------------------------------------------------------------------------
# Panels C/D: the conductance schedule that bazh_net.py plays into the model
# ---------------------------------------------------------------------------
# bazh_net.py builds a single time vector
#     t = [awake_to_s2_start, awake_to_s2_end,
#          s2_to_s3_start,    s2_to_s3_end,
#          s3_to_rem_start,   s3_to_rem_end,
#          rem_to_s2_start,   rem_to_s2_end, rem_to_s2_end + dt]
# and plays every parameter with the value pattern
#     [awake, s2, s2, s3, s3, rem, rem, s2, s2]
# using Vector.play(..., continuous=True), i.e. linear interpolation between
# the listed points and a hold outside the range.
TRANSITIONS_MS = [
    config.awake_to_s2_start,
    config.awake_to_s2_end,
    config.s2_to_s3_start,
    config.s2_to_s3_end,
    config.s3_to_rem_start,
    config.s3_to_rem_end,
    config.rem_to_s2_start,
    config.rem_to_s2_end,
]
T_END_S = 360.0  # full Figure-2 cycle; last transition ends at 307.5 s
STATE_EDGES_S = [0.0] + [TRANSITIONS_MS[i] / 1000.0 for i in (0, 2, 4, 6)] + [T_END_S]
STATE_NAMES = ["Awake", "N2", "N3", "REM", "N2"]


def schedule(awake, s2, s3, rem):
    """Return (t_sec, percent-of-awake) for one played parameter."""
    t = np.array([0.0] + [ms / 1000.0 for ms in TRANSITIONS_MS] + [T_END_S])
    y = np.array([awake, awake, s2, s2, s3, s3, rem, rem, s2, s2], dtype=float)
    return t, 100.0 * y / awake


# ---------------------------------------------------------------------------
# Panel A helpers: connectors drawn in display space so curvature and
# arrow/ball ends look right regardless of the axes aspect ratio.
# ---------------------------------------------------------------------------
def connector(ax, p0, p1, color="k", rad=0.0, head="arrow", shrink_a=11.0,
              shrink_b=11.0, lw=1.15, zorder=4, ball_ms=5.2, head_scale=10.0):
    tr = ax.transData
    a = np.asarray(tr.transform(p0), dtype=float)
    b = np.asarray(tr.transform(p1), dtype=float)
    d = b - a
    length = float(np.hypot(*d))
    if length < 1e-9:
        return
    perp = np.array([-d[1], d[0]]) / length
    ctrl = 0.5 * (a + b) + rad * length * perp

    t = np.linspace(0.0, 1.0, 400)[:, None]
    pts = (1 - t) ** 2 * a + 2 * (1 - t) * t * ctrl + t**2 * b
    seg = np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1]))
    s = np.concatenate([[0.0], np.cumsum(seg)])
    keep = (s >= shrink_a) & (s <= s[-1] - shrink_b)
    pts = pts[keep]
    if len(pts) < 3:
        return

    xy = tr.inverted().transform(pts)
    ax.plot(xy[:, 0], xy[:, 1], color=color, lw=lw, solid_capstyle="round",
            zorder=zorder)
    if head == "arrow":
        ax.annotate(
            "",
            xy=tuple(xy[-1]),
            xytext=tuple(xy[-4]),
            arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, shrinkA=0,
                            shrinkB=0, mutation_scale=head_scale),
            zorder=zorder,
        )
    elif head == "ball":
        ax.plot(xy[-1, 0], xy[-1, 1], "o", color=color, ms=ball_ms,
                markeredgewidth=0, zorder=zorder + 1)


def draw_schematic(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    box_x0, box_w = 0.165, 0.815
    cortex_y0, thal_y0, box_h = 0.585, 0.175, 0.345
    for y0, label in ((cortex_y0, "Cortex"), (thal_y0, "Thalamus")):
        ax.add_patch(
            FancyBboxPatch(
                (box_x0, y0), box_w, box_h,
                boxstyle="round,pad=0,rounding_size=0.012",
                facecolor=BOX_FACE, edgecolor="none", zorder=0,
            )
        )
        ax.text(box_x0 - 0.018, y0 + box_h / 2.0, label, ha="right", va="center",
                fontsize=9.5)

    rows = {
        "PY": cortex_y0 + 0.245,
        "IN": cortex_y0 + 0.085,
        "TC": thal_y0 + 0.245,
        "RE": thal_y0 + 0.085,
    }
    for name, y in rows.items():
        ax.text(box_x0 + 0.035, y, name, ha="left", va="center", fontsize=9.5,
                zorder=2)

    # TC and RE are staggered relative to each other, as in the printed panel,
    # so that the cortex->RE arrow does not have to pass through a TC cell.
    py = np.linspace(0.30, 0.955, 10)
    inn = np.linspace(0.33, 0.90, 4)
    tc = np.linspace(0.320, 0.880, 4)
    re = np.linspace(0.355, 0.915, 4)
    xs = {"PY": py, "IN": inn, "TC": tc, "RE": re}
    for name, x in xs.items():
        ax.plot(x, np.full_like(x, rows[name]), "o", color=CELL_COLORS[name],
                ms=9.5, markeredgewidth=0, zorder=6, clip_on=False)

    def P(name, i):
        return (xs[name][i], rows[name])

    black, blue, red = "#000000", "#2b3a9c", "#e03030"

    # --- intracortical -----------------------------------------------------
    # PY->PY: AMPA_D2 + NMDA_D1 (radius 5 each); arcs above the PY row
    connector(ax, P("PY", 0), P("PY", 2), blue, rad=0.42)
    connector(ax, P("PY", 1), P("PY", 3), blue, rad=0.42)
    # PY->IN: AMPA_D2 + NMDA_D1 (radius 1 each)
    connector(ax, P("PY", 1), P("IN", 0), blue, rad=0.16)
    # IN->PY: GABA_A_D2 (radius 5)
    connector(ax, P("IN", 0), P("PY", 0), red, rad=0.16, head="ball")

    # --- thalamus -> cortex (AMPA) ----------------------------------------
    connector(ax, P("TC", 1), P("PY", 4), black)   # TC->PY, radius 10
    connector(ax, P("TC", 1), P("PY", 6), black)
    connector(ax, P("TC", 1), P("IN", 1), black)   # TC->IN, radius 2

    # --- cortex -> thalamus (AMPA) ----------------------------------------
    connector(ax, P("PY", 7), P("TC", 2), black)   # PY->TC, radius 10
    connector(ax, P("PY", 8), P("RE", 3), black, rad=0.10)  # PY->RE, radius 8

    # --- intrathalamic -----------------------------------------------------
    connector(ax, P("TC", 1), P("RE", 0), black)   # TC->RE AMPA, radius 8
    connector(ax, P("TC", 1), P("RE", 2), black)
    connector(ax, P("RE", 1), P("TC", 1), red)     # RE->TC GABA-A + GABA-B
    connector(ax, P("RE", 3), P("TC", 3), red)
    # RE->RE GABA-A (radius 5), drawn as arcs below the RE row
    connector(ax, P("RE", 0), P("RE", 1), red, rad=-0.32, head="ball")
    connector(ax, P("RE", 2), P("RE", 3), red, rad=-0.32, head="ball")

    # --- legend ------------------------------------------------------------
    seg = 0.075
    entries = [
        (0.055, 0.085, black, "arrow", "AMPA"),
        (0.505, 0.085, blue, "arrow", "AMPA + NMDA"),
        (0.055, 0.015, red, "ball", "GABA-A"),
        (0.505, 0.015, red, "arrow", "GABA-A + GABA-B"),
    ]
    for x0, y0, color, head, label in entries:
        connector(ax, (x0, y0), (x0 + seg, y0), color, head=head, shrink_a=0,
                  shrink_b=0, lw=1.3, ball_ms=5.0)
        ax.text(x0 + seg + 0.022, y0, label, ha="left", va="center", fontsize=8.5)


# ---------------------------------------------------------------------------
def draw_state_bands(ax, ymax):
    """Vertical separators + Awake/N2/N3/REM/N2 labels, as in the paper."""
    for edge in STATE_EDGES_S[1:-1]:
        ax.axvline(edge, color="black", lw=1.2, zorder=1)
    for name, lo, hi in zip(STATE_NAMES, STATE_EDGES_S[:-1], STATE_EDGES_S[1:]):
        ax.text(0.5 * (lo + hi), ymax * 1.035, name, ha="center", va="bottom",
                fontsize=9.5, fontweight="bold")


def style_schedule_axis(ax, ymax):
    ax.set_xlim(0, T_END_S)
    ax.set_ylim(0, ymax)
    ax.set_xticks([0, 100, 200, 300])
    ax.set_yticks([0, 100, 200])
    ax.set_xlabel("Time (sec)", fontweight="bold")
    ax.set_ylabel("Normalized\nConductance (%)", fontweight="bold")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(direction="out", length=4, width=1.3)


# ---------------------------------------------------------------------------
def main():
    verify_ih_formula()

    fig = plt.figure(figsize=(11.0, 10.4))
    gs = GridSpec(
        3, 2, figure=fig,
        height_ratios=[1.08, 1.0, 1.0], width_ratios=[1.12, 1.0],
        left=0.075, right=0.985, top=0.955, bottom=0.105,
        hspace=0.62, wspace=0.26,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])
    ax_d = fig.add_subplot(gs[2, :])

    for ax, label, dx in ((ax_a, "A", -0.06), (ax_b, "B", -0.10),
                          (ax_c, "C", -0.055), (ax_d, "D", -0.055)):
        ax.text(dx, 1.10, label, transform=ax.transAxes, fontsize=15,
                fontweight="bold", va="top", ha="left")

    # ---------------- A ----------------
    draw_schematic(ax_a)

    # ---------------- B ----------------
    v = np.linspace(-100.0, -40.0, 1201)
    shifts = [
        ("Awake", config.gh_TC_awake),
        ("N2", config.gh_TC_s2),
        ("N3", config.gh_TC_s3),
        ("REM", config.gh_TC_rem),
    ]
    for name, fac in shifts:
        ax_b.plot(v, ih_activation(v, fac), color=STATE_COLORS[name], lw=1.9,
                  label=name, zorder=3)

    v_half_rem = -(75.0 + config.gh_TC_rem)
    v_half_awake = -(75.0 + config.gh_TC_awake)
    ax_b.annotate(
        "", xy=(v_half_awake, 0.5), xytext=(v_half_rem, 0.5),
        arrowprops=dict(arrowstyle="-|>", color="black", lw=1.3,
                        mutation_scale=11), zorder=5,
    )
    ax_b.text(v_half_rem - 1.5, 0.5, "Shift$_{\\mathrm{HA}}$", ha="right",
              va="center", fontsize=9.5)
    ax_b.text(
        0.03, 0.06,
        "fac_gh_TC = %g / %g / %g / %g mV" % tuple(f for _, f in shifts),
        transform=ax_b.transAxes, fontsize=7.5, color="#555555",
    )

    ax_b.set_xlim(-100, -40)
    ax_b.set_ylim(-0.02, 1.06)
    ax_b.set_xticks([-100, -80, -60, -40])
    ax_b.set_yticks([0, 0.5, 1.0])
    ax_b.set_yticklabels(["0", "0.5", "1.0"])
    ax_b.set_xlabel("V (mV)", fontweight="bold")
    ax_b.set_ylabel("H-Current Activation", fontweight="bold")
    for side in ("top", "right"):
        ax_b.spines[side].set_visible(False)
    ax_b.tick_params(direction="out", length=4, width=1.3)
    ax_b.legend(loc="upper right", frameon=False, handlelength=1.6,
                labelspacing=0.35, borderpad=0.2)

    # ---------------- C: cortical ----------------
    # gkl_pyr_* and gkl_inh_* have identical awake-normalized profiles
    # (both are awake * [1, s2_scale, s3_scale, 0.9]), so one K-leak trace
    # covers both, exactly as the paper draws it.
    t, kleak_pyr = schedule(config.gkl_pyr_awake, config.gkl_pyr_s2,
                            config.gkl_pyr_s3, config.gkl_pyr_rem)
    _, kleak_inh = schedule(config.gkl_inh_awake, config.gkl_inh_s2,
                            config.gkl_inh_s3, config.gkl_inh_rem)
    assert np.allclose(kleak_pyr, kleak_inh)
    _, ampa_cort = schedule(config.awake_AMPA_cort, config.s2_AMPA_cort,
                            config.s3_AMPA_cort, config.rem_AMPA_cort)
    _, ampa_pyrpyr = schedule(config.awake_AMPA_pyrpyr, config.s2_AMPA_pyrpyr,
                              config.s3_AMPA_pyrpyr, config.rem_AMPA_pyrpyr)
    _, gaba_cort = schedule(config.awake_GABA_D2, config.s2_GABA_D2,
                            config.s3_GABA_D2, config.rem_GABA_D2)

    ymax_c = 285.0
    draw_state_bands(ax_c, ymax_c)
    ax_c.plot(t, ampa_pyrpyr, color=C_AMPA, lw=1.6, ls=(0, (4, 2.4)),
              zorder=2, label="AMPA (PY$\\rightarrow$PY, this port)")
    ax_c.plot(t, kleak_pyr, color=C_KLEAK, lw=2.1, zorder=3, label="K-leak")
    ax_c.plot(t, ampa_cort, color=C_AMPA, lw=2.1, zorder=3, label="AMPA")
    ax_c.plot(t, gaba_cort, color=C_GABA, lw=2.1, zorder=3, label="GABA")
    style_schedule_axis(ax_c, ymax_c)
    ax_c.set_title("Cortical Manipulations", fontweight="bold", pad=22)
    handles, labels = ax_c.get_legend_handles_labels()
    order = [1, 2, 3, 0]
    ax_c.legend([handles[i] for i in order], [labels[i] for i in order],
                loc="upper right", bbox_to_anchor=(1.0, 1.02), frameon=False,
                handlelength=1.9, labelspacing=0.3, borderpad=0.2)

    # ---------------- D: thalamic ----------------
    _, kleak_tc = schedule(config.gkl_TC_awake, config.gkl_TC_s2,
                           config.gkl_TC_s3, config.gkl_TC_rem)
    _, kleak_re = schedule(config.gkl_RE_awake, config.gkl_RE_s2,
                           config.gkl_RE_s3, config.gkl_RE_rem)
    _, gaba_thal = schedule(config.awake_GABA_thal, config.s2_GABA_thal,
                            config.s3_GABA_thal, config.rem_GABA_thal)
    _, ampa_thal = schedule(config.awake_AMPA_thal, config.s2_AMPA_thal,
                            config.s3_AMPA_thal, config.rem_AMPA_thal)

    ymax_d = 285.0
    draw_state_bands(ax_d, ymax_d)
    ax_d.plot(t, ampa_thal, color=C_PORT, lw=1.5, ls=(0, (4, 2.4)), zorder=2,
              label="AMPA (thalamic, flat in this port)")
    ax_d.plot(t, kleak_tc, color=C_KLEAK, lw=2.1, zorder=3, label="K-leak-TC")
    ax_d.plot(t, kleak_re, color=C_AMPA, lw=2.1, zorder=3, label="K-leak-RE")
    ax_d.plot(t, gaba_thal, color=C_GABA, lw=2.1, zorder=3, label="GABA")
    style_schedule_axis(ax_d, ymax_d)
    ax_d.set_title("Thalamic Manipulations", fontweight="bold", pad=22)
    handles, labels = ax_d.get_legend_handles_labels()
    order = [1, 2, 3, 0]
    ax_d.legend([handles[i] for i in order], [labels[i] for i in order],
                loc="upper right", bbox_to_anchor=(1.0, 1.02), frameon=False,
                handlelength=1.9, labelspacing=0.3, borderpad=0.2)

    fig.text(
        0.5, 0.010,
        "Krishnan et al. 2016 eLife 5:e18607, Figure 1, redrawn from this NEURON port: "
        "per-state values from config.py, schedule shape from the h.Vector.play calls in bazh_net.py. "
        "No simulation is involved.\n"
        "Dashed traces are port-specific and do not appear in the published panel. "
        "PY$\\rightarrow$PY AMPA follows its own multiplier (100 / 124 / 270 / 106 %) rather than the shared\n"
        "cortical AMPA multiplier (100 / 124 / 196 / 88 %), and the thalamic AMPA multiplier is held at 100 % in every state.",
        ha="center", va="bottom", fontsize=7.6, color="#444444", linespacing=1.5,
    )

    os.makedirs(os.path.dirname(OUTFILE), exist_ok=True)
    fig.savefig(OUTFILE, dpi=180, facecolor="white")
    plt.close(fig)

    # Print the table the panels encode, so the numbers are checkable by eye.
    print("wrote %s" % OUTFILE)
    print("\nnormalized conductance (% of awake)   Awake    N2     N3    REM")
    rows = [
        ("cortex K-leak (PY, IN)", config.gkl_pyr_awake, config.gkl_pyr_s2,
         config.gkl_pyr_s3, config.gkl_pyr_rem),
        ("cortex AMPA (TC->PY/IN, PY->IN)", config.awake_AMPA_cort,
         config.s2_AMPA_cort, config.s3_AMPA_cort, config.rem_AMPA_cort),
        ("cortex AMPA (PY->PY)", config.awake_AMPA_pyrpyr,
         config.s2_AMPA_pyrpyr, config.s3_AMPA_pyrpyr, config.rem_AMPA_pyrpyr),
        ("cortex GABA (IN->PY)", config.awake_GABA_D2, config.s2_GABA_D2,
         config.s3_GABA_D2, config.rem_GABA_D2),
        ("thalamus K-leak-TC", config.gkl_TC_awake, config.gkl_TC_s2,
         config.gkl_TC_s3, config.gkl_TC_rem),
        ("thalamus K-leak-RE", config.gkl_RE_awake, config.gkl_RE_s2,
         config.gkl_RE_s3, config.gkl_RE_rem),
        ("thalamus GABA (RE->TC A/B, RE->RE)", config.awake_GABA_thal,
         config.s2_GABA_thal, config.s3_GABA_thal, config.rem_GABA_thal),
        ("thalamus AMPA (TC->RE, PY->TC/RE)", config.awake_AMPA_thal,
         config.s2_AMPA_thal, config.s3_AMPA_thal, config.rem_AMPA_thal),
    ]
    for name, a, s2, s3, rem in rows:
        print("%-36s %6.1f %6.1f %6.1f %6.1f"
              % (name, 100.0, 100.0 * s2 / a, 100.0 * s3 / a, 100.0 * rem / a))
    print("\nIh half-activation V (mV):  " + "  ".join(
        "%s %.1f" % (n, -(75.0 + f)) for n, f in shifts))
    print("Shift_HA (awake - REM) = %.1f mV" % (v_half_awake - v_half_rem))


if __name__ == "__main__":
    main()
