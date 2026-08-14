#!/usr/bin/env python
"""Reproduce Figure 2 of Krishnan et al. 2016 eLife from this port's output.

Figure 2 of "Cellular and neurochemical basis of sleep stages in the
thalamocortical network" (eLife 2016;5:e18607) shows one Awake -> N2 -> N3 ->
REM -> N2 ultradian cycle in four stacked panels:

  A  activity of all 800 neurons, y = neuron index, x = time
  B  somatic voltage of one representative neuron per population
  C  the LFP
  D  a sliding-window FFT spectrogram of the LFP, 0-100 Hz

This script rebuilds those four panels from the text files written by
``bazh_net.py``.  One deliberate departure from the paper is called out on the
figure itself: panel A there is an *image of membrane voltage* for every cell,
which this port does not save (800 cells x 14.4M samples), so panel A here is
the spike raster the port does save, coloured by population.  Everything the
paper's panel A is used for -- when each population fires, how synchronous it
is, and where the Up/Down structure of N3 sits -- survives that substitution;
the subthreshold hyperpolarisation across stages does not, and panel B is where
to read that instead.

Preprocessing follows the conventions already declared in
``analyze_lfp_states.py``: median centering, polyphase anti-alias resampling,
and zero-phase filtering.  The LFP is resampled to 200 Hz, whose Nyquist is
exactly the 100 Hz top of the paper's spectrogram axis.

Usage
-----
Full 360 s sleep-cycle run (the paper figure)::

    .venv/bin/python plot_paper_fig02.py

Concatenated static-state runs, for validating the plotting path before the
sleep-cycle run finishes::

    .venv/bin/python plot_paper_fig02.py --validation

Neither mode runs a simulation; both only read files from ``--root``.
"""

from __future__ import annotations

import argparse
import ast
import re
import sys
import textwrap
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from scipy import signal

# --------------------------------------------------------------------------
# Populations.
#
# Display order is the paper's: PY, IN, TC, RE, top to bottom.  Note that TC is
# drawn above RE even though the gids run PY, IN, RE, TC -- the two thalamic
# blocks are swapped for display so the rows line up with Figure 2.
#
# Colours are categorical slots 1, 2, 3 and 7 of the project's validated
# palette.  On a white surface the four pass the lightness band, the chroma
# floor, all-pairs CVD separation (worst dE 9.2, deuteranopia) and the
# normal-vision floor (worst dE 16.3); aqua sits below 3:1 contrast, which the
# per-block direct labels relieve.  Identity is never carried by colour alone:
# every block and every trace row is labelled in text as well.
# --------------------------------------------------------------------------
POPULATIONS = [
    # name, first gid, count, colour
    ("PY", 0, 500, "#2a78d6"),
    ("IN", 500, 100, "#eb6834"),
    ("TC", 700, 100, "#1baf7a"),
    ("RE", 600, 100, "#4a3aa7"),
]

INK = "#0b0b0b"
INK_SOFT = "#52514e"
INK_MUTED = "#8a8880"
GRID = "#dedcd5"
TRANSITION_FILL = "#e7e5dd"
INIT_FILL = "#f0eee6"

RAW_DT_MS = 0.025  # h.dt in bazh_net.py; the LFP is written every time step


# --------------------------------------------------------------------------
# Configuration: stage transition times
# --------------------------------------------------------------------------
TRANSITION_NAMES = (
    "awake_to_s2_start",
    "awake_to_s2_end",
    "s2_to_s3_start",
    "s2_to_s3_end",
    "s3_to_rem_start",
    "s3_to_rem_end",
    "rem_to_s2_start",
    "rem_to_s2_end",
)


def read_transition_times(config_path: Path) -> dict[str, float]:
    """Pull the sleep-stage transition times out of ``config.py``.

    ``config.py`` only defines these names when ``do_sleepstates`` is true, and
    importing it pulls in NEURON, so the file is parsed rather than executed.
    The values are plain integer literals in milliseconds.
    """
    tree = ast.parse(config_path.read_text(encoding="utf-8"), str(config_path))
    found: dict[str, float] = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign) or not isinstance(node.value, ast.Constant):
            continue
        if not isinstance(node.value.value, (int, float)):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id in TRANSITION_NAMES:
                found[target.id] = float(node.value.value)
    missing = [name for name in TRANSITION_NAMES if name not in found]
    if missing:
        raise KeyError(f"{config_path} defines no {', '.join(missing)}")
    return found


def build_stage_plan(
    transitions_ms: dict[str, float], duration_s: float
) -> tuple[list[tuple[str, float, float]], list[tuple[float, float]]]:
    """Turn the transition times into steady stages and ramp windows, in seconds.

    Both lists are clipped to ``duration_s`` so a partially written run still
    plots, showing only the stages it actually reached.
    """
    t = {name: value / 1000.0 for name, value in transitions_ms.items()}
    steady = [
        ("Awake", 0.0, t["awake_to_s2_start"]),
        ("N2", t["awake_to_s2_end"], t["s2_to_s3_start"]),
        ("N3", t["s2_to_s3_end"], t["s3_to_rem_start"]),
        ("REM", t["s3_to_rem_end"], t["rem_to_s2_start"]),
        ("N2", t["rem_to_s2_end"], duration_s),
    ]
    ramps = [
        (t["awake_to_s2_start"], t["awake_to_s2_end"]),
        (t["s2_to_s3_start"], t["s2_to_s3_end"]),
        (t["s3_to_rem_start"], t["s3_to_rem_end"]),
        (t["rem_to_s2_start"], t["rem_to_s2_end"]),
    ]
    steady = [
        (name, start, min(stop, duration_s))
        for name, start, stop in steady
        if start < duration_s
    ]
    ramps = [
        (start, min(stop, duration_s)) for start, stop in ramps if start < duration_s
    ]
    return steady, ramps


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def load_scalar_column(path: Path) -> np.ndarray:
    """Read a one-column text file, tolerating a partly written final line."""
    try:
        values = np.loadtxt(path, dtype=np.float64)
    except ValueError:
        # A run still in progress can leave the last line half flushed.
        text = path.read_text(encoding="utf-8", errors="ignore")
        lines = text.splitlines()
        values = np.array(
            [float(line) for line in lines[:-1] if line.strip()], dtype=np.float64
        )
    values = np.atleast_1d(values)
    if not np.isfinite(values).all():
        raise ValueError(f"Non-finite samples in {path}")
    return values


def resample_lfp(
    raw: np.ndarray, raw_fs: float, view_fs: float
) -> np.ndarray:
    """Median-centre and polyphase-resample, as ``analyze_lfp_states`` does.

    ``resample_poly`` applies the anti-alias filter for us, so the 200 Hz view
    is band-limited to 100 Hz -- the paper's spectrogram range -- without a
    separate low-pass.  No band-pass is applied here: unlike the metrics script
    this figure has to show everything up to 100 Hz.
    """
    centred = raw - np.median(raw)
    ratio = Fraction(view_fs / raw_fs).limit_denominator(10000)
    return signal.resample_poly(
        centred, ratio.numerator, ratio.denominator, padtype="line"
    )


def highpass(x: np.ndarray, fs: float, cutoff: float) -> np.ndarray:
    """Zero-phase Butterworth high-pass, matching the project's filter style."""
    if cutoff <= 0:
        return x
    sos = signal.butter(4, cutoff, btype="highpass", fs=fs, output="sos")
    padlen = min(len(x) - 1, int(3 * fs / cutoff))
    return signal.sosfiltfilt(sos, x, padtype="even", padlen=padlen)


def load_raster(path: Path, duration_s: float) -> tuple[np.ndarray, np.ndarray]:
    """Read ``time_ms  gid`` pairs.  Parallel runs are unsorted; that is fine."""
    data = np.loadtxt(path, dtype=np.float64, ndmin=2)
    if data.size == 0:
        return np.empty(0), np.empty(0, dtype=np.int64)
    t_s = data[:, 0] / 1000.0
    gid = data[:, 1].astype(np.int64)
    keep = t_s <= duration_s
    return t_s[keep], gid[keep]


def load_traces(
    root: Path, runtag: str, nhost: int, duration_s: float
) -> tuple[dict[str, tuple[np.ndarray, np.ndarray]], float]:
    """Collect the representative somatic traces, keyed by population name.

    The files are globbed rather than named so a run that recorded a different
    set of ``FINK_TRACE_GIDS`` still plots.  The longest trace's end time is
    returned alongside: it is the only file in the set that carries explicit
    timestamps, which makes it the check on the sample-rate assumption used to
    put the LFP on a time axis.
    """
    tail = re.escape(f"_{runtag}_nhost={nhost:g}.txt")
    pattern = re.compile(r"^trace_(PY|IN|RE|TC)_gid(\d+)" + tail + "$")
    traces: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    end_s = 0.0
    for path in sorted(root.glob(f"trace_*_{runtag}_nhost={nhost:g}.txt")):
        match = pattern.match(path.name)
        if match is None:
            continue
        species = match.group(1)
        if species in traces:
            continue
        data = np.loadtxt(path, dtype=np.float64, comments="#", ndmin=2)
        if data.size == 0:
            continue
        t_s = data[:, 0] / 1000.0
        end_s = max(end_s, float(t_s[-1]))
        keep = t_s <= duration_s
        traces[species] = (t_s[keep], data[keep, 1])
    return traces, end_s


# --------------------------------------------------------------------------
# Segments
#
# A segment is one contiguous stretch of simulated time placed on the figure's
# display timeline.  The sleep-cycle run is a single segment; the validation
# mode is one segment per static-state run, laid end to end.  Everything
# downstream of here works on a list of segments and so is shared by the two
# modes.
# --------------------------------------------------------------------------
@dataclass
class Segment:
    label: str
    t0_s: float
    duration_s: float
    view_fs: float
    lfp: np.ndarray
    raster_t: np.ndarray
    raster_gid: np.ndarray
    traces: dict[str, tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
    source: str = ""

    @property
    def t1_s(self) -> float:
        return self.t0_s + self.duration_s

    @property
    def lfp_time(self) -> np.ndarray:
        return self.t0_s + np.arange(len(self.lfp)) / self.view_fs


def load_segment(
    root: Path,
    runtag: str,
    nhost: int,
    view_fs: float,
    t0_s: float,
    label: str,
    lfp_source: str,
    duration_s: float | None = None,
) -> Segment:
    prefix = "lfp" if lfp_source == "biophysical" else "vcort"
    lfp_path = root / f"{prefix}_{runtag}_nhost={nhost:g}.txt"
    if not lfp_path.exists():
        raise FileNotFoundError(lfp_path)

    raw = load_scalar_column(lfp_path)
    raw_fs = 1000.0 / RAW_DT_MS
    if duration_s is None:
        duration_s = len(raw) / raw_fs
    else:
        raw = raw[: int(round(duration_s * raw_fs))]
    lfp = resample_lfp(raw, raw_fs, view_fs)

    raster_path = root / f"raster_{runtag}_nhost={nhost:g}.txt"
    if raster_path.exists():
        raster_t, raster_gid = load_raster(raster_path, duration_s)
    else:
        raster_t, raster_gid = np.empty(0), np.empty(0, dtype=np.int64)

    traces, trace_end_s = load_traces(root, runtag, nhost, duration_s)
    # The LFP file is a bare column of numbers, so its time axis rests entirely
    # on RAW_DT_MS.  The trace files do carry timestamps; if the two disagree the
    # simulation used a different h.dt and every time on this figure is wrong.
    if trace_end_s > 0 and abs(trace_end_s - duration_s) > 0.01 * duration_s:
        print(
            f"warning: {lfp_path.name} implies {duration_s:.1f} s at "
            f"{raw_fs / 1000:g} kHz, but the trace files end at "
            f"{trace_end_s:.1f} s. Check h.dt against RAW_DT_MS.",
            file=sys.stderr,
        )

    return Segment(
        label=label,
        t0_s=t0_s,
        duration_s=duration_s,
        view_fs=view_fs,
        lfp=lfp,
        raster_t=raster_t + t0_s,
        raster_gid=raster_gid,
        traces={
            name: (t + t0_s, v) for name, (t, v) in traces.items()
        },
        source=lfp_path.name,
    )


# --------------------------------------------------------------------------
# Decimation
# --------------------------------------------------------------------------
def minmax_decimate(
    t: np.ndarray, y: np.ndarray, n_cols: int
) -> tuple[np.ndarray, np.ndarray]:
    """Reduce a trace to two points per pixel column, keeping the extremes.

    Plain striding would throw away spike peaks, which is most of what panel B
    is for.  Emitting the min and the max of each bin in their original time
    order preserves the envelope exactly at a fraction of the vertices.
    """
    n = len(y)
    if n <= 2 * n_cols or n_cols < 1:
        return t, y
    step = n // n_cols
    used = step * n_cols
    yb = y[:used].reshape(n_cols, step)
    tb = t[:used].reshape(n_cols, step)
    rows = np.arange(n_cols)
    i_min = yb.argmin(axis=1)
    i_max = yb.argmax(axis=1)
    first = np.minimum(i_min, i_max)
    second = np.maximum(i_min, i_max)
    t_out = np.empty(2 * n_cols)
    y_out = np.empty(2 * n_cols)
    t_out[0::2] = tb[rows, first]
    y_out[0::2] = yb[rows, first]
    t_out[1::2] = tb[rows, second]
    y_out[1::2] = yb[rows, second]
    if used < n:
        t_out = np.concatenate([t_out, t[used:]])
        y_out = np.concatenate([y_out, y[used:]])
    return t_out, y_out


def subsample_spikes(
    t: np.ndarray, gid: np.ndarray, cap: int, seed: int = 0
) -> tuple[np.ndarray, np.ndarray, float]:
    """Thin a spike train to at most ``cap`` markers, without bias in time.

    A 360 s run produces far more spikes than a raster has pixels, so drawing
    them all costs minutes and changes nothing visible.  A uniform random
    subsample keeps the local density proportional to the true rate, which is
    what the eye reads off a raster.
    """
    n = len(t)
    if n <= cap:
        return t, gid, 1.0
    rng = np.random.default_rng(seed)
    idx = rng.choice(n, cap, replace=False)
    return t[idx], gid[idx], cap / n


# --------------------------------------------------------------------------
# Panel A geometry
# --------------------------------------------------------------------------
def population_layout(gap: float = 24.0) -> dict[str, dict[str, float]]:
    """Map gids onto display rows in the paper's PY / IN / TC / RE order."""
    layout: dict[str, dict[str, float]] = {}
    cursor = 0.0
    for name, gid_lo, count, colour in POPULATIONS:
        layout[name] = {
            "gid_lo": gid_lo,
            "count": count,
            "colour": colour,
            "y_lo": cursor,
            "y_hi": cursor + count,
            "centre": cursor + count / 2.0,
        }
        cursor += count + gap
    layout["_total"] = {"y_hi": cursor - gap}
    return layout


# --------------------------------------------------------------------------
# Drawing
# --------------------------------------------------------------------------
def draw_stage_furniture(
    ax,
    ramps: list[tuple[float, float]],
    boundaries: list[float],
    init_windows: list[tuple[float, float]],
    xlim: tuple[float, float],
    label_row: bool = False,
    steady: list[tuple[str, float, float]] | None = None,
    shade: bool = True,
) -> None:
    """Shade transition ramps and the initialisation epochs on one axes.

    There is one initialisation window per segment, not one per figure: each
    concatenated static-state run starts from its own highly synchronous
    initial condition and carries its own transient.

    ``shade=False`` is for the spectrogram, whose pcolormesh would sit on top of
    the fills anyway; there the boundaries are carried by the dashed lines alone
    so no data is obscured.
    """
    if shade:
        for start, stop in ramps:
            if stop > start:
                ax.axvspan(start, stop, color=TRANSITION_FILL, lw=0, zorder=0)
        for start, stop in init_windows:
            if stop > start:
                ax.axvspan(start, stop, color=INIT_FILL, lw=0, zorder=0)
    for x in boundaries:
        if xlim[0] < x < xlim[1]:
            ax.axvline(x, color=INK_MUTED, lw=0.9, ls=(0, (4, 3)), zorder=6)
    for _, stop in init_windows:
        if xlim[0] < stop < xlim[1]:
            ax.axvline(stop, color=INK_MUTED, lw=0.8, ls=(0, (1, 2)), zorder=6)
    if label_row and steady:
        for name, start, stop in steady:
            if stop <= start:
                continue
            ax.text(
                (start + stop) / 2.0,
                1.13,
                name,
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=12.5,
                fontweight="bold",
                color=INK,
            )


def style_axes(ax, xlim: tuple[float, float]) -> None:
    ax.set_xlim(xlim)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=INK_SOFT, labelsize=8.5, length=3, width=0.8)


def panel_label(fig, ax, letter: str) -> None:
    """Place the panel letter in the left margin, in figure coordinates.

    Axes coordinates would put it a fixed fraction of each panel's width to the
    left, which lands it on top of the y tick labels of the short panels.
    """
    fig.text(
        0.012,
        ax.get_position().y1,
        letter,
        ha="left",
        va="top",
        fontsize=16,
        fontweight="bold",
        color=INK,
    )


def draw_raster(ax, segments: list[Segment], cap: int, marker_size: float) -> dict:
    layout = population_layout()
    stats: dict[str, dict] = {}
    for name, gid_lo, count, colour in POPULATIONS:
        block = layout[name]
        ax.axhspan(
            block["y_lo"], block["y_hi"], color=colour, alpha=0.055, lw=0, zorder=0
        )
        times = []
        rows = []
        for seg in segments:
            mask = (seg.raster_gid >= gid_lo) & (seg.raster_gid < gid_lo + count)
            times.append(seg.raster_t[mask])
            rows.append(
                block["y_lo"] + (seg.raster_gid[mask] - gid_lo).astype(np.float64)
            )
        t_all = np.concatenate(times) if times else np.empty(0)
        y_all = np.concatenate(rows) if rows else np.empty(0)
        n_total = len(t_all)
        t_plot, y_plot, kept = subsample_spikes(t_all, y_all, cap)
        stats[name] = {"n": n_total, "kept": kept}
        if n_total:
            ax.plot(
                t_plot,
                y_plot,
                linestyle="none",
                marker=".",
                markersize=marker_size,
                markeredgewidth=0,
                color=colour,
                alpha=0.9,
                rasterized=True,
                zorder=3,
            )
    total_y = layout["_total"]["y_hi"]
    ax.set_ylim(total_y + 12, -12)
    ax.set_yticks([layout[name]["centre"] for name, _, _, _ in POPULATIONS])
    ax.set_yticklabels([name for name, _, _, _ in POPULATIONS])
    for tick, (name, _, _, colour) in zip(ax.get_yticklabels(), POPULATIONS):
        tick.set_color(colour)
        tick.set_fontweight("bold")
        tick.set_fontsize(11)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    return stats


def draw_trace_row(
    ax,
    segments: list[Segment],
    name: str,
    colour: str,
    n_cols: int,
    init_skip_s: float,
):
    drawn = False
    lo, hi = np.inf, -np.inf
    for seg in segments:
        if name not in seg.traces:
            continue
        t, v = seg.traces[name]
        if len(t) == 0:
            continue
        td, vd = minmax_decimate(t, v, n_cols)
        ax.plot(td, vd, color=INK, lw=0.45, solid_joinstyle="round", rasterized=True)
        # The network starts highly synchronous, so the first seconds carry a
        # transient far larger than anything physiological.  It is drawn (and
        # shaded) but must not be allowed to set the scale for the rest.
        settled = t >= seg.t0_s + init_skip_s
        if settled.any():
            lo = min(lo, float(np.min(v[settled])))
            hi = max(hi, float(np.max(v[settled])))
        drawn = True
    if drawn:
        if np.isfinite(lo):
            span = max(hi - lo, 1.0)
            ax.set_ylim(lo - 0.06 * span, hi + 0.06 * span)
            # A silent hyperpolarised cell can span only a millivolt or two, so
            # integer ticks would collapse to a single value.
            digits = 0 if hi - lo >= 10.0 else 1
            ax.set_yticks([lo, hi])
            ax.set_yticklabels([f"{lo:.{digits}f}", f"{hi:.{digits}f}"])
        else:
            ax.set_yticks([])
    else:
        ax.text(
            0.5,
            0.5,
            f"no {name} trace file",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=8.5,
            color=INK_MUTED,
        )
        ax.set_yticks([])
    ax.set_ylabel(name, color=colour, fontweight="bold", fontsize=11.5, rotation=0,
                  ha="right", va="center", labelpad=16)
    return drawn


def compute_spectrogram(
    segments: list[Segment], window_s: float, overlap: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sliding-window FFT, computed per segment and stitched afterwards.

    Concatenating segments before the FFT would put a step discontinuity at
    every join and paint a spurious broadband stripe there.  Transforming each
    segment separately and then placing its columns on the display timeline
    avoids that entirely, and is identical to the single-segment case when
    there is only one segment.
    """
    freqs = None
    times: list[np.ndarray] = []
    powers: list[np.ndarray] = []
    for seg in segments:
        fs = seg.view_fs
        nperseg = int(round(window_s * fs))
        # The window length has to be identical across segments: it sets the
        # frequency grid, and columns computed on different grids cannot be
        # placed on one axes.  A segment too short to hold a whole window is
        # dropped from the spectrogram rather than transformed at a different
        # resolution.
        if len(seg.lfp) < nperseg:
            print(
                f"  note: segment at t={seg.t0_s:.0f} s is shorter than the "
                f"{window_s:g} s window; omitted from the spectrogram",
                file=sys.stderr,
            )
            continue
        noverlap = int(round(nperseg * overlap))
        f, t, sxx = signal.spectrogram(
            seg.lfp,
            fs=fs,
            window="hann",
            nperseg=nperseg,
            noverlap=min(noverlap, nperseg - 1),
            detrend="constant",
            scaling="density",
            mode="psd",
        )
        freqs = f
        times.append(t + seg.t0_s)
        powers.append(sxx)
    if freqs is None:
        raise ValueError("No segment was long enough to transform")
    return freqs, np.concatenate(times), np.hstack(powers)


def build_figure(
    segments: list[Segment],
    steady: list[tuple[str, float, float]],
    ramps: list[tuple[float, float]],
    boundaries: list[float],
    args,
    title: str,
    subtitle: str,
    note: str,
) -> tuple[plt.Figure, dict]:
    xlim = (segments[0].t0_s, segments[-1].t1_s)
    n_cols = max(800, int(args.width * args.dpi))

    # The header is laid out first so the panels can be pushed down by exactly
    # as much as the wrapped text needs.  Otherwise a subtitle that happens to
    # wrap onto a second line lands on top of the stage labels.
    title_y = 0.985
    subtitle_y = 0.951
    subtitle_lines = textwrap.wrap(subtitle, 150)
    note_y = subtitle_y - 0.021 * len(subtitle_lines)
    note_lines = textwrap.wrap(note, 150)
    header_bottom = note_y - 0.019 * len(note_lines)
    top = min(0.86, max(0.70, header_bottom - 0.060))

    # One initialisation window per segment: every run starts synchronous.
    init_windows = [
        (seg.t0_s, min(seg.t0_s + args.init_skip_s, seg.t1_s)) for seg in segments
    ]

    fig = plt.figure(figsize=(args.width, args.height), facecolor="white")
    outer = GridSpec(
        4,
        2,
        figure=fig,
        height_ratios=[4.6, 2.5, 1.15, 2.15],
        width_ratios=[1.0, 0.016],
        left=0.098,
        right=0.935,
        top=top,
        bottom=0.062,
        hspace=0.32,
        wspace=0.012,
    )

    ax_a = fig.add_subplot(outer[0, 0])
    inner_b = GridSpecFromSubplotSpec(
        4, 1, subplot_spec=outer[1, 0], hspace=0.16
    )
    ax_b = [fig.add_subplot(inner_b[i, 0], sharex=ax_a) for i in range(4)]
    ax_c = fig.add_subplot(outer[2, 0], sharex=ax_a)
    ax_d = fig.add_subplot(outer[3, 0], sharex=ax_a)
    ax_cbar = fig.add_subplot(outer[3, 1])

    # ---- A: spike raster -------------------------------------------------
    raster_stats = draw_raster(ax_a, segments, args.spike_cap, args.marker_size)
    style_axes(ax_a, xlim)
    draw_stage_furniture(
        ax_a, ramps, boundaries, init_windows, xlim, label_row=True, steady=steady
    )
    panel_label(fig, ax_a, "A")
    ax_a.set_title(
        "Network activity of all 800 neurons: spike raster "
        "(the paper shows membrane voltage as an image)",
        fontsize=10.5,
        color=INK_SOFT,
        pad=7,
        loc="left",
    )
    plt.setp(ax_a.get_xticklabels(), visible=False)

    # ---- B: representative somatic voltages ------------------------------
    drawn_any = False
    for ax, (name, _, _, colour) in zip(ax_b, POPULATIONS):
        drawn = draw_trace_row(ax, segments, name, colour, n_cols, args.init_skip_s)
        drawn_any = drawn_any or drawn
        style_axes(ax, xlim)
        draw_stage_furniture(ax, ramps, boundaries, init_windows, xlim)
        plt.setp(ax.get_xticklabels(), visible=False)
    panel_label(fig, ax_b[0], "B")
    ax_b[0].set_title(
        "Membrane voltage of one representative neuron per population (mV)",
        fontsize=10.5,
        color=INK_SOFT,
        pad=7,
        loc="left",
    )

    # ---- C: LFP ----------------------------------------------------------
    lfp_lo, lfp_hi = np.inf, -np.inf
    for seg in segments:
        x = highpass(seg.lfp, seg.view_fs, args.lfp_highpass_hz)
        t = seg.lfp_time
        keep = t >= seg.t0_s + args.init_skip_s
        if keep.any():
            lfp_lo = min(lfp_lo, float(np.min(x[keep])))
            lfp_hi = max(lfp_hi, float(np.max(x[keep])))
        td, xd = minmax_decimate(t, x, n_cols)
        ax_c.plot(td, xd, color=INK, lw=0.45, rasterized=True)
    if np.isfinite(lfp_lo):
        span = max(lfp_hi - lfp_lo, 1e-6)
        ax_c.set_ylim(lfp_lo - 0.10 * span, lfp_hi + 0.10 * span)
    style_axes(ax_c, xlim)
    draw_stage_furniture(ax_c, ramps, boundaries, init_windows, xlim)
    ax_c.set_ylabel("LFP (mV)", color=INK_SOFT, fontsize=9.5)
    panel_label(fig, ax_c, "C")
    ax_c.set_title(
        f"LFP ({args.lfp_source}), median-centred, resampled to {args.view_fs:g} Hz, "
        f"{args.lfp_highpass_hz:g} Hz high-pass",
        fontsize=10.5,
        color=INK_SOFT,
        pad=7,
        loc="left",
    )
    plt.setp(ax_c.get_xticklabels(), visible=False)

    # ---- D: spectrogram --------------------------------------------------
    freqs, times, power = compute_spectrogram(
        segments, args.spec_window_s, args.spec_overlap
    )
    floor = np.max(power) * 1e-12 if np.max(power) > 0 else 1e-30
    db = 10.0 * np.log10(np.maximum(power, floor))
    settled = np.ones(len(times), dtype=bool)
    for seg in segments:
        settled &= ~(
            (times >= seg.t0_s) & (times < seg.t0_s + args.init_skip_s)
        )
    reference = db[:, settled] if settled.any() else db
    vmin, vmax = np.percentile(reference, [2.0, 99.6])
    mesh = ax_d.pcolormesh(
        times,
        freqs,
        db,
        cmap=args.cmap,
        vmin=vmin,
        vmax=vmax,
        shading="nearest",
        rasterized=True,
    )
    ax_d.set_ylim(args.fmax, 0)  # 0 Hz at the top, as in the paper
    ax_d.set_yticks([0, 25, 50, 75, 100][: int(args.fmax // 25) + 1])
    style_axes(ax_d, xlim)
    draw_stage_furniture(ax_d, ramps, boundaries, init_windows, xlim, shade=False)
    ax_d.set_ylabel("Frequency (Hz)", color=INK_SOFT, fontsize=9.5)
    ax_d.set_xlabel("Time (s)", color=INK_SOFT, fontsize=10)
    panel_label(fig, ax_d, "D")
    ax_d.set_title(
        f"Spectrogram of the LFP ({args.spec_window_s:g} s Hann window, "
        f"{args.spec_overlap * 100:g}% overlap)",
        fontsize=10.5,
        color=INK_SOFT,
        pad=7,
        loc="left",
    )
    cbar = fig.colorbar(mesh, cax=ax_cbar)
    cbar.set_label("Power (dB re. arb. units)", color=INK_SOFT, fontsize=8.5)
    cbar.ax.tick_params(colors=INK_SOFT, labelsize=8)
    cbar.outline.set_visible(False)

    # ---- header ----------------------------------------------------------
    fig.text(0.012, title_y, title, fontsize=16.5, fontweight="bold", color=INK,
             va="top")
    fig.text(
        0.012,
        subtitle_y,
        "\n".join(subtitle_lines),
        fontsize=10,
        color=INK_SOFT,
        va="top",
        linespacing=1.5,
    )
    fig.text(
        0.012,
        note_y,
        "\n".join(note_lines),
        fontsize=9,
        color=INK_MUTED,
        va="top",
        linespacing=1.55,
    )
    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=6, color=colour,
               label=name)
        for name, _, _, colour in POPULATIONS
    ]
    if ramps:
        handles.append(
            plt.Rectangle((0, 0), 1, 1, color=TRANSITION_FILL, label="stage transition")
        )
    handles.append(
        plt.Rectangle((0, 0), 1, 1, color=INIT_FILL, label="initialisation")
    )
    fig.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.935, 0.995),
        ncol=3,
        frameon=False,
        fontsize=9.5,
        handletextpad=0.5,
        columnspacing=1.6,
        labelcolor=INK_SOFT,
    )
    return fig, {"raster": raster_stats, "traces_drawn": drawn_any}


# --------------------------------------------------------------------------
# Modes
# --------------------------------------------------------------------------
STATIC_STATES = [(0, "Awake"), (1, "N2"), (2, "N3"), (42, "State 42")]


def sweep_mode(args) -> tuple[list[Segment], list, list, list, str, str, str]:
    runtag = f"sweep_seed{args.seed}"
    segment = load_segment(
        args.root, runtag, args.nhost, args.view_fs, 0.0, "sleep cycle",
        args.lfp_source,
    )
    transitions = read_transition_times(args.root / "config.py")
    steady, ramps = build_stage_plan(transitions, segment.duration_s)
    boundaries = [start for start, _ in ramps] + [stop for _, stop in ramps]
    ramp_s = transitions["awake_to_s2_end"] - transitions["awake_to_s2_start"]
    # The run is only the whole figure once it has passed the last ramp.
    complete = segment.duration_s >= transitions["rem_to_s2_end"] / 1000.0
    reached = steady[-1][0] if steady else "nothing"
    title = "Krishnan et al. 2016 (eLife 5:e18607) Figure 2, reproduced"
    subtitle = (
        (
            "Full Awake / N2 / N3 / REM / N2 ultradian cycle"
            if complete
            else f"Sleep-cycle run still in progress, reached {reached}"
        )
        + f", {segment.duration_s:.0f} s  |  500 PY, 100 IN, 100 RE, 100 TC  |  "
        f"seed {args.seed}, {args.nhost} MPI ranks"
    )
    note = (
        "Panel A is a spike raster coloured by population, not the paper's "
        "per-cell voltage image: this port does not save voltage for all 800 "
        "cells, so subthreshold hyperpolarisation is read from panel B instead. "
        f"Shaded bands are the {ramp_s / 1000:g} s neuromodulator ramps between "
        f"stages and the first {args.init_skip_s:g} s of initialisation. "
        f"Source: {segment.source}"
    )
    return [segment], steady, ramps, boundaries, title, subtitle, note


def validation_mode(args) -> tuple[list[Segment], list, list, list, str, str, str]:
    segments: list[Segment] = []
    steady: list[tuple[str, float, float]] = []
    boundaries: list[float] = []
    cursor = 0.0
    for state, label in STATIC_STATES:
        runtag = f"state{state}_seed{args.seed}"
        try:
            seg = load_segment(
                args.root, runtag, args.nhost, args.view_fs, cursor, label,
                args.lfp_source,
            )
        except FileNotFoundError:
            print(f"  skipping state {state}: no LFP file yet")
            continue
        # A run that has only just started is all initialisation transient, and
        # contributes a sliver too short to carry a spectrogram window.
        if seg.duration_s < args.min_segment_s:
            print(
                f"  skipping state {state}: only {seg.duration_s:.0f} s written "
                f"so far (need {args.min_segment_s:g} s)"
            )
            continue
        segments.append(seg)
        steady.append((f"{label}  (state {state})", cursor, cursor + seg.duration_s))
        cursor += seg.duration_s
        boundaries.append(cursor)
        print(
            f"  state {state} ({label}): {seg.duration_s:.0f} s, "
            f"{len(seg.raster_t)} spikes, traces {sorted(seg.traces)}"
        )
    if not segments:
        raise SystemExit(
            "No static-state runs found. Expected files like "
            f"lfp_state0_seed{args.seed}_nhost={args.nhost}.txt in {args.root}"
        )
    boundaries = boundaries[:-1]
    title = "Figure 2 plotting path, validated on static-state runs"
    subtitle = (
        f"{len(segments)} constant-neuromodulator runs concatenated end to end "
        f"({cursor:.0f} s total)  |  500 PY, 100 IN, 100 RE, 100 TC  |  "
        f"seed {args.seed}, {args.nhost} MPI ranks"
    )
    note = (
        "This is NOT the paper figure: each block is a separate simulation held "
        "at one state of vigilance, so there are no transitions and no shared "
        "network history across the joins. It exercises every panel of "
        "plot_paper_fig02.py on real output while the 360 s sleep-cycle run "
        "completes. Panel A is a spike raster, not the paper's voltage image."
    )
    return segments, steady, [], boundaries, title, subtitle, note


# --------------------------------------------------------------------------
def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--validation",
        action="store_true",
        help="build the concatenated static-state figure instead of the sweep",
    )
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--nhost", type=int, default=8)
    parser.add_argument(
        "--lfp-source",
        choices=("biophysical", "vcort"),
        default="biophysical",
        help="biophysical LFP (this port's addition) or summed cortical voltage "
        "(what the paper plots)",
    )
    parser.add_argument("--view-fs", type=float, default=200.0,
                        help="LFP resampling rate; Nyquist sets the spectrogram top")
    parser.add_argument("--fmax", type=float, default=100.0)
    parser.add_argument("--spec-window-s", type=float, default=2.0)
    parser.add_argument("--spec-overlap", type=float, default=0.875)
    parser.add_argument("--lfp-highpass-hz", type=float, default=0.3)
    parser.add_argument("--init-skip-s", type=float, default=10.0,
                        help="initialisation epoch, shaded and excluded from scaling")
    parser.add_argument("--min-segment-s", type=float, default=20.0,
                        help="validation mode: ignore static-state runs shorter "
                        "than this, which are still all initialisation transient")
    parser.add_argument("--spike-cap", type=int, default=250_000,
                        help="max spike markers drawn per population")
    parser.add_argument("--marker-size", type=float, default=0.8)
    parser.add_argument("--cmap", default="turbo")
    parser.add_argument("--width", type=float, default=14.5)
    parser.add_argument("--height", type=float, default=12.5)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.fmax > args.view_fs / 2:
        print(
            f"fmax {args.fmax:g} Hz exceeds the {args.view_fs / 2:g} Hz Nyquist of "
            f"the resampled LFP; clamping.",
            file=sys.stderr,
        )
        args.fmax = args.view_fs / 2

    if args.output is None:
        args.output = args.root / "docs" / "figures" / (
            "paper_fig02_validation_static_states.png"
            if args.validation
            else "paper_fig02_sleep_cycle.png"
        )

    print(f"mode: {'validation (static states)' if args.validation else 'sleep cycle'}")
    builder = validation_mode if args.validation else sweep_mode
    try:
        segments, steady, ramps, boundaries, title, subtitle, note = builder(args)
    except FileNotFoundError as missing:
        print(
            f"error: {missing} does not exist.\n"
            "The 360 s sleep-cycle run writes it as FINK_SLEEPSTATES=1; until it "
            "finishes, build the static-state validation figure instead:\n"
            "    .venv/bin/python plot_paper_fig02.py --validation",
            file=sys.stderr,
        )
        return 1

    fig, stats = build_figure(
        segments, steady, ramps, boundaries, args, title, subtitle, note
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=args.dpi, facecolor="white")
    plt.close(fig)

    for name, info in stats["raster"].items():
        shown = "all" if info["kept"] >= 1.0 else f"{info['kept'] * 100:.0f}% shown"
        print(f"  {name}: {info['n']} spikes ({shown})")
    if not stats["traces_drawn"]:
        print("  warning: no trace files matched; panel B is empty")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
