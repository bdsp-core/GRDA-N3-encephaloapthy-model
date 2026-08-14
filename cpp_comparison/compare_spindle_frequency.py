"""Compare the N2 oscillation frequency of this port against the original C++.

Both implementations are measured on the same quantity: the mean membrane
potential of the cortical population, which is the LFP definition used in
Krishnan et al. The C++ writes per-cell voltage, so the mean is taken across
cells; this port writes the summed cortical voltage, which is the same signal
up to a scale factor that cannot affect frequency.

The comparison matters because the port's N2 oscillation sits near 6.9 Hz,
below the 7-15 Hz the paper states for the model and well below the 11-16 Hz
of human sleep spindles. Whether that is a porting artifact or a property of
the published model determines whether it can be fixed at all.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy import signal

BANDS = [
    ("delta", 0.5, 4.0),
    ("theta", 4.0, 8.0),
    ("alpha", 8.0, 11.0),
    ("sigma (human spindle)", 11.0, 16.0),
    ("beta", 16.0, 30.0),
]


def load_cpp(path: Path) -> tuple[np.ndarray, float]:
    """Mean cortical voltage from the C++ time_cx file (col 0 = time in ms)."""
    data = np.loadtxt(path)
    t_ms = data[:, 0]
    fs = 1000.0 / float(np.median(np.diff(t_ms)))
    return data[:, 1:].mean(axis=1), fs


def load_port(path: Path, duration_s: float) -> tuple[np.ndarray, float]:
    """Summed cortical voltage from this port, decimated toward ~1 kHz."""
    raw = np.loadtxt(path)
    fs = len(raw) / duration_s
    factor = max(1, int(round(fs / 1000.0)))
    return signal.resample_poly(raw, 1, factor), fs / factor


def describe(x: np.ndarray, fs: float, drop_s: float, label: str) -> dict[str, float]:
    """Peak frequency and band composition on the steady-state epoch."""
    x = x[int(drop_s * fs) :]
    x = x - np.mean(x)
    nperseg = int(min(len(x), 8.0 * fs))
    freq, power = signal.welch(x, fs=fs, nperseg=nperseg, noverlap=nperseg // 2)

    total_mask = (freq >= 0.5) & (freq <= 40.0)
    total = np.trapezoid(power[total_mask], freq[total_mask])
    row = {"peak_hz": float(freq[total_mask][np.argmax(power[total_mask])])}
    for name, low, high in BANDS:
        mask = (freq >= low) & (freq <= high)
        row[name] = 100.0 * float(np.trapezoid(power[mask], freq[mask])) / total

    print(f"\n{label}")
    print(f"  epoch {len(x)/fs:.1f} s at {fs:.0f} Hz, first {drop_s:.0f} s dropped")
    print(f"  dominant peak: {row['peak_hz']:.2f} Hz")
    for name, _, _ in BANDS:
        print(f"    {name:<22s} {row[name]:5.1f} %")
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpp-time-cx", type=Path, required=True)
    parser.add_argument("--port-vcort", type=Path, required=True)
    parser.add_argument("--port-duration-s", type=float, default=120.0)
    parser.add_argument("--drop-s", type=float, default=10.0)
    args = parser.parse_args()

    cpp, cpp_fs = load_cpp(args.cpp_time_cx)
    port, port_fs = load_port(args.port_vcort, args.port_duration_s)

    cpp_row = describe(cpp, cpp_fs, args.drop_s, "Original C++ (mean cortical Vm)")
    port_row = describe(port, port_fs, args.drop_s, "NEURON port (summed cortical Vm)")

    delta = abs(cpp_row["peak_hz"] - port_row["peak_hz"])
    print(f"\nPeak frequency difference: {delta:.2f} Hz")
    print(
        "Both implementations agree, so the low N2 frequency is a property of the\n"
        "published model rather than of this port."
        if delta < 1.0
        else "The implementations disagree; the port deviates from the original."
    )
    print(
        "\nFor reference, the paper states 7-15 Hz for the model (citing cat and\n"
        "ferret work), 8-15 Hz for its cat and mouse recordings, and 9-17 Hz for\n"
        "its human ECoG, noting that the human band is 'conventional in human\n"
        "studies'. Human sleep spindles are 11-16 Hz."
    )


if __name__ == "__main__":
    main()
