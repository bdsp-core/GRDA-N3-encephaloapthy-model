"""Recompute saved-LFP metrics; this does not run any simulations.

Uses the accompanying unchanged analyze_lfp_states.py. All results use 10–119 s
of 120 s runs; states 99/103/106 use the named −0.30 uA/cm² intervention.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json

import numpy as np
import scipy
import analyze_lfp_states as analysis
import raw_package_io


STATES = (2, 35, 38, 39, 40, 42, 44, 49, 50, 57, 74, 99, 103, 106)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()

    tasks = []
    for state in STATES:
        seeds = range(1, 6) if state in (39, 44) else range(1, 5) if state == 2 else range(1, 4)
        for seed in seeds:
            tag = 'tcIm0p3_r5000ms' if state in (99, 103, 106) else ''
            tasks.append((state, seed, tag))

    def run(task):
        state, seed, tag = task
        row, _, _ = raw_package_io.analyze_run(args.root, state, seed, 8, 120.0,
                                         200.0, (0.5, 30.0), 10.0, 119.0, tag)
        row['nhost'] = 8
        row['tc_bias_ua_cm2'] = -0.30 if tag else 0.0
        row['tc_ramp_ms'] = 5000 if tag else 0
        return row

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(run, tasks))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    provenance = {
        'scope': 'Recomputed from existing raw LFP; no simulation reruns',
        'analysis_script_sha256': hashlib.sha256(Path(analysis.__file__).read_bytes()).hexdigest(),
        'numpy': np.__version__, 'scipy': scipy.__version__,
        'runs': len(rows), 'epoch_s': [10.0, 119.0],
        'duration_s': 120.0, 'band_hz': [0.5, 30.0], 'view_fs_hz': 200,
        'welch_segment_s': 5.12,
        'source_sha256_scope': 'Original uncompressed text bytes; .gz is decompressed as a stream',
    }
    args.output.with_suffix('.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(f'Wrote {len(rows)} runs to {args.output}')


if __name__ == '__main__':
    main()
