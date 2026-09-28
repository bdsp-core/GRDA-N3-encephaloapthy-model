"""Check the legacy README values against both Welch integration conventions.

The README's five oldest rows numerically match simple bin summation. The saved
analyzer and committed CSV instead use trapezoidal integration; this reproduces
the discrepancy from the original raw files without guessing different seeds.
"""
import argparse
import csv
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import statistics
from scipy import signal
import analyze_lfp_states as analysis
import raw_package_io


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    def run(task):
        state, seed = task
        row, eeg, fs = raw_package_io.analyze_run(args.root, state, seed, 8, 120., 200., (.5,30.), 10., 119.)
        f, power = signal.welch(eeg[2000:23800], fs=fs, nperseg=1024, noverlap=512, detrend='constant')
        use=(f>=.5)&(f<=30.);f=f[use];power=power[use]
        return {'state':state, 'seed':seed, 'delta_sum':float(100*power[(f>=.5)&(f<=4.)].sum()/power.sum()),
                'centroid_sum':float((f*power).sum()/power.sum()),
                'delta_trapezoid':row['delta_pct'], 'centroid_trapezoid':row['spectral_centroid_hz']}
    tasks=[(state,seed) for state in (2,35,38,39,40) for seed in (1,2,3)]
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows=list(pool.map(run,tasks))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    for state in (2,35,38,39,40):
        group=[r for r in rows if r['state']==state]
        for method in ('sum','trapezoid'):
            d=[r[f'delta_{method}'] for r in group];c=[r[f'centroid_{method}'] for r in group]
            print(f'{state} {method}: {statistics.mean(d):.1f} ± {statistics.stdev(d):.1f}%; '
                  f'{statistics.mean(c):.2f} ± {statistics.stdev(c):.2f} Hz')


if __name__=='__main__':
    main()
