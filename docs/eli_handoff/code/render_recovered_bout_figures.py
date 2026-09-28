"""Render the historical bout-screen data without remeasuring any LFP.

This renderer was reconstructed for Eli's handoff. The original plotting-only
commands were not recovered. It uses the exact saved condition/run CSV values
behind the archived figures, and deliberately uses new output filenames.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CV = 'median_compatible_bout_period_cv_robust'
SHAPE = 'median_compatible_bout_template_corr_median'
N3 = 'template_corr_to_continuous_n3'
DURATION = 'median_compatible_bout_duration_s'


def number(row, field):
    try:
        return float(row[field])
    except (ValueError, KeyError):
        return math.nan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'data')
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).resolve().parents[1] / 'reconstructed_figures')
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader((args.data_dir / 'rda_condition_summary.csv').open()))
    runs = list(csv.DictReader((args.data_dir / 'rda_run_metrics.csv').open()))
    control = next(r for r in runs if r['condition_id'] == 'canonical_continuous_n3')
    purple, green, blue = '#785895', '#31866e', '#5884b1'

    def style(row):
        if row['condition_id'] == 'state_2':
            return purple, 'D'
        return (green if row['replicated'] == 'True' else blue,
                'X' if row['sharp_fast_warning'] == 'True' else 'o')

    def draw(ax, row, field):
        color, marker = style(row)
        ax.scatter(number(row, field), number(row, SHAPE),
                   s=45+2*number(row, 'rda_compatible_burden_pct'),
                   c=color, marker=marker, alpha=.8, zorder=4)

    legend = [Line2D([], [], color=c, marker=m, lw=0, markersize=9, label=l) for c,m,l in
              [(green,'o','Replicated across seeds'), (blue,'o','One-seed / not replicated'),
               (purple,'*','Continuous N3'), (purple,'D','State 2 N3'),
               ('#555555','X','Sharp/fast warning in at least one run')]]
    for field, filename, title, candidates in [
        (CV, 'reconstructed_rda_bout_screen_all_completed.png', 'RDA bout screen across completed simulations', None),
        (N3, 'reconstructed_rda_candidate_n3_tradeoff.png', 'Repeated cycle shape versus N3 resemblance',
         {'state_2', 'state_39', 'state_44', 'state_74', 'state_99_tcIm0p3_r5000ms'}),
    ]:
        fig, ax = plt.subplots(figsize=(11, 7), layout='constrained')
        plotted=[]
        for row in rows:
            if candidates is not None and row['condition_id'] not in candidates:
                continue
            if not all(math.isfinite(number(row, k)) for k in [field, SHAPE]):
                continue
            draw(ax, row, field)
            plotted.append(row)
            if candidates is not None or row['condition_id'] in {'state_2','state_39','state_44','state_74','state_99_tcIm0p3_r5000ms'}:
                label = 'State 2 N3' if row['condition_id']=='state_2' else row['state'] + (' + TC bias' if row['run_tag'] else '')
                if candidates is not None:
                    label += f' ({number(row, DURATION):.0f} s)'
                ax.annotate(label,(number(row,field),number(row,SHAPE)),xytext=(5,6),textcoords='offset points',fontsize=10)
        cx = number(control, CV) if field == CV else 1.0
        cy = number(control, SHAPE)
        ax.scatter(cx, cy, c=purple, marker='*', s=240, zorder=5)
        ax.axvline(cx,color=purple,alpha=.4,lw=1)
        ax.axhline(cy,color=purple,alpha=.4,lw=1)
        ax.annotate('Continuous N3',(cx,cy),xytext=(-90,8),textcoords='offset points',fontsize=10)
        ax.set(title=title, ylabel='Median cycle-template correlation')
        ax.set_xlabel('Robust cycle-period CV (lower = more regular)' if field==CV else 'Template correlation with continuous N3')
        if field == CV:
            ax.invert_xaxis()
        ax.grid(color='#dce1e8',lw=.7)
        ax.legend(handles=legend,frameon=False,fontsize=9,loc='lower left')
        fig.savefig(args.output_dir / filename,dpi=170)
        plt.close(fig)
        fields=['condition_id','state','seeds_available','n_seeds',CV,SHAPE,N3,DURATION,
                'rda_compatible_burden_pct','replicated','sharp_fast_warning']
        control_point = {**control, 'seeds_available': control['seed'], 'n_seeds': 1,
                         N3: 1.0, 'replicated': False,
                         'sharp_fast_warning': float(control.get('sharp_fast_compatible_bout_fraction') or 0) > 0}
        with (args.output_dir / filename.replace('.png','_plot_values.csv')).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows([*plotted, control_point])
    print(f'Wrote two reconstructed figures to {args.output_dir}')


if __name__ == '__main__':
    main()
