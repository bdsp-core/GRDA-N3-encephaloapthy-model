"""Run the unchanged historical screen directly on the nested gzip package.

Default scope is exactly the saved historical 126 static runs plus continuous
N3 control. --scope all-legacy-120s is a new expanded analysis; it will not have
the same cohort as the archived figures. No raw files are copied or extracted.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import analyze_rda_bouts as screen
from raw_package_io import read_manifest, legacy_row, locate


def resolve_specs(root, historical_csv, scope='historical', nhost=8):
    manifest = read_manifest(root)
    if scope == 'historical':
        with historical_csv.open(newline='') as stream:
            historical = list(csv.DictReader(stream))
        targets = [row for row in historical if row['condition_id'] != 'canonical_continuous_n3']
        entries = [(legacy_row(manifest,row['state'],row['seed'],nhost,row['run_tag']),row)
                   for row in targets]
    else:
        keys = {(row['state'],row['seed'],row['run_tag']) for row in manifest
                if row['collection'].startswith('legacy_') and int(row['nhost']) == nhost
                and abs(float(row['duration_ms'])-120000) < .5}
        entries = [(legacy_row(manifest,state,seed,nhost,tag),None) for state,seed,tag in sorted(keys)]
    specs=[]
    for entry,historical in entries:
        state,seed,tag=int(entry['state']),int(entry['seed']),entry['run_tag']
        condition=f'state_{state}' + (f'_{tag}' if tag else '')
        specs.append(screen.RunSpec(state=state,seed=seed,run_tag=tag,nhost=nhost,
                                   lfp_path=locate(root,entry,'lfp'),
                                   raster_path=locate(root,entry,'raster'),
                                   duration_s=float(entry['duration_ms'])/1000,
                                   analysis_start_s=float(historical['analysis_start_s']) if historical else 10.,
                                   analysis_end_s=float(historical['analysis_end_s']) if historical else 119.,
                                   condition_id=condition,
                                   label=f'State {state}, seed {seed}' + (f', {tag}' if tag else '')))
    controls=[row for row in manifest if row['collection']=='canonical_dynamic' and int(row['nhost'])==nhost]
    if len(controls)!=1:
        raise ValueError(f'Expected one canonical dynamic control; found {len(controls)}')
    control=controls[0]
    n3=screen.RunSpec(state=None,seed=int(control['seed']),run_tag='',nhost=nhost,
                      lfp_path=locate(root,control,'lfp'),raster_path=locate(root,control,'raster'),
                      duration_s=float(control['duration_ms'])/1000,
                      analysis_start_s=180.,analysis_end_s=210.,
                      condition_id='canonical_continuous_n3',label='Canonical continuous N3, 180–210 s',is_control=True)
    return sorted(specs,key=lambda s:(s.state,s.seed,s.run_tag)),n3


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-package',type=Path,required=True)
    parser.add_argument('--historical-csv',type=Path,default=Path(__file__).resolve().parents[1]/'data/rda_run_metrics.csv')
    parser.add_argument('--scope',choices=['historical','all-legacy-120s'],default='historical')
    parser.add_argument('--nhost',type=int,default=8)
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--review-pdf',type=Path)
    parser.add_argument('--check-layout',action='store_true',help='Check paths and cohort only; no EEG measurements')
    args=parser.parse_args()
    root=args.raw_package.resolve()
    specs,n3=resolve_specs(root,args.historical_csv,args.scope,args.nhost)
    if args.check_layout:
        print(json.dumps({'scope':args.scope,'static_runs':len(specs),
                          'conditions':len({s.condition_id for s in specs}),
                          'control':n3.condition_id,'files_resolved':2*(len(specs)+1),
                          'raw_files_extracted':0},indent=2))
        return
    if args.output_dir is None or args.review_pdf is None:
        parser.error('--output-dir and --review-pdf are required for a measurement run')
    # Replace path discovery only. Original filtering, measurements and reporting
    # functions remain unchanged, including gzip-aware trace/raster readers.
    screen.discover_runs=lambda root,nhost:specs
    screen.canonical_n3_spec=lambda root,nhost:n3
    sys.argv=[str(Path(screen.__file__)), '--root',str(root),'--nhost',str(args.nhost),
              '--output-dir',str(args.output_dir.resolve()),'--review-pdf',str(args.review_pdf.resolve())]
    screen.main()


if __name__=='__main__':
    main()
