"""Resolve nested raw-package files; gzip streams are read without extraction."""
from __future__ import annotations

import csv
import gzip
import hashlib
from pathlib import Path

import analyze_lfp_states as analysis


def read_manifest(root: Path):
    with (root / 'manifest.csv').open(newline='') as stream:
        return list(csv.DictReader(stream))


def locate(root: Path, row, kind='lfp') -> Path:
    path = (root / row[f'{kind}_file']).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError('Manifest path escapes the raw package')
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def legacy_row(rows, state, seed, nhost=8, run_tag=''):
    matches = [row for row in rows
               if row['collection'].startswith('legacy_')
               and row['state'] == str(state) and row['seed'] == str(seed)
               and row['nhost'] == str(nhost) and row['run_tag'] == run_tag]
    if len(matches) > 1:
        # State-2 seed-1 legacy and explicit filenames contain identical LFP
        # and raster bytes. The original detector prefers the explicit name.
        signatures = {(row['lfp_sha256_txt'], row['raster_sha256_txt']) for row in matches}
        if len(signatures) == 1:
            matches = [sorted(matches, key=lambda r: (f'_seed{seed}_' not in Path(r['lfp_file']).name,
                                                      r['run_id']))[0]]
    if len(matches) != 1:
        raise ValueError(f'Expected one manifest run for {(state,seed,run_tag,nhost)}, got {len(matches)}')
    return matches[0]


def text_sha256(path: Path) -> str:
    """Hash the original text bytes, independent of gzip encoding metadata."""
    digest = hashlib.sha256()
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rb') as stream:
        for block in iter(lambda: stream.read(8*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def analyze_run(root, state, seed, nhost, duration_s, view_fs, band, start_s,
                end_s, run_tag='', allow_seed1_legacy=True):
    """Use the original analyzer with either its flat files or a raw package."""
    root = Path(root)
    if not (root / 'manifest.csv').is_file():
        row, eeg, fs = analysis.analyze_run(root,state,seed,nhost,duration_s,
                                          view_fs,band,start_s,end_s,run_tag,allow_seed1_legacy)
        source = root / row['source_file']
    else:
        entry = legacy_row(read_manifest(root), state, seed, nhost, run_tag)
        if abs(float(entry['duration_ms']) / 1000 - duration_s) > 1e-6:
            raise ValueError(f'Unexpected duration for {entry["run_id"]}')
        source = locate(root, entry)
        # The unchanged loader uses numpy.loadtxt, which supports .gz directly.
        eeg, fs, raw_fs = analysis.load_eegband(source,duration_s,view_fs,band)
        row = {'state':state,'seed':seed,'source_file':source.name.removesuffix('.gz'),
               'run_tag':run_tag,'samples_raw':int(round(raw_fs*duration_s)),
               'duration_s':duration_s,'raw_fs_hz':raw_fs,'view_fs_hz':view_fs,
               'band_low_hz':band[0],'band_high_hz':band[1]}
        row.update(analysis.compute_metrics(eeg,fs,start_s,end_s))
        row['source_archive_path'] = entry['lfp_file']
    row['source_sha256'] = text_sha256(source)
    if (root / 'manifest.csv').is_file() and row['source_sha256'] != entry['lfp_sha256_txt']:
        raise ValueError(f'Original-text hash mismatch for {entry["run_id"]}')
    return row,eeg,fs
