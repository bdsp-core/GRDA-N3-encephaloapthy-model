# Model analysis handoff for Eli

[Raw simulation downloads](https://github.com/graham-mcleod/GRDA-N3-encephaloapthy-model/releases/tag/fink-raw-output-handoff-20260927) · [Restored model code and launch instructions](MODEL_CODE.md) · [Raw-run manifest](manifest.csv)

## References and candidates

- **N3 reference: state 2.** The historical bout screen also uses the original continuous sleep simulation's **180–210 s** N3 interval as its template reference. These are distinct controls.
- **State 39** is the most repetitive model candidate; **state 74** provides the clearest separation from the continuous-N3 waveform among the four plotted candidates. **44** and **99 with −0.30 µA/cm² TC bias** are the other candidate-tradeoff examples. All remain exploratory model patterns.
- States **49/50** are intermediate slowing variants; **57** shifts toward delta–theta activity; **103/106** are thalamocortical pathway tests. The requested per-seed tables include every saved complete 120 s run of 49, 50, 57, 74, 99, 103 and 106: **three seeds each**.
- In the human comparison, **NESI supplies GRDA; MASS and LENS supply N3**. LENS is the selected adult clinical EEG N3 group, not a GRDA group.

[One-line definitions of all states 49–106](data/states_49_to_106_descriptions.csv) are copied from the local state matrix. The TC bias is a separate run intervention: the state number alone does not specify it.

## Why the README and replication CSV differ

The README mixes two numerical conventions. Its five oldest rows (N3, 35, 38, 39 and 40) match **simple sums of Welch frequency-bin powers**. The committed CSV and `analyze_lfp_states.py` use **trapezoidal integration**. Recomputing the same 15 raw traces with bin summation reproduces all five README rows at their displayed precision; trapezoidal integration reproduces the committed CSV exactly. Frequency-band endpoints receive different weights in these calculations.

Both calculations below use seeds 1–3, 10–119 s of the 120 s traces, the same 200 Hz resampling, 0.5–30 Hz filter and 5.12 s Welch windows. Values are mean ± sample SD.

| State | Original README: delta % / centroid Hz | CSV and current analyzer: delta % / centroid Hz |
|---|---|---|
| N3 (2) | 96.5 ± 0.9 / 1.28 ± 0.13 | 96.0 ± 0.9 / 1.35 ± 0.13 |
| 35 | 78.7 ± 0.5 / 3.21 ± 0.03 | 77.2 ± 0.5 / 3.38 ± 0.02 |
| 38 | 97.0 ± 0.8 / 1.22 ± 0.13 | 96.6 ± 0.8 / 1.28 ± 0.13 |
| 39 | 92.2 ± 4.2 / 1.84 ± 0.56 | 91.4 ± 4.6 / 1.93 ± 0.59 |
| 40 | 84.9 ± 1.3 / 2.57 ± 0.16 | 83.5 ± 1.2 / 2.74 ± 0.15 |

The later README rows (42, 44, 49, 50, 57, 74 and 99) already match trapezoidal integration at their displayed precision. Use that convention consistently for the historical Welch summary. The original README numbers are retained in [this excerpt](original_current_results_excerpt.md); [per-seed reconciliation](data/readme_integration_reconciliation_per_seed.csv), [reproduction code](code/reconcile_readme_metrics.py) and [validation](validation.json) document the check. This establishes the numerical explanation; the original ad hoc command that produced the bin-sum table was not recovered.

## Recovered bout screen and figures

The original, previously uncommitted [analyze_rda_bouts.py](code/analyze_rda_bouts.py) and its saved outputs are included unchanged:

- [Per-run measurements](data/rda_run_metrics.csv): 126 complete static 120 s runs plus one continuous N3 control.
- [Per-bout measurements](data/rda_bout_metrics.csv).
- [Per-condition summary](data/rda_condition_summary.csv): 91 conditions, including separately tagged TC-current interventions.
- [Original provenance and thresholds](data/rda_screen_provenance.json).
- [Requested states, per-seed bout measurements](data/requested_states_bout_metrics.csv).

The **95-state atlas** has a different scope: it includes short causal runs that the 120 s bout screen excludes. The delivered raw folder also includes later legacy seeds (states 0 and 2, seed 4) and the separate 25-run N3 batch. The archived bout figures remain based on their saved **126 static runs plus continuous control**. The historical screen's state-2 summary uses seeds 1–3; later state-2 seed 4 is included in the newly exported spectral metrics but does not silently alter the historical screen.

This historical model screen defines cycles by **Hilbert phase crossings in the 0.5–4 Hz signal**. It measures cycle shapes in 0.5–30 Hz, interpolates each to 128 points, median-centers it and divides by its 95th–5th percentile amplitude. It is a separate method from the later trough-to-trough human EEG analysis.

For a six-cycle window to pass, its median frequency must be 0.5–4 Hz; at least 80% of cycle periods must be 0.25–2 s; more than half the neighboring period ratios must be <1.5; delta power fraction must be ≥0.35; median leave-one-out template correlation must be ≥0.65, with at least two-thirds of cycles reaching r≥0.50; and median neighboring-cycle correlation must be ≥0.50. Overlapping windows are merged only while their union still passes. The separate continuity check classifies a bout as RDA-compatible when its quiet-cycle fraction is ≤0.20. There is no absolute voltage gate.

Key exported measures:

| Measure | Saved definition |
|---|---|
| Cycle-period CV | Robust CV = 1.4826 × median absolute deviation / median period. |
| Cycle-template correlation | Median leave-one-out correlation between normalized cycles and the median template of the remaining cycles. |
| Bout burden | Union of qualifying bout intervals divided by the analyzed time; overlapping intervals count once. |
| N3 similarity | Correlation of a condition's median run template with the continuous-N3 template. |

Condition summaries take medians across runs. A “replicated” label requires at least two seeds with compatible bouts, across-seed frequency CV≤0.20 and across-seed template correlation≥0.60. The plotted cross marks a sharp/fast warning in at least one run; it is not an automatic exclusion. Read the original script for its explicitly separate sharpness, fast-power and phase-locking criteria.

Both [original archived PNGs](original_figures/) are included. Their plotting-only command was not found. A [new renderer](code/render_recovered_bout_figures.py) reproduces their axes and data from the recovered CSVs, writes explicit plot-value CSVs, and labels its outputs `reconstructed_*`; it does not claim pixel-identical reproduction. The candidate tradeoff plot uses states 2, 39, 44, 74, 99+TC bias and the continuous N3 control. The underlying analysis code is the unchanged original.

## Spectral exports and reproduction

[Requested states, per-seed spectral metrics](data/requested_states_spectral_metrics.csv) contains 21 rows. [Requested and reference summaries](data/requested_and_reference_summary_seeds1_to_3.csv) contains 14 states, all summarized over seeds 1–3. [The full recomputed table](data/recomputed_requested_and_reference_metrics.csv) contains 47 runs, including additional saved reference/candidate seeds. It records source-file SHA-256 hashes. These are historical Welch metrics provided for traceability; Eli's multitaper analysis can use the original unfiltered traces independently.

Run from this handoff directory, substituting the delivered raw folder containing `manifest.csv`. The adapter resolves the nested collection/run folders and reads `.txt.gz` streams directly; it creates no uncompressed raw-file copies. `source_sha256` refers to the **original text bytes after decompression**, so it matches the manifest's `lfp_sha256_txt`, not the compressed-file checksum.

```bash
RAW="/path/to/Fink model raw outputs for Eli"
python code/recompute_requested_metrics.py --root "$RAW" --output rerun/recomputed_requested_and_reference_metrics.csv --workers 4
python code/reconcile_readme_metrics.py --root "$RAW" --output rerun/readme_integration_reconciliation_per_seed.csv
python code/render_recovered_bout_figures.py
```

The [bout-screen wrapper](code/run_bout_screen_from_package.py) resolves the canonical control through the manifest and selects the exact historical scope from `data/rda_run_metrics.csv`. It replaces file discovery only; the original `analyze_rda_bouts.py` remains unchanged. Check the layout first, then run:

```bash
python code/run_bout_screen_from_package.py --raw-package "$RAW" --check-layout
python code/run_bout_screen_from_package.py --raw-package "$RAW" --output-dir rerun/bout_metrics --review-pdf rerun/bout_review.pdf
```

The optional `--scope all-legacy-120s` adds the two later 120 s legacy runs; its results are a new 128-run analysis. It still excludes the separate longer N3 batch and short causal pilots. Legacy and explicit state-2 seed-1 aliases have identical LFP/raster hashes and count once, as in the original screen.

Checks completed: original synthetic detector self-tests passed; all 21 committed replication rows reproduced exactly for delta percent and centroid; all five discrepant README rows reproduced at displayed precision using bin summation. The gzip adapter reproduced state-49 seed-1 delta power, centroid and original-text hash exactly. Source copies and their hashes are listed in [source_manifest.json](source_manifest.json). No simulations or full bout-screen reruns were run for this handoff.
