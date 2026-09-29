# Continuous 60-second GRDA search

The search uses the same frozen random sample and the GRDA rules described in the [main README](../../README.md#search-for-a-full-minute-of-grda).

A minute qualifies only if all 51 contained 10-second MORGOTH predictions pass the cutoff and a continuous chain of passing six-cycle windows covers the entire minute on one fixed left and one fixed right electrode. Every adjacent shape correlation must be at least 0.90; at least three of five neighboring duration changes must be below 50%. No voltage threshold applies.

Candidate starts are one second apart, with ten seconds of real EEG context on either side. Complete trough-to-trough cycles may extend across the minute edges. One failing six-cycle window breaks a chain. Bilateral coverage does not establish phase synchrony.

| MORGOTH cutoff | Files with candidate minutes | Non-overlapping candidate minutes |
|---|---:|---:|
| At least 0.85 | 638 | 2,007 |
| At least 0.95 | 276 | 669 |

These are classifier-only counts. Waveform screening examined 1,515 overlapping candidate minutes from fourteen locally available files / thirteen adults. None qualified for a full minute. The longest qualifying bilateral stretch was 15.18 seconds. Seven tested files came from the frozen sample and seven from earlier pilots. The full sampled population has not undergone minute-long waveform screening.

Local code: `outputs/grda_60s_search/code/scan_scores.py` and `outputs/grda_60s_search/code/search_waveforms.py`. Original manifests, source EEG and participant-level results remain in the local study workspace and require the original data access.
