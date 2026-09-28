## Current results

Values below are mean +/- SD across three runs with different random seeds,
analyzed over 10--119 seconds. See `analyze_lfp_states.py` for details.

| State | Delta power | Spectral centroid | Interpretation |
|---:|---:|---:|---|
| N3 | 96.5 +/- 0.9% | 1.28 +/- 0.13 Hz | Slow-sleep reference. |
| 35 | 78.7 +/- 0.5% | 3.21 +/- 0.03 Hz | Consistent irregular slowing. |
| 38 | 97.0 +/- 0.8% | 1.22 +/- 0.13 Hz | Consistently slow but essentially N3-like. |
| 39 | 92.2 +/- 4.2% | 1.84 +/- 0.56 Hz | Regularity varied across runs; not robust RDA. |
| 40 | 84.9 +/- 1.3% | 2.57 +/- 0.16 Hz | Consistent slowing, less organized than N3. |
| 42 | 90.1 +/- 1.2% | 2.11 +/- 0.16 Hz | Consistent irregular slowing; strongest current diffuse-slowing candidate. |
| 44 | 91.9 +/- 0.6% | 1.86 +/- 0.07 Hz | Consistent stronger slowing, but more organized and N3-like. |
| 49 | 91.1 +/- 0.6% | 1.98 +/- 0.06 Hz | Consistent midpoint between states 42 and 44. |
| 50 | 90.5 +/- 1.6% | 2.02 +/- 0.18 Hz | Moderate slowing, but the waveform varied across runs. |
| 57 | 72.7 +/- 1.4% | 3.71 +/- 0.12 Hz | Shifted away from delta toward mixed delta-theta activity. |
| 74 | 95.8 +/- 0.4% | 1.40 +/- 0.04 Hz | Spectrally N3-like, but least N3-like among the four candidates in the bout-shape comparison. |
| 99 | 96.3 +/- 0.5% | 1.54 +/- 0.10 Hz | Shorter bouts, but irregular and highly N3-like in shape. |

**No state yet provides a robust, validated RDA phenotype.** The clearest result is a reproducible continuum from irregular slowing toward increasingly organized, N3-like slow activity.

State 39 is the most repetitive candidate. State 74 currently provides the clearest separation from N3. Neither has yet been validated against real GRDA EEG.

