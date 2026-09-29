New readme, slightly updated since Friday but still pending redo of spectral analyses, and specific hypothesis-driven probing of the Fink model 

## Table of contents

### 1. [Comparison of N3 and GRDA](#comparison-of-n3-and-grda)

### 2. [Tinkering with the Fink model](#tinkering-with-the-fink-model)

# Comparison of N3 and GRDA

We compare generalized rhythmic delta activity (GRDA) from **NESI** with N3 slow-wave sleep from **MASS and LENS**. How do they differ in waveform uniformity and rhythmicity? We test this with autocorrelation and cycle measurements. We also compare GRDA with **cyclic N3**, and human N3 with **25 Fink N3 model runs**.

## Selecting representative GRDA from NESI

### Initial selection and review

Arka supplied files with **GRDA labels in ≥80% of 591 overlapping 10-second predictions**. → **8,903 distinct 10-minute NESI filenames** (10,503 entries across overlapping cohorts). [Original selection code](https://github.com/rsarka34/Graham-s-Project/blob/main/SelectionProcessOfPatternsForAnalysis/GRDA_FilenameFinding.py).

We checked five eligibility criteria:

- **Age:** confirmed adult age; we excluded unresolved, conflicting or implausible ages.
- **Breach:** we excluded documented cases and held possible cases. Missing reports leave status unknown.
- **Scores:** valid MORGOTH data; GRDA uniquely highest in ≥80% of predictions.
- **Reference:** verified original EEG reference.
- **Usable excerpt:** at least one 10-second GRDA score ≥0.95, with 10 seconds of EEG on each side for filtering.

After retaining one entry per filename, we had **8,421 files**. We used Python’s random generator to select **1,000 files without replacement**, seed **20260923**. These represented **741 sessions, 553 adults**.

In the first pass, we selected **one 10-second excerpt per file**, with **MORGOTH GRDA ≥0.95**. This gave **1,000 excerpts**. We visually reviewed a sample.

**MORGOTH alone was not enough.** Some high-scoring excerpts looked too irregular. We added waveform rules to find more prototypical GRDA.

[Initial waveform review PDF](docs/grda_n3/MORGOTH_GRDA_initial_waveform_review.pdf). Ten 10-second excerpts. The first four preceded the waveform rules; six additional examples come from later reviews.

### Defining an individual wave

We define **one wave = trough to trough**, lasting **0.25–2 seconds** (0.5–4 Hz). [Trough-to-trough definition PDF](outputs/grda_n3_publication_all_excerpts/reports/GRDA_N3_cycle_examples.pdf).

The detector finds troughs in the **0.5–4 Hz** signal, with spacing ≥0.25 seconds and prominence ≥0.5 × signal SD. The code measures shape and peak-to-peak amplitude between those troughs in the **0.5–30 Hz** signal. The GRDA rules below use these waves.

### Requiring monomorphic, rhythmic runs

GRDA repeats similar waves for at least **six consecutive cycles**. [ACNS 2021 RDA definition, Figure 21, p. 12](https://www.acns.org/UserFiles/file/ACNSStandardizedCriticalCareEEGTerminology_rev2021.pdf#page=12). ACNS describes 0.5–4 Hz activity with similar morphology and duration, no intervening interval, and <50% duration variation in most neighboring pairs.

**Shape:** all five neighboring pairs must have **correlation ≥0.90**. [Shape-rule review PDF](outputs/grda_n3_publication_all_excerpts/grda_run_review/GRDA_six_cycle_rule_review.pdf). The code resamples each wave to 101 points, including both troughs, subtracts its mean and divides by its population SD. This compares shape while allowing amplitude and duration to differ. A run must contain six original neighbors; the code cannot skip a failed wave.

We tried **0.70, 0.80 and 0.90**. In our review, RUN-02/04 lacked uniformity; RUN-05/06 looked excellent. We chose **0.90**. This later review refined the cutoff; the pilot already required six cycles and the timing rule. **We chose the correlation cutoff; ACNS does not specify one.** It rejects some good-looking runs too.

**Timing:** at least **three of five neighboring duration changes must be <50%** of the earlier duration. [Timing-rule review PDF](outputs/grda_n3_publication_all_excerpts/timing_rule_review/GRDA_timing_rule_review.pdf). The code checks `abs(next − previous) / previous < 0.50`. After a 1-second wave, the next wave passes if it lasts >0.5 and <1.5 seconds. Three of five follows ACNS’s majority rule; we use the earlier duration as denominator. Every wave must last 0.25–2 seconds.

**Both sides:** we require a passing run on one left electrode (**Fp1, F7, F3 or C3**) and one right electrode (**Fp2, F8, F4 or C4**). The runs must have **>0 seconds of overlap**. We measure Fz too, but do not count it as either side. This confirms qualifying rhythmic delta on both sides; MORGOTH also favors GRDA over LRDA. We do not measure symmetry or phase synchrony.

**No GRDA voltage cutoff:** we reviewed lower-amplitude candidates, kept the shape/timing rules and dropped the ≥75 µV rule for GRDA. [Review without a GRDA amplitude cutoff](outputs/grda_n3_publication_all_excerpts/grda_run_review_no_amplitude/GRDA_uniformity_without_amplitude_cutoff_review.pdf). We retained ≥75 µV for human N3 waves.

### Applying the rules to the comparison

We searched the sampled files for **non-overlapping 20-second excerpts**. Each excerpt had to pass the waveform rules, with every contained MORGOTH prediction **≥0.85 or ≥0.95**. This left **228 excerpts / 51 adults at ≥0.85**, including **148 excerpts / 42 adults at ≥0.95**.

A qualifying run can fill only part of the excerpt. **We use all 20 seconds for autocorrelation and waves in qualifying runs for GRDA cycle measurements.**

We select GRDA for similar wave shapes and timing. High cycle correlations are therefore partly expected.

## Selecting N3 from MASS and LENS

**MASS** contains sleep-lab PSG with expert staging. We require **expert N3 + MORGOTH N3 ≥0.85**. [MASS sleep annotations](https://doi.org/10.5683/SP3/YD8AYI).

**LENS** contains clinical EEG with sleep-stage predictions on BDSP. We use adult recordings with reviewed normal reports and no identified breach, and require **N3 score ≥0.85**. Expert N3 labels are unavailable here. [LENS dataset](https://www.bdsp.io/content/q8qpxsk3sgq57vkm5abp/1.0.0/).

**Every prediction in the 20 seconds must pass:** both MASS windows; all three LENS windows (0–10, 5–15, 10–20 seconds). We retained every qualifying excerpt from the selected recordings. This gave **MASS: 21 excerpts / eight adults. LENS: 332 excerpts / nine adults.**

### Selecting clearly formed N3 waves

We use the same **trough-to-trough definition** and retain N3 waves with **peak-to-peak amplitude ≥75 µV**. [75 µV rule-review PDF](outputs/grda_n3_publication_all_excerpts/rule_review_75uv/GRDA_N3_cycle_75uV_rule_review.pdf).

Some waves drop out. The code preserves these gaps:

- **Successive-wave correlation:** the code compares original neighbors only when both pass.
- **Template correlation:** the code compares each wave with the median shape of the other retained waves. Gaps are allowed; the calculation requires at least three waves.

This amplitude rule selects individual waves. We use the full 20 seconds for autocorrelation.

### Selecting cyclic N3

Our cyclic N3 rule requires **six consecutive waves**, each ≥75 µV and 0.25–2 seconds long, with **duration SD ÷ mean ≤0.20** and overlapping left/right runs. Their shapes can differ; we apply no shape-correlation cutoff. [N3 cyclicity rule-review PDF](outputs/cyclic_n3_pilot/reports/N3_cyclicity_rule_review.pdf).

This left **MASS: 6 excerpts / four adults. LENS: 112 excerpts / six adults.** We use the full excerpt for autocorrelation and qualifying-run waves for cycle measurements. Selecting regular timing partly determines the timing result.

## Datasets and inclusion

We included **adults only**, excluded known breach and held possible NESI cases. We found no breach in reviewed LENS reports. We cannot confirm breach absence in MASS or NESI, especially when reports are missing.

| Dataset | Cutoff and selection | 20-s excerpts | EEG sessions | Adults |
|---|---|---:|---:|---:|
| NESI GRDA | MORGOTH ≥0.95 + waveform rule | 148 | 46 | 42 |
| NESI GRDA | MORGOTH ≥0.85 + waveform rule | 228 | 60 | 51 |
| MASS N3 | MORGOTH ≥0.85 + expert N3 | 21 | 8 | 8 |
| LENS N3 | N3 score ≥0.85 + normal EEG report | 332 | 10 | 9 |
| MASS cyclic N3 | MASS criteria above + cyclicity rule | 6 | 4 | 4 |
| LENS cyclic N3 | LENS criteria above + cyclicity rule | 112 | 7 | 6 |
| Fink N3 LFP | Canonical N3 preset, state 2 | 225 | 25 model runs | — |

## Fink N3

We ran the canonical N3 model **25 times**, seeds `1001–1025`, for 220 seconds each. We dropped the first 30 seconds and retained nine 20-second excerpts per run from 30–210 seconds. This gave **225 excerpts**.

We measure autocorrelation and cycles in the local LFP. We apply **no voltage cutoff** because we lack a conversion from model LFP amplitude to scalp EEG amplitude.

## Electrodes and filtering

We use **Fp1, Fp2, F7, F8, F3, F4, C3, C4 and Fz**, all present in every included human recording, with the same 15-electrode common-average reference. MASS’s PSG includes these electrodes. [Electrode/reference details](docs/grda_n3/Electrodes.md).

The code filters to **0.5–30 Hz** for autocorrelation and wave shapes, and **0.5–4 Hz** to find troughs. We assume a shared original reference in LENS. NESI absolute voltage units remain unconfirmed.

## Measurements

### Autocorrelation

**Autocorrelation compares a signal with a time-shifted copy of itself.** Repeated peaks show recurring wave structure. The code subtracts the signal mean and normalizes the curve to 1 at zero lag.

We plot the full curve. The code finds its first positive local peak after the zero crossing, at a lag of **0.25–2 seconds**. **Height:** strength of repetition. **Lag:** time until the signal resembles itself again. [ACF rule-review PDF](outputs/grda_n3_publication_all_excerpts/rule_review_75uv/GRDA_N3_ACF_rule_review.pdf).

### Cycle measurements

- **Successive-cycle correlation:** the code compares normalized shapes of original neighboring waves. Both must pass; the code never bridges a rejected wave.
- **Template correlation:** the code compares each retained wave with the median shape of the other retained waves. Gaps are allowed; the calculation requires at least three waves.
- **Period variability:** duration SD ÷ mean duration. Lower values mean more regular timing. The calculation requires at least two waves.

### Power spectrum (PSD)

Deferred.

### Spectral peaks

Deferred.

### Bandpower

Deferred.

## Results

| PDF | Groups compared |
|---|---|
| [Human GRDA versus N3](outputs/grda_n3_report_suite/reports/Human_GRDA_vs_N3.pdf) | Both GRDA cutoffs versus original and cyclic MASS/LENS N3 |
| [Human N3 versus Fink N3](outputs/grda_n3_report_suite/reports/N3_vs_Fink.pdf) | Original and cyclic MASS/LENS N3 versus Fink N3 |
| [All groups](outputs/grda_n3_report_suite/reports/GRDA_vs_all_N3.pdf) | Both GRDA cutoffs versus each human or model N3 group |

The human PDF shows GRDA beside N3, then GRDA beside cyclic N3. The model PDF includes Fink examples too. The plots combine electrodes and show one measurement per slide. Tables give median [25th–75th percentile], n. Original N3 includes its cyclic subset.

<details>
<summary>Detailed selection counts, the 60-second search and MASS assumptions</summary>

## Detailed selection counts and reproducibility

### Eligibility and random sampling

For the sampling pool, we require:

- Verified adult age; we hold missing, conflicting or implausible ages for review.
- No documented breach rhythm or unresolved possible-breach flag. Unavailable reports leave breach status unknown and do not automatically exclude a file.
- Verified source reference and valid activation data.
- GRDA as the uniquely highest-scoring class in at least 80% of predictions.
- At least one 10-second GRDA prediction **≥0.95**, with 10 seconds of available EEG on each side for filtering.

The code retains one eligible entry per filename, preferring RASS, then GCS, CAMS and ICANS. This deduplication does not establish that alternative cohort entries contain identical EEG samples.

We recalculated the highest-scoring class from each prediction’s six scores. Eighteen listed ICANS files fell below Arka’s 80% requirement and were excluded; their presence in the original list remains unexplained.

The resulting pool contains **8,421 files from 2,395 sessions and 1,346 adults**. We excluded 482 filenames: documented breach rhythm (365), possible breach rhythm (55), unverified reference (35), <80% GRDA predictions (18) and implausible age (9). These counts include each filename once.

The code sorts the eligible files and selects **1,000 without replacement** with `random.Random(20260923).sample(...)`. → **741 sessions / 553 adults**. Each file has equal selection probability. Adults may contribute several files. Sampling precedes waveform screening.

We use the same files for both cutoffs. Each sampled file already contains a ≥0.95 window. The ≥0.85 group therefore covers this selected pool; it does not cover all files eligible only at ≥0.85. Its excerpts include the ≥0.95 subset.

### Classifier selection of 20-second excerpts

The code tests excerpts beginning at **10, 30, 50, …, 570 seconds** within each file. Excerpts do not overlap within a file, and each has 10 seconds of real source context on either side for filtering.

At excerpt start `t`, the code checks **11 predictions**: `t, t+1, …, t+10`. GRDA must be uniquely highest and ≥0.85 or ≥0.95 in every prediction. The code records the **lowest of the 11 scores**.

### Selection counts

A **file** is a named 10-minute NESI extract. A **session** is its parent participant/session identifier; one session may contribute multiple files.

| Selection stage | 20-s excerpts | Files | EEG sessions | Adults |
|---|---:|---:|---:|---:|
| Eligible population | — | 8,421 | 2,395 | 1,346 |
| Random sample | — | 1,000 | 741 | 553 |
| All 11 predictions ≥0.85 | 12,319 | 996 | 738 | 551 |
| Above + bilateral waveform rule | **228** | **67** | **60** | **51** |
| All 11 predictions ≥0.95 | 6,074 | 894 | 675 | 506 |
| Above + bilateral waveform rule | **148** | **51** | **46** | **42** |

The ≥0.95 rows are contained within their corresponding ≥0.85 rows.

We excluded one reviewed NESI session: waveform review raised concern for status epilepticus. The diagnosis remains unconfirmed. Its sampled file already failed the 20-second classifier check, so this removed no accepted excerpts.

Our selection prevents overlap within each file. Files from the same session may still overlap. We lack verified session-time offsets; `seg` numbers alone do not give them.

### Search for a full minute of GRDA

The 60-second search requires **all 60 seconds to pass on the same left and right electrodes**. Every sliding six-cycle window must pass. One failed window breaks the stretch. All **51 contained MORGOTH predictions** must also pass. We apply no voltage cutoff.

The code tests starts one second apart and includes 10 seconds of real EEG on each side for filtering. Complete cycles may cross the minute’s edges. The code sets trough prominence from the SD of the full minute.

Classifier screening of the sampled files, after the session exclusion, found:

| Cutoff | Files with a candidate minute | Non-overlapping candidate minutes |
|---|---:|---:|
| ≥0.85 | 638 | 2,007 |
| ≥0.95 | 276 | 669 |

These counts pass the classifier only. We checked waveforms in **1,515 overlapping minutes from 14 local files / 13 adults**. **None passed for a full minute.** Longest passing bilateral stretch: **15.18 seconds**.

Seven tested files come from the random sample; seven come from earlier pilots. The windows overlap heavily. Most sampled files still need waveform screening. We cannot yet say whether the full sample contains a passing minute. This search does not change the existing 20-second results.

[60-second search methods and code](docs/grda_n3/Sixty_second_search.md).

## MASS timing and coverage

We screened SS1 and SS3. SS1 gave three passing 10-second windows, with no adjacent pair. All accepted 20-second excerpts therefore come from **eight SS3 adults**.

Two MASS questions remain:

- **Timing:** we assume each row covers 10 seconds, starting at `filename seconds + 10 × suffix`. `_0`, `_1`, `_2` split a 30-second epoch. Example: `120sec_1` means 130–140 seconds. We still lack confirmation from the export script.
- **Coverage:** scores cover **115 recordings in SS1/SS3**; EEG covers **200 in SS1–SS5**. We do not know why SS2/SS4/SS5 and many SS1/SS3 intervals lack scores.

</details>

## Data and code

Detailed participant-level records remain in the local study workspace.

Numerical summaries appear in the results PDFs. Measurement CSVs, column definitions and comparison scripts remain in the local study workspace.

| Code | Purpose |
|---|---|
| NESI eligibility (`outputs/grda_n3/prepare_nesi_full.py`) | Establish the eligible file population |
| Random sampling (`outputs/grda_n3/full_scale.py`) | Select the frozen sample with Python’s random generator |
| 20-second excerpt selection (`outputs/grda_n3/prepare_strict20.py`) | Apply classifier cutoffs to every contained prediction |
| GRDA waveform rules (`outputs/grda_n3/strict20_cycles.py`) | Detect qualifying six-cycle bilateral runs |
| GRDA–N3 comparison (`outputs/grda_n3/compare_strict20.py`) | Run the human comparison |
| Fink N3 batch (`outputs/fink_n3/run_fink_n3.sh`) | Generate the model runs |
| Report builder (`outputs/grda_n3_report_suite/code/run_report_suite.sh`) | Rebuild PDFs from completed measurements |
| Continuous 60-second GRDA search (`outputs/grda_60s_search/code/search_waveforms.py`) | Test uninterrupted bilateral waveform qualification |

These local scripts use the original environments, data access and saved intermediates. The command below rebuilds the reports in that workspace.

```bash
bash outputs/grda_n3_report_suite/code/run_report_suite.sh
```

---
# Tinkering with the Fink model

# Exploratory thalamocortical modeling of encephalopathic slowing and RDA

This repository extends [Fink et al.'s NEURON port](https://pubmed.ncbi.nlm.nih.gov/39345726/) of the [Krishnan et al. thalamocortical sleep model](https://doi.org/10.7554/eLife.18607). It tests parameter combinations that may produce diffuse background slowing and/or brief 0.5--4 Hz GRDA-like trains.

No sepsis-specific biology or new oscillators have been added. The parameters are not derived from real patient data nor any validated
model of sepsis-associated encephalopathy.

- [State atlas (completed states 0-106)](docs/figures/SAE_state_atlas_completed_states_0_to_106.pdf)
- [State definitions and hypotheses](SAE_STATE_MATRIX.md)
- Results: [all states](docs/figures/state_atlas_metrics_seed1.csv) and [repeated runs](docs/figures/lead_replication_metrics_seeds1_to_3.csv)
- [Fink et al.'s original README](README_FINK_ET_AL.md)
- [Recovered states 49–106, bout-screen code, per-seed results and raw-output handoff](docs/eli_handoff/README.md)

## What changed from Fink et al.

The underlying receptor and channel equations are unchanged (the mechanism files were only made thread-safe; see [Performance](#performance)). Synaptic, intrinsic, and GABA-A timing parameters were modified (preserving original defaults). GABA-A and GABA-B receptor mechanisms existed but scaled multiple pathway strengths shared a single multiplier; these were separated.

| File | Change |
|---|---|
| `config.py` | Adds experimental states and configurable synaptic and intrinsic parameters. |
| `bazh_net.py` | Applies settings, records parameters, and labels outputs by state and run. |
| `network_class.py` | Separates six thalamic connections that previously shared gains. |
| `cell_classes.py` | Makes selected pyramidal-cell conductances configurable. |
| `lfp_recorder.py` | Records the summed cortical voltage and LFP without per-time-step Python work. |
| `run_sweep.py` | Launches single runs and sweeps, compiling the mechanisms when needed. |
| `parallel_layout.py` | Detects cores and MPI, chooses the parallel layout, and calibrates it per machine. |
| `mod/` | Makes seven synapse mechanisms and `xtra.mod` thread-safe. |

Fink et al.'s README is preserved unchanged.

## Experimental path

Completed simulations now cover 95 numbered states: the physiological reference states 0--3 and experimental states through state 106. Most states were run once; selected states were repeated with different random seeds. Testing was iterative, with each group of states informed by the preceding results.

**Recurrent cortical drive** refers to the PYR->PYR AMPA-D2 gain. It controls excitation between pyramidal neurons and also scales small, randomly timed background excitatory inputs generated by the model. TC and RE are thalamocortical relay and reticular neurons.

### Approach 1: start near wake and induce slowing (states 5--33)

| States | Hypothesis family | Result |
|---:|---|---|
| 5-6 | Combined reduced cortical excitation, stronger inhibition, and altered K leak to test a low-excitability SAE hypothesis. | Both states suppressed network activity too strongly (~flat output). |
| 7-8 | Restored only recurrent cortical drive on the state-5/6 backgrounds, to test whether inadequate self-excitation caused suppression. | Activity returned as sharp packets rather than a continuous slow background. |
| 9-12 | Built a milder N1-like wake-to-N2 bridge and separated cortical recurrence, cortical E/I, and thalamic recruitment tests. | Stable active backgrounds with only modest slowing; state 11 became the background anchor. |
| 13-15 | Prolonged cortical, thalamic, or both GABA-A decay constants by 1.5x. | Propofol-inspired kinetic sensitivity tests, not SAE calibration; no compelling target phenotype. |
| 16-18 | Increased recurrent drive and cortical K leak together along the wake-to-N2 axis. | State 17 produced repeated theta packets. State 18 entered a low-output, fast-dominated regime (about 24% delta; centroid 11.2 Hz), not the target. |
| 19-24 | Perturbed separated thalamic pathways. | Packet incidence/duration changed, but frequency remained approximately 6.4--6.9 Hz: theta, not delta. |
| 25-29 | Tested TC Ih/K-leak timing and PYR dendritic NaP/Km/KCa modules. | Changed packet recruitment without reaching the target phenotype. |
| 30-33 | Fine-bracketed recurrence alone versus recurrence with matched cortical K leak. | Located a transition into theta packets; packet-free states showed only modest slowing. |

These wake-based tests produced suppression, modest slowing, or theta packets rather than the target. 
This was not an exhaustive search, and other combinations of synaptic strengths, membrane properties, receptor kinetics, and pathway interactions remain unexplored. 
The results therefore do not rule out a successful wake-based configuration.
Nevertheless, a pivot was made: start from N3 and weaken sleep organization while retaining slow activity.

### Approach 2: start from N3 and remove sleep organization (states 34--48)

| States | Hypothesis family | Result |
|---:|---|---|
| 34-40 | Reduced cortical recurrence/leak or altered thalamic inhibition and membrane properties. | States 35 and 40 produced consistent irregular slowing. State 38 remained N3-like; state 39 varied across runs. |
| 41-44 | Fine-tuned recurrent drive, with or without matched cortical K leak. | States 42 and 44 reproduced consistently across three random seeds. State 42 produced irregular slowing; state 44 was slower but more organized and N3-like. |
| 45-48 | Transferred the state-39 thalamic module or state-38 fast RE->TC GABA-A change onto states 35/40. | State 45 developed intermittent ~6.5-Hz theta packets; state 46 had smaller fast bouts. States 47/48 gave little morphological improvement. |

### Later tests (states 49--106)

| States | What was tested | Main result |
|---:|---|---|
| 49-59 | Further separated cortical recurrence, PYR/INH K leak, and TC/RE membrane settings. | State 49 was a stable midpoint. State 57 shifted toward mixed delta-theta activity. State 59 became the starting point for later tests, but retained broad N3-like plateaus. |
| 60-89 | Tested individual cortical, thalamocortical, GABA-A/GABA-B, K-leak, and TC Ih changes using N3 or state 59. | These changes altered timing and waveform sharpness, but none clearly separated the slow activity from N3. |
| 90-97 | Removed the TC T-current or selected cortex-thalamus connections in short runs with a direct negative TC current. | These tests examined which currents and connections affected the pattern rather than searching for a final candidate. |
| 98-106 | Varied PYR-to-TC and TC-to-cortex excitation, separately and together, using the same TC current. | Some states produced shorter bouts. State 99 was repeated, but remained irregular and highly N3-like. |

Not every defined state in these ranges was run.

## Current results

Values below are mean +/- SD across three runs with different random seeds,
analyzed over 10--119 seconds. All rows use trapezoidal integration, matching `analyze_lfp_states.py` and the saved CSVs. Five earlier README rows used summed Welch bins; the [reconciliation](docs/eli_handoff/README.md#why-the-readme-and-replication-csv-differ) preserves both calculations.

| State | Delta power | Spectral centroid | Interpretation |
|---:|---:|---:|---|
| N3 | 96.0 +/- 0.9% | 1.35 +/- 0.13 Hz | Slow-sleep reference. |
| 35 | 77.2 +/- 0.5% | 3.38 +/- 0.02 Hz | Consistent irregular slowing. |
| 38 | 96.6 +/- 0.8% | 1.28 +/- 0.13 Hz | Consistently slow but essentially N3-like. |
| 39 | 91.4 +/- 4.6% | 1.93 +/- 0.59 Hz | Regularity varied across runs; not robust RDA. |
| 40 | 83.5 +/- 1.2% | 2.74 +/- 0.15 Hz | Consistent slowing, less organized than N3. |
| 42 | 90.1 +/- 1.2% | 2.11 +/- 0.16 Hz | Consistent irregular slowing; strongest current diffuse-slowing candidate. |
| 44 | 91.9 +/- 0.6% | 1.86 +/- 0.07 Hz | Consistent stronger slowing, but more organized and N3-like. |
| 49 | 91.1 +/- 0.6% | 1.98 +/- 0.06 Hz | Consistent midpoint between states 42 and 44. |
| 50 | 90.5 +/- 1.6% | 2.02 +/- 0.18 Hz | Moderate slowing, but the waveform varied across runs. |
| 57 | 72.7 +/- 1.4% | 3.71 +/- 0.12 Hz | Shifted away from delta toward mixed delta-theta activity. |
| 74 | 95.8 +/- 0.4% | 1.40 +/- 0.04 Hz | Spectrally N3-like, but least N3-like among the four candidates in the bout-shape comparison. |
| 99 | 96.3 +/- 0.5% | 1.54 +/- 0.10 Hz | Shorter bouts, but irregular and highly N3-like in shape. |

**No state yet provides a robust, validated RDA phenotype.** The clearest result is a reproducible continuum from irregular slowing toward increasingly organized, N3-like slow activity.

State 39 is the most repetitive candidate. State 74 currently provides the clearest separation from N3. Neither has yet been validated against real GRDA EEG.

## Analysis and atlas

![Wake and N3 LFP, autocorrelation, and power-spectrum comparison](docs/figures/wake_vs_n3_analysis_proof_of_concept.png)

Figure 1. Proof of concept for LFP output analysis: autocorrelation and power spectrum

The atlas contains one page for each of the 95 completed states from 0 through 106, labeled `state_n`, showing the full available trace and final 30 seconds. Traces are filtered to 0.5--30 Hz and scaled separately for readability; amplitude and spectral measures appear on each page. The first 10 seconds are initialization.

Rebuild the atlas and metrics with:

```bash
python build_state_atlas.py
```

## Running a static state

Install the Python packages (`pip install -r requirements.txt`) and, for parallel runs, an MPI library. Without MPI, simulations run on NEURON threads instead, which is slower.

| System | MPI |
|---|---|
| macOS | `brew install open-mpi` |
| Debian/Ubuntu | `sudo apt install openmpi-bin libopenmpi-dev` |
| Fedora/RHEL | `sudo dnf install openmpi`, then `module load mpi/openmpi-$(uname -m)` |
| conda | `conda install -c conda-forge openmpi` |
| Cluster | load the site's MPI module |

`run_sweep.py` runs one or many (state, seed) simulations and writes a `run_state<S>_seed<N>.log` for each. It detects the machine's physical (performance) cores and its MPI installation, respects CPU affinity and container limits, compiles the mechanisms when `mod/` is new or has changed, and checks MPI with a two-rank test before using it. Run `--calibrate` once per machine: it times short simulations in candidate layouts (1-3 minutes) and saves the fastest in `.parallel_layout.json`, which later runs on the same hardware use automatically. Without a calibration, a conservative default is used.

```bash
python run_sweep.py --calibrate                    # once per machine
python run_sweep.py --states 35                    # one run
python run_sweep.py --states 35 40 42 --seeds 1-3  # 9 runs
python run_sweep.py --states 42 --ranks 8          # reproduce an earlier 8-rank run exactly
python run_sweep.py --states 41-44 --dry-run       # show the detected machine, layout, and commands
```

`--dry-run` prints the exact command for the current machine. Manually, one run is:

```bash
FINK_STATE=35 FINK_SEED=1 FINK_DURATION_MS=120000 \
  mpiexec -n <cores> nrniv -mpi -python bazh_net.py 2>&1 | tee run_state35_seed1.log
```

If NEURON cannot find the MPI library (common on macOS), set `MPI_LIB_NRN_PATH`, for example `export MPI_LIB_NRN_PATH=$(brew --prefix open-mpi)/lib/libmpi.dylib`; `run_sweep.py` does this automatically. Without MPI, `FINK_NTHREAD=<cores> nrniv -python bazh_net.py` runs the same simulation on NEURON threads. `--cores` (or `FINK_CORES`) limits the cores used, `--mpiexec` (or `FINK_MPIEXEC`) selects the launcher, e.g. `srun`, `FINK_MPIEXEC_ARGS` adds launcher options, and `FINK_VERBOSE=1` restores the original per-cell setup printout. `python -m unittest discover tests` checks the layout logic, including simulated Linux machines.

Outputs are named by state, seed, and MPI rank count. The analysis scripts use the `--nhost` (default 8) file when present and otherwise the single available rank count. See `requirements.txt` for tested Python packages.

### Performance

Measured on an Apple M3 Max (12 performance and 4 efficiency cores), state 35, 120 simulated seconds:

| Code and layout | Wall time |
|---|---:|
| Original code, 8 MPI ranks | 9.4 min |
| Current code, 8 MPI ranks | 5.5 min |
| Current code, 12 MPI ranks | 4.7 min |

- **LFP recording.** The original code summed the voltage and LFP of every cortical segment in a Python callback on every 0.025-ms step, which cost about as much as the simulation itself. `lfp_recorder.py` records the same values in C (`Vector.record`) and sums them with NumPy once per 50-ms chunk. The output files are byte-identical.
- **Parallel layout.** By default only physical (performance) cores are used: every rank waits for the slowest at each spike exchange, so hyperthreads and efficiency cores usually slow the whole job (`--calibrate` tests them). On the M3 Max, 16 ranks including the efficiency cores were twice as slow as 12, and MPI ranks outperformed NEURON threads (2 simulated seconds: 3.7 s on 12 ranks, 6.3 s on 12 threads).
- **Sweeps.** On the M3 Max, concurrent narrower MPI jobs gave the same throughput as full-width runs one after another, whereas without MPI, single-thread runs side by side gave about 40% more throughput than one 12-thread run at a time. `--calibrate` measures this on each machine.
- **Thread safety.** Seven synapse mechanisms and `xtra.mod` are `THREADSAFE`. `Rinf` and `Rtau` became per-instance `RANGE` variables because NEURON cannot give `NET_RECEIVE` blocks thread-specific globals; every instance computes the same values, so the output files are byte-identical.
- **Analysis.** `analyze_lfp_states.py` and `build_state_atlas.py` process states in parallel (`--jobs`).

A run with the same MPI rank count as before reproduces the original output files byte for byte (checked for a full 120-s run at 8 ranks and for 2-s runs at 1, 2, and 12 ranks). A different rank count changes the order in which simultaneous synaptic events are summed, so membrane voltages differ at floating-point rounding level (median ~1e-13 mV), which action potentials can briefly amplify to ~1e-3 mV. In the full 120-s run at 8 versus 12 ranks, all 511,921 spikes were identical, and 32 of 4.8 million LFP samples differed, each by 0.001. Identical spikes across rank counts are therefore expected but not guaranteed: use `--ranks 8` to reproduce an earlier 8-rank run exactly, and `compare_runs.py` to compare two runs. NEURON threads gave identical spikes and byte-identical LFP and vcort in 2-s tests with 1 to 8 threads.

## Limitations
- Most states have one run. Selected states have additional runs with different random seeds.
- One local cortical LFP cannot establish the bilaterally synchronous, symmetric distribution required for clinical generalized RDA.
- Parameter changes are mechanistic probes, not measurements from septic human or animal brains.
- Raw outputs are excluded from GitHub because they total several gigabytes.
