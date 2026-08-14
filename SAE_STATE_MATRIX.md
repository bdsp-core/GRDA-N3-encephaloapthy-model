# SAE / GRDA exploratory state matrix

These presets are mechanistic experiments, not patient-derived estimates.
State 2 is Fink et al.'s N3 state. The initial experiments started near wake:
state 11 provided a stable active background, and state 17 produced repeated
theta-frequency packets. States 34--48 instead start from N3 and change
selected parameters to determine which ones maintain its organized sleep
waves and which changes preserve slowing while making the activity less
sleep-like.

## Route-specific thalamic refactor

In Fink et al.'s model, one `GABA_thal` setting controlled three inhibitory
pathways, and one `AMPA_thal` setting controlled three excitatory pathways.
For states 0--18, those original shared settings are retained so that their
behavior is unchanged. For later experiments, the six pathway strengths can
be adjusted separately:

- `init_RE_TC_GABA_A`
- `init_RE_TC_GABA_B`
- `init_RE_RE_GABA_A`
- `init_TC_RE_AMPA`
- `init_PYR_TC_AMPA`
- `init_PYR_RE_AMPA`

States 19--24 override individual routes. In the N3-centered branch, states
38, 47, and 48 use the same route-specific controls to weaken fast RE->TC
GABA-A while retaining N3-strength slow RE->TC GABA-B.

## New states

States 19--29 start with all state-17 settings and change only the parameters
listed in the table. States 30--33 similarly start with state 11. States
34--48 start with N3. States 49--54 are the Krishnan et al. Figure 11
comparison states and are described in their own section below. Parameters
not mentioned for a state are therefore unchanged from that state's stated
starting point.

| State | Explicit change | Hypothesis |
|---:|---|---|
| 19 | RE->TC GABA-B factor `0.59125 -> 0.70950` | Stronger slow TC inhibition/rebound |
| 20 | RE->RE GABA-A factor `0.59125 -> 0.4434375` | Reticular disinhibition / stronger coherent RE bursts |
| 21 | RE->TC GABA-A factor `0.59125 -> 0.47300` | Shift relative weight from fast GABA-A toward slow GABA-B |
| 22 | PYR->RE AMPA factor `1.00 -> 1.10` | Selectively strengthen corticoreticular recruitment |
| 23 | State 19 + state 22 | Corticoreticular recruitment plus slow TC inhibition |
| 24 | State 19 + state 20 | Reticular disinhibition plus slow TC inhibition |
| 25 | TC Ih activation shift `-6 -> -5 mV` | Sleepward TC recovery-timing probe |
| 26 | TC Ih activation shift `-6 -> -7 mV` | Opposite-direction sign control |
| 27 | TC K leak `2.0856e-5 -> 2.1804e-5 S/cm2` | 75% wake-to-N2 TC hyperpolarization probe |
| 28 | PYR dendritic Km and KCa `x1.15625` | Outward-adaptation-only test |
| 29 | PYR dendritic NaP, Km, and KCa `x1.15625` | Balanced Figure-11 intrinsic-current module |
| 30 | PYR->PYR AMPA factor `1.06 -> 1.09`; state-11 leaks | Recurrence-only lower bracket |
| 31 | PYR->PYR AMPA factor `1.06 -> 1.105`; state-11 leaks | Recurrence-only upper bracket |
| 32 | PYR->PYR `1.09`, PYR leak `2.24675e-6`, INH leak `1.83825e-6` | Balanced cortical-ACh axis at 37.5% |
| 33 | PYR->PYR `1.105`, PYR leak `2.272875e-6`, INH leak `1.859625e-6` | Balanced cortical-ACh axis at 43.75% |
| 34 | N3; PYR->PYR AMPA `2.7048 -> 2.35` | Mildly weaken cortical Up/Down-state synchronization while preserving the complete N3 substrate |
| 35 | N3; PYR->PYR AMPA `2.7048 -> 2.00` | Moderate recurrence subtraction |
| 36 | N3; PYR->PYR AMPA `2.7048 -> 1.60` | Strong recurrence subtraction and lower boundary for slow-wave survival |
| 37 | N3; PYR and INH K leak each `x0.90` | Test whether N3 Down-state polarization depends on cortical leak independently of synaptic recurrence |
| 38 | N3; RE->TC GABA-A factor `0.715 -> 0.572` | Weaken fast reticulothalamic inhibition while preserving N3 GABA-B |
| 39 | N3; thalamic intrinsic module set to N2: RE leak `9.72e-6`, TC leak `2.2752e-5 S/cm2`, TC Ih shift `-4 mV` | Test whether N3 thalamic membrane polarization is necessary while cortex and all thalamic synapses remain N3 |
| 40 | N3; PYR->PYR `2.00`, PYR/INH K leak each `x0.90` | Test whether modest leak relief preserves activity after recurrence subtraction |
| 41 | N3; PYR->PYR `2.10` | Fine bracket between diffuse-slow state 35 and N3-like state 34 |
| 42 | N3; PYR->PYR `2.20` | Upper fine recurrence bracket |
| 43 | N3; PYR->PYR `2.10`, PYR/INH K leak each `x0.95` | Test modest leak relief at the lower recurrence bracket |
| 44 | N3; PYR->PYR `2.20`, PYR/INH K leak each `x0.95` | Test modest leak relief at the upper recurrence bracket |
| 45 | State 35 plus the state-39 N2 thalamic intrinsic module | Test whether the regularizing thalamic module can organize an irregular slow background |
| 46 | State 40 plus the state-39 N2 thalamic intrinsic module | Same transfer onto the recurrence-plus-leak background |
| 47 | State 35 plus RE->TC GABA-A `0.715 -> 0.572`; GABA-B remains `0.715` | Transfer state 38's selective fast-inhibition reduction onto diffuse slowing |
| 48 | State 40 plus RE->TC GABA-A `0.715 -> 0.572`; GABA-B remains `0.715` | Same route transfer onto the recurrence-plus-leak background |
| 49 | N2 plus PYR dendritic NaP, Km, and KCa `x1.2` | Krishnan Figure 11A extended model, held at N2 |
| 50 | N3 plus PYR dendritic NaP, Km, and KCa `x1.8` | Figure 11A extended model, held at N3 |
| 51 | REM plus PYR dendritic NaP, Km, and KCa `x0.9` | Figure 11A extended model, held at REM |
| 52 | Wake K-leak/Ih/GABA/AMPA plus PYR dendritic NaP, Km, and KCa `x1.2` | Krishnan Figure 11B reduced model at the N2 ACh level |
| 53 | Wake K-leak/Ih/GABA/AMPA plus PYR dendritic NaP, Km, and KCa `x1.8` | Figure 11B reduced model at the N3 ACh level |
| 54 | Wake K-leak/Ih/GABA/AMPA plus PYR dendritic NaP, Km, and KCa `x0.9` | Figure 11B reduced model at the REM ACh level |

States 34--38 each test one change, while state 39 tests a calibrated set of
three thalamic membrane changes. State 40 tests two changes together: reduced
recurrent PYR->PYR AMPA and reduced cortical K leak. These experiments change
selected parts of the N3 model while leaving its other settings unchanged.

Recurrent cortical drive refers to the PYR->PYR AMPA-D2 gain. It scales both
spike-triggered excitation between pyramidal neurons and small, randomly timed
PYR->PYR AMPA events generated by the model. The tested changes are exploratory
model perturbations, not measurements from SAE brains.

States 41--44 test two nearby levels of PYR->PYR AMPA, each with and without a
5% reduction in cortical K leak. This lets us compare states 41 vs 43 and 42
vs 44 to isolate the effect of the K-leak change. States 35 vs 45 and 40 vs 46
isolate the effect of replacing the N3 thalamic membrane settings with N2
settings. States 35 vs 47 and 40 vs 48 isolate the effect of reducing only
RE->TC GABA-A. Repeated runs use the same state number with a different random
seed.

## Results through state 48

- States 19--29 changed packet incidence and duration, but detected packet
  frequencies remained approximately 6.4--6.9 Hz rather than delta.
- State 24 was the most repeatable route-specific packet scaffold; state 29
  supplied the strongest mixed intrinsic-current/background lead.
- State 31 showed that a small further increase in PYR->PYR AMPA changed the
  output from a packet-free background to theta packets.
- States 32--33 remained packet-free and shifted background spectral measures
  modestly toward slowing, but neither was yet convincingly encephalopathic.
- State 34 remained close to N3. State 35 produced irregular predominantly
  delta slowing; state 40 retained more delta but remained moderately
  rhythmic. These are the two diffuse-slowing backgrounds advanced here.
- State 36 crossed into a roughly 5-Hz theta regime. State 37 produced an
  activated, desynchronized cortical field with near-silent thalamus rather
  than cortical suppression.
- Across three seeds, state 35 was a reproducible irregular diffuse-slow
  phenotype and state 40 was a reproducible stronger-slowing compromise.
  State 38 remained essentially N3-like. State 39's unusually regular seed-1
  trace did not reproduce, so state 39 is considered sensitive to the random
  seed rather than a robust RDA candidate.
- States 41--44 strengthened delta as recurrence and paired leak relief
  increased. Across three seeds, state 42 consistently produced irregular
  slowing, while state 44 consistently produced stronger but more organized,
  N3-like slowing.
- State 45 developed intermittent approximately 6.5-Hz theta packets; state
  46 had smaller fast bouts. States 47--48 did not clearly improve morphology.
  States 41, 43, and 45--48 currently have one seed each; states 42 and 44
  have three.

For state 28, dendritic NaP remains `4.2e-5 S/cm2`, while Km and KCa are
`2.3125e-5` and `5.78125e-5 S/cm2`. State 29 additionally sets dendritic NaP
to `4.85625e-5 S/cm2`. Soma NaP and inhibitory-cell intrinsic currents do not
change.

## Figure 11 of Krishnan et al. 2016: extended intrinsic-current model

The baseline model in this port lets ACh and HA act on four things only:
the K leak, the TC Ih activation shift, the GABA conductances, and the AMPA
conductances. Figure 11 of Krishnan et al. reports an extended model in which
ACh also scales three pyramidal *dendritic* intrinsic conductances. The
Materials and methods section "Intrinsic currents" states it directly: "In the
extended model (Figure 11), the conductances of I_d^KCa, I_d^Km and I_d^Nap
were also scaled by AChPY variable." AChPY is the same variable that already
sets the pyramidal K-leak strength, so the multipliers are the ones this port
already uses for K leak: 1 for wake, 1.2 for N2, 1.8 for N3, and 0.9 for REM
(C++ `main.cpp` lines 556--558; the paper's prose quotes 1/1.25/1.8/0.85, and
this port follows the released code so that one ACh axis drives everything).

The resulting dendritic conductances are:

| Current | Wake | N2 | N3 | REM |
|---|---:|---:|---:|---:|
| NaP (`gnabar_nap`) | 4.2e-5 | 5.04e-5 | 7.56e-5 | 3.78e-5 |
| Km (`gkbar_km`) | 2.0e-5 | 2.4e-5 | 3.6e-5 | 1.8e-5 |
| KCa (`gkcabar_kca`) | 5.0e-5 | 6.0e-5 | 9.0e-5 | 4.5e-5 |

All values are S/cm2. The Figure 11 caption also lists the K-leak current and
the H-current, but those are the baseline ACh and HA effects that this port
already implements (`gkl_*` and `gh_TC_*`); no thalamic conductance needs to
change for Figure 11.

Two environment variables select the variant, and both default to off, so with
neither set every state and the baseline cycle behave exactly as before:

- `FINK_FIG11=1` (or `extended`) is panel A: the normal Wake->N2->N3->REM->N2
  cycle plus the ACh scaling of dendritic NaP, Km, and KCa.
- `FINK_FIG11=2` (or `reduced`) is panel B: only the NaP/Km/KCa scaling
  follows the schedule, while K leak, Ih, GABA, and AMPA stay pinned at their
  awake values. In the paper this model "failed to simulate transitions
  between sleep stages."
- `FINK_FIG11_INH=1` additionally scales the interneuron dendritic Km and KCa.
  The Methods sentence quoted above names only the pyramidal dendritic
  currents (I_d^Nap does not exist in the IN model), while the Figure 11
  caption speaks more loosely of "all K+ channels", so the paper does not
  settle whether IN cells were included. The default follows the literal
  Methods reading (pyramidal cells only) rather than inventing a value.

Production commands (each needs the full 360 s of the cycle):

```
FINK_SLEEPSTATES=1 FINK_FIG11=1 FINK_SEED=1 FINK_DURATION_MS=360000 \
    mpiexec -n 8 ./.venv/bin/nrniv -mpi -python bazh_net.py
FINK_SLEEPSTATES=1 FINK_FIG11=2 FINK_SEED=1 FINK_DURATION_MS=360000 \
    mpiexec -n 8 ./.venv/bin/nrniv -mpi -python bazh_net.py
```

Both write `sweep_seed1` output files, so run them in separate directories or
rename the outputs between runs; `plot_paper_fig11.py` expects the panel A and
panel B files to be given explicitly. States 49--54 are the corresponding
static states, for holding one stage indefinitely instead of running the
cycle. Extended-model wake needs no state ID because its ACh factor is 1,
which makes it state 0 exactly.

Newly exposed conductances: `config.init_inh_dend_km` and
`config.init_inh_dend_kca` (consumed in `InhCell.defineBiophysics`). The three
pyramidal dendritic conductances used by Figure 11 were already exposed as
`config.init_pyr_dend_nap`, `config.init_pyr_dend_km`, and
`config.init_pyr_dend_kca` for states 28--29.

## Predeclared advancement criteria

Analyze the edge-free physiological epoch rather than initialization or
filter boundaries. Advance a state if it does either of the following:

1. Produces at least six reasonably uniform cycles with median frequency in
   the 0.5--4 Hz RDA range.
2. Improves event-excluded delta/theta power, slow/fast ratio, and spectral
   centroid relative to state 11 without cortical suppression.

A low-frequency alternating sharp-discharge/silent pattern should be retained
as mechanistically informative but labeled spike-and-wave-like, not GRDA.
One local model LFP cannot establish the bilateral synchrony and symmetry
required to call a pattern clinically "generalized."

Mechanistic background:

- [Bal, Debay & Destexhe 2000](https://pmc.ncbi.nlm.nih.gov/articles/PMC6772790/)
- [Destexhe 1998](https://pmc.ncbi.nlm.nih.gov/articles/PMC6793559/)
- [Krishnan et al. 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5111887/)
- [ACNS critical-care EEG terminology 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8135051/)
