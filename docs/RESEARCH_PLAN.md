# Research plan: finding a GRDA regime

This plan incorporates external expert review of the first round of results.
Where the review corrected us, the correction is recorded rather than quietly
absorbed, because several of the corrections change what we should measure.

Companion documents: [FINDINGS.md](FINDINGS.md) for results to date,
[../cpp_comparison/README.md](../cpp_comparison/README.md) for the port
fidelity check.

---

## What the review changed

**1. Our headline conclusion was overstated.** We concluded from the isolated-cell
null that GRDA *must* be a slowed RE-TC loop. It does not follow. The sweep
covered two parameters, TC K-leak and the I_h activation shift; intrinsic TC
delta also depends on g_T, the voltage dependence of I_T, g_h, I_h kinetics,
calcium entry and decay, calcium-dependent regulation of I_h, and tonic input.
The null holds at these conductances, over that slice. Conclusion withdrawn.

**2. Damped ringing is a lead, not a dead end.** Monotonic decay to a fixed
point is the signature of a resonant, near-Hopf system. A weakly damped
resonator under common noise or fluctuating afferent drive produces
intermittent rhythmic runs with no deterministic limit cycle anywhere. This is
attractive precisely because it explains regularity *within* runs and
intermittency *across* runs with one mechanism, instead of requiring a
separate intermittency generator.

**3. GRDA is not "N3 made worse", but neither is it obviously a slowed
spindle.** A third possibility we never considered is a **cortical near-Hopf
E-I oscillator**. Supporting this: a neural-mass study reproduced FIRDA as a
noisy limit cycle and polymorphic delta as filtered noise around a stable
point, and the FIRDA-like regime arose from *increased* effective excitation.
Our entire search reduced excitation. We may have been walking away from the
answer in the cortical dimension as well as the thalamic one.

**4. Our target definition was too narrow.** The ACNS category permits GRDA
with superimposed fast activity, sharply contoured components, and
stimulus-induced forms. "No faster activity" and "suppressed by arousal" are
not universal properties. We scope the initial target explicitly as:

> Spontaneous, non-evolving, frontally predominant GRDA without sharp or spike
> modifiers, occurring in runs, distinct from physiological N3 spindle-slow-wave
> coupling.

Retain absence of spindle *events* and absence of sigma-delta phase-amplitude
coupling as the N3 discriminators. Drop the requirement of no fast activity;
treat GRDA+F as a neighbouring regime. Test arousal reactivity only for the
spontaneous subtype.

**5. Two clinical claims in our briefing were wrong and are corrected here.**
The classical association of FIRDA with deep midline lesions, hydrocephalus
and raised intracranial pressure is *not* supported as a specific association
in modern series; structural lesions and encephalopathy are independently
associated but no etiology is privileged. And in a 4,772-patient critical-care
cohort, GRDA was **not** associated with seizures, unlike LPDs, LRDA and GPDs
— which supports treating canonical GRDA as a distinct target from generalized
spike-and-wave rather than a point on the same continuum.

**6. Anti-NMDAR encephalitis is an etiologic family we omitted.** GRDA
occurred in 50% of patients in one cohort. NMDAR antibodies alter cortical
microcircuit connectivity rather than simply reducing excitation, so NMDA
manipulations must be pathway-specific, not global.

**7. Our hyperammonemia model was too isolated.** Ammonia impairs astrocytic
potassium buffering, raising extracellular K+, which overactivates NKCC1 and
depolarizes E_GABA. So a depolarizing E_Cl shift is a good *dynamical probe*
but not a model of hyperammonemia. Worse for our reasoning: raised [K+]o
shifts E_K, which changes both the K-leak and the potassium current carrying
most of the GABA-B response. Our claim that an E_Cl shift leaves GABA-B and
K-leak "intact" is therefore false in any realistic ammonia model.

**8. Prolonging GABA-A decay alone is demoted.** In the paper's own propofol
simulations, increasing the GABA-A time constant by itself did not produce the
relevant slow and spindle-band activity; coordinated ACh and HA changes were
required.

---

## Phase 0: port timing, in parallel with everything else

The proposed mechanisms all depend on thalamic rebound and synaptic timing, so
the timing discrepancy is not something to defer. But it need not block
qualitative phase mapping either. Two tracks run together.

**Track A, immediate mitigation.** Report every frequency **normalized to the
port's own spindle frequency**. A manipulation that slows the loop by a factor
of three or four is meaningful qualitatively even while the absolute scale is
in doubt. No absolute claim of clinical 2-3 Hz correspondence until Track B
closes.

**Track B, a component-level validation suite.** The order-of-magnitude error
in *event duration* is more concerning than the spindle peak sitting just
below 7 Hz, and the pattern — N2 and REM too slow and too prolonged, N3
roughly right — points at the thalamic subsystem rather than a global
time-base error. Compare against the original C++ and published traces:

| Component | Measurement |
|---|---|
| TC rebound | Latency after matched hyperpolarizing steps |
| TC LTS burst | Spike count and intraburst interval |
| RE burst | Duration and intraburst frequency |
| RE->TC IPSPs | Single-spike and burst-evoked GABA-A and GABA-B amplitude and decay |
| I_h | Activation time course and calcium-dependent shift after a burst |
| I_T | De-inactivation and recovery time course |
| Synaptic depression | Recovery time constant |
| Housekeeping | Synaptic delays, Q10 and temperature, integration timestep, unit conversions |

Suspects for the prolonged events specifically: I_h-mediated termination,
GABA-B kinetics, synaptic depression, and miniature-event statistics.

---

## Phase 1: the first three parameter maps

Run these before any high-dimensional automated search. All three use a
deliberately non-physiological **split thalamocortical state** rather than
starting from N3, following the logic of the paper's own Figure 5: cortical
K-leak alone makes cortex quiescent, recurrent AMPA alone activates it, and
only the combination produces the N3 slow oscillation. For GRDA we want the
quiescent half without the N3 AMPA increase, and then a rhythmic thalamic
drive strong enough to make that damped cortex follow.

**Thalamus:** N2/N3-like TC hyperpolarization so I_T is available; sweep g_h
and g_T; increased CT->RE recruitment; sweep GABA-B dominance over GABA-A;
optionally a true tonic TC inhibitory conductance.

**Cortex:** raise K-leak or tonic inhibition enough to suppress autonomous
cortical oscillation; hold recurrent AMPA at wake/N2 level or below rather
than applying the N3 increase; raise TC->PY and TC->IN gain; retain enough
feedforward inhibition to prevent a sharp population spike.

| # | Map | Axes | Question |
|---|---|---|---|
| 1 | Thalamic loop | GABA-B/GABA-A balance × CT->RE strength, cortex held as damped follower | Is there a smooth rhythmic-delta window before spike-and-wave? |
| 2 | Cortical follower | TC->PY/IN gain × PY->PY recurrence, at 2-3 cortical K-leak levels | Can cortex be made to follow without generating its own rhythm? |
| 3 | Intrinsic TC | g_h × g_T, repeated across several tonic TC conductance or K-leak levels | Does an intrinsic or resonant delta regime exist in the wider parameter set? |

Map 1 must classify outcomes, not just report frequency: spindle-like rhythm,
smooth rhythmic delta, rhythmic delta with a weak sharp component, generalized
spike-and-wave, and silence or depolarization block. A wave-only window is
plausible — the slow thalamic cycle can exist while cortical recruitment stays
below the threshold for a sharp synchronous discharge — but it is not
guaranteed, and what matters is a broad contiguous region, not a lucky point.

---

## Phase 2: characterize resonance, not just oscillation

This is the methodological change that follows from correction 2, and it
applies to every map above. At each parameter point, measure more than whether
spontaneous bursting occurs:

- dominant frequency **and damping constant** of the response to a brief
  perturbation;
- leading eigenvalues of the linearized system where practical;
- response to **independent versus common** stochastic drive;
- whether coupling several heterogeneous TC cells yields coherent quasi-cycles;
- whether calcium-dependent I_h regulation produces natural waxing and waning.

A stable focus with a resonant frequency near 2 Hz and a small negative real
eigenvalue may be as useful as a true limit cycle. `tc_delta_map.py` already
performs the perturbation and records the decay envelope; it needs to fit a
damped sinusoid and report frequency and damping per grid point instead of
classifying points as silent.

For the cortical subsystem, the analogous exercise is to linearize around the
fixed point and look for a complex eigenpair with imaginary part corresponding
to 1.5-3 Hz and real part just below zero.

---

## Phase 3: the cortical E-I oscillator

Untouched so far, and the review rates it a serious alternative. The
experiment is the *opposite* of what has been run: maintain or modestly
increase effective recurrent excitation while balancing it with slow
restorative mechanisms — increased K-leak, increased feedforward or recurrent
inhibition, tonic inhibition, M-current or calcium-dependent potassium
adaptation, slow synaptic depression, cortical GABA-B if added, or cell-type
specific AMPA versus NMDA changes.

The target is a small- to moderate-amplitude oscillation around a **nonzero**
firing state with both PY and IN populations active throughout the cycle — not
an Up/Down cycle.

**Predicted signature:** the rhythm persists when thalamus is clamped or
disconnected; PY-IN phase relationships set the frequency; modest reductions
in cortical excitation abolish it into irregular slowing; firing-rate and
membrane-potential distributions are unimodal rather than bimodal.

A weak prior from a five-subject dipole study places the FIRDA source in
medial frontal or anterior cingulate cortex, which is consistent with a
cortical oscillator amplified through thalamocortical connections. Treat it as
a hint, not a localization.

---

## Phase 4: ionic, pump, and astrocytic mechanisms

The highest-value branch if the objective is to explain **acute metabolic**
GRDA rather than merely to reproduce the waveform. A later model in the same
lineage added extracellular potassium, Na/K pumps and glial buffering, and
raising extracellular potassium drove the cortical network into highly
synchronized 2-3 Hz oscillations — a striking frequency match, though those
authors described it as paroxysmal rather than as clinical GRDA.

Staged so that each stage is useful on its own:

- **Cheap:** clamp [K+]o at several levels and recompute every potassium
  reversal potential; independently shift E_GABA in TC cells, cortical PY
  cells and cortical IN cells. Cell-type specificity matters and these must be
  separate experiments — shifting TC E_GABA may suppress rebound and spindle
  dynamics, while shifting cortical E_GABA may raise effective cortical
  excitation and generate a cortical limit cycle.
- **Intermediate:** dynamic intracellular chloride, or an NKCC1/KCC2 proxy.
- **Full:** import extracellular potassium, Na/K pump and glial buffering
  dynamics from the later model.

For sepsis, hypoxia, hypoglycemia and pump failure this branch is far more
defensible than representing energy failure as an arbitrary leak change.

---

## Phase 5: etiology-specific and structural

**NMDA.** Test pathway-specific rather than global: PY->PY NMDA, PY->IN NMDA,
corticothalamic NMDA where represented, and NMDA changes combined with altered
inhibitory gain.

**Structural.** Do not begin with random TC cell deletion, which mainly
reduces network size and has no specific clinical warrant. More informative:
reduced tonic ascending input to TC cells; loss of independent afferents
leaving a stronger common drive; selective weakening of thalamocortical or
corticothalamic pathways; unilateral disruption to model asymmetric GRDA.
These test deafferentation and synchronization. Run after an oscillator is
established.

---

## Spatial organization and the simulated signal

Two separate problems that we had conflated.

**Neural synchrony.** Introduce several cortical and thalamic modules with
modest parameter heterogeneity, and compare independent local oscillators,
common thalamic drive, corticothalamically synchronized oscillators, and
direct corticocortical synchronization. Measure phase locking, phase delays,
travelling-wave velocity, and robustness of synchrony to heterogeneity.

**Scalp generalization.** A deep or broadly projecting medial-frontal source
can produce a widespread symmetric scalp field without zero phase lag between
cortical patches. Synchrony and generalization are not the same measurement.

**A better output signal.** The paper estimates LFP by averaging PY membrane
voltage; this port already improves on that with a biophysical extracellular
LFP. The next step is a dipole-like output from soma-dendrite axial currents
or summed transmembrane synaptic currents, passed through a simple volume
conduction or scalp forward model. Waveform morphology, harmonic content and
apparent sinusoidality can all change substantially under a realistic
measurement transform, and those are exactly the features we score on.

---

## Evaluation: what would make a simulated pattern credibly GRDA

Five layers. A narrow-band smooth delta rhythm satisfies only the first.

**1. Waveform.** Center frequency and instantaneous-frequency distribution;
spectral peak width or Q; cycle-length CV; autocorrelation decay and side
lobes; rise-decay symmetry; peak and trough sharpness; harmonic power and an
explicit spike-component score; run duration and inter-run interval.

**2. Distinct from N3.** Fraction of time with near-complete cortical silence;
bimodality of PY membrane potential and population rate; Up- and Down-state
duration distributions; physiological spindle count; sigma-delta
phase-amplitude coupling; spatial propagation versus near-synchronous
modulation.

**3. Distinct from generalized spike-and-wave.** No brief high-slope
population discharge preceding each slow wave; low high-frequency burst
fraction; no progressive evolution; lower cycle-locked firing synchrony;
robust separation on waveform sharpness and harmonic structure.

**4. Spatial and behavioural.** Interregional coherence and symmetry; a
plausible frontal or generalized scalp field after forward modeling;
intermittency on the clinical timescale; appropriate response to a simulated
arousal pulse, tested only for the spontaneous subtype.

**5. Robustness and mechanism.** Five to ten seeds; local perturbation of all
parameters by 5-10%; a broad contiguous parameter region rather than a point;
a characteristic phase sequence and a mechanism-specific pharmacological
perturbation.

**Acceptance ranges must be learned from patient data, not chosen by
intuition.** The strongest available validation is to process real adult GRDA,
N3, polymorphic delta and generalized spike-and-wave through the same filters
and montage, train a classifier to separate those human patterns, and apply it
blind to forward-modeled simulation output.

---

## Code to build

Roughly in dependency order.

1. **Damped-sinusoid fitting** in `tc_delta_map.py`: report resonant frequency
   and damping constant per grid point; replace the silent/tonic/bursting
   classification as the primary output.
2. **Expose the intrinsic thalamic conductances**: `gcabar_it_tc`,
   `ghbar_iar`, calcium decay `taur_cad` and depth, and the calcium-I_h
   coupling constants, all currently hard-coded in `cell_classes.py`.
3. **A true tonic inhibitory conductance** on TC cells, distinct from the
   miniature-IPSP mechanism, which is not equivalent.
4. **Split cortical and thalamic state selection**, so a split thalamocortical
   state can be specified without adding a numbered preset per combination.
5. **Independent E_GABA per cell type** — already available for the three
   GABA-A pathways; extend to cortical PY versus IN separately.
6. **[K+]o clamp** that recomputes every potassium reversal potential
   consistently, including the GABA-B current.
7. **Sub-population LFPs and phase-locking metrics**; then a dipole output and
   a simple forward model.
8. **Waveform-morphology metrics**: rise-decay symmetry, peak and trough
   sharpness, harmonic power, spike-component score, sigma-delta PAC.
9. **Multi-module network** with heterogeneity, for the synchrony experiments.

Items 1-3 unblock Phase 1 map 3, item 4 unblocks maps 1 and 2, and items 7-8
are needed before any candidate can be scored honestly.
