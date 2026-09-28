'''
Define baseline parameters for model. References to C++ code are to the code 
originally used in the paper "Cellular and neurochemical basis of sleep stages 
in the thalamocortical network" 
(Krishnan et. al., eLife, 2016, https://doi.org/10.7554/eLife.18607) 
 
The C++ code may be accessed here: 
https://github.com/bazhlab-ucsd/sleep-stage-transition/blob/main
'''

import math
import os

from neuron import h
#from neuron import gui # enable when not running on cluster

#### New ParallelContext object 
pc = h.ParallelContext()
# see https://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/network/parcon.html#ParallelContext.set_maxstep
# as well as section 2.4 of "Simulation Neurotechnologies for Advancing Brain 
# Research: Parallelizing Large Networks in NEURON" (Lytton et. al, 2016)
pc.set_maxstep(10) 
idhost = int(pc.id())
nhost = int(pc.nhost())

# duration of simulation. These environment variables let batch runs select
# their state, seed, and duration without editing this file between runs.
duration = float(os.environ.get("FINK_DURATION_MS", "120000.0")) # ms
t_seg=50.0 #(ms) simulation time between each data dump to node 0

# NEURON threads per MPI rank (all mechanisms are THREADSAFE). Threads and MPI
# ranks can be combined; see run_sweep.py. Output files are reproduced byte for
# byte only with the same MPI rank count; other layouts differ at floating-point
# rounding level (see "Performance" in README.md).
nthread = int(os.environ.get("FINK_NTHREAD", "1"))
if nthread < 1:
    raise ValueError("FINK_NTHREAD must be at least 1")

# Set FINK_VERBOSE=1 to restore the per-cell and per-section setup printout
# and the per-chunk gather messages.
verbose = os.environ.get("FINK_VERBOSE", "0") not in ("", "0")

# set randomizer seed
randSeed = int(os.environ.get("FINK_SEED", "1")) # global RNG seed
h.Random().Random123_globalindex(randSeed) #this changes ALL Random123 streams

# this is True if you want to run through sleep states, according to
# Wake->N2->N3->REM->N2 make this False if you want to just simulate one state
# of vigilance (in which case, select that state by setting the appropriate
# value for 'sleep_state'). FINK_SLEEPSTATES=1 selects the sweep without
# editing this file.

# if do_sleepstates is set to True, then sleep_state will be ignored.
# State IDs: 0=wake, 1=N2, 2=N3, 3=REM, 4=reserved for the later GRDA
# hypothesis, 5=SAE-theta starter, 6=SAE-delta starter,
# 7=SAE-theta rescue R1, 8=SAE-delta rescue R1,
# 9=SAE wakeful-slow/N1-like rescue R2,
# 10=SAE recurrence-edge probe, 11=SAE cortical E/I rebalance,
# 12=SAE conservative thalamic recruitment,
# 13=cortical GABA-A decay probe, 14=thalamic GABA-A decay probe,
# 15=combined cortical+thalamic GABA-A decay probe,
# 16--18=balanced cortical wake-to-N2 trajectory at 50%, 62.5%, and 75%,
# 19--24=route-specific thalamic packet-frequency probes,
# 25--27=TC intrinsic-membrane probes,
# 28--29=PYR dendritic intrinsic-current probes,
# 30--33=stable-background cortical probes,
# 34--40=N3-centered subtraction and decoupling probes,
# 41--48=N3-derived fine brackets and mechanism-transfer probes,
# 49--59=N3-centered cortical fine mapping and gentle thalamic decomposition,
# 60--67=state-59 cortical pathway and random-event amplitude probes,
# 68--89=matched K-leak, thalamocortical-input, GABA, and TC-Ih probes,
# 90--106=causal TC-rebound, pathway-necessity, and loop-gain tests.
do_sleepstates = os.environ.get("FINK_SLEEPSTATES", "0") == "1"
sleep_state = int(os.environ.get("FINK_STATE", "18"))

# Optional conductance-independent hyperpolarizing-current probe for every
# thalamocortical relay (TC) cell.  This is disabled unless a non-zero current
# is supplied at launch, so ordinary states and their filenames are unchanged.
# Current density is in uA/cm2; negative values hyperpolarize.  A short linear
# ramp avoids imposing a simultaneous step transient on all TC cells.
tc_iclamp_uA_cm2 = float(os.environ.get("FINK_TC_ICLAMP_UA_CM2", "0.0"))
tc_iclamp_ramp_ms = float(os.environ.get("FINK_TC_ICLAMP_RAMP_MS", "5000.0"))
tc_iclamp_allow_overwrite = os.environ.get(
    "FINK_TC_ICLAMP_ALLOW_OVERWRITE", "0"
).strip().lower() in {"1", "true", "yes", "on"}
if not math.isfinite(tc_iclamp_uA_cm2):
    raise ValueError("FINK_TC_ICLAMP_UA_CM2 must be finite")
if not -1.0 <= tc_iclamp_uA_cm2 <= 0.0:
    raise ValueError(
        "FINK_TC_ICLAMP_UA_CM2 must be between -1.0 and 0.0 uA/cm2"
    )
tc_iclamp_enabled = tc_iclamp_uA_cm2 != 0.0
if not math.isfinite(tc_iclamp_ramp_ms):
    raise ValueError("FINK_TC_ICLAMP_RAMP_MS must be finite")
if tc_iclamp_ramp_ms < 0.0:
    raise ValueError("FINK_TC_ICLAMP_RAMP_MS must be nonnegative")
if tc_iclamp_enabled and tc_iclamp_ramp_ms > duration:
    raise ValueError(
        "FINK_TC_ICLAMP_RAMP_MS must be between 0 and the run duration"
    )

# Every TC compartment has exactly 29,000 um2 of membrane area, so a current
# density in uA/cm2 converts to nA by multiplying by 0.29.
tc_iclamp_target_nA = 0.29 * tc_iclamp_uA_cm2


def _filename_float(value):
    """Return a deterministic, filename-safe representation of a float."""
    if value == 0.0:
        value = 0.0
    return (
        format(value, ".12g")
        .replace("-", "m")
        .replace("+", "p")
        .replace(".", "p")
    )


if tc_iclamp_enabled:
    tc_iclamp_suffix = "_tcI%s_r%sms" % (
        _filename_float(tc_iclamp_uA_cm2),
        _filename_float(tc_iclamp_ramp_ms),
    )
else:
    tc_iclamp_suffix = ""

# determine whether or not to record LFP
doextra = True

# determines cell density (micrometers^2 per cell) for when cells are placed in
# concentric rings (only need this if doextra==True)
area_cell=100
# x coordinate of recording electrode (in micrometers); see setCellLocations
# method in Net class
XE=2000.0
YE=0.0 #y coordinate of recording electrode (in micrometers)
ZE=0.0 #z coordinate of recording electrode (in micrometers) 

# The summed cortical voltage and biophysical LFP are recorded on every time
# step by lfp_recorder.CorticalFieldRecorder (C-level Vector.record plus one
# NumPy reduction per data dump), which replaced a per-step Python callback
# and produces bit-identical sums.

# set numbers of each cell type (see C++ code network.cfg); note 
# that the method 'connectCells' assumes Npyr=500, Ninh=100, Nre=100, and 
# Ntc=100

Npyr = 500
Ninh = 100
Nre = 100
Ntc = 100

# Optional, behavior-neutral diagnostic for testing whether representative
# thalamocortical relay (TC) cells actually occupy the hyperpolarized voltage
# range associated with intrinsic delta oscillations.  Recording is disabled
# by default so ordinary runs retain their original output and overhead.
record_tc_voltage = os.environ.get(
    "FINK_RECORD_TC_VOLTAGE", "0"
).strip().lower() in {"1", "true", "yes", "on"}
tc_voltage_dt_ms = float(os.environ.get("FINK_TC_VOLTAGE_DT_MS", "1.0"))
tc_voltage_count = min(
    Ntc,
    max(1, int(os.environ.get("FINK_TC_VOLTAGE_COUNT", "8"))),
)
if tc_voltage_dt_ms <= 0:
    raise ValueError("FINK_TC_VOLTAGE_DT_MS must be positive")
if tc_iclamp_enabled and not record_tc_voltage:
    raise ValueError(
        "TC current calibration requires FINK_RECORD_TC_VOLTAGE=1"
    )

_tc_gid_start = Npyr + Ninh + Nre
if tc_voltage_count == 1:
    tc_voltage_gids = (_tc_gid_start + (Ntc - 1) // 2,)
else:
    tc_voltage_gids = tuple(
        _tc_gid_start + (i * Ntc) // tc_voltage_count
        for i in range(tc_voltage_count)
    )

# threshold for detecting spikes (for recording) and for initiating NetCon
# events (mV); see associatedGid and createNetcon and connect2Source in Cell
# class, and connectCells in Net class; this is equivalent to 'Prethresh' in
# C++ currents.cpp
thresh=0

# Synaptic connectivity parameters (see C++ network.cfg)
# note that all strengths are the TOTAL synaptic strength impinging on each
# post-synaptic cell; the actual weight of any one synapse will be the total
# syn strength DIVIDED BY the number of presynaptic connections received by the
# particular neuron species in question
#
# note that all strengths are by default prescribed in NEURON in units of
# microSiemens, while C++ currents.cpp uses milliSiemens. Also,
# currents.cpp prescribes a delay of 0 for all synapses, but this is difficult 
# to implement in NEURON, so they all delays have been set to 0.1 ms
#
# There are many factors of "0.75" because we found that the NEURON model gives
# similar to results to the C++ code when the synaptic strength is
# 0.75 that of the C++ code for pyramidal and inhibitory post-synatpic neurons
# (because we used a full model for these cells, while Bazhenov et. al.
# used a reduced model)

# These factors are equal to 2.0 - ach_level, with ach_level ranging 
# from 0.2 (for s3) to 1.1 (for REM) (It is equal to 1.0 for wake)

s2_scale=1.2  #from C++ main.cpp line 556
s3_scale = 1.8  #from C++ main.cpp line 557
rem_scale=0.9  #from C++ main.cpp line 558

# following factors apply to connections terminating in thalamus (RE->TC GABA-A, 
# RE->TC GABA-B, and RE->RE GABA-A connections, defined below)
# originally referred to as awake_GABA_TC, s2_GABA_TC, s3_GABA_TC, and 
# rem_GABA_TC in lines 597-601 of main.cpp of C++ code
awake_GABA_thal     =0.55
s2_GABA_thal        =awake_GABA_thal*1.15
s3_GABA_thal        =awake_GABA_thal*1.3
rem_GABA_thal       =awake_GABA_thal*0.7

re2tc_gaba_a_rad = 8 #C++ line 23 of network.cfg
re2tc_gaba_a_str = 0.05 #uS
re2tc_gaba_a_del = 0.1 #ms 

re2tc_gaba_b_rad = 8 #C++ line 24 of network.cfg
re2tc_gaba_b_str = 0.002 #uS
re2tc_gaba_b_del = 0.1 #ms

re2re_gaba_a_rad = 5 #C++ line 25 of network.cfg
re2re_gaba_a_str = 0.1 #uS
re2re_gaba_a_del = 0.1 #ms 

# following factors apply to AMPA connections terminating in the thalamus
# (PYR->TC AMPA, PYR->RE AMPA, and TC->RE AMPA connections, defined below)
# note the following parameters were originally named awake_AMPA_TC, 
# s2_AMPA_TC, s3_AMPA_TC, and rem_AMPA_TC in lines 605-608 of main.cpp
# of the original C++ code. They were all originally set to 0.5, but I here set
# them to 1.0, and then multiply the originals synaptic strength values by 0.5
awake_AMPA_thal     = 1.0
s2_AMPA_thal        = 1.0
s3_AMPA_thal        = 1.0
rem_AMPA_thal       = 1.0

pyr2tc_ampa_rad = 10 #C++ line 35 of network.cfg
pyr2tc_ampa_str = 0.5*0.05 #uS; original value of 0.05 is multiplied by 0.5, as noted above
pyr2tc_ampa_del = 0.1 #ms

pyr2re_ampa_rad = 8 #C++ line 36 of network.cfg
pyr2re_ampa_str = 0.5*0.15 #uS; original value of 0.05 is multiplied by 0.5, as noted above
pyr2re_ampa_del = 0.1 #ms

tc2re_ampa_rad = 8 #C++ line 26 of network.cfg
tc2re_ampa_str = 0.5*0.05 #uS; original value of 0.05 is multiplied by 0.5, as noted above
tc2re_ampa_del = 0.1 #ms

# following factors apply to all AMPA connections termining in cortex, other
# than PYR->PYR connections (so this includes TC->PYR, TC->INH, and PYR->INH
# connections, as defined below)
awake_AMPA_cort     =1.0 #this factor is named 'awake_AMPAd2' in C++ main.cpp line 580 and was set to 0.2; here we make this 1.0 instead of 0.2, to emphasize that this is the baseline value
s2_AMPA_cort        =awake_AMPA_cort*(s2_scale + (s2_scale-1)*0.2) # see C++ main.cpp line 582
s3_AMPA_cort        =awake_AMPA_cort*(s3_scale + (s3_scale-1)*0.2) # see C++ main.cpp line 583
rem_AMPA_cort       =awake_AMPA_cort*(rem_scale + (rem_scale-1)*0.2) # see C++ main.cpp line 584

# these are connected with "normal" AMPA synapses (using ampa.mod, not ampa_D2.mod)
tc2pyr_ampa_rad = 10 #C++ line 28 of network.cfg
tc2pyr_ampa_str = 0.75*0.2/5.0 #uS; divide by 5 because we increased awake_AMPA_cort from 0.2 ('awake_AMPAd2' in the C++ code) to 1.0
tc2pyr_ampa_del = 0.1 #ms

# these are connected with "normal" AMPA synapses (using ampa.mod, not ampa_D2.mod)
tc2inh_ampa_rad = 2 #C++ line 29 of network.cfg
tc2inh_ampa_str = 0.75*0.2/5.0 #uS, divide by 5 because we increased awake_AMPA_cort from 0.2 ('awake_AMPAd2' in the C++ code) to 1.0
tc2inh_ampa_del = 0.1 #ms

# These are AMPA_D2 synapses, which have both short-term depression and 
# stochastic EPSP's
# NOTE: if pyr2inh_ampa_d2_str is set to zero, then the program will force
# pyr2inh_ampa_d2_mini_str to zero as well (see createSynapses methods in
# cell_classes.py)
pyr2inh_ampa_d2_rad = 1  #C++ line 29 of network.cfg
pyr2inh_ampa_d2_str = 0.75*0.12/5.0 #uS, divide by 5 because we increased awake_AMPA_cort from 0.2 ('awake_AMPAd2' in the C++ code) to 1.0
pyr2inh_ampa_d2_del = 0.1 #ms
# strength of stochastic EPSP's
pyr2inh_ampa_d2_mini_str = 0.75*0.20/5.0  #uS; divide by 5 because we increased awake_AMPA_cort from 0.2 ('awake_AMPAd2' in the C++ code) to 1.0
# this corresponds to 'mini_fre' in the C++ code (see line 520 of currents.cpp)
# despite the name implying this is a frequency, it is really value of time. The
# largeer this value is, the less frequent the stochastic EPSP's (see function 
# 'gen_nextpsp' in ampa_D2.mod)
pyr2inh_ampa_d2_mini_f = 20.0 #ms; 

# scaling of PYR->PYR AMPA D2 synapses for each sleep stage. Note that this is 
# treated differently than in the C++ code, and these values were selected 
# in order to get obtain better spindles in N2
awake_AMPA_pyrpyr = 1.0
s2_AMPA_pyrpyr = 1.24
s3_AMPA_pyrpyr = 2.7048
rem_AMPA_pyrpyr = 1.056 

# These are AMPA_D2 synapses, which have both short-term depression and 
# stochastic EPSP's
# NOTE: if pyr2pyr_ampa_d2_str is set to zero, then the program will force
# pyr2pyr2_ampa_d2_mini_str to zero as well (see createSynapses methods in
# cell_classes.py)
pyr2pyr_ampa_d2_rad = 5 #C++ line 31 of network.cfg
pyr2pyr_ampa_d2_str = 0.03 #this value deviates from the C++ code in order to obtain better spindles in N2
pyr2pyr_ampa_d2_del = 0.1 #ms
pyr2pyr2_ampa_d2_mini_str = (0.33/0.24) * 0.03 #this value deviates from the C++ code in order to obtain better spindles in N2; the ratio 0.33/0.24 comes from the original scaling in the C++ code (see line 31 of network.cfg)
# note that ratio pyr2pyr2_ampa_d2_mini_str/pyr2pyr_ampa_d2_str determines PSP
# weight in AMPA_D2.mod, so that gmax associated with stochastic EPSP's changes 
# when gmax associated with pyr2pyr_ampa_d2_str changes (as it does in the 
# full sleep states simulation)
# this corresponds to 'mini_fre' in the C++ code (see line 520 of currents.cpp)
# despite the name implying this is a frequency, it is really value of time. The
# largeer this value is, the less frequent the stochastic EPSP's (see function 
# 'gen_nextpsp' in ampa_D2.mod)
pyr2pyr2_ampa_d2_mini_f = 20.0 #ms

# D1 synapses have short-term depression, but no stochastic EPSP's (these
# values do not vary with sleep stage)
pyr2pyr_nmda_d1_rad = 5 #C++ line 32 of network.cfg
pyr2pyr_nmda_d1_str = 0.75*0.01 #uS
pyr2pyr_nmda_d1_del = 0.1 #ms
# unitless factor (between 0 and 1) which determines the degree of short-term
# depression experienced with each presynaptic spike; C++ code sets this to
# 0.0, which means there is no depression
pyr2pyr_nmda_d1_Use = 0.0

# D1 synapses have short-term depression, but no stochastic EPSP's (these
# values do not vary with sleep stage)
pyr2inh_nmda_d1_rad = 1 #C++ line 34 of network.cfg
pyr2inh_nmda_d1_str = 0.75*0.01 #uS
pyr2inh_nmda_d1_del = 0.1 #ms
# unitless factor (between 0 and 1) which determines the degree of short-term
# depression experienced with each presynaptic spike; Krishnan sets this to
# 0.0, which means there is no depression
pyr2inh_nmda_d1_Use = 0.0

# following factors apply to INH->PYR GABA-A connections. See lines 591-595 
# of C++ main.cpp
awake_GABA_D2      =0.22
s2_GABA_D2         =awake_GABA_D2*1.15
s3_GABA_D2         =awake_GABA_D2*1.3
rem_GABA_D2        =awake_GABA_D2*0.7

# These are GABA_A_D2 synapses, which have both short-term depression and 
# stochastic EPSP's
# NOTE: if inh2pyr_gaba_a_d2_str is set to zero, then the program will force
# inh2pyr_gaba_a_d2_mini_str to zero as well (see createSynapses methods in
# cell_classes.py)
inh2pyr_gaba_a_d2_rad = 5 #C++ line 38 of network.cfg
inh2pyr_gaba_a_d2_str = 0.75*0.24 #uS
inh2pyr_gaba_a_d2_del = 0.1 #ms
# strength of stochastic IPSP's
inh2pyr_gaba_a_d2_mini_str = 0.75*0.20  #uS
# this corresponds to 'mini_fre' in the C++ code (see line 520 of currents.cpp)
# despite the name implying this is a frequency, it is really value of time. The
# largeer this value is, the less frequent the stochastic EPSP's (see function 
# 'gen_nextpsp' in gaba_A_D2.mod)
inh2pyr_gaba_a_d2_mini_f = 20.0 #ms; this parameter is involved in calculating the stochastic EPSP times 


# cellular properties that vary with sleep stage
gkl_pyr_awake         = 0.19 * 0.000011 #S/cm2; factor of 0.19 from C++ main.cpp line 561, and 0.000011 from C++ CellSyn.h line 362
gkl_pyr_s2            = gkl_pyr_awake*s2_scale
gkl_pyr_s3            = gkl_pyr_awake*s3_scale
gkl_pyr_rem           = gkl_pyr_awake*.9

gkl_inh_awake         = 0.19 * 0.000009 #S/cm2; factor of 0.19 from C++ main.cpp line 561, and 0.000011 from C++ CellSyn.h line 525
gkl_inh_s2            = gkl_inh_awake*s2_scale
gkl_inh_s3            = gkl_inh_awake*s3_scale
gkl_inh_rem           = gkl_inh_awake*.9

gkl_TC_awake      = 0.79 * 0.000024 # S/cm2; factor of 0.79 from C++ main.cpp line 567, and 0.000024 from C++ CellSyn.h line 241
gkl_TC_s2         = gkl_TC_awake*s2_scale
gkl_TC_s3         = gkl_TC_awake*s3_scale
gkl_TC_rem        = gkl_TC_awake*.9

gkl_RE_awake      = 0.9 * 0.000012 # S/cm2; factor of 0.9 from C++ main.cpp line 573, and 0.000012 from C++ CellSyn.h line 177
gkl_RE_s2         = gkl_RE_awake*((2-s2_scale/2)-0.5)
gkl_RE_s3         = gkl_RE_awake*((2-s3_scale/2)-0.5)
gkl_RE_rem        = gkl_RE_awake*1.1

gh_TC_awake       =-8.0 #mV; see lines 586-589 from C++ main.cpp
gh_TC_s2          =-4.0 #C++ code uses -3.0, but -4.0 was found to give better spindles in N2
gh_TC_s3          =-2.0
gh_TC_rem         = 0.0


# Exploratory sepsis-associated encephalopathy (SAE) starter presets.
# These are mechanistically informed hypotheses for parameter exploration,
# not direct patient-derived measurements.  They intentionally remain within
# or close to the model's wake-to-N3 parameter envelope.  State 4 remains
# reserved for the later GRDA-focused hypothesis.

# SAE-theta: mild/moderate encephalopathic slowing with an N2-like increase in
# cortical leak/inhibition, slightly reduced cortical excitation, and TC Ih
# kept closer to wake than sleep.
sae_theta_GABA_thal   = 0.55
sae_theta_AMPA_thal   = 1.00
sae_theta_AMPA_cort   = 0.95
sae_theta_AMPA_pyrpyr = 0.95
sae_theta_GABA_D2     = 0.253
sae_theta_gkl_pyr     = 2.508e-6   # S/cm2
sae_theta_gkl_inh     = 2.052e-6   # S/cm2
sae_theta_gkl_RE      = 9.720e-6   # S/cm2
sae_theta_gkl_TC      = 2.2752e-5  # S/cm2
sae_theta_gh_TC       = -6.0       # mV shift of Ih kinetics; not conductance

# SAE-delta: deeper diffuse slowing, with stronger membrane
# hyperpolarization/inhibition and reduced excitatory gains, but without the
# very strong PYR->PYR recurrence used to generate N3 Up/Down oscillations.
sae_delta_GABA_thal   = 0.62
sae_delta_AMPA_thal   = 0.95
sae_delta_AMPA_cort   = 0.90
sae_delta_AMPA_pyrpyr = 0.90
sae_delta_GABA_D2     = 0.264
sae_delta_gkl_pyr     = 3.135e-6   # S/cm2
sae_delta_gkl_inh     = 2.565e-6   # S/cm2
sae_delta_gkl_RE      = 9.720e-6   # S/cm2
sae_delta_gkl_TC      = 2.844e-5   # S/cm2
sae_delta_gh_TC       = -4.0       # mV shift of Ih kinetics; not conductance

# Rescue R1 after states 5 and 6 produced genuine suppression.  Change only
# recurrent PYR->PYR AMPA_D2 gain, restoring it to the calibrated wake value;
# all other SAE background settings remain fixed for a clean causal test.
sae_theta_r1_AMPA_pyrpyr = 1.00
sae_delta_r1_AMPA_pyrpyr = 1.00

# Rescue R2: conservative N1-like calibration bridge between wake and N2.
# These values are the midpoint of the calibrated wake-to-N2 trajectory,
# except recurrent PYR->PYR AMPA is held at the wake value to avoid promoting
# the spindle-like population packets observed in the R1 experiments.
sae_n1_r2_GABA_thal   = 0.59125
sae_n1_r2_AMPA_thal   = 1.00
sae_n1_r2_AMPA_cort   = 1.12
sae_n1_r2_AMPA_pyrpyr = 1.00
sae_n1_r2_GABA_D2     = 0.2365
sae_n1_r2_gkl_pyr     = 2.299e-6   # S/cm2
sae_n1_r2_gkl_inh     = 1.881e-6   # S/cm2
sae_n1_r2_gkl_RE      = 1.026e-5   # S/cm2
sae_n1_r2_gkl_TC      = 2.0856e-5  # S/cm2
sae_n1_r2_gh_TC       = -6.0       # mV shift of Ih kinetics; not conductance

# State 10: recurrence-edge hypothesis.  This changes only recurrent PYR->PYR
# AMPA_D2 gain relative to state 9, using one quarter of the remaining
# state-9-to-N2 gain difference (1.00 -> 1.24).
sae_rec_edge_GABA_thal   = 0.59125
sae_rec_edge_AMPA_thal   = 1.00
sae_rec_edge_AMPA_cort   = 1.12
sae_rec_edge_AMPA_pyrpyr = 1.06
sae_rec_edge_GABA_D2     = 0.2365
sae_rec_edge_gkl_pyr     = 2.299e-6   # S/cm2
sae_rec_edge_gkl_inh     = 1.881e-6   # S/cm2
sae_rec_edge_gkl_RE      = 1.026e-5   # S/cm2
sae_rec_edge_gkl_TC      = 2.0856e-5  # S/cm2
sae_rec_edge_gh_TC       = -6.0       # mV shift of Ih kinetics

# State 11: cortical E/I-rebalance hypothesis.  Keep state 9's membrane and
# thalamic background, reduce PYR->INH recruitment and INH->PYR inhibition,
# and add the same mild recurrent gain tested independently in state 10.
sae_ei_GABA_thal   = 0.59125
sae_ei_AMPA_thal   = 1.00
sae_ei_AMPA_cort   = 1.06
sae_ei_AMPA_pyrpyr = 1.06
sae_ei_GABA_D2     = 0.22825
sae_ei_gkl_pyr     = 2.299e-6   # S/cm2
sae_ei_gkl_inh     = 1.881e-6   # S/cm2
sae_ei_gkl_RE      = 1.026e-5   # S/cm2
sae_ei_gkl_TC      = 2.0856e-5  # S/cm2
sae_ei_gh_TC       = -6.0       # mV shift of Ih kinetics

# State 12: conservative thalamic-recruitment hypothesis.  Keep the proven
# state 9 cortical background, return RE leak to wake, and move TC leak and
# the Ih activation shift only one cautious step toward N2.
sae_thal_recruit_GABA_thal   = 0.59125
sae_thal_recruit_AMPA_thal   = 1.00
sae_thal_recruit_AMPA_cort   = 1.12
sae_thal_recruit_AMPA_pyrpyr = 1.00
sae_thal_recruit_GABA_D2     = 0.2365
sae_thal_recruit_gkl_pyr     = 2.299e-6   # S/cm2
sae_thal_recruit_gkl_inh     = 1.881e-6   # S/cm2
sae_thal_recruit_gkl_RE      = 1.080e-5   # S/cm2
sae_thal_recruit_gkl_TC      = 2.2752e-5  # S/cm2
sae_thal_recruit_gh_TC       = -5.0       # mV shift of Ih kinetics

# States 16--18: balanced cortical-ACh-axis sweep.  State 11 placed cortical
# K-leak halfway from wake to N2, but recurrent PYR->PYR AMPA only one quarter
# of the way.  These states first repair that mismatch and then advance leak
# and recurrence together.  Every other cortical and thalamic parameter is
# held exactly at state 11.
sae_bal50_AMPA_pyrpyr = 1.12
sae_bal50_gkl_pyr     = 2.299e-6
sae_bal50_gkl_inh     = 1.881e-6

sae_bal625_AMPA_pyrpyr = 1.15
sae_bal625_gkl_pyr     = 2.35125e-6
sae_bal625_gkl_inh     = 1.92375e-6

sae_bal75_AMPA_pyrpyr = 1.18
sae_bal75_gkl_pyr     = 2.4035e-6
sae_bal75_gkl_inh     = 1.9665e-6

# States 19--24: route-specific thalamic packet-frequency probes on the
# complete state-17 background.  The legacy model uses one GABA multiplier
# for RE->TC GABA_A, RE->TC GABA_B, and RE->RE GABA_A, and one AMPA
# multiplier for TC->RE, PYR->TC, and PYR->RE.  These values are applied only
# after those shared controls have been expanded into pathway-specific
# initializers below, leaving states 0--18 unchanged.
sae_re_tc_gaba_b_1p20 = sae_ei_GABA_thal * 1.20
sae_re_re_gaba_a_0p75 = sae_ei_GABA_thal * 0.75
sae_re_tc_gaba_a_0p80 = sae_ei_GABA_thal * 0.80
sae_pyr_re_ampa_1p10 = sae_ei_AMPA_thal * 1.10

# States 25--27: conservative TC intrinsic-membrane probes on state 17.
# The leak value is the 75% point between calibrated wake and N2, while
# state 17 itself is at the 50% point.
sae_tc_ih_sleepward = -5.0
sae_tc_ih_wakeward = -7.0
sae_tc_gkl_75 = gkl_TC_awake + 0.75 * (gkl_TC_s2 - gkl_TC_awake)

# States 28--29: Figure-11-inspired PYR dendritic current module.  The
# 1.15625 factor is 62.5% of the calibrated wake-to-N2 change (1.0 -> 1.25),
# matching the cortical depth of state 17.  State 28 changes the two outward
# adaptation currents only; state 29 adds persistent sodium for the balanced
# module.  Soma NaP and all INH intrinsic currents remain fixed.
pyr_dend_nap_baseline = 0.000042
pyr_dend_km_baseline = 0.000020
pyr_dend_kca_baseline = 0.000050
pyr_intrinsic_62p5_scale = 1.15625

init_pyr_dend_nap = pyr_dend_nap_baseline
init_pyr_dend_km = pyr_dend_km_baseline
init_pyr_dend_kca = pyr_dend_kca_baseline

# States 30--33: two recurrence-only brackets between state 11 (1.06) and
# state 16 (1.12), plus paired AMPA/K-leak versions at the corresponding
# 37.5% and 43.75% positions along the calibrated cortical ACh axis.
sae_bg_rec37_AMPA_pyrpyr = 1.090
sae_bg_rec44_AMPA_pyrpyr = 1.105

sae_bg_bal37_AMPA_pyrpyr = 1.090
sae_bg_bal37_gkl_pyr = gkl_pyr_awake + 0.3750 * (gkl_pyr_s2 - gkl_pyr_awake)
sae_bg_bal37_gkl_inh = gkl_inh_awake + 0.3750 * (gkl_inh_s2 - gkl_inh_awake)

sae_bg_bal44_AMPA_pyrpyr = 1.105
sae_bg_bal44_gkl_pyr = gkl_pyr_awake + 0.4375 * (gkl_pyr_s2 - gkl_pyr_awake)
sae_bg_bal44_gkl_inh = gkl_inh_awake + 0.4375 * (gkl_inh_s2 - gkl_inh_awake)

# States 34--40: N3-centered subtraction and decoupling screen.  These are
# not points on a presumed N3-to-wake trajectory.  Each state starts from the
# complete calibrated N3 parameter vector and selectively weakens one
# candidate organizer of normal N3 activity. State 40 tests two changes
# together: reduced recurrent PYR->PYR AMPA and reduced cortical K leak.
# Slow RE->TC GABA-B remains at N3 in all seven states.
n3_sub_recur_mild = 2.35
n3_sub_recur_moderate = 2.00
n3_sub_recur_strong = 1.60
n3_sub_cortical_leak_scale = 0.90
n3_sub_gkl_pyr = gkl_pyr_s3 * n3_sub_cortical_leak_scale
n3_sub_gkl_inh = gkl_inh_s3 * n3_sub_cortical_leak_scale
n3_sub_RE_TC_GABA_A = s3_GABA_thal * 0.80

# States 41--48: follow-up around diffuse-slowing states 35 and 40.  States
# 41--44 form a small recurrence/leak grid. States 45--48 transfer the two
# regularizing mechanisms from states 39 and 38 separately onto states 35
# and 40.  Replicate seeds are repeated runs, not new state IDs.
n3_follow_recur_low = 2.10
n3_follow_recur_high = 2.20
n3_follow_leak_scale = 0.95
n3_follow_gkl_pyr = gkl_pyr_s3 * n3_follow_leak_scale
n3_follow_gkl_inh = gkl_inh_s3 * n3_follow_leak_scale

# States 49--51: start from state 42 (PYR->PYR AMPA 2.20). State 49 tests a
# midpoint between states 42 and 44 by reducing both cortical K leaks 2.5%.
# States 50 and 51 apply the full 5% reduction to PYR or INH K leak alone.
n3_leak_mid_scale = 0.975
n3_leak_mid_gkl_pyr = gkl_pyr_s3 * n3_leak_mid_scale
n3_leak_mid_gkl_inh = gkl_inh_s3 * n3_leak_mid_scale

# States 52--59: complete interpretable local comparisons around states 44
# and 49. States 52--53 isolate the 2.5% leak change by cortical cell type;
# states 54--55 and 59 refine the recurrence/leak grid; states 56--58 move
# one thalamic intrinsic control halfway from N3 toward N2 while every other
# state-49 setting is held fixed.
n3_recur_mid = 0.5 * (n3_follow_recur_low + n3_follow_recur_high)
n3_gh_TC_mid_to_N2 = 0.5 * (gh_TC_s3 + gh_TC_s2)
n3_gkl_TC_mid_to_N2 = 0.5 * (gkl_TC_s3 + gkl_TC_s2)
n3_gkl_RE_mid_to_N2 = 0.5 * (gkl_RE_s3 + gkl_RE_s2)

# States 60--67: state-59 cortical-route decomposition. States 60--65 vary
# one of the three pathways formerly grouped under AMPA_cort. States 66--67
# vary only the amplitude of randomly timed PYR->PYR AMPA-D2 events while
# spike-triggered PYR->PYR recurrence remains fixed at 2.15.
cortical_route_probe_low = 0.90
cortical_route_probe_high = 1.10
init_PYR_PYR_AMPA_D2_random_event_scale = 1.0

# States 68--89: matched tests requested after the initial N3-centered search.
# The TC K-leak series is crossed with either full N3 recurrence or the
# state-59 value (2.15). Cortical K-leak tests start from full N3. Functional
# thalamocortical-input and GABA tests use both full N3 and/or state 59 as
# stated below. These are sensitivity probes, not patient-derived estimates.
n3_tc_gkl_low = gkl_TC_s3 * 0.90
n3_tc_gkl_high = gkl_TC_s3 * 1.10
n3_tc_gkl_high2 = gkl_TC_s3 * 1.20

n3_pyr_gkl_high = gkl_pyr_s3 * 1.05
n3_pyr_gkl_high2 = gkl_pyr_s3 * 1.10
n3_inh_gkl_high = gkl_inh_s3 * 1.05

n3_tc_input_low = 0.90

# Coarse route-gain brackets used by the causal closed-loop screen.  These
# scale one already-separated pathway (or the paired ascending TC routes)
# without changing connection counts, receptor kinetics, or other N3/state-59
# controls.
causal_gain_0p25 = 0.25
causal_gain_0p50 = 0.50
causal_gain_0p75 = 0.75
causal_gain_1p25 = 1.25
n3_thal_gaba_high = s3_GABA_thal * 1.10
n3_thal_gaba_high2 = s3_GABA_thal * 1.20
n3_cort_gaba_high = s3_GABA_D2 * 1.05
n3_cort_gaba_high2 = s3_GABA_D2 * 1.10

# Fink's TC Ih maximal conductance was previously hard-coded in TCCell.
# Exposing it with the exact old value as the default preserves states 0--87.
TC_Ih_gbar_baseline = 1.6e-5
init_TC_Ih_gbar = TC_Ih_gbar_baseline
TC_Ih_gbar_low = TC_Ih_gbar_baseline * 0.90

# Fink's TC T-current maximal conductance was also hard-coded in TCCell.
# Expose it with the exact original value as the default so states 0--89 are
# unchanged. State 90 turns it off as a causal rebound-current test.
TC_T_gbar_baseline = 0.0025
init_TC_T_gbar = TC_T_gbar_baseline

# Figure-9-inspired GABA_A kinetic sensitivity branches.  States 13--15 use
# the complete state-11 conductance background and alter only GABA_A decay.
# Since the mechanisms decay as exp(-Beta*t), dividing Beta by 1.5 lengthens
# the nominal off-decay time constant by 1.5.  GABA_B is not changed.
Beta_GABA_A_baseline = 0.166       # /ms; RE->TC and RE->RE
Beta_GABA_A_D2_baseline = 0.18     # /ms; cortical INH->PYR
Beta_GABA_A_1p5x_decay = Beta_GABA_A_baseline / 1.5
Beta_GABA_A_D2_1p5x_decay = Beta_GABA_A_D2_baseline / 1.5

# Explicit defaults prevent kinetic settings from carrying into other states.
init_Beta_GABA_A = Beta_GABA_A_baseline
init_Beta_GABA_A_D2 = Beta_GABA_A_D2_baseline

if do_sleepstates:
    # this is where you specify the initial state of vigilance; these values
    # are used to instantiate the network in the 'connectCells' method
    init_GABA_thal = awake_GABA_thal
    init_AMPA_thal = awake_AMPA_thal
    init_AMPA_cort = awake_AMPA_cort
    init_AMPA_pyrpyr = awake_AMPA_pyrpyr
    init_GABA_D2 = awake_GABA_D2
    init_gkl_pyr = gkl_pyr_awake
    init_gkl_inh  = gkl_inh_awake
    init_gkl_RE  = gkl_RE_awake
    init_gkl_TC  = gkl_TC_awake
    init_gh_TC   = gh_TC_awake
    # specify transition times between sleep states (in order to replicate
    # Figs. 1 and 2 in Bazhenov 2016). this assumes all the 'init' variables in
    # the block above are set to the 'awake' state  
    awake_to_s2_start = 80000
    awake_to_s2_end = 97500
    s2_to_s3_start = 150000
    s2_to_s3_end = 167500
    s3_to_rem_start = 220000
    s3_to_rem_end = 237500
    rem_to_s2_start = 290000
    rem_to_s2_end = 307500
    
else: #if do_sleepstates != True, then just simulate one sleep state
    #these values are used to instantiate the network in the 'connectCells' method
    
    if sleep_state == 0:
        init_GABA_thal = awake_GABA_thal
        init_AMPA_thal = awake_AMPA_thal
        init_AMPA_cort = awake_AMPA_cort
        init_AMPA_pyrpyr = awake_AMPA_pyrpyr
        init_GABA_D2 = awake_GABA_D2
        init_gkl_pyr = gkl_pyr_awake
        init_gkl_inh = gkl_inh_awake
        init_gkl_RE  = gkl_RE_awake
        init_gkl_TC  = gkl_TC_awake
        init_gh_TC   = gh_TC_awake
    elif sleep_state == 1: 
        init_GABA_thal = s2_GABA_thal
        init_AMPA_thal = s2_AMPA_thal
        init_AMPA_cort = s2_AMPA_cort
        init_AMPA_pyrpyr = s2_AMPA_pyrpyr
        init_GABA_D2 = s2_GABA_D2  
        init_gkl_pyr = gkl_pyr_s2
        init_gkl_inh = gkl_inh_s2
        init_gkl_RE  = gkl_RE_s2
        init_gkl_TC  = gkl_TC_s2
        init_gh_TC   = gh_TC_s2
    elif sleep_state == 2: 
        init_GABA_thal = s3_GABA_thal
        init_AMPA_thal = s3_AMPA_thal
        init_AMPA_cort = s3_AMPA_cort
        init_AMPA_pyrpyr = s3_AMPA_pyrpyr
        init_GABA_D2 = s3_GABA_D2
        init_gkl_pyr = gkl_pyr_s3
        init_gkl_inh = gkl_inh_s3
        init_gkl_RE  = gkl_RE_s3
        init_gkl_TC  = gkl_TC_s3
        init_gh_TC   = gh_TC_s3
    elif sleep_state == 3:
        init_GABA_thal = rem_GABA_thal
        init_AMPA_thal = rem_AMPA_thal
        init_AMPA_cort = rem_AMPA_cort
        init_AMPA_pyrpyr = rem_AMPA_pyrpyr
        init_GABA_D2 = rem_GABA_D2
        init_gkl_pyr = gkl_pyr_rem
        init_gkl_inh = gkl_inh_rem
        init_gkl_RE  = gkl_RE_rem
        init_gkl_TC  = gkl_TC_rem
        init_gh_TC   = gh_TC_rem
    elif sleep_state == 5:
        init_GABA_thal = sae_theta_GABA_thal
        init_AMPA_thal = sae_theta_AMPA_thal
        init_AMPA_cort = sae_theta_AMPA_cort
        init_AMPA_pyrpyr = sae_theta_AMPA_pyrpyr
        init_GABA_D2 = sae_theta_GABA_D2
        init_gkl_pyr = sae_theta_gkl_pyr
        init_gkl_inh = sae_theta_gkl_inh
        init_gkl_RE  = sae_theta_gkl_RE
        init_gkl_TC  = sae_theta_gkl_TC
        init_gh_TC   = sae_theta_gh_TC
    elif sleep_state == 6:
        init_GABA_thal = sae_delta_GABA_thal
        init_AMPA_thal = sae_delta_AMPA_thal
        init_AMPA_cort = sae_delta_AMPA_cort
        init_AMPA_pyrpyr = sae_delta_AMPA_pyrpyr
        init_GABA_D2 = sae_delta_GABA_D2
        init_gkl_pyr = sae_delta_gkl_pyr
        init_gkl_inh = sae_delta_gkl_inh
        init_gkl_RE  = sae_delta_gkl_RE
        init_gkl_TC  = sae_delta_gkl_TC
        init_gh_TC   = sae_delta_gh_TC
    elif sleep_state == 7:
        init_GABA_thal = sae_theta_GABA_thal
        init_AMPA_thal = sae_theta_AMPA_thal
        init_AMPA_cort = sae_theta_AMPA_cort
        init_AMPA_pyrpyr = sae_theta_r1_AMPA_pyrpyr
        init_GABA_D2 = sae_theta_GABA_D2
        init_gkl_pyr = sae_theta_gkl_pyr
        init_gkl_inh = sae_theta_gkl_inh
        init_gkl_RE  = sae_theta_gkl_RE
        init_gkl_TC  = sae_theta_gkl_TC
        init_gh_TC   = sae_theta_gh_TC
    elif sleep_state == 8:
        init_GABA_thal = sae_delta_GABA_thal
        init_AMPA_thal = sae_delta_AMPA_thal
        init_AMPA_cort = sae_delta_AMPA_cort
        init_AMPA_pyrpyr = sae_delta_r1_AMPA_pyrpyr
        init_GABA_D2 = sae_delta_GABA_D2
        init_gkl_pyr = sae_delta_gkl_pyr
        init_gkl_inh = sae_delta_gkl_inh
        init_gkl_RE  = sae_delta_gkl_RE
        init_gkl_TC  = sae_delta_gkl_TC
        init_gh_TC   = sae_delta_gh_TC
    elif sleep_state == 9:
        init_GABA_thal = sae_n1_r2_GABA_thal
        init_AMPA_thal = sae_n1_r2_AMPA_thal
        init_AMPA_cort = sae_n1_r2_AMPA_cort
        init_AMPA_pyrpyr = sae_n1_r2_AMPA_pyrpyr
        init_GABA_D2 = sae_n1_r2_GABA_D2
        init_gkl_pyr = sae_n1_r2_gkl_pyr
        init_gkl_inh = sae_n1_r2_gkl_inh
        init_gkl_RE  = sae_n1_r2_gkl_RE
        init_gkl_TC  = sae_n1_r2_gkl_TC
        init_gh_TC   = sae_n1_r2_gh_TC
    elif sleep_state == 10:
        init_GABA_thal = sae_rec_edge_GABA_thal
        init_AMPA_thal = sae_rec_edge_AMPA_thal
        init_AMPA_cort = sae_rec_edge_AMPA_cort
        init_AMPA_pyrpyr = sae_rec_edge_AMPA_pyrpyr
        init_GABA_D2 = sae_rec_edge_GABA_D2
        init_gkl_pyr = sae_rec_edge_gkl_pyr
        init_gkl_inh = sae_rec_edge_gkl_inh
        init_gkl_RE  = sae_rec_edge_gkl_RE
        init_gkl_TC  = sae_rec_edge_gkl_TC
        init_gh_TC   = sae_rec_edge_gh_TC
    elif sleep_state in (11, 13, 14, 15):
        # States 13--15 inherit state 11 exactly, allowing the cortical and
        # thalamic GABA_A decay hypotheses to be tested independently.
        init_GABA_thal = sae_ei_GABA_thal
        init_AMPA_thal = sae_ei_AMPA_thal
        init_AMPA_cort = sae_ei_AMPA_cort
        init_AMPA_pyrpyr = sae_ei_AMPA_pyrpyr
        init_GABA_D2 = sae_ei_GABA_D2
        init_gkl_pyr = sae_ei_gkl_pyr
        init_gkl_inh = sae_ei_gkl_inh
        init_gkl_RE  = sae_ei_gkl_RE
        init_gkl_TC  = sae_ei_gkl_TC
        init_gh_TC   = sae_ei_gh_TC
        if sleep_state in (14, 15):
            init_Beta_GABA_A = Beta_GABA_A_1p5x_decay
        if sleep_state in (13, 15):
            init_Beta_GABA_A_D2 = Beta_GABA_A_D2_1p5x_decay
    elif sleep_state == 12:
        init_GABA_thal = sae_thal_recruit_GABA_thal
        init_AMPA_thal = sae_thal_recruit_AMPA_thal
        init_AMPA_cort = sae_thal_recruit_AMPA_cort
        init_AMPA_pyrpyr = sae_thal_recruit_AMPA_pyrpyr
        init_GABA_D2 = sae_thal_recruit_GABA_D2
        init_gkl_pyr = sae_thal_recruit_gkl_pyr
        init_gkl_inh = sae_thal_recruit_gkl_inh
        init_gkl_RE  = sae_thal_recruit_gkl_RE
        init_gkl_TC  = sae_thal_recruit_gkl_TC
        init_gh_TC   = sae_thal_recruit_gh_TC
    elif 16 <= sleep_state <= 29:
        # States 16--29 share the state-11 cortical E/I and thalamic
        # background. States 19--29 inherit the complete state-17 packet
        # scaffold before applying documented changes.
        init_GABA_thal = sae_ei_GABA_thal
        init_AMPA_thal = sae_ei_AMPA_thal
        init_AMPA_cort = sae_ei_AMPA_cort
        init_GABA_D2 = sae_ei_GABA_D2
        init_gkl_RE = sae_ei_gkl_RE
        init_gkl_TC = sae_ei_gkl_TC
        init_gh_TC = sae_ei_gh_TC

        if sleep_state == 16:
            init_AMPA_pyrpyr = sae_bal50_AMPA_pyrpyr
            init_gkl_pyr = sae_bal50_gkl_pyr
            init_gkl_inh = sae_bal50_gkl_inh
        elif sleep_state == 18:
            init_AMPA_pyrpyr = sae_bal75_AMPA_pyrpyr
            init_gkl_pyr = sae_bal75_gkl_pyr
            init_gkl_inh = sae_bal75_gkl_inh
        else:
            # State 17 itself, plus every state-17-derived probe (19--29).
            init_AMPA_pyrpyr = sae_bal625_AMPA_pyrpyr
            init_gkl_pyr = sae_bal625_gkl_pyr
            init_gkl_inh = sae_bal625_gkl_inh

        if sleep_state == 25:
            init_gh_TC = sae_tc_ih_sleepward
        elif sleep_state == 26:
            init_gh_TC = sae_tc_ih_wakeward
        elif sleep_state == 27:
            init_gkl_TC = sae_tc_gkl_75
        elif sleep_state == 28:
            init_pyr_dend_km = pyr_dend_km_baseline * pyr_intrinsic_62p5_scale
            init_pyr_dend_kca = pyr_dend_kca_baseline * pyr_intrinsic_62p5_scale
        elif sleep_state == 29:
            init_pyr_dend_nap = pyr_dend_nap_baseline * pyr_intrinsic_62p5_scale
            init_pyr_dend_km = pyr_dend_km_baseline * pyr_intrinsic_62p5_scale
            init_pyr_dend_kca = pyr_dend_kca_baseline * pyr_intrinsic_62p5_scale
    elif sleep_state in (30, 31, 32, 33):
        # Stable-background screen on the complete state-11 E/I and thalamic
        # background. States 30--31 change recurrence only; states 32--33
        # pair recurrence with aligned cortical K leak.
        init_GABA_thal = sae_ei_GABA_thal
        init_AMPA_thal = sae_ei_AMPA_thal
        init_AMPA_cort = sae_ei_AMPA_cort
        init_GABA_D2 = sae_ei_GABA_D2
        init_gkl_RE = sae_ei_gkl_RE
        init_gkl_TC = sae_ei_gkl_TC
        init_gh_TC = sae_ei_gh_TC

        if sleep_state == 30:
            init_AMPA_pyrpyr = sae_bg_rec37_AMPA_pyrpyr
            init_gkl_pyr = sae_ei_gkl_pyr
            init_gkl_inh = sae_ei_gkl_inh
        elif sleep_state == 31:
            init_AMPA_pyrpyr = sae_bg_rec44_AMPA_pyrpyr
            init_gkl_pyr = sae_ei_gkl_pyr
            init_gkl_inh = sae_ei_gkl_inh
        elif sleep_state == 32:
            init_AMPA_pyrpyr = sae_bg_bal37_AMPA_pyrpyr
            init_gkl_pyr = sae_bg_bal37_gkl_pyr
            init_gkl_inh = sae_bg_bal37_gkl_inh
        elif sleep_state == 33:
            init_AMPA_pyrpyr = sae_bg_bal44_AMPA_pyrpyr
            init_gkl_pyr = sae_bg_bal44_gkl_pyr
            init_gkl_inh = sae_bg_bal44_gkl_inh
    elif 34 <= sleep_state <= 106:
        # N3-centered screen: initialize every legacy synaptic and membrane
        # control at the calibrated N3 value, then make only the explicitly
        # documented subtraction. Route-specific changes for states 38,
        # 47, and 48 are applied after the legacy thalamic bundles are
        # expanded below.
        init_GABA_thal = s3_GABA_thal
        init_AMPA_thal = s3_AMPA_thal
        init_AMPA_cort = s3_AMPA_cort
        init_AMPA_pyrpyr = s3_AMPA_pyrpyr
        init_GABA_D2 = s3_GABA_D2
        init_gkl_pyr = gkl_pyr_s3
        init_gkl_inh = gkl_inh_s3
        init_gkl_RE = gkl_RE_s3
        init_gkl_TC = gkl_TC_s3
        init_gh_TC = gh_TC_s3

        if sleep_state == 34:
            init_AMPA_pyrpyr = n3_sub_recur_mild
        elif sleep_state == 35:
            init_AMPA_pyrpyr = n3_sub_recur_moderate
        elif sleep_state == 36:
            init_AMPA_pyrpyr = n3_sub_recur_strong
        elif sleep_state == 37:
            init_gkl_pyr = n3_sub_gkl_pyr
            init_gkl_inh = n3_sub_gkl_inh
        elif sleep_state == 39:
            # Figure-5-style necessity test: retain complete N3 cortex and
            # N3 thalamic synapses, changing only the calibrated thalamic
            # intrinsic membrane module from N3 to N2.
            init_gkl_RE = gkl_RE_s2
            init_gkl_TC = gkl_TC_s2
            init_gh_TC = gh_TC_s2
        elif sleep_state == 40:
            init_AMPA_pyrpyr = n3_sub_recur_moderate
            init_gkl_pyr = n3_sub_gkl_pyr
            init_gkl_inh = n3_sub_gkl_inh
        elif sleep_state == 41:
            init_AMPA_pyrpyr = n3_follow_recur_low
        elif sleep_state == 42:
            init_AMPA_pyrpyr = n3_follow_recur_high
        elif sleep_state == 43:
            init_AMPA_pyrpyr = n3_follow_recur_low
            init_gkl_pyr = n3_follow_gkl_pyr
            init_gkl_inh = n3_follow_gkl_inh
        elif sleep_state == 44:
            init_AMPA_pyrpyr = n3_follow_recur_high
            init_gkl_pyr = n3_follow_gkl_pyr
            init_gkl_inh = n3_follow_gkl_inh
        elif sleep_state == 45:
            init_AMPA_pyrpyr = n3_sub_recur_moderate
            init_gkl_RE = gkl_RE_s2
            init_gkl_TC = gkl_TC_s2
            init_gh_TC = gh_TC_s2
        elif sleep_state == 46:
            init_AMPA_pyrpyr = n3_sub_recur_moderate
            init_gkl_pyr = n3_sub_gkl_pyr
            init_gkl_inh = n3_sub_gkl_inh
            init_gkl_RE = gkl_RE_s2
            init_gkl_TC = gkl_TC_s2
            init_gh_TC = gh_TC_s2
        elif sleep_state == 47:
            init_AMPA_pyrpyr = n3_sub_recur_moderate
        elif sleep_state == 48:
            init_AMPA_pyrpyr = n3_sub_recur_moderate
            init_gkl_pyr = n3_sub_gkl_pyr
            init_gkl_inh = n3_sub_gkl_inh
        elif sleep_state == 49:
            init_AMPA_pyrpyr = n3_follow_recur_high
            init_gkl_pyr = n3_leak_mid_gkl_pyr
            init_gkl_inh = n3_leak_mid_gkl_inh
        elif sleep_state == 50:
            init_AMPA_pyrpyr = n3_follow_recur_high
            init_gkl_pyr = n3_follow_gkl_pyr
        elif sleep_state == 51:
            init_AMPA_pyrpyr = n3_follow_recur_high
            init_gkl_inh = n3_follow_gkl_inh
        elif sleep_state == 52:
            # Mild PYR-leak effect alone: completes a 2x2 comparison with
            # states 42 (neither), 53 (INH only), and 49 (both).
            init_AMPA_pyrpyr = n3_follow_recur_high
            init_gkl_pyr = n3_leak_mid_gkl_pyr
        elif sleep_state == 53:
            # Mild INH-leak effect alone.
            init_AMPA_pyrpyr = n3_follow_recur_high
            init_gkl_inh = n3_leak_mid_gkl_inh
        elif sleep_state == 54:
            # Completes the recurrence x leak grid at recurrence 2.10.
            init_AMPA_pyrpyr = n3_follow_recur_low
            init_gkl_pyr = n3_leak_mid_gkl_pyr
            init_gkl_inh = n3_leak_mid_gkl_inh
        elif sleep_state == 55:
            init_AMPA_pyrpyr = n3_recur_mid
            init_gkl_pyr = n3_leak_mid_gkl_pyr
            init_gkl_inh = n3_leak_mid_gkl_inh
        elif sleep_state == 56:
            init_AMPA_pyrpyr = n3_follow_recur_high
            init_gkl_pyr = n3_leak_mid_gkl_pyr
            init_gkl_inh = n3_leak_mid_gkl_inh
            init_gh_TC = n3_gh_TC_mid_to_N2
        elif sleep_state == 57:
            init_AMPA_pyrpyr = n3_follow_recur_high
            init_gkl_pyr = n3_leak_mid_gkl_pyr
            init_gkl_inh = n3_leak_mid_gkl_inh
            init_gkl_TC = n3_gkl_TC_mid_to_N2
        elif sleep_state == 58:
            init_AMPA_pyrpyr = n3_follow_recur_high
            init_gkl_pyr = n3_leak_mid_gkl_pyr
            init_gkl_inh = n3_leak_mid_gkl_inh
            init_gkl_RE = n3_gkl_RE_mid_to_N2
        elif sleep_state == 59:
            # Matched no-leak control for state 55.
            init_AMPA_pyrpyr = n3_recur_mid
        elif 60 <= sleep_state <= 67:
            # State-59 background. Route-specific changes are applied after
            # the legacy AMPA_cort bundle is expanded below.
            init_AMPA_pyrpyr = n3_recur_mid
            if sleep_state == 66:
                init_PYR_PYR_AMPA_D2_random_event_scale = 0.80
            elif sleep_state == 67:
                init_PYR_PYR_AMPA_D2_random_event_scale = 1.20
        elif sleep_state == 68:
            # Full N3 recurrence; TC K leak 10% below N3.
            init_gkl_TC = n3_tc_gkl_low
        elif sleep_state == 69:
            # Full N3 recurrence; TC K leak 10% above N3.
            init_gkl_TC = n3_tc_gkl_high
        elif sleep_state == 70:
            # State-59 recurrence; TC K leak 10% below N3.
            init_AMPA_pyrpyr = n3_recur_mid
            init_gkl_TC = n3_tc_gkl_low
        elif sleep_state == 71:
            # State-59 recurrence; TC K leak 10% above N3.
            init_AMPA_pyrpyr = n3_recur_mid
            init_gkl_TC = n3_tc_gkl_high
        elif sleep_state == 72:
            # Second upward TC K-leak dose with full N3 recurrence.
            init_gkl_TC = n3_tc_gkl_high2
        elif sleep_state == 73:
            # Second upward TC K-leak dose with state-59 recurrence.
            init_AMPA_pyrpyr = n3_recur_mid
            init_gkl_TC = n3_tc_gkl_high2
        elif sleep_state == 74:
            init_gkl_pyr = n3_pyr_gkl_high
        elif sleep_state == 75:
            init_gkl_pyr = n3_pyr_gkl_high2
        elif sleep_state == 76:
            init_gkl_inh = n3_inh_gkl_high
        elif sleep_state == 77:
            init_gkl_pyr = n3_pyr_gkl_high
            init_gkl_inh = n3_inh_gkl_high
        elif sleep_state == 79:
            # State 78 keeps full N3 recurrence; state 79 is its matched
            # low-recurrence counterpart. Route changes occur below.
            init_AMPA_pyrpyr = n3_recur_mid
        elif 80 <= sleep_state <= 87:
            # Route-specific GABA dose series on the state-59 scaffold.
            init_AMPA_pyrpyr = n3_recur_mid
            if sleep_state == 86:
                init_GABA_D2 = n3_cort_gaba_high
            elif sleep_state == 87:
                init_GABA_D2 = n3_cort_gaba_high2
        elif sleep_state == 88:
            # Reduced TC Ih maximal conductance on state 59.
            init_AMPA_pyrpyr = n3_recur_mid
            init_TC_Ih_gbar = TC_Ih_gbar_low
        elif sleep_state == 89:
            # Cytokine-inspired interaction: modestly greater TC K leak plus
            # modestly smaller TC Ih maximal conductance on state 59.
            init_AMPA_pyrpyr = n3_recur_mid
            init_gkl_TC = n3_tc_gkl_high
            init_TC_Ih_gbar = TC_Ih_gbar_low
        elif 90 <= sleep_state <= 106:
            # Causal circuit-dissection states use the exact state-59
            # background. Route lesions and gain changes for states 91--106
            # are applied after the legacy bundled pathway factors are
            # expanded below.
            init_AMPA_pyrpyr = n3_recur_mid
            if sleep_state == 90:
                init_TC_T_gbar = 0.0
    else:
        raise ValueError(
            "Unsupported sleep_state. Implemented states: 0=wake, 1=N2, "
            "2=N3, 3=REM, 5=SAE-theta, 6=SAE-delta, "
            "7=SAE-theta rescue R1, 8=SAE-delta rescue R1, "
            "9=SAE wakeful-slow/N1-like rescue R2, "
            "10=SAE recurrence-edge probe, 11=SAE cortical E/I "
            "rebalance, 12=SAE conservative thalamic recruitment, "
            "13=cortical GABA-A decay, 14=thalamic GABA-A decay, "
            "15=combined GABA-A decay, 16--18=balanced cortical "
            "wake-to-N2 trajectory, 19--24=route-specific thalamic "
            "packet probes, 25--27=TC intrinsic probes, 28--29=PYR "
            "intrinsic-current probes, 30--33=stable-background "
            "cortical probes, 34--40=N3-centered subtraction and "
            "decoupling probes, 41--48=N3-derived fine brackets, "
            "49--59=cortical fine mapping and thalamic decomposition, "
            "60--67=cortical pathway and random-event probes, "
            "68--89=K-leak, thalamocortical-input, GABA and TC-Ih probes, "
            "and 90--106=causal TC-rebound, pathway and loop-gain tests. "
            "State 4 is reserved for the later GRDA-focused hypothesis."
        )

# Expand the two legacy bundled thalamic multipliers into route-specific
# controls. For all legacy states these are exact aliases, so the refactor is
# behavior-preserving. States 19--24 override a single route or a predeclared
# two-route combination. States 38, 47, and 48 selectively weaken fast
# RE->TC GABA-A from its N3 value while preserving slow RE->TC GABA-B.
init_RE_TC_GABA_A = init_GABA_thal
init_RE_TC_GABA_B = init_GABA_thal
init_RE_RE_GABA_A = init_GABA_thal
init_TC_RE_AMPA = init_AMPA_thal
init_PYR_TC_AMPA = init_AMPA_thal
init_PYR_RE_AMPA = init_AMPA_thal

# Expand the remaining AMPA_cort bundle into route-specific controls. These
# aliases preserve every older state exactly unless a state below overrides
# one route explicitly.
init_TC_PYR_AMPA = init_AMPA_cort
init_TC_INH_AMPA = init_AMPA_cort
init_PYR_INH_AMPA_D2 = init_AMPA_cort


if not do_sleepstates:
    if sleep_state in (19, 23, 24):
        init_RE_TC_GABA_B = sae_re_tc_gaba_b_1p20
    if sleep_state in (20, 24):
        init_RE_RE_GABA_A = sae_re_re_gaba_a_0p75
    if sleep_state == 21:
        init_RE_TC_GABA_A = sae_re_tc_gaba_a_0p80
    if sleep_state in (22, 23):
        init_PYR_RE_AMPA = sae_pyr_re_ampa_1p10

    if sleep_state in (38, 47, 48):
        init_RE_TC_GABA_A = n3_sub_RE_TC_GABA_A

    if sleep_state == 60:
        init_TC_PYR_AMPA *= cortical_route_probe_low
    elif sleep_state == 61:
        init_TC_PYR_AMPA *= cortical_route_probe_high
    elif sleep_state == 62:
        init_TC_INH_AMPA *= cortical_route_probe_low
    elif sleep_state == 63:
        init_TC_INH_AMPA *= cortical_route_probe_high
    elif sleep_state == 64:
        init_PYR_INH_AMPA_D2 *= cortical_route_probe_low
    elif sleep_state == 65:
        init_PYR_INH_AMPA_D2 *= cortical_route_probe_high

    # States 78--79 reduce functional thalamocortical AMPA strength while
    # leaving connection counts and the cortical PYR->INH route unchanged.
    if sleep_state in (78, 79):
        init_TC_PYR_AMPA *= n3_tc_input_low
        init_TC_INH_AMPA *= n3_tc_input_low

    # N3-strength route-specific GABA dose series on state 59.
    if sleep_state == 80:
        init_RE_TC_GABA_A = n3_thal_gaba_high
    elif sleep_state == 81:
        init_RE_TC_GABA_A = n3_thal_gaba_high2
    elif sleep_state == 82:
        init_RE_TC_GABA_B = n3_thal_gaba_high
    elif sleep_state == 83:
        init_RE_TC_GABA_B = n3_thal_gaba_high2
    elif sleep_state == 84:
        init_RE_RE_GABA_A = n3_thal_gaba_high
    elif sleep_state == 85:
        init_RE_RE_GABA_A = n3_thal_gaba_high2

    # Causal lesions on the state-59 background. States 91 and 93--97 block
    # ascending thalamic output so descending recruitment and the internal
    # TC-RE feedback loop can be tested without closing the loop through
    # cortex. State 93 blocks both directions and supplies the neither-route
    # control.
    if sleep_state in (91, 93, 94, 95, 96, 97):
        init_TC_PYR_AMPA = 0.0
        init_TC_INH_AMPA = 0.0
    if sleep_state in (92, 93):
        init_PYR_TC_AMPA = 0.0
        init_PYR_RE_AMPA = 0.0
    if sleep_state in (94, 96, 97):
        # Direct cortical input to TC remains; the cortical RE route is off.
        init_PYR_RE_AMPA = 0.0
    elif sleep_state == 95:
        # Cortical input to RE remains; direct cortical input to TC is off.
        init_PYR_TC_AMPA = 0.0

    if sleep_state == 96:
        # Direct PYR->TC input remains, but TC cannot recruit RE.  This tests
        # whether the forward arm of the internal TC->RE->TC reset loop is
        # required for repeated/synchronized low-threshold bursts.
        init_TC_RE_AMPA = 0.0
    elif sleep_state == 97:
        # TC can recruit RE, but neither fast nor slow RE inhibition returns
        # to TC.  Receptor-specific follow-up is deferred until this combined
        # return-arm test shows an effect.
        init_RE_TC_GABA_A = 0.0
        init_RE_TC_GABA_B = 0.0

    # Closed-loop gain screen at the state-59 operating point.  The tonic TC
    # current used with these states is supplied at launch and is deliberately
    # not encoded in the state itself.
    if sleep_state == 98:
        init_PYR_TC_AMPA *= causal_gain_0p75
    elif sleep_state == 99:
        init_PYR_TC_AMPA *= causal_gain_0p50
    elif sleep_state == 100:
        init_PYR_TC_AMPA *= causal_gain_0p25
    elif sleep_state == 101:
        init_TC_PYR_AMPA *= causal_gain_0p75
        init_TC_INH_AMPA *= causal_gain_0p75
    elif sleep_state == 102:
        init_TC_PYR_AMPA *= causal_gain_0p50
        init_TC_INH_AMPA *= causal_gain_0p50
    elif sleep_state == 103:
        init_TC_PYR_AMPA *= causal_gain_0p25
        init_TC_INH_AMPA *= causal_gain_0p25
    elif sleep_state == 104:
        init_PYR_TC_AMPA *= causal_gain_0p50
        init_TC_PYR_AMPA *= causal_gain_0p50
        init_TC_INH_AMPA *= causal_gain_0p50
    elif sleep_state == 105:
        # Exploratory balance: recruit TC bursts somewhat more often while
        # reducing how strongly each burst is expressed in cortex.
        init_PYR_TC_AMPA *= causal_gain_1p25
        init_TC_PYR_AMPA *= causal_gain_0p50
        init_TC_INH_AMPA *= causal_gain_0p50
    elif sleep_state == 106:
        init_PYR_TC_AMPA *= causal_gain_1p25
        init_TC_PYR_AMPA *= causal_gain_0p25
        init_TC_INH_AMPA *= causal_gain_0p25
