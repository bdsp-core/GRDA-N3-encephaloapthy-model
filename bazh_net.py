"""
Main simulation file for replicating 
"Cellular and neurochemical basis of sleep stages in the thalamocortical network," eLife, 2016.
"""

import config
from datetime import datetime
from neuron import h
import numpy as np
import os
runtag = (
    "sweep_seed%d" % config.randSeed
    if config.do_sleepstates
    else "state%d_seed%d" % (config.sleep_state, config.randSeed)
) + config.tc_iclamp_suffix

def preflight_tc_current_outputs():
    """Refuse to overwrite an earlier calibration run by accident."""
    if not config.tc_iclamp_enabled or config.tc_iclamp_allow_overwrite:
        return
    output_names = [
        "%s_%s_nhost=%g.txt" % (prefix, runtag, config.nhost)
        for prefix in ("raster", "lfp", "vcort")
    ]
    output_names.append(
        "tc_voltage_%s_nhost=%g.npz" % (runtag, config.nhost)
    )
    existing = [name for name in output_names if os.path.exists(name)]
    if existing:
        raise FileExistsError(
            "TC-current calibration output already exists; choose a new "
            "state/seed/current or set FINK_TC_ICLAMP_ALLOW_OVERWRITE=1: %s"
            % ", ".join(existing)
        )


preflight_tc_current_outputs()
h('load_file("stdgui.hoc")') # need this instead of import gui to get the simulation to be reproducible and not give an LFP flatline
from network_class import Net
from lfp_recorder import CorticalFieldRecorder

h('CVode[0].use_fast_imem(1)') #see use_fast_imem() at https://www.neuron.yale.edu/neuron/static/new_doc/simctrl/cvode.html  

def setup_tc_voltage_recording(net):
    """Create low-rate, local diagnostic vectors for selected TC cells."""
    if not config.record_tc_voltage:
        return None, {}, {}, {}, {}

    time_vec = None
    if config.idhost == 0:
        time_vec = h.Vector()
        time_vec.record(net.cells[0].synlist[0], h._ref_t, config.tc_voltage_dt_ms)

    voltage_vecs = {}
    t_h_gate_vecs = {}
    t_current_vecs = {}
    iclamp_current_vecs = {}
    selected = set(config.tc_voltage_gids)
    for gid in net.tc_gidList:
        if gid not in selected:
            continue
        cell = config.pc.gid2cell(gid)
        segment = cell.soma(0.5)
        owner = cell.synlist[0]  # associate scheduled recording with this cell's NEURON thread
        voltage_vecs[gid] = h.Vector()
        voltage_vecs[gid].record(owner, segment._ref_v, config.tc_voltage_dt_ms)
        # h_it_tc is the T-current availability/de-inactivation gate: larger
        # values mean more rebound current is available after hyperpolarizing.
        t_h_gate_vecs[gid] = h.Vector()
        t_h_gate_vecs[gid].record(
            owner, segment._ref_h_it_tc, config.tc_voltage_dt_ms
        )
        # it_tc is the only calcium-current mechanism in a TC cell, so the
        # segment's calcium current reports the actual inward T-current.
        t_current_vecs[gid] = h.Vector()
        t_current_vecs[gid].record(owner, segment._ref_ica, config.tc_voltage_dt_ms)
        if config.tc_iclamp_enabled:
            cell = config.pc.gid2cell(gid)
            iclamp_current_vecs[gid] = h.Vector()
            iclamp_current_vecs[gid].record(
                cell.Iext, cell.Iext._ref_i, config.tc_voltage_dt_ms
            )
    return (
        time_vec,
        voltage_vecs,
        t_h_gate_vecs,
        t_current_vecs,
        iclamp_current_vecs,
    )


def save_tc_voltage_recording(
    time_vec,
    voltage_vecs,
    t_h_gate_vecs,
    t_current_vecs,
    iclamp_current_vecs,
):
    """Gather selected TC traces on rank 0 and save one compressed NPZ."""
    payload = [None] * config.nhost
    payload[0] = {
        "voltage": {
            int(gid): np.asarray(vec, dtype=np.float32)
            for gid, vec in voltage_vecs.items()
        },
        "t_h_gate": {
            int(gid): np.asarray(vec, dtype=np.float32)
            for gid, vec in t_h_gate_vecs.items()
        },
        "t_current": {
            int(gid): np.asarray(vec, dtype=np.float32)
            for gid, vec in t_current_vecs.items()
        },
        "iclamp_current": {
            int(gid): np.asarray(vec, dtype=np.float32)
            for gid, vec in iclamp_current_vecs.items()
        },
    }
    config.pc.barrier()
    gathered = config.pc.py_alltoall(payload)
    config.pc.barrier()

    if config.idhost != 0:
        return

    voltage_by_gid = {}
    t_h_gate_by_gid = {}
    t_current_by_gid = {}
    iclamp_current_by_gid = {}
    for rank_payload in gathered:
        voltage_by_gid.update(rank_payload["voltage"])
        t_h_gate_by_gid.update(rank_payload["t_h_gate"])
        t_current_by_gid.update(rank_payload["t_current"])
        iclamp_current_by_gid.update(rank_payload["iclamp_current"])

    gids = np.asarray(sorted(voltage_by_gid), dtype=np.int32)
    expected_gids = np.asarray(sorted(config.tc_voltage_gids), dtype=np.int32)
    if not np.array_equal(gids, expected_gids):
        raise RuntimeError(
            "TC voltage recording gathered gids %s; expected %s"
            % (gids.tolist(), expected_gids.tolist())
        )

    time_ms = np.asarray(time_vec, dtype=np.float32)
    voltage_mv = np.vstack([voltage_by_gid[int(gid)] for gid in gids])
    t_h_gate = np.vstack(
        [t_h_gate_by_gid[int(gid)] for gid in gids]
    )
    t_current = np.vstack(
        [t_current_by_gid[int(gid)] for gid in gids]
    )
    if voltage_mv.shape[1] != len(time_ms):
        raise RuntimeError(
            "TC voltage/time length mismatch: %d versus %d"
            % (voltage_mv.shape[1], len(time_ms))
        )
    if t_h_gate.shape != voltage_mv.shape:
        raise RuntimeError("TC T-current gate and voltage shapes do not match")
    if t_current.shape != voltage_mv.shape:
        raise RuntimeError("TC T-current and voltage shapes do not match")

    output_data = dict(
        time_ms=time_ms,
        gids=gids,
        voltage_mv=voltage_mv,
        t_h_gate=t_h_gate,
        t_current_mA_cm2=t_current,
        sample_dt_ms=np.asarray(config.tc_voltage_dt_ms, dtype=np.float32),
        tc_t_gbar_s_cm2=np.asarray(
            config.init_TC_T_gbar, dtype=np.float32
        ),
    )
    if config.tc_iclamp_enabled:
        iclamp_current_nA = np.vstack(
            [iclamp_current_by_gid[int(gid)] for gid in gids]
        )
        if iclamp_current_nA.shape != voltage_mv.shape:
            raise RuntimeError(
                "TC clamp-current and voltage shapes do not match"
            )
        output_data.update(
            tc_iclamp_current_nA=iclamp_current_nA,
            tc_iclamp_uA_cm2=np.asarray(
                config.tc_iclamp_uA_cm2, dtype=np.float64
            ),
            tc_iclamp_target_nA=np.asarray(
                config.tc_iclamp_target_nA, dtype=np.float64
            ),
            tc_iclamp_ramp_ms=np.asarray(
                config.tc_iclamp_ramp_ms, dtype=np.float64
            ),
        )

    output = "tc_voltage_%s_nhost=%g.npz" % (runtag, config.nhost)
    np.savez_compressed(output, **output_data)
    print(
        "Saved TC voltage diagnostic %s: %d cells x %d samples"
        % (output, len(gids), len(time_ms))
    )

def write_trace(stream, values):
    '''write one value per line in the legacy "%.3f " format, formatting a whole chunk in one operation'''
    if len(values):
        stream.write(("%.3f \n" * len(values)) % tuple(values.tolist()))

def write_raster(stream, times, ids):
    '''write "time  gid" lines in the legacy format, formatting a whole chunk in one operation'''
    if len(times):
        stream.write(("%.3f  %g\n" * len(times)) % tuple(np.column_stack((times, ids)).ravel().tolist()))
  
def onerun(randSeed,Npyr,Ninh,Nre,Ntc):
       
    h.Random().Random123_globalindex(randSeed) #this changes ALL Random123 streams

    # These GABA_A kinetic parameters are GLOBAL within their mechanisms.
    # Set them on every MPI rank before cells/synapses are created and before
    # finitialize, so kinetic probes do not alter conductance strength.
    h.Beta_GABA_A = config.init_Beta_GABA_A
    h.Beta_GABA_A_D2 = config.init_Beta_GABA_A_D2
    
    # create network
    net = Net(config.Npyr,config.Ninh,config.Nre,config.Ntc)
    if config.tc_iclamp_enabled:
        local_tc_clamps = 0
        for gid in net.tc_gidList:
            cell = config.pc.gid2cell(gid)
            if not hasattr(cell, "Iext"):
                raise RuntimeError("TC gid %d has no current clamp" % gid)
            expected_nA = (
                config.tc_iclamp_uA_cm2
                * cell.tc_iclamp_area_um2
                * 1e-5
            )
            if not np.isclose(
                cell.tc_iclamp_target_nA,
                expected_nA,
                rtol=0.0,
                atol=1e-15,
            ):
                raise RuntimeError(
                    "TC gid %d current-density conversion mismatch" % gid
                )
            local_tc_clamps += 1
        total_tc_clamps = int(config.pc.allreduce(local_tc_clamps, 1))
        if total_tc_clamps != config.Ntc:
            raise RuntimeError(
                "TC current probe attached to %d cells; expected %d"
                % (total_tc_clamps, config.Ntc)
            )
    else:
        total_tc_clamps = 0
    (
        tc_voltage_time_vec,
        tc_voltage_vecs,
        tc_t_h_gate_vecs,
        tc_t_current_vecs,
        tc_iclamp_current_vecs,
    ) = setup_tc_voltage_recording(net)

    
    #build list of sections in cortical cells, because only these are the ones that will generate the LFP
    cort_secs = []
    for pyr_gid in net.pyr_gidList:
        secs_to_add = config.pc.gid2cell(pyr_gid).all
        for sec in secs_to_add:
            cort_secs.append(sec)
        
    for inh_gid in net.inh_gidList:
        secs_to_add = config.pc.gid2cell(inh_gid).all
        for sec in secs_to_add:
            cort_secs.append(sec)
         
    if config.doextra:
        recorder = CorticalFieldRecorder(cort_secs) #records summed cortical voltage and LFP on every time step (see lfp_recorder.py)
    
    '''if do_sleepstates==True, specify how parameters should change to induce different sleep states (see lines 640-777 of C++ main.cpp)'''
    if config.do_sleepstates:
        # see https://www.neuron.yale.edu/neuron/static/py_doc/programming/math/vector.html?highlight=vector#Vector.play
        t=h.Vector([config.awake_to_s2_start,config.awake_to_s2_end,config.s2_to_s3_start,config.s2_to_s3_end,config.s3_to_rem_start,config.s3_to_rem_end,config.rem_to_s2_start,config.rem_to_s2_end,config.rem_to_s2_end+h.dt]) #last entry (with "+h.dt") ensures that last parameter values remain constant to end of simulation ("If a constant outside the range is desired, make sure the last two points have the same y value and have different t values")
               
        #the next three entries are for changing GABA_A and GABA_B strengths in thalamus
        RE_TC_GABA_A_vec=h.Vector(config.re2tc_gaba_a_str * np.array([config.awake_GABA_thal,config.s2_GABA_thal,config.s2_GABA_thal,config.s3_GABA_thal,config.s3_GABA_thal,config.rem_GABA_thal,config.rem_GABA_thal,config.s2_GABA_thal,config.s2_GABA_thal]))
        RE_TC_A_veclist = []
        for tc_gid in net.tc_gidList:
            norm_vec = RE_TC_GABA_A_vec / config.pc.gid2cell(tc_gid).k_RE_TC_GABA_A
            RE_TC_A_veclist.append(norm_vec)
            RE_TC_A_veclist[-1].play(config.pc.gid2cell(tc_gid).synlist[0]._ref_gmax, t, True)
        
        RE_TC_GABA_B_vec=h.Vector(config.re2tc_gaba_b_str * np.array([config.awake_GABA_thal,config.s2_GABA_thal,config.s2_GABA_thal,config.s3_GABA_thal,config.s3_GABA_thal,config.rem_GABA_thal,config.rem_GABA_thal,config.s2_GABA_thal,config.s2_GABA_thal]))
        RE_TC_B_veclist = []
        for tc_gid in net.tc_gidList:
            norm_vec = RE_TC_GABA_B_vec / config.pc.gid2cell(tc_gid).k_RE_TC_GABA_B
            RE_TC_B_veclist.append(norm_vec)
            RE_TC_B_veclist[-1].play(config.pc.gid2cell(tc_gid).synlist[1]._ref_gmax, t, True)
        
        RE_RE_GABA_A_vec=h.Vector(config.re2re_gaba_a_str * np.array([config.awake_GABA_thal,config.s2_GABA_thal,config.s2_GABA_thal,config.s3_GABA_thal,config.s3_GABA_thal,config.rem_GABA_thal,config.rem_GABA_thal,config.s2_GABA_thal,config.s2_GABA_thal]))
        RE_RE_veclist = []
        for re_gid in net.re_gidList:
            norm_vec = RE_RE_GABA_A_vec / config.pc.gid2cell(re_gid).k_RE_RE
            RE_RE_veclist.append(norm_vec)
            RE_RE_veclist[-1].play(config.pc.gid2cell(re_gid).synlist[2]._ref_gmax, t, True)
            
        #next three entries address changing AMPA strengths for connections terminating in thalamus
        TC_RE_AMPA_vec=h.Vector(config.tc2re_ampa_str * np.array([config.awake_AMPA_thal,config.s2_AMPA_thal,config.s2_AMPA_thal,config.s3_AMPA_thal,config.s3_AMPA_thal,config.rem_AMPA_thal,config.rem_AMPA_thal,config.s2_AMPA_thal,config.s2_AMPA_thal]))
        TC_RE_veclist =[]
        for re_gid in net.re_gidList:
            norm_vec = TC_RE_AMPA_vec / config.pc.gid2cell(re_gid).k_TC_RE
            TC_RE_veclist.append(norm_vec)
            TC_RE_veclist[-1].play(config.pc.gid2cell(re_gid).synlist[0]._ref_gmax, t, True)
        
        PYR_TC_AMPA_vec=h.Vector(config.pyr2tc_ampa_str * np.array([config.awake_AMPA_thal,config.s2_AMPA_thal,config.s2_AMPA_thal,config.s3_AMPA_thal,config.s3_AMPA_thal,config.rem_AMPA_thal,config.rem_AMPA_thal,config.s2_AMPA_thal,config.s2_AMPA_thal]))
        PYR_TC_veclist = []
        for tc_gid in net.tc_gidList:
            norm_vec = PYR_TC_AMPA_vec / config.pc.gid2cell(tc_gid).k_PY_TC
            PYR_TC_veclist.append(norm_vec)
            PYR_TC_veclist[-1].play(config.pc.gid2cell(tc_gid).synlist[2]._ref_gmax, t, True)
        
        PYR_RE_AMPA_vec=h.Vector(config.pyr2re_ampa_str * np.array([config.awake_AMPA_thal,config.s2_AMPA_thal,config.s2_AMPA_thal,config.s3_AMPA_thal,config.s3_AMPA_thal,config.rem_AMPA_thal,config.rem_AMPA_thal,config.s2_AMPA_thal,config.s2_AMPA_thal]))
        PYR_RE_veclist = []
        for re_gid in net.re_gidList:
            norm_vec = PYR_RE_AMPA_vec / config.pc.gid2cell(re_gid).k_PY_RE
            PYR_RE_veclist.append(norm_vec)
            PYR_RE_veclist[-1].play(config.pc.gid2cell(re_gid).synlist[1]._ref_gmax, t, True)
            
        #next three entries address all AMPA connections termining in cortex, other than PYR->PYR connections
        TC_PYR_AMPA_D2_vec=h.Vector(config.tc2pyr_ampa_str * np.array([config.awake_AMPA_cort,config.s2_AMPA_cort,config.s2_AMPA_cort,config.s3_AMPA_cort,config.s3_AMPA_cort,config.rem_AMPA_cort,config.rem_AMPA_cort,config.s2_AMPA_cort,config.s2_AMPA_cort]))
        TC_PYR_veclist = [] #pretty sure I need to make a list of all the play vectors, bc. according to documentation, "The system maintains a set of play vectors and the vector will be removed from the list if the vector or var is destroyed." 
        for pyr_gid in net.pyr_gidList:
            norm_vec =  TC_PYR_AMPA_D2_vec / config.pc.gid2cell(pyr_gid).k_TC_PY  #need to noramlize each connection by in-degree of post-synaptic cell
            TC_PYR_veclist.append(norm_vec)
            TC_PYR_veclist[-1].play( config.pc.gid2cell(pyr_gid).synlist[0]._ref_gmax, t, True)
            
        TC_INH_AMPA_D2_vec=h.Vector(config.tc2inh_ampa_str * np.array([config.awake_AMPA_cort,config.s2_AMPA_cort,config.s2_AMPA_cort,config.s3_AMPA_cort,config.s3_AMPA_cort,config.rem_AMPA_cort,config.rem_AMPA_cort,config.s2_AMPA_cort,config.s2_AMPA_cort]))
        TC_INH_veclist = [] #pretty sure I need to make a list of all the play vectors, bc. according to documentation, "The system maintains a set of play vectors and the vector will be removed from the list if the vector or var is destroyed." 
        for inh_gid in net.inh_gidList:
            norm_vec = TC_INH_AMPA_D2_vec / config.pc.gid2cell(inh_gid).k_TC_IN
            TC_INH_veclist.append(norm_vec)
            TC_INH_veclist[-1].play(config.pc.gid2cell(inh_gid).synlist[0]._ref_gmax, t, True)
            
        PYR_INH_AMPA_D2_vec=h.Vector(config.pyr2inh_ampa_d2_str * np.array([config.awake_AMPA_cort,config.s2_AMPA_cort,config.s2_AMPA_cort,config.s3_AMPA_cort,config.s3_AMPA_cort,config.rem_AMPA_cort,config.rem_AMPA_cort,config.s2_AMPA_cort,config.s2_AMPA_cort]))
        PYR_INH_veclist = []
        for inh_gid in net.inh_gidList:
            norm_vec = PYR_INH_AMPA_D2_vec / config.pc.gid2cell(inh_gid).k_PY_IN_AMPA
            PYR_INH_veclist.append(norm_vec)
            PYR_INH_veclist[-1].play(config.pc.gid2cell(inh_gid).synlist[1]._ref_gmax, t, True)
        
        #next one addresses scaling of PYR->PYR connections only
        PYR_PYR_AMPA_D2_vec=h.Vector( config.pyr2pyr_ampa_d2_str * np.array([config.awake_AMPA_pyrpyr,config.s2_AMPA_pyrpyr,config.s2_AMPA_pyrpyr,config.s3_AMPA_pyrpyr,config.s3_AMPA_pyrpyr,config.rem_AMPA_pyrpyr,config.rem_AMPA_pyrpyr,config.s2_AMPA_pyrpyr,config.s2_AMPA_pyrpyr]))
        PYR_PYR_veclist = []
        for pyr_gid in net.pyr_gidList:
            norm_vec = PYR_PYR_AMPA_D2_vec / config.pc.gid2cell(pyr_gid).k_PY_PY_AMPA
            PYR_PYR_veclist.append(norm_vec)
            PYR_PYR_veclist[-1].play(config.pc.gid2cell(pyr_gid).synlist[1]._ref_gmax, t, True)
        
        #this last entry is for changing GABA_A_D2 strength for INH->PYR connections alone
        INH_PYR_GABA_D2_vec=h.Vector(config.inh2pyr_gaba_a_d2_str * np.array([config.awake_GABA_D2,config.s2_GABA_D2,config.s2_GABA_D2,config.s3_GABA_D2,config.s3_GABA_D2,config.rem_GABA_D2,config.rem_GABA_D2,config.s2_GABA_D2,config.s2_GABA_D2]))
        INH_PYR_veclist = []
        for pyr_gid in net.pyr_gidList:
            norm_vec = INH_PYR_GABA_D2_vec / config.pc.gid2cell(pyr_gid).k_IN_PY
            INH_PYR_veclist.append(norm_vec)
            INH_PYR_veclist[-1].play(config.pc.gid2cell(pyr_gid).synlist[3]._ref_gmax, t, True)
        
        #cellular properties
        gkl_pyr_vec=h.Vector( np.array([config.gkl_pyr_awake,config.gkl_pyr_s2,config.gkl_pyr_s2,config.gkl_pyr_s3,config.gkl_pyr_s3,config.gkl_pyr_rem,config.gkl_pyr_rem,config.gkl_pyr_s2,config.gkl_pyr_s2]))
        for pyr_gid in net.pyr_gidList:
            gkl_pyr_vec.play(config.pc.gid2cell(pyr_gid).dend(0.5)._ref_gkL_kL, t, True) #'True' makes it so values are linearly interpolated between points specified in vectors
            
        gkl_inh_vec=h.Vector( np.array([config.gkl_inh_awake,config.gkl_inh_s2,config.gkl_inh_s2,config.gkl_inh_s3,config.gkl_inh_s3,config.gkl_inh_rem,config.gkl_inh_rem,config.gkl_inh_s2,config.gkl_inh_s2]))
        for inh_gid in net.inh_gidList:
            gkl_inh_vec.play(config.pc.gid2cell(inh_gid).dend(0.5)._ref_gkL_kL, t, True)
        
        gkl_TC_vec=h.Vector( np.array([config.gkl_TC_awake,config.gkl_TC_s2,config.gkl_TC_s2,config.gkl_TC_s3,config.gkl_TC_s3,config.gkl_TC_rem,config.gkl_TC_rem,config.gkl_TC_s2,config.gkl_TC_s2]))
        for tc_gid in net.tc_gidList:
            gkl_TC_vec.play(config.pc.gid2cell(tc_gid).soma(0.5)._ref_gkL_kL, t, True)
        
        gkl_RE_vec=h.Vector( np.array([config.gkl_RE_awake,config.gkl_RE_s2,config.gkl_RE_s2,config.gkl_RE_s3,config.gkl_RE_s3,config.gkl_RE_rem,config.gkl_RE_rem,config.gkl_RE_s2,config.gkl_RE_s2]))
        for re_gid in net.re_gidList:
            gkl_RE_vec.play(config.pc.gid2cell(re_gid).soma(0.5)._ref_gkL_kL, t, True)
        
        gh_TC_vec=h.Vector([config.gh_TC_awake,config.gh_TC_s2,config.gh_TC_s2,config.gh_TC_s3,config.gh_TC_s3,config.gh_TC_rem,config.gh_TC_rem,config.gh_TC_s2,config.gh_TC_s2])
        for tc_gid in net.tc_gidList:
            gh_TC_vec.play(config.pc.gid2cell(tc_gid).soma(0.5)._ref_fac_gh_TC_iar, t, True)
    
    '''set up custom initialization to set RE voltage values'''
    def set_RE_voltages():
        #loop through all RE cells and set their initial voltages to -65 mV
        for re_gid in net.re_gidList:
            config.pc.gid2cell(re_gid).soma.v=-65 
        if(h.cvode.active()):
            h.cvode.re_init()
        else:
            h.fcurrent()
        h.frecord_init() 
       
    # run sim and gather spikes
    if config.nthread > 1:
        config.pc.nthread(config.nthread) #split this rank's cells among NEURON threads
    config.pc.set_maxstep(10) #see https://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/network/parcon.html#ParallelContext.set_maxstep, as well as section 2.4 of the Lytton/Salvador paper
    h.dt = 0.025
    
    fih = h.FInitializeHandler(set_RE_voltages)
    h.finitialize(-68) #set initial voltages of all cells *except RE cells* to -68 mV
    #h.stdinit()
    
    if config.idhost==0: 
        print(
            "Configuration: state=%d seed=%d duration_ms=%g nhost=%d nthread=%d "
            "Beta_GABA_A=%g Beta_GABA_A_D2=%g"
            % (
                config.sleep_state,
                config.randSeed,
                config.duration,
                config.nhost,
                config.nthread,
                config.init_Beta_GABA_A,
                config.init_Beta_GABA_A_D2,
            )
        )
        print(
            "Thalamic pathway factors: RE_TC_GABA_A=%g RE_TC_GABA_B=%g "
            "RE_RE_GABA_A=%g TC_RE_AMPA=%g PYR_TC_AMPA=%g "
            "PYR_RE_AMPA=%g"
            % (
                config.init_RE_TC_GABA_A,
                config.init_RE_TC_GABA_B,
                config.init_RE_RE_GABA_A,
                config.init_TC_RE_AMPA,
                config.init_PYR_TC_AMPA,
                config.init_PYR_RE_AMPA,
            )
        )
        print(
            "Cortical synaptic factors: AMPA_cort=%g AMPA_pyrpyr=%g "
            "GABA_D2=%g"
            % (
                config.init_AMPA_cort,
                config.init_AMPA_pyrpyr,
                config.init_GABA_D2,
            )
        )
        print(
            "Membrane factors: gkl_pyr=%g gkl_inh=%g gkl_RE=%g "
            "gkl_TC=%g gh_TC=%g PYR_dend_NaP=%g PYR_dend_Km=%g "
            "PYR_dend_KCa=%g"
            % (
                config.init_gkl_pyr,
                config.init_gkl_inh,
                config.init_gkl_RE,
                config.init_gkl_TC,
                config.init_gh_TC,
                config.init_pyr_dend_nap,
                config.init_pyr_dend_km,
                config.init_pyr_dend_kca,
            )
        )
        if config.tc_iclamp_enabled:
            print(
                "TC holding-current probe: density_uA_cm2=%g "
                "target_nA_per_cell=%g ramp_ms=%g cells=%d runtag=%s"
                % (
                    config.tc_iclamp_uA_cm2,
                    config.tc_iclamp_target_nA,
                    config.tc_iclamp_ramp_ms,
                    total_tc_clamps,
                    runtag,
                )
            )
        if config.record_tc_voltage:
            print(
                "TC voltage diagnostic: dt_ms=%g gids=%s"
                % (config.tc_voltage_dt_ms, list(config.tc_voltage_gids))
            )
        print('Running sim...')
        startTime = datetime.now() # store sim start time
        
        raster_file = open("raster_%s_nhost=%g.txt"%(runtag, config.nhost), 'w') # prepare file to print raster data to file
        lfp_file = open("lfp_%s_nhost=%g.txt"%(runtag, config.nhost),'w') #prepare to print biophysical lfp data to file
        vcort_file = open("vcort_%s_nhost=%g.txt"%(runtag, config.nhost), 'w') #prepare to print summed intracellular voltage trace data to file
    

    #actually run the simulation on all nodes
    t_curr = 0
    while (t_curr < config.duration-h.dt): # include the '-dt' to account for rounding error; otherwise, may get error in writeVoltages
        #step forward in periods of 'config.t_seg', and dump data to file after each step forward. Then resize vectors, so program does not run out of memory
        if(t_curr + config.t_seg < config.duration):
            config.pc.psolve(t_curr+config.t_seg)  
            if config.idhost==0: print("Numerically integrated through %.2f ms"%(t_curr+config.t_seg))
        else:
            config.pc.psolve(config.duration)
            
        if config.doextra:
            v_part, lfp_part = recorder.collect() #this host's summed cortical voltage and LFP for every step of this chunk
        else:
            v_part = lfp_part = None
        net.gatherChunk(v_part, lfp_part)  # gather spikes (and LFP data) from all nodes onto master node
        
        if(config.idhost==0):
            write_raster(raster_file, net.tVecAll, net.idVecAll) # use the bash command 'sort -k 1n,1n -k 2n,2n raster_nhost=4 > raster_nhost=4_sorted' to sort the raster plots when nhost>1
            if config.doextra:
                write_trace(vcort_file, net.v_sum) #print cortical voltage data to file
                write_trace(lfp_file, net.lfp_sum) #print LFP data to file
                    
        net.tVec.resize(0) #reset the raster vectors on each individual node, so they do not grow too large as the simulation progresses
        net.idVec.resize(0)
            
        t_curr = t_curr + config.t_seg
    
    if config.idhost==0: 
        raster_file.close() #close raster file
        lfp_file.close()
        vcort_file.close()
        
        runTime = (datetime.now() - startTime).total_seconds()  # calculate run time
        print("Run time for %d sec sim = %.2f sec"%(int(config.duration/1000.0), runTime) )
        
    
    
    # plot net raster, save net data and plot cell 0 traces 
    '''net.gatherSpikes()  # gather spikes from all nodes onto master node
    if config.idhost==0: #if statement because we don't want every host to plot and save data (that would be redundant)
        #to plot raster data, we need to load the file, then define tVecAll and idVecAll lists so that they exist when the plotRaster method is called
        rasterdat=np.loadtxt("raster_%s_nhost=%g.txt"%(runtag, config.nhost))
        if(len(rasterdat)>0):
            net.tVecAll=list(rasterdat[:,0])
            net.idVecAll=list(rasterdat[:,1])
            #net.plotRaster()
        else:
            print("No cells spiked, so there is no raster plot to display.")
        
        net.saveData()
        #net.cells[plot_cell].plotTraces()'''
        
    if config.record_tc_voltage:
        save_tc_voltage_recording(
            tc_voltage_time_vec,
            tc_voltage_vecs,
            tc_t_h_gate_vecs,
            tc_t_current_vecs,
            tc_iclamp_current_vecs,
        )

    del net

onerun(config.randSeed,config.Npyr,config.Ninh,config.Nre,config.Ntc)

config.pc.barrier()
h.quit()
