"""Optional checks against the actual compiled NEURON cells and synapses.

After ``nrnivmodl mod``:
    FINK_RUN_NEURON_TESTS=1 python -m unittest discover -s tests

Each state runs in its own process, so NEURON cell/global state is isolated.
"""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
PROBE = r'''
import json
from neuron import h
import config
from network_class import Net
h.CVode().use_fast_imem(1)
net = Net(config.Npyr, config.Ninh, config.Nre, config.Ntc)
tc = config.pc.gid2cell(net.tc_gidList[0])
re = config.pc.gid2cell(net.re_gidList[0])
pyr = config.pc.gid2cell(net.pyr_gidList[0])
inh = config.pc.gid2cell(net.inh_gidList[0])
values = {
    "tc_t_gbar": tc.soma.gcabar_it_tc,
    "tc_ih_gbar": tc.soma.ghbar_iar,
    "tc_to_pyr": pyr.synlist[0].gmax,
    "tc_to_inh": inh.synlist[0].gmax,
    "pyr_to_tc": tc.synlist[2].gmax,
    "pyr_to_re": re.synlist[1].gmax,
    "tc_to_re": re.synlist[0].gmax,
    "re_to_tc_fast": tc.synlist[0].gmax,
    "re_to_tc_slow": tc.synlist[1].gmax,
    "clamp_count": sum(hasattr(config.pc.gid2cell(gid), "Iext") for gid in net.tc_gidList),
}
if config.tc_iclamp_enabled:
    config.pc.nthread(2)
    config.pc.set_maxstep(10)
    h.dt = 0.025
    h.finitialize(-68)
    config.pc.psolve(20)
    values["clamp_targets_nA"] = [config.pc.gid2cell(gid).tc_iclamp_target_nA for gid in net.tc_gidList]
    values["clamp_currents_nA"] = [config.pc.gid2cell(gid).Iext.i for gid in net.tc_gidList]
print("PROBE_JSON=" + json.dumps(values))
config.pc.gid_clear()
h.quit()
'''


@unittest.skipUnless(os.environ.get('FINK_RUN_NEURON_TESTS') == '1',
                     'set FINK_RUN_NEURON_TESTS=1 after compiling mechanisms')
class NeuronInterventions(unittest.TestCase):
    def probe(self, state, bias=False):
        env = {k: v for k, v in os.environ.items() if not k.startswith('FINK_')}
        env.update(NEURON_MODULE_OPTIONS='-nogui', FINK_STATE=str(state), FINK_DURATION_MS='100')
        if bias:
            env.update(FINK_TC_ICLAMP_UA_CM2='-0.3', FINK_TC_ICLAMP_RAMP_MS='10', FINK_RECORD_TC_VOLTAGE='1')
        result = subprocess.run([sys.executable, '-c', PROBE], cwd=ROOT, env=env,
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        line = next(line for line in result.stdout.splitlines() if line.startswith('PROBE_JSON='))
        return json.loads(line.removeprefix('PROBE_JSON='))

    def test_baseline_and_t_current_removal(self):
        base = self.probe(59)
        lesion = self.probe(90)
        self.assertEqual(base['tc_t_gbar'], 0.0025)
        self.assertEqual(lesion['tc_t_gbar'], 0)
        self.assertEqual(base['tc_ih_gbar'], lesion['tc_ih_gbar'])
        self.assertEqual(base['clamp_count'], 0)
        for route in ('tc_to_pyr', 'tc_to_inh', 'pyr_to_tc', 'pyr_to_re', 'tc_to_re'):
            self.assertGreater(base[route], 0)
            self.assertEqual(lesion[route], base[route])

    def test_cortex_thalamus_disconnection_preserves_internal_thalamic_loop(self):
        lesion = self.probe(93)
        for route in ('tc_to_pyr', 'tc_to_inh', 'pyr_to_tc', 'pyr_to_re'):
            self.assertEqual(lesion[route], 0)
        for route in ('tc_to_re', 're_to_tc_fast', 're_to_tc_slow'):
            self.assertGreater(lesion[route], 0)

    def test_biased_run_records_complete_traces_with_two_threads(self):
        import numpy as np
        env = {k: v for k, v in os.environ.items() if not k.startswith('FINK_')}
        env.update(NEURON_MODULE_OPTIONS='-nogui', FINK_STATE='99', FINK_SEED='9001',
                   FINK_DURATION_MS='100', FINK_NTHREAD='2', FINK_TC_ICLAMP_UA_CM2='-0.3',
                   FINK_TC_ICLAMP_RAMP_MS='10', FINK_RECORD_TC_VOLTAGE='1')
        with tempfile.TemporaryDirectory(prefix='fink-current-test-') as tmp:
            directory = Path(tmp)
            for path in list(ROOT.glob('*.py')) + list(ROOT.glob('*.hoc')):
                (directory / path.name).symlink_to(path)
            for architecture in ('arm64', 'x86_64'):
                if (ROOT / architecture).exists():
                    (directory / architecture).symlink_to(ROOT / architecture, target_is_directory=True)
            result = subprocess.run([sys.executable, 'bazh_net.py'], cwd=directory, env=env,
                                    capture_output=True, text=True, timeout=60)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            tag = 'state99_seed9001_tcIm0p3_r10ms_nhost=1'
            for prefix in ('lfp', 'vcort'):
                trace = np.loadtxt(directory / f'{prefix}_{tag}.txt')
                self.assertEqual(trace.shape, (4000,))
                self.assertTrue(np.isfinite(trace).all())
            with np.load(directory / f'tc_voltage_{tag}.npz') as diagnostic:
                self.assertEqual(diagnostic['voltage_mv'].shape, (8, 101))
                np.testing.assert_array_equal(diagnostic['time_ms'], np.arange(101))
                np.testing.assert_allclose(diagnostic['tc_iclamp_current_nA'][:, -1], -0.087,
                                           rtol=0, atol=1e-8)
            # The diagnostic name includes the bias, and rerunning is protected.
            repeated = subprocess.run([sys.executable, 'bazh_net.py'], cwd=directory, env=env,
                                      capture_output=True, text=True, timeout=60)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertIn('output already exists', repeated.stderr)

    def test_current_reaches_all_tc_cells_with_two_threads(self):
        values = self.probe(99, bias=True)
        self.assertEqual(values['clamp_count'], 100)
        self.assertEqual(len(values['clamp_targets_nA']), 100)
        for value in values['clamp_targets_nA'] + values['clamp_currents_nA']:
            self.assertAlmostEqual(value, -0.087, places=12)


if __name__ == '__main__':
    unittest.main()
