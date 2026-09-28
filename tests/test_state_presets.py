"""Frozen preset and TC-intervention checks; these tests do not need NEURON.

The fixture was extracted from the historical uncommitted configuration, not
from the restored implementation. Its SHA256 is recorded alongside the values.
Run ``python -m unittest discover -s tests`` from the repository root.
"""

import json
import os
from pathlib import Path
import runpy
import sys
import types
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


class FakeParallelContext:
    def set_maxstep(self, value):
        pass

    def id(self):
        return 0

    def nhost(self):
        return 1


class FakeRandom:
    def Random123_globalindex(self, value):
        pass


def load_config(state, path=None, **env):
    neuron = types.ModuleType('neuron')
    neuron.h = types.SimpleNamespace(ParallelContext=FakeParallelContext, Random=FakeRandom)
    settings = {'FINK_STATE': str(state), **{k: str(v) for k, v in env.items()}}
    with patch.dict(os.environ, settings, clear=True), patch.dict(sys.modules, {'neuron': neuron}):
        return runpy.run_path(str(path or ROOT / 'config.py'))


class StatePresets(unittest.TestCase):
    def test_all_presets_match_historical_configuration(self):
        fixture = json.loads((ROOT / 'tests/fixtures/historical_state_parameters.json').read_text())
        for state, expected in fixture['states'].items():
            with self.subTest(state=state):
                config = load_config(state)
                actual = {k: v for k, v in config.items() if k.startswith('init_')}
                self.assertEqual(actual, expected)

    def test_reserved_and_unrestored_states_are_rejected(self):
        for state in (-1, 4, 107, 119):
            with self.subTest(state=state), self.assertRaises(ValueError):
                load_config(state)

    def test_current_is_explicit_and_does_not_change_presets(self):
        ordinary = load_config(99)
        biased = load_config(99, FINK_TC_ICLAMP_UA_CM2=-0.30, FINK_RECORD_TC_VOLTAGE=1)
        self.assertFalse(ordinary['tc_iclamp_enabled'])
        self.assertEqual(ordinary['tc_iclamp_suffix'], '')
        self.assertEqual(biased['tc_iclamp_suffix'], '_tcIm0p3_r5000ms')
        self.assertAlmostEqual(biased['tc_iclamp_target_nA'], -0.087)
        self.assertEqual(biased['tc_voltage_gids'], (700, 712, 725, 737, 750, 762, 775, 787))
        self.assertEqual({k: v for k, v in ordinary.items() if k.startswith('init_')},
                         {k: v for k, v in biased.items() if k.startswith('init_')})

    def test_invalid_current_settings_are_rejected(self):
        for env in (
            {'FINK_TC_ICLAMP_UA_CM2': 'nan'},
            {'FINK_TC_ICLAMP_UA_CM2': -1.1},
            {'FINK_TC_ICLAMP_UA_CM2': 0.1},
            {'FINK_TC_ICLAMP_UA_CM2': -0.3},  # diagnostics required
            {'FINK_TC_ICLAMP_RAMP_MS': -1},
            {'FINK_TC_ICLAMP_UA_CM2': -0.3, 'FINK_RECORD_TC_VOLTAGE': 1, 'FINK_DURATION_MS': 100},
        ):
            with self.subTest(env=env), self.assertRaises(ValueError):
                load_config(99, **env)

    def test_causal_lesions(self):
        base = load_config(59)
        tc_t_off = load_config(90)
        self.assertEqual(tc_t_off['init_TC_T_gbar'], 0)
        self.assertEqual(base['init_TC_T_gbar'], 0.0025)
        self.assertEqual(tc_t_off['init_TC_Ih_gbar'], base['init_TC_Ih_gbar'])
        isolated = load_config(93)
        for route in ('init_TC_PYR_AMPA', 'init_TC_INH_AMPA', 'init_PYR_TC_AMPA', 'init_PYR_RE_AMPA'):
            self.assertEqual(isolated[route], 0)
            self.assertGreater(base[route], 0)
        # Removing cortex--thalamus drive retains the internal TC--RE loop.
        self.assertEqual(isolated['init_TC_RE_AMPA'], base['init_TC_RE_AMPA'])
        self.assertEqual(isolated['init_RE_TC_GABA_A'], base['init_RE_TC_GABA_A'])

    def test_main_runtime_controls_are_preserved(self):
        config = load_config(99, FINK_NTHREAD=3, FINK_VERBOSE=1, FINK_SLEEPSTATES=1)
        self.assertEqual(config['nthread'], 3)
        self.assertTrue(config['verbose'])
        self.assertTrue(config['do_sleepstates'])
        self.assertEqual(config['init_AMPA_pyrpyr'], config['awake_AMPA_pyrpyr'])


if __name__ == '__main__':
    unittest.main()
