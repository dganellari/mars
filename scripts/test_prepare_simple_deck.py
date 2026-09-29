#!/usr/bin/env python3
"""Exercise the deck bridge using only the synthetic public channel."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import prepare_simple_deck as bridge

ROOT = Path(__file__).resolve().parents[1]
PUBLIC = ROOT / 'tests/data/public_openaccel_reference/channel_smoke.i'


class DeckTests(unittest.TestCase):
    def setUp(self):
        self.doc = bridge.load_deck(PUBLIC.read_bytes())
        self.control = self.doc['simulation']['solver']['solver_control']
        self.domain = self.doc['simulation']['physical_analysis']['domains'][0]

    def reject(self, doc=None):
        with self.assertRaises(bridge.Unsupported):
            bridge.translate(self.doc if doc is None else doc)

    def test_public_and_high_resolution(self):
        args, mesh = bridge.translate(self.doc)
        values = dict(zip(args[::2], args[1::2]))
        self.assertEqual(mesh, 'channel.exo')
        self.assertEqual(values['--rho'], '1.0')
        self.assertEqual(values['--mu'], '0.1')
        self.assertEqual(values['--inlet-velocity'], '0.1')
        self.assertEqual(values['--outlet-beta'], '0.05')
        self.assertEqual(values['--pseudo-dt'], '0.01')
        self.assertEqual(values['--relax-mass'], '0.75')
        self.control['basic_settings']['advection_scheme'] = 'high_resolution'
        self.control['expert_parameters']['blend_factor_max'] = 1
        self.assertIn('high-resolution', bridge.translate(self.doc)[0])

    def test_each_expert_mismatch(self):
        for key, value in dict(consistent=True, fractional_step_method=True,
                              incremental_gradient_change=False, limit_gradients=True,
                              relax_gradients=True, correct_gradients=True,
                              false_mass_accumulation=False, nonlinear_stabilisation=True,
                              disable_momentum_predictor=True, high_speed_blend_damping=True,
                              blend_factor_max=2).items():
            with self.subTest(key=key):
                doc = copy.deepcopy(self.doc)
                doc['simulation']['solver']['solver_control']['expert_parameters'][key] = value
                self.reject(doc)
        del self.control['expert_parameters']['relax_gradients']
        self.reject()  # OpenAccel defaults to relaxed gradients; MARS does not.

    def test_decomposition_is_not_a_physical_control(self):
        expected = bridge.translate(self.doc)
        for method in ('RCB', 'RIB', 'KWAY', ''):
            with self.subTest(method=method):
                self.doc['mesh']['automatic_decomposition_type'] = method
                self.assertEqual(bridge.translate(self.doc), expected)
        for method in (None, True, 1, [], {}):
            with self.subTest(method=method):
                self.doc['mesh']['automatic_decomposition_type'] = method
                with self.assertRaisesRegex(bridge.Unsupported, 'automatic_decomposition_type: expected a string'):
                    bridge.translate(self.doc)
        self.doc['mesh']['automatic_decomposition_type'] = 'RCB'
        for key in ('transformation', 'decomposition_properties', 'unknown'):
            with self.subTest(key=key):
                doc = copy.deepcopy(self.doc)
                doc['mesh'][key] = {}
                self.reject(doc)

    def test_unsupported_physics(self):
        for section, key, value in [
                ('domain', 'motion', {'option': 'rotating'}),
                ('domain_models', 'gravity', [0, 0, -9.81]),
                ('fluid_models', 'energy', {'option': 'total_energy'}),
                ('domain_models', 'reference_pressure', 101325)]:
            with self.subTest(key=key):
                doc = copy.deepcopy(self.doc)
                domain = doc['simulation']['physical_analysis']['domains'][0]
                target = domain if section == 'domain' else domain[section]
                target[key] = value
                self.reject(doc)
        self.domain['fluid_models']['turbulence']['option'] = 'sst'
        self.reject()

    def test_initial_and_time_state(self):
        self.domain['initialization']['velocity']['velocity'] = [0, 0, 1]
        self.reject()
        self.domain['initialization']['velocity']['velocity'] = [0, 0, 0]
        self.doc['simulation']['physical_analysis']['analysis_type']['option'] = 'transient'
        self.reject()

    def test_global_initialization(self):
        self.doc['simulation']['physical_analysis']['initialization'] = self.domain.pop('initialization')
        bridge.translate(self.doc)

    def test_boundary_contract(self):
        for option in ('static_pressure', 'mass_flow_rate'):
            doc = copy.deepcopy(self.doc)
            doc['simulation']['physical_analysis']['domains'][0]['boundaries'][2]['boundary_details']['mass_and_momentum']['option'] = option
            self.reject(doc)
        self.domain['boundaries'][0]['boundary_details'] = {'mass_and_momentum': {'option': 'no_slip_wall', 'wall_velocity': [1, 0, 0]}}
        self.reject()

    def test_duplicate_sets_and_names(self):
        self.domain['boundaries'][2]['location'] = ['inlet']
        self.reject()
        for value in [['with,comma'], ['with\nnewline'], [3], []]:
            with self.subTest(value=value), self.assertRaises(bridge.Unsupported):
                bridge.names(value, 'location')

    def test_numeric_fail_closed(self):
        for value in (True, 'NaN', 'inf', -1, 0, {'value': 1}):
            with self.subTest(value=value), self.assertRaises(bridge.Unsupported):
                bridge.number(value, 'density', positive=True)

    def test_duplicate_or_tagged_yaml(self):
        for text in ('mesh: {}\nmesh: {}', '!!python/object:object {}', 'x: &x {a: 1}\ny: {<<: *x}'):
            with self.subTest(text=text), self.assertRaises(bridge.Unsupported):
                bridge.load_deck(text)

    def test_interpolation_and_subiterations(self):
        self.control['basic_settings']['interpolation_scheme']['velocity_interpolation_type'] = 'linear_linear'
        self.reject()
        self.control['basic_settings']['interpolation_scheme']['velocity_interpolation_type'] = 'trilinear'
        self.control['advanced_options']['equation_controls']['sub_iterations']['pressure_correction'] = 2
        self.reject()

    def test_no_mesh_read_and_literal_arguments(self):
        scratch = ROOT / '.local-worktrees/simple-deck-tests'
        scratch.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=str(scratch)) as tmp:
            tmp = Path(tmp)
            deck, mesh, out = tmp/'input.i', tmp/'channel.exo', tmp/'prepared'
            # Not a real mesh: preparation must not try to parse it.
            mesh.write_bytes(b'not a mesh')
            text = PUBLIC.read_text().replace('location: [walls]', 'location: ["wall $(false); name"]')
            deck.write_text(text)
            argv = [str(ROOT/'scripts/prepare_simple_deck.py'), '--deck', str(deck), '--mesh', str(mesh), '--output', str(out)]
            original = Path.read_bytes
            def protected(path):
                self.assertNotEqual(path.resolve(), mesh.resolve())
                return original(path)
            with patch.object(sys, 'argv', argv), patch.object(Path, 'read_bytes', protected):
                bridge.main()
            args = (out/'args.nul').read_bytes().decode().split('\0')[:-1]
            self.assertIn('wall $(false); name', args)
            self.assertEqual(args, json.loads((out/'case.json').read_text())['arguments'])
            result = subprocess.run([sys.executable] + argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertNotEqual(result.returncode, 0)  # Never overwrite an existing preparation.
            self.assertIn(b'choose a fresh directory', result.stderr)
            argv[-1] = str(tmp/'new')
            argv[4] = str(tmp/'other.exo')
            result = subprocess.run([sys.executable] + argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(b'supplied mesh differs', result.stderr)
            self.assertFalse((tmp/'new').exists())


if __name__ == '__main__':
    unittest.main()
