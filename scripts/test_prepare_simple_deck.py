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
PUBLIC = ROOT / 'tests/reference/openaccel/distributed_simple/preparation_fixture.i'


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

    def test_named_linear_solvers(self):
        expected = bridge.translate(self.doc)
        solver = self.doc['simulation']['solver']
        settings = self.control['advanced_options']['linear_solver_settings']
        for family in ('PETSc', 'HYPRE', 'Trilinos', 'AMGsolver', 'gmres'):
            with self.subTest(family=family):
                solver['user named solver'] = {'family': family, 'options': {'backend setting': 3}}
                settings['default'] = {'lookup': 'user named solver'}
                self.assertEqual(bridge.translate(self.doc), expected)
        # Unused definitions are also part of OpenAccel's solver library.
        solver['unused solver'] = {'family': 'hypre'}
        solver['pressure solver'] = {'family': 'Hypre'}
        settings['pressure_correction'] = {'lookup': 'pressure solver'}
        settings['coupled_navier_stokes'] = {'family': 'PETSc'}
        self.assertEqual(bridge.translate(self.doc), expected)

    def test_invalid_linear_solver_references(self):
        solver = self.doc['simulation']['solver']
        settings = self.control['advanced_options']['linear_solver_settings']
        for name in ('missing private name', 'solver_control', 'output_control', '', None, True, [], {}):
            with self.subTest(name=name):
                settings['default'] = {'lookup': name}
                with self.assertRaisesRegex(bridge.Unsupported, 'linear_solver_settings.entry.lookup') as error:
                    bridge.translate(self.doc)
                self.assertNotIn('missing private name', str(error.exception))
        settings['default'] = {'lookup': 'named'}
        for config in (None, [], 'private value', {'family': 'unknown'}, {'lookup': 'named'}):
            with self.subTest(config=config):
                solver['named'] = config
                self.reject()
        del solver['named']
        for config in ([], None, 'private value', {}, {'family': []}):
            with self.subTest(config=config):
                settings['default'] = config
                self.reject()
        self.control['advanced_options']['linear_solver_settings'] = []
        self.reject()

    def test_restart_is_never_ignored_as_solver_metadata(self):
        for config in ({'file_path': '/private/saved-fields.e'}, {}, None, {'family': 'Hypre'}):
            with self.subTest(config=config):
                self.doc['simulation']['solver']['restart_control'] = config
                with self.assertRaisesRegex(bridge.Unsupported, 'solver.restart_control: loads saved fields') as error:
                    bridge.translate(self.doc)
                self.assertNotIn('/private/saved-fields.e', str(error.exception))

    def test_explicit_subsonic_boundary(self):
        expected = bridge.translate(self.doc)
        for b in self.domain['boundaries'][1:]:
            b['boundary_details']['option'] = 'subsonic'
        self.assertEqual(bridge.translate(self.doc), expected)
        for index in (1, 2):
            doc = copy.deepcopy(self.doc)
            doc['simulation']['physical_analysis']['domains'][0]['boundaries'][index]['boundary_details']['option'] = 'supersonic'
            self.reject(doc)

    def test_independent_errors_reported_together_without_private_values(self):
        self.doc['mesh']['private mesh key'] = 'secret-mesh-value'
        self.domain['initialization']['pressure']['pressure'] = 42
        self.domain['boundaries'][2]['boundary_details']['mass_and_momentum']['option'] = 'secret-outlet-mode'
        self.doc['simulation']['solver']['private solver key'] = {'private setting': 'secret-solver-value'}
        self.control['expert_parameters']['relax_gradients'] = True
        with self.assertRaises(bridge.Unsupported) as error:
            bridge.translate(self.doc)
        message = str(error.exception)
        for path in ('mesh:', 'initialization.pressure:', 'outlet.option:', 'solver:', 'expert_parameters.relax_gradients:'):
            self.assertIn(path, message)
        for private in ('private', 'secret', '42'):
            self.assertNotIn(private, message)

    def test_malformed_sections_are_rejected_without_crashing(self):
        paths = [('mesh',), ('simulation',), ('simulation', 'physical_analysis'),
                 ('simulation', 'physical_analysis', 'domains'),
                 ('simulation', 'physical_analysis', 'domains', 0),
                 ('simulation', 'physical_analysis', 'domains', 0, 'boundaries'),
                 ('simulation', 'physical_analysis', 'domains', 0, 'boundaries', 0),
                 ('simulation', 'solver'), ('simulation', 'solver', 'solver_control'),
                 ('simulation', 'solver', 'solver_control', 'basic_settings'),
                 ('simulation', 'solver', 'solver_control', 'basic_settings', 'advection_scheme')]
        for path in paths:
            for value in (None, [], True, 'private text'):
                with self.subTest(path=path, value=value):
                    doc = copy.deepcopy(self.doc)
                    target = doc
                    for key in path[:-1]:
                        target = target[key]
                    target[path[-1]] = value
                    self.reject(doc)

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
        for option in ('total_pressure', 'mass_flow_rate', 'normal_speed', None):
            doc = copy.deepcopy(self.doc)
            doc['simulation']['physical_analysis']['domains'][0]['boundaries'][2]['boundary_details']['mass_and_momentum']['option'] = option
            self.reject(doc)
        self.domain['boundaries'][0]['boundary_details'] = {'mass_and_momentum': {'option': 'no_slip_wall', 'wall_velocity': [1, 0, 0]}}
        self.reject()

    def test_constant_static_pressure(self):
        mm = self.domain['boundaries'][2]['boundary_details']['mass_and_momentum']
        mm['option'] = 'static_pressure'
        # An average-pressure blend left in the deck does not affect this mode.
        for pressure in (0, 7.5, -3.5):
            with self.subTest(pressure=pressure):
                mm['relative_pressure'] = pressure
                with_blend = bridge.translate(self.doc)
                mm.pop('pressure_profile_blend', None)
                self.assertEqual(bridge.translate(self.doc), with_blend)
                args = dict(zip(with_blend[0][::2], with_blend[0][1::2]))
                self.assertEqual(float(args['--outlet-pressure']), pressure)
                self.assertEqual(float(args['--outlet-beta']), 1)
                mm['pressure_profile_blend'] = 0.05
        for pressure in (None, True, 'nan', 'inf', [0], {'value': 0}):
            with self.subTest(pressure=pressure):
                mm['relative_pressure'] = pressure
                self.reject()
        mm['relative_pressure'] = 0
        mm['option'] = 'average_static_pressure'
        del mm['pressure_profile_blend']
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

    def test_cli_rejection_has_no_partial_arguments_or_private_text(self):
        scratch = ROOT / '.local-worktrees/simple-deck-tests'
        scratch.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=str(scratch)) as tmp:
            tmp = Path(tmp)
            deck, out = tmp/'private-deck.i', tmp/'prepared'
            argv = [sys.executable, str(ROOT/'scripts/prepare_simple_deck.py'),
                    '--deck', str(deck), '--mesh', str(tmp/'private-mesh.exo'), '--output', str(out)]
            # A file-access failure must not print its private path.
            result = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(result.returncode, 1)
            self.assertIn(b'file access', result.stderr)
            self.assertNotIn(b'private-', result.stderr)
            self.doc['mesh']['private-key'] = 'secret-value'
            self.doc['simulation']['solver']['restart_control'] = {'file_path': '/private/saved.e'}
            self.control['expert_parameters']['relax_gradients'] = True
            # JSON is a YAML subset; no additional fixture or private input is needed.
            deck.write_text(json.dumps(self.doc))
            result = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(result.returncode, 1)
            for path in (b'ERROR: mesh:', b'ERROR: solver.restart_control:', b'ERROR: expert_parameters.relax_gradients:'):
                self.assertIn(path, result.stderr)
            for private in (b'private', b'secret', str(tmp).encode()):
                self.assertNotIn(private, result.stderr)
            self.assertFalse(out.exists())
            self.assertFalse(result.stdout)

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
