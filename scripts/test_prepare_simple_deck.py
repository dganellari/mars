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

    def reference_values(self):
        args, _ = bridge.translate(self.doc, pressure_linear_policy='reference')
        return dict(zip(args[::2], args[1::2]))

    def test_reference_pressure_precedence_and_defaults(self):
        settings = self.control['advanced_options']['linear_solver_settings']
        settings['default'] = {'family': 'Hypre', 'rtol': 1e-4, 'atol': 1e-10}
        settings['segregated_flow'] = {'family': 'Hypre', 'rtol': 1e-5, 'atol': 0}
        settings['pressure_correction'] = {'family': 'Hypre'}
        expected_default = bridge.translate(self.doc)
        self.assertEqual(expected_default, bridge.translate(self.doc, pressure_linear_policy='mars'))
        self.assertNotIn('--pressure-linear-rtol', expected_default[0])
        for key, rtol, atol in [('pressure_correction', 1e-6, 1e-16),
                                ('segregated_flow', 1e-5, 0), ('default', 1e-4, 1e-10)]:
            values = self.reference_values()
            self.assertEqual(float(values['--pressure-linear-rtol']), rtol)
            self.assertEqual(float(values['--pressure-linear-atol']), atol)
            normal = dict(zip(expected_default[0][::2], expected_default[0][1::2]))
            for flag, value in normal.items():
                self.assertEqual(values[flag], value)
            del settings[key]
        with self.assertRaisesRegex(bridge.Unsupported, 'missing or malformed pressure solver definition'):
            self.reference_values()

    def test_reference_pressure_named_solver_and_case_sensitive_lookup(self):
        solver = self.doc['simulation']['solver']
        settings = self.control['advanced_options']['linear_solver_settings']
        solver['Private Solver Name'] = {'family': 'hYpRe', 'rtol': '2e-5', 'atol': '0',
                                         'options': {'solver': {'type': 'fLeXgMrEs', 'kdim': 17},
                                                     'preconditioner': {'type': 'boomeramg'}}}
        settings['pressure_correction'] = {'lookup': 'Private Solver Name', 'rtol': 0.25}
        self.assertEqual(float(self.reference_values()['--pressure-linear-rtol']), 2e-5)
        self.assertEqual(float(self.reference_values()['--pressure-linear-atol']), 0)
        settings['pressure_correction']['lookup'] = 'private solver name'
        with self.assertRaises(bridge.Unsupported) as error:
            self.reference_values()
        self.assertNotIn('private solver name', str(error.exception))

    def test_reference_pressure_scaling_rejected(self):
        settings = self.control['advanced_options']['linear_solver_settings']
        for key in ('normalize_matrix', 'diagonal_scaling'):
            settings['pressure_correction'] = {'family': 'Hypre', key: False}
            self.reference_values()
            for value in (True, None, 0, 1, 'false', 'private scaling', [], {}):
                with self.subTest(key=key, value=value):
                    settings['pressure_correction'][key] = value
                    with self.assertRaisesRegex(bridge.Unsupported, key) as error:
                        self.reference_values()
                    self.assertNotIn('private scaling', str(error.exception))

    def test_reference_pressure_tolerances_rejected(self):
        settings = self.control['advanced_options']['linear_solver_settings']
        for key, values in [('rtol', (0, -1, 1, 2, True, None, 'nan', 'inf', [], {}, 'private value')),
                            ('atol', (-1, True, None, 'nan', 'inf', [], {}, 'private value'))]:
            for value in values:
                with self.subTest(key=key, value=value):
                    settings['pressure_correction'] = {'family': 'Hypre', key: value}
                    with self.assertRaisesRegex(bridge.Unsupported, 'reference.' + key) as error:
                        self.reference_values()
                    self.assertNotIn('private value', str(error.exception))

    def test_reference_pressure_unsupported_or_malformed_solver(self):
        settings = self.control['advanced_options']['linear_solver_settings']
        configs = [None, [], True, 'private value', {}, {'family': 'Trilinos'}, {'family': []},
                   {'lookup': []}, {'lookup': 'private missing name'}]
        configs += [{'family': 'Hypre', 'options': value}
                    for value in (None, [], True, 'private value', {}, {'solver': None},
                                  {'solver': []}, {'solver': {}}, {'solver': {'type': []}})]
        configs += [{'family': 'Hypre', 'options': {'solver': {'type': value}}}
                    for value in ('boomeramg', 'mgr', 'private value')]
        for config in configs:
            with self.subTest(config=config):
                settings['pressure_correction'] = config
                with self.assertRaises(bridge.Unsupported) as error:
                    self.reference_values()
                self.assertNotIn('private', str(error.exception))
        for value in (None, [], True, 'private value'):
            self.control['advanced_options']['linear_solver_settings'] = value
            with self.assertRaises(bridge.Unsupported):
                self.reference_values()

    def test_reference_pressure_cli_and_redaction(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            deck, mesh, out = tmp/'private-deck.i', tmp/'channel.exo', tmp/'prepared'
            mesh.write_bytes(b'not a mesh')
            settings = self.control['advanced_options']['linear_solver_settings']
            settings['pressure_correction'] = {'family': 'Hypre', 'options': {'solver': {'type': 'GMRES'}}}
            deck.write_text(json.dumps(self.doc))
            argv = [sys.executable, str(ROOT/'scripts/prepare_simple_deck.py'), '--deck', str(deck),
                    '--mesh', str(mesh), '--output', str(out), '--pressure-linear-policy', 'reference']
            result = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(result.returncode, 0, result.stderr)
            record = json.loads((out/'case.json').read_text())
            self.assertEqual(record['pressure_linear_policy'], 'reference')
            self.assertIn('--pressure-linear-rtol', record['arguments'])
            self.assertIn('--pressure-linear-atol', record['arguments'])
            self.assertTrue(any('full solver configuration are not translated' in note for note in record['notes']))
            self.assertEqual((out/'args.nul').read_bytes().decode().split('\0')[:-1], record['arguments'])
            argv[argv.index('--output') + 1] = str(tmp/'rejected')
            settings['pressure_correction'] = {'lookup': 'private missing solver'}
            deck.write_text(json.dumps(self.doc))
            result = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(result.returncode, 1)
            self.assertNotIn(b'private', result.stderr)
            self.assertNotIn(str(tmp).encode(), result.stderr)
            self.assertFalse(result.stdout)
            self.assertFalse((tmp/'rejected').exists())
            argv[-1] = 'private invalid policy'
            result = subprocess.run(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(result.returncode, 2)
            self.assertNotIn(b'private', result.stderr)
            self.assertFalse((tmp/'rejected').exists())

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

    def test_unused_coupling_key_does_not_change_arguments(self):
        expected = bridge.translate(self.doc)
        expert = self.control['expert_parameters']
        for value in (False, True):
            expert['coupled_pressure_velocity'] = value
            self.assertEqual(bridge.translate(self.doc), expected)
        for value in (None, 0, 1, 'false', 'private-value', [], {}):
            with self.subTest(value=value):
                expert['coupled_pressure_velocity'] = value
                with self.assertRaisesRegex(bridge.Unsupported, 'coupled_pressure_velocity: expected a boolean') as error:
                    bridge.translate(self.doc)
                self.assertNotIn('private-value', str(error.exception))
        expert['coupled_pressure_velocity'] = False
        expert['private-unknown-key'] = 'private-value'
        with self.assertRaisesRegex(bridge.Unsupported, 'expert_parameters: unrecognized keys') as error:
            bridge.translate(self.doc)
        self.assertNotIn('private-', str(error.exception))

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
        interp = self.control['basic_settings']['interpolation_scheme']
        for scheme in ('trilinear', 'linear_linear'):
            interp['velocity_interpolation_type'] = scheme
            args, _ = bridge.translate(self.doc)
            self.assertEqual(dict(zip(args[::2], args[1::2]))['--velocity-interpolation'], scheme.replace('_', '-'))
        for invalid in ('unsupported private setting', None, [], True):
            interp['velocity_interpolation_type'] = invalid
            self.reject()
        interp['velocity_interpolation_type'] = 'trilinear'
        for key in ('pressure_interpolation_type', 'velocity_gradient_interpolation_type', 'pressure_gradient_interpolation_type'):
            interp[key] = 'trilinear'
            self.reject()
            interp[key] = 'linear_linear'
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
