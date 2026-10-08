"""Only synthetic decks and capture records are used; mesh and field files need not exist."""

import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import yaml

import simple_pressure_settings as settings
import simple_startup_probe as startup


class SourceSettingsTests(unittest.TestCase):
    def setUp(self):
        self.config = dict(family='Hypre', rtol=1e-10, atol=0., max_iterations=2000)
        self.arguments = {'--pressure-linear-rtol': '1e-10', '--pressure-linear-atol': '0'}

    def compare(self, env=None):
        return settings.compare_settings(self.config, self.arguments, env or {})

    def amg(self):
        self.config['options'] = dict(solver=dict(type='gmres', kdim=100),
                                      precond=dict(type='boomeramg'))
        return self.config['options']['precond']

    def test_absent_options_means_unpreconditioned_reference(self):
        self.config.pop('max_iterations')
        result = self.compare()
        self.assertEqual(result['known_differences'], ['maximum_iterations', 'preconditioner'])
        self.assertIn('krylov_method', result['known_matches'])
        self.assertIn('restart_dimension', result['unresolved_settings'])
        self.assertFalse(result['amg_comparison_applicable'])
        self.assertFalse(any(k.startswith('amg_') for k in result['known_matches']))

    def test_explicit_amg_settings_and_enforced_overrides(self):
        precond = self.amg()
        precond.update(coarsentype=8, interptype=6, relaxtype=18, relaxorder=0,
                       strongthreshold=.25, aggnumlevels=0, numsweeps=2,
                       maxlevels=25, maxcoarsesize=128, maxiter=70, tol=.03)
        result = self.compare()
        self.assertEqual(result['known_differences'], [])
        for name in ('amg_max_iterations', 'amg_tolerance', 'amg_coarsen_type',
                     'amg_num_sweeps', 'restart_dimension'):
            self.assertIn(name, result['known_matches'])
        # Even two omissions are unresolved across differently built libraries.
        for name in ('minimum_iterations', 'amg_numpaths', 'amg_coarse_relax_type'):
            self.assertIn(name, result['unresolved_settings'])

    def test_captured_environment_overrides_mars_source_defaults(self):
        precond = self.amg()
        precond.update(coarsentype=10, strongthreshold=.5, numsweeps=1)
        env = dict(MARS_AMG_COARSEN='10', MARS_AMG_STRONG='.5', MARS_AMG_SWEEPS='1',
                   MARS_HYPRE_FLEXGMRES='1', MARS_HYPRE_ABSTOL='999')
        result = self.compare(env)
        self.assertEqual(result['known_differences'], ['krylov_method'])
        self.assertIn('amg_coarsen_type', result['known_matches'])
        self.assertIn('pressure_atol', result['known_matches'])
        self.assertEqual(self.compare()['known_differences'],
                         ['amg_coarsen_type', 'amg_num_sweeps', 'amg_strong_threshold'])

    def test_flex_flag_uses_cpp_string_semantics(self):
        self.amg()
        self.config['options']['solver']['type'] = 'FlexGMRES'
        for flag in ('1', '', 'false', 'PRIVATE_STRING'):
            self.assertIn('krylov_method', self.compare({'MARS_HYPRE_FLEXGMRES': flag})['known_matches'])
        self.assertIn('krylov_method', self.compare({'MARS_HYPRE_FLEXGMRES': '0'})['known_differences'])

    def test_malformed_numeric_environment_stays_unresolved(self):
        self.amg().update(coarsentype=8, strongthreshold=.25)
        for value in ('8suffix', 'garbage', '8.0', '1e3', str(2**90)):
            result = self.compare({'MARS_AMG_COARSEN': value})
            self.assertIn('amg_coarsen_type', result['unresolved_settings'])
        for value in ('NaN', 'inf', '.25suffix'):
            result = self.compare({'MARS_AMG_STRONG': value})
            self.assertIn('amg_strong_threshold', result['unresolved_settings'])
        self.assertIn('amg_coarsen_type', self.compare({'MARS_AMG_COARSEN': ''})['known_matches'])

    def test_ignored_reference_options_do_not_become_settings(self):
        self.amg().update(pmaxelmts=4, PRIVATE_OPTION='PRIVATE_VALUE')
        self.config['min_iterations'] = 3
        result = self.compare()
        self.assertTrue(result['reference_ignored_solver_option_keys_present'])
        self.assertIn('amg_p_max_elmts', result['unresolved_settings'])
        self.assertIn('minimum_iterations', result['unresolved_settings'])
        self.assertFalse(result['reference_min_iterations_forwarded'])
        self.assertNotIn('PRIVATE', json.dumps(result))

    def test_case_insensitive_options_but_required_type_key_is_exact(self):
        precond = self.amg()
        precond.update(CoarsenType=8)
        self.assertIn('amg_coarsen_type', self.compare()['known_matches'])
        precond['coarsentype'] = 9
        with self.assertRaises(ValueError):
            self.compare()
        del precond['coarsentype']
        precond['TYPE'] = precond.pop('type')
        with self.assertRaises(ValueError):
            self.compare()

    def test_unmodelled_preconditioner_is_not_amg_parity(self):
        self.amg()['type'] = 'mgr'
        result = self.compare()
        self.assertIn('preconditioner', result['known_differences'])
        self.assertFalse(result['reference_preconditioner_modelled'])
        self.assertFalse(result['amg_comparison_applicable'])

    def test_tolerance_and_scaling_mismatches_are_reported(self):
        self.config.update(rtol=1e-6, atol=1e-8, diagonal_scaling=True, normalize_matrix=True)
        result = self.compare()
        for label in ('pressure_rtol', 'pressure_atol', 'diagonal_scaling', 'normalize_matrix'):
            self.assertIn(label, result['known_differences'])


class CaptureTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.pair = self.root / 'PRIVATE_PAIR'
        (self.pair / 'reference').mkdir(parents=True)
        (self.pair / 'mars').mkdir()
        config = dict(family='Hypre', rtol=1e-10, atol=0., max_iterations=20,
                      options=dict(solver=dict(type='gmres', kdim=30),
                                   precond=dict(type='boomeramg', numsweeps=1)))
        self.deck = dict(simulation=dict(solver=dict(solver_control=dict(advanced_options=dict(
            linear_solver_settings=dict(pressure_correction=config))))))
        self.case = dict(format='mars-simple-deck-v1', pressure_linear_policy='reference',
                         arguments=['--mesh', str(self.pair / 'PRIVATE_MISSING_MESH.exo'),
                                    '--mesh-format', 'exodus', '--pressure-linear-rtol', '1e-10',
                                    '--pressure-linear-atol', '0'])
        self.record = dict(schema=startup.SCHEMA, steps=1, first_step_audit=True,
                           mesh=str(self.pair / 'PRIVATE_MISSING_MESH.exo'), mesh_sha256='a'*64)
        for directory in ('reference', 'mars'):
            (self.pair / directory / 'run.log').write_text('PRIVATE LOG DATA\n')
            (self.pair / directory / 'run.exit').write_text('0\n')
        self.output = self.root / 'public.json'
        self.env = {}
        self.bind()

    def bind(self):
        (self.pair / 'reference/input.i').write_text(yaml.safe_dump(self.deck))
        self.case['deck_sha256'] = startup.digest(self.pair / 'reference/input.i')
        (self.pair / 'case.json').write_text(json.dumps(self.case))
        self.record.update(deck_sha256=self.case['deck_sha256'], case_sha256=startup.digest(self.pair / 'case.json'))
        (self.pair / 'pair.json').write_text(json.dumps(self.record))
        for solver, directory in (('openaccel', 'reference'), ('mars', 'mars')):
            launch = dict(schema=startup.SCHEMA, solver=solver, ranks=4, status='started',
                          pair_sha256=startup.digest(self.pair / 'pair.json'),
                          executable='/PRIVATE_EXECUTABLE', executable_sha256='c'*64,
                          libraries={'/PRIVATE_LIBRARY': 'd'*64}, environment=self.env if solver == 'mars' else {},
                          command=['PRIVATE_LAUNCHER', '/PRIVATE_EXECUTABLE'] + startup.solver_arguments(self.pair, solver))
            path = self.pair / directory
            (path / 'launch-start.json').write_text(json.dumps(launch))
            files = {name: startup.digest(path / name) for name in ('run.log', 'run.exit')}
            if solver == 'openaccel':
                files['input.i'] = startup.digest(path / 'input.i')
            files['PRIVATE_MISSING_FIELDS.csv'] = 'e'*64
            launch.update(status='finished', exit_code=0, files=files)
            (path / 'launch.json').write_text(json.dumps(launch))

    def invoke(self):
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
            code = settings.main(['--pair', str(self.pair), '--output', str(self.output)])
        text = self.output.read_text()
        self.assertNotIn('PRIVATE', stream.getvalue() + text)
        self.assertNotIn(str(self.root), stream.getvalue() + text)
        return code, json.loads(text)

    def test_read_only_saved_settings_without_mesh_fields_or_launch(self):
        before = {p: p.read_bytes() for p in self.pair.rglob('*') if p.is_file()}
        with patch('subprocess.Popen', side_effect=AssertionError('must not launch')), \
                patch('subprocess.run', side_effect=AssertionError('must not inspect runtime')):
            code, public = self.invoke()
        self.assertEqual(code, 0)
        self.assertEqual(public['known_differences'], ['amg_num_sweeps', 'maximum_iterations', 'restart_dimension'])
        for key in ('binary_source_identity_verified', 'effective_library_defaults_verified',
                    'launcher_child_environment_verified', 'field_provenance_verified',
                    'identical_linear_solvers_verified', 'reference_runtime_convergence_verified'):
            self.assertFalse(public[key])
        self.assertEqual(before, {p: p.read_bytes() for p in self.pair.rglob('*') if p.is_file()})

    def test_saved_environment_used_not_current_shell(self):
        self.env = {'MARS_AMG_SWEEPS': '1'}
        self.bind()
        with patch.dict(os.environ, {'MARS_AMG_SWEEPS': '7'}):
            _, result = self.invoke()
        self.assertIn('amg_num_sweeps', result['known_matches'])

    def test_named_and_fallback_reference_resolution(self):
        solver = self.deck['simulation']['solver']
        config = solver['solver_control']['advanced_options']['linear_solver_settings']['pressure_correction']
        solver['PRIVATE_LOOKUP'] = config
        for key in ('pressure_correction', 'segregated_flow', 'default'):
            solver['solver_control']['advanced_options']['linear_solver_settings'] = {key: {'lookup': 'PRIVATE_LOOKUP'}}
            self.bind()
            self.assertEqual(settings.inspect(self.pair)['comparison_status'], 'completed')

    def test_missing_or_changed_capture_is_rejected(self):
        for name in ('case.json', 'reference/input.i', 'mars/run.log', 'reference/run.exit',
                     'mars/launch.json', 'reference/launch-start.json'):
            path = self.pair / name
            saved = path.read_bytes()
            for replacement in (None, saved + b'PRIVATE_TAMPERING'):
                with self.subTest(name=name, replacement=replacement is None):
                    if replacement is None:
                        path.unlink()
                    else:
                        path.write_bytes(replacement)
                    result = settings.inspect(self.pair)
                    self.assertEqual(result['comparison_status'], 'invalid_evidence')
                    self.assertNotIn('PRIVATE', json.dumps(result))
                    path.write_bytes(saved)

    def test_changed_command_or_environment_cannot_pass_metadata(self):
        path = self.pair / 'mars/launch.json'
        original = json.loads(path.read_text())
        for key, value in (('command', ['PRIVATE_WRONG']), ('environment', {'MARS_HYPRE_FLEXGMRES': '1'}),
                           ('pair_sha256', 'f'*64), ('exit_code', 2), ('executable_sha256', 'bad')):
            path.write_text(json.dumps(dict(original, **{key: value})))
            self.assertEqual(settings.inspect(self.pair)['failed_check'], 'saved_launch_metadata')
        path.write_text(json.dumps(original))

    def test_failed_history_still_has_auditable_controls(self):
        path = self.pair / 'mars'
        (path / 'run.exit').write_text('255\n')
        self.bind()
        record = json.loads((path / 'launch.json').read_text())
        record.update(status='failed', exit_code=255)
        (path / 'launch.json').write_text(json.dumps(record))
        self.assertEqual(settings.inspect(self.pair)['comparison_status'], 'completed')

    def test_recovery_enabled_or_attempted_rejected(self):
        path = self.pair / 'mars/run.log'
        for text in ('pressure_refinement=1\n', '[simple-pressure-refinement] recovered=1\n',
                     'pressure_refinement=0 pressure_refinement=0\n'):
            path.write_text(text)
            self.bind()
            self.assertEqual(settings.inspect(self.pair)['failed_check'], 'experimental_recovery')
        path.write_text('pressure_refinement=0\n')
        self.bind()
        self.assertEqual(settings.inspect(self.pair)['comparison_status'], 'completed')

    def test_recovery_argument_rejected_without_a_log_marker(self):
        self.case['arguments'] += ['--pressure-refinement', '1']
        self.bind()
        self.assertEqual(settings.inspect(self.pair)['failed_check'], 'experimental_recovery')

    def test_parser_errors_do_not_echo_arguments(self):
        stream = io.StringIO()
        with contextlib.redirect_stderr(stream), self.assertRaises(SystemExit):
            settings.main(['--pair', 'PRIVATE_PAIR', '--output', 'PRIVATE_OUTPUT', '--PRIVATE_OPTION'])
        self.assertNotIn('PRIVATE', stream.getvalue())

    def test_invalid_config_and_existing_output_do_not_leak(self):
        config = self.deck['simulation']['solver']['solver_control']['advanced_options']['linear_solver_settings']['pressure_correction']
        config['family'] = 'PRIVATE_INVALID'
        self.bind()
        code, public = self.invoke()
        self.assertEqual(code, 1)
        self.assertEqual(public['failed_check'], 'pressure_configuration')
        saved = self.output.read_bytes()
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(settings.main(['--pair', str(self.pair), '--output', str(self.output)]), 1)
        self.assertEqual(self.output.read_bytes(), saved)


if __name__ == '__main__':
    unittest.main()
