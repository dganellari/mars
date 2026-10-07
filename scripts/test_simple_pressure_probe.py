"""Pressure-only experiments use synthetic decks, matrices and launch records."""
import contextlib
import copy
import io
import json
import os
from pathlib import Path
import shutil
import unittest
from unittest.mock import patch

import yaml

import simple_pressure_probe as pressure
import simple_startup_probe as startup
import test_simple_first_step_audit as fixtures


class PressureProbeTests(unittest.TestCase):
    def setUp(self):
        self.capture = fixtures.FirstStepTests()
        self.capture.setUp()
        self.addCleanup(self.capture.doCleanups)
        self.baseline = self.capture.pair
        self.root = self.capture.fixture.root
        self.pair = self.root / 'tightened'
        self.output = self.root / 'public.json'
        self.details = self.root / 'details'

    def deck(self):
        return startup.load_deck((self.baseline / 'reference/input.i').read_bytes())

    def prepare(self):
        pressure.prepare(self.baseline, self.pair)

    def hashes(self, path):
        return {str(p.relative_to(path)): startup.digest(p) for p in path.rglob('*') if p.is_file()}

    def invoke(self, *args):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return pressure.main(list(args))

    def finish(self):
        self.prepare()
        (self.pair / 'mars').mkdir()
        for directory in ('reference', 'mars'):
            for path in (self.baseline / directory).iterdir():
                if path.name not in ('input.i', 'launch-start.json', 'launch.json'):
                    shutil.copyfile(path, self.pair / directory / path.name)
        log = self.pair / 'mars/run.log'
        log.write_text(log.read_text().replace('pressure_linear_rtol=0.0001 pressure_linear_atol=1e-10',
                                              'pressure_linear_rtol=1e-10 pressure_linear_atol=0.0'))
        self.capture.pair = self.pair
        for solver in ('mars', 'openaccel'):
            self.capture.record(solver)
        self.capture.pair = self.baseline
        self.capture.mars = self.pair / 'mars'

    def comparison(self):
        code = self.invoke('compare', '--pair', str(self.pair), '--output', str(self.output),
                           '--detail-dir', str(self.details))
        text = self.output.read_text()
        self.assertNotIn(str(self.root), text)
        self.assertNotIn('PRIVATE', text)
        return code, json.loads(text)

    def test_prepare_preserves_baseline_and_changes_only_pressure(self):
        before = self.hashes(self.baseline)
        self.prepare()
        pressure.check_inputs(self.pair)
        original, expected = self.deck(), self.deck()
        expected['simulation']['solver']['solver_control']['advanced_options']['linear_solver_settings'][
            'pressure_correction'].update(rtol=1e-10, atol=0.)
        actual = startup.load_deck((self.pair / 'reference/input.i').read_bytes())
        self.assertEqual(actual, expected)
        self.assertEqual(self.deck(), original)
        self.assertEqual(self.hashes(self.baseline), before)
        args = startup.options(startup.read_json(self.pair / 'case.json')['arguments'])
        self.assertEqual(float(args['--pressure-linear-rtol']), 1e-10)
        self.assertEqual(float(args['--pressure-linear-atol']), 0.)
        self.assertEqual(startup.read_json(self.pair / 'case.json')['pressure_linear_policy'], 'reference')

    def test_shared_fallback_and_named_solver_are_cloned(self):
        for name in ('default', 'segregated_flow'):
            for named in (False, True):
                with self.subTest(fallback=name, lookup=named):
                    original = self.deck()
                    solver = original['simulation']['solver']
                    config = dict(family='Hypre', rtol=1e-5, atol=1e-8, max_iterations=432,
                                  options={'solver': {'type': 'FlexGMRES', 'kdim': 37},
                                           'precond': {'type': 'boomeramg'}}, write_system=True)
                    if named:
                        solver['PRIVATE'] = copy.deepcopy(config)
                    solver['solver_control']['advanced_options']['linear_solver_settings'] = {
                        name: {'lookup': 'PRIVATE'} if named else copy.deepcopy(config)}
                    saved = copy.deepcopy(original)
                    changed, _ = pressure.tightened_deck(original)
                    self.assertEqual(original, saved)
                    self.assertEqual(pressure.resolved_solver(changed, 'coupled_navier_stokes'), config)
                    self.assertEqual(pressure.resolved_solver(changed, 'pressure_correction'),
                                     dict(config, rtol=1e-10, atol=0.))
                    self.assertEqual(changed['simulation']['solver']['solver_control']['advanced_options'][
                        'linear_solver_settings'][name], solver['solver_control']['advanced_options']['linear_solver_settings'][name])

    def test_tighter_existing_relative_target_is_not_relaxed(self):
        original = self.deck()
        config = original['simulation']['solver']['solver_control']['advanced_options']['linear_solver_settings']['pressure_correction']
        config.update(rtol=1e-12, atol=1e-8)
        changed, _ = pressure.tightened_deck(original)
        self.assertEqual(pressure.resolved_solver(changed, 'pressure_correction')['rtol'], 1e-12)
        config['atol'] = 0.
        with self.assertRaisesRegex(startup.EvidenceError, 'pressure_target_unchanged'):
            pressure.tightened_deck(original)

    def test_unsupported_reference_configuration_rejected(self):
        for changes in ({'family': 'PETSc'}, {'normalize_matrix': True}, {'diagonal_scaling': True},
                        {'rtol': float('nan')}, {'options': {'solver': {'type': 'cg'}}}):
            with self.subTest(changes=changes):
                original = self.deck()
                original['simulation']['solver']['solver_control']['advanced_options']['linear_solver_settings'][
                    'pressure_correction'].update(changes)
                with self.assertRaises(ValueError):
                    pressure.tightened_deck(original)

    def test_missing_or_changed_baseline_capture_rejected_before_output(self):
        path = self.baseline / 'mars/flow-audit-rank000001.bin'
        with path.open('ab') as stream:
            stream.write(b'PRIVATE')
        with self.assertRaises(ValueError):
            self.prepare()
        self.assertFalse(self.pair.exists())

    def test_changed_fresh_controls_rejected_even_with_updated_hashes(self):
        self.prepare()
        path = self.pair / 'reference/input.i'
        deck = startup.load_deck(path.read_bytes())
        config = deck['simulation']['solver']['solver_control']['advanced_options']['linear_solver_settings']['pressure_correction']
        config['max_iterations'] = 999
        path.write_text(yaml.safe_dump(deck))
        case = startup.read_json(self.pair / 'case.json')
        case['deck_sha256'] = startup.digest(path)
        (self.pair / 'case.json').write_text(json.dumps(case))
        pair = startup.read_json(self.pair / 'pair.json')
        pair.update(deck_sha256=case['deck_sha256'], case_sha256=startup.digest(self.pair / 'case.json'))
        (self.pair / 'pair.json').write_text(json.dumps(pair))
        with self.assertRaisesRegex(startup.EvidenceError, 'nonpressure_controls_changed'):
            pressure.check_inputs(self.pair)

    def test_complete_comparison_is_bounded_and_preserves_baseline(self):
        self.finish()
        before = self.hashes(self.baseline)
        code, result = self.comparison()
        self.assertEqual(code, 0, result)
        self.assertEqual(result['outcome'], 'first_step_matches')
        self.assertTrue(result['tighter_pressure_residuals_pass'])
        self.assertFalse(result['root_cause_proven'])
        self.assertFalse(result['nonlinear_convergence_verified'])
        self.assertTrue(result['pressure_only_controls_verified'])
        self.assertEqual(self.hashes(self.baseline), before)

    def test_changed_binary_library_rank_launcher_or_environment_rejected(self):
        self.finish()
        record = startup.verified_launch(self.pair, 'mars')
        old = startup.verified_launch(self.baseline, 'mars')
        mutations = ({'executable_sha256': '0'*64}, {'libraries': {'changed': 'b'*64}}, {'ranks': 4},
                     {'command': ['changed-launcher'] + record['command'][1:]},
                     {'environment': {'MARS_HYPRE_ABSTOL': '1'}})
        for changes in mutations:
            with self.subTest(changes=changes), self.assertRaises(startup.EvidenceError):
                pressure.runtime_matches(old, dict(record, **changes), self.pair, 'mars')

    def test_matching_fields_with_failed_tight_residual_cannot_pass(self):
        self.finish()
        self.capture.pair = self.pair
        self.capture.phi = self.capture.phi + 1e-9
        self.capture.solutions['pressure'] = self.capture.phi
        self.capture.write_mars_parts()
        self.capture.change_mars_final()
        self.capture.record('mars')
        self.capture.pair = self.baseline
        code, result = self.comparison()
        self.assertEqual(code, 0, result)
        self.assertTrue(all(result['tightened_first_step']['first_step_stage_matches'].values()))
        self.assertFalse(result['tighter_pressure_residuals_pass'])
        self.assertEqual(result['outcome'], 'pressure_target_not_met')

    def test_upstream_change_is_not_classified_as_pressure_disagreement(self):
        self.finish()
        self.capture.pair = self.pair
        self.capture.write_mars_parts(fault='momentum_matrix')
        self.capture.record('mars')
        self.capture.pair = self.baseline
        code, result = self.comparison()
        self.assertEqual(code, 0, result)
        self.assertEqual(result['outcome'], 'upstream_stage_mismatch')

    def test_wrong_baseline_policy_and_history_rejected(self):
        record = startup.read_json(self.baseline / 'pair.json')
        with patch.object(startup, 'pair_inputs', return_value=dict(record, steps=20, first_step_audit=False)):
            with self.assertRaisesRegex(startup.EvidenceError, 'baseline_first_step_required'):
                self.prepare()
        original = startup.read_json
        def without_reference_policy(path):
            value = original(path)
            if path.name == 'case.json':
                value['pressure_linear_policy'] = 'mars'
            return value
        with patch.object(startup, 'pair_inputs', return_value=record), \
             patch.object(startup, 'read_json', side_effect=without_reference_policy):
            with self.assertRaisesRegex(startup.EvidenceError, 'baseline_reference_policy_required'):
                self.prepare()

    def test_current_solver_environment_override_prevents_launch(self):
        self.prepare()
        record = startup.verified_launch(self.baseline, 'mars')
        with patch.object(startup, 'runtime_libraries', return_value=record['libraries']), \
             patch.dict(os.environ, {'MARS_HYPRE_ABSTOL': '999'}, clear=True), \
             patch.object(startup, 'launch') as launch:
            with self.assertRaisesRegex(startup.EvidenceError, 'solver_environment_changed'):
                pressure.run(self.pair, 'mars')
        launch.assert_not_called()

    def test_launch_reuses_captured_launcher_and_rank_count(self):
        self.prepare()
        record = startup.verified_launch(self.baseline, 'openaccel')
        environment = {'OMP_NUM_THREADS': '1', 'OMP_PROC_BIND': 'close', 'OMP_PLACES': 'cores'}
        # Synthetic records do not contain reference OMP defaults; supply them for this launch test.
        record['environment'] = environment
        with patch.object(pressure, 'check_inputs', return_value=(self.baseline, {'openaccel': record})), \
             patch.object(startup, 'runtime_libraries', return_value=record['libraries']), \
             patch.dict(os.environ, {}, clear=True), patch.object(startup, 'launch') as launch, \
             patch.object(startup, 'verified_launch', return_value=record), patch.object(pressure, 'runtime_matches'):
            pressure.run(self.pair, 'openaccel')
        launch.assert_called_once_with(self.pair, 'openaccel', Path(record['executable']), 2, ['launcher'],
                                       environment=environment)

    def test_changed_runtime_prevents_launch(self):
        self.prepare()
        with patch.object(startup, 'runtime_libraries', return_value={'changed': 'b'*64}), \
             patch.object(startup, 'launch') as launch:
            with self.assertRaisesRegex(startup.EvidenceError, 'libraries_changed'):
                pressure.run(self.pair, 'mars')
        launch.assert_not_called()

    def test_environment_report_never_exports_unknown_names_or_values(self):
        current = {'MARS_HYPRE_ABSTOL': 'PRIVATE', 'MARS_PRIVATE_NAME': 'PRIVATE',
                   'CUDA_VISIBLE_DEVICES': 'PRIVATE', 'PATH': 'PRIVATE'}
        result = pressure.environment_report(pressure.environment_changes(current, {}))
        self.assertEqual(result['changed_known_names'], ['CUDA_VISIBLE_DEVICES', 'MARS_HYPRE_ABSTOL'])
        self.assertTrue(result['other_changed_names_present'])
        self.assertTrue(result['protected_or_unknown_changes_present'])
        self.assertNotIn('PRIVATE', json.dumps(result))

    def test_restore_sets_removes_and_preserves_protected_environment(self):
        current = {'MARS_HYPRE_ABSTOL': 'changed', 'MARS_AMG_RELAX': '18',
                   'CUDA_VISIBLE_DEVICES': '3', 'MARS_SIGNAL_DIR': '/current', 'LD_LIBRARY_PATH': '/current/lib'}
        recorded = {'MARS_HYPRE_ABSTOL': '0', 'MARS_HYPRE_FLEXGMRES': '1',
                    'CUDA_VISIBLE_DEVICES': '0', 'MARS_SIGNAL_DIR': '/old', 'LD_LIBRARY_PATH': '/old/lib'}
        expected = dict(current, MARS_HYPRE_ABSTOL='0', MARS_HYPRE_FLEXGMRES='1')
        expected.pop('MARS_AMG_RELAX')
        saved = dict(current)
        self.assertEqual(pressure.restore_environment(current, recorded), expected)
        self.assertEqual(current, saved)

    def test_restore_launch_passes_only_restored_child_environment(self):
        self.prepare()
        record = startup.verified_launch(self.baseline, 'mars')
        record['environment'] = {'MARS_HYPRE_ABSTOL': '0', 'MARS_HYPRE_FLEXGMRES': '1'}
        parent = {'MARS_HYPRE_ABSTOL': '999', 'MARS_AMG_RELAX': '18', 'PATH': '/parent'}
        public = {}
        with patch.object(pressure, 'check_inputs', return_value=(self.baseline, {'mars': record})), \
             patch.object(startup, 'runtime_libraries', return_value=record['libraries']), \
             patch.dict(os.environ, parent, clear=True), patch.object(startup, 'launch') as launch, \
             patch.object(startup, 'verified_launch', return_value=record), patch.object(pressure, 'runtime_matches'):
            pressure.run(self.pair, 'mars', True, public)
            self.assertEqual(dict(os.environ), parent)
        self.assertEqual(launch.call_args[1]['environment'], dict(record['environment'], PATH='/parent'))
        self.assertEqual(public['restored_solver_controls'],
                         ['MARS_AMG_RELAX', 'MARS_HYPRE_ABSTOL', 'MARS_HYPRE_FLEXGMRES'])
        self.assertNotIn('999', json.dumps(public))

    def test_restore_still_rejects_gpu_binding_unknown_names_and_options(self):
        self.prepare()
        record = startup.verified_launch(self.baseline, 'mars')
        for name in ('CUDA_VISIBLE_DEVICES', 'MARS_PRIVATE_NAME', 'PETSC_OPTIONS', 'MARS_VERBOSE_MESH'):
            with self.subTest(name=name), patch.object(startup, 'runtime_libraries', return_value=record['libraries']), \
                 patch.dict(os.environ, {name: 'PRIVATE'}, clear=True), patch.object(startup, 'launch') as launch:
                public = {}
                with self.assertRaisesRegex(startup.EvidenceError, 'solver_environment_changed'):
                    pressure.run(self.pair, 'mars', True, public)
                self.assertNotIn('PRIVATE', json.dumps(public))
                launch.assert_not_called()

    def test_restore_is_explicit_and_mars_only(self):
        self.prepare()
        record = startup.verified_launch(self.baseline, 'openaccel')
        with patch.object(startup, 'runtime_libraries', return_value=record['libraries']), \
             patch.dict(os.environ, {}, clear=True), patch.object(startup, 'launch') as launch:
            with self.assertRaisesRegex(startup.EvidenceError, 'environment_restore_requires_mars'):
                pressure.run(self.pair, 'openaccel', True)
        launch.assert_not_called()

    def test_cli_mismatch_identifies_safe_names_without_launch(self):
        self.prepare()
        record = startup.verified_launch(self.baseline, 'mars')
        with patch.object(startup, 'runtime_libraries', return_value=record['libraries']), \
             patch.dict(os.environ, {'MARS_HYPRE_ABSTOL': 'PRIVATE', 'MARS_PRIVATE_NAME': 'PRIVATE'}, clear=True), \
             patch.object(startup, 'launch') as launch:
            code = self.invoke('run', '--pair', str(self.pair), '--solver', 'mars', '--output', str(self.output))
        self.assertEqual(code, 1)
        text = self.output.read_text()
        self.assertNotIn('PRIVATE', text)
        report = json.loads(text)
        self.assertEqual(report['failed_check'], 'solver_environment_changed')
        self.assertEqual(report['environment_check']['changed_known_names'], ['MARS_HYPRE_ABSTOL'])
        self.assertTrue(report['environment_check']['other_changed_names_present'])
        launch.assert_not_called()

    def test_solver_rejection_has_public_diagnostics_and_cannot_pass(self):
        self.prepare()
        (self.pair / 'mars').mkdir()
        (self.pair / 'mars/run.log').write_text('ERROR: pressure correction at SIMPLE iteration 1: linear solve failed PRIVATE\n')
        (self.pair / 'mars/run.exit').write_text('1\n')
        with patch.object(pressure, 'run', side_effect=startup.EvidenceError('launcher_exit')):
            code = self.invoke('run', '--pair', str(self.pair), '--solver', 'mars', '--output', str(self.output))
        self.assertEqual(code, 1)
        text = self.output.read_text()
        self.assertNotIn('PRIVATE', text)
        result = json.loads(text)
        self.assertEqual(result['outcome'], 'inconclusive')
        self.assertEqual(result['launch_diagnostics']['process_exit_code'], 1)
        self.assertEqual(self.invoke('compare', '--pair', str(self.pair), '--detail-dir', str(self.details),
                                    '--output', str(self.root / 'comparison.json')), 1)

    def test_repeated_launch_command_cannot_overwrite_report_or_launch(self):
        self.output.write_text('preserve')
        with patch.object(pressure, 'run') as run:
            self.assertEqual(self.invoke('run', '--pair', str(self.pair), '--solver', 'mars',
                                        '--output', str(self.output)), 1)
        run.assert_not_called()
        self.assertEqual(self.output.read_text(), 'preserve')


if __name__ == '__main__':
    unittest.main()
