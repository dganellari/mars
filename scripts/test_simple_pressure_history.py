"""Twenty-step comparisons on synthetic fields; no private case data."""
import copy
import csv
import io
import json
import os
import shutil
import unittest
from unittest.mock import patch

import numpy as np
from netCDF4 import Dataset
import yaml

import simple_pressure_profile as profile
import simple_startup_probe as startup
from simple_snapshot_compare import options
import test_simple_pressure_profile as fixtures
from test_simple_pressure_log import report


class HistoryTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.ProfileTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        f = self.fixture
        self.pair, self.root = f.pair, f.root

    def prepare(self):
        f = self.fixture
        with patch.object(startup, 'runtime_libraries', return_value=f.libs):
            profile.prepare(f.baseline, f.probe, f.exe, self.pair, history=True, both_linear_solvers=f.both_linear_solvers)

    def record(self, solver):
        f = self.fixture
        old_pair = f.capture.pair
        f.capture.pair = self.pair
        f.capture.record(solver)
        f.capture.pair = old_pair
        old = startup.verified_launch(f.baseline, solver)
        for name in ('launch-start.json', 'launch.json'):
            path = self.pair / ('mars' if solver == 'mars' else 'reference') / name
            record = startup.read_json(path)
            exe = f.exe if solver == 'mars' else old['executable']
            record.update(executable=str(exe), executable_sha256=startup.digest(exe), libraries=f.libs,
                          command=['launcher', str(exe)] + startup.solver_arguments(self.pair, solver))
            path.write_text(json.dumps(record))

    def finish(self):
        self.prepare()
        f = self.fixture
        mars = self.pair / 'mars'
        mars.mkdir()
        for step in range(21):
            shutil.copyfile(f.baseline / ('mars/flow-step-{}-fields.csv'.format(min(step, 1))),
                            mars / ('flow-step-{}-fields.csv'.format(step)))
        shutil.copyfile(f.baseline / 'mars/flow-fields.csv', mars / 'flow-fields.csv')
        text = (f.baseline / 'mars/run.log').read_text().replace('iterations=1 ', 'iterations=20 ')
        (mars / 'run.log').write_text(text + '\npressure_solver_profile=explicit_gpu\npressure_refinement=0 pressure_expansion=0\n')
        (mars / 'run.exit').write_text('2\n')
        with (f.baseline / 'mars/flow-metrics.csv').open() as stream:
            rows = list(csv.DictReader(stream))
        with (mars / 'flow-metrics.csv').open('w') as stream:
            writer = csv.DictWriter(stream, list(rows[0]))
            writer.writeheader()
            for step in range(21):
                writer.writerow(dict(rows[min(step, 1)], iteration=step))
        reference = self.pair / 'reference'
        (reference / 'run.log').write_text('\n'.join('Iter = ' + str(i) for i in range(1, 21)) + '\nSimulation is complete\n')
        (reference / 'run.exit').write_text('0\n')
        for source in (f.baseline / 'reference').glob('results.e.*'):
            with Dataset(str(source)) as old, Dataset(str(reference / source.name), 'w') as new:
                for name, dimension in old.dimensions.items():
                    new.createDimension(name, 21 if name == 'time_step' else len(dimension))
                for name, variable in old.variables.items():
                    target = new.createVariable(name, variable.dtype, variable.dimensions)
                    if name == 'time_whole':
                        target[:] = np.arange(21)
                    elif variable.dimensions[0] == 'time_step':
                        target[0] = variable[0]
                        for step in range(1, 21): target[step] = variable[1]
                    else:
                        target[:] = variable[:]
        _, _, values, _ = profile.check_inputs(self.pair)
        actual = dict(values, hypre_release=23300, effective_levels=3, effective_relax_3=18)
        actual['effective_relax_1'] = actual.pop('relax_down')
        actual['effective_relax_2'] = actual.pop('relax_up')
        for step in range(1, 21):
            for rank in range(2):
                (mars / ('flow-pressure-step-{}-rank-{}.settings'.format(step, rank))).write_text(profile.text_profile(actual))
        if f.both_linear_solvers:
            momentum, _, _ = profile.momentum_data(f.baseline, f.probe, startup.verified_launch(f.baseline, 'openaccel'))
            actual = dict(momentum, hypre_release=23300, effective_levels=3, effective_relax_3=18)
            actual['effective_relax_1'] = actual.pop('relax_down')
            actual['effective_relax_2'] = actual.pop('relax_up')
            for step in range(1,21):
                for rank in range(2):
                    (mars / ('flow-momentum-step-{}-rank-{}.settings'.format(step,rank))).write_text(profile.text_profile(actual))
            log=mars/'run.log'
            log.write_text(log.read_text() + 'momentum_solver_profile=explicit_gpu\n'
                + 'momentum_linear_rtol={} momentum_linear_atol={}\n'.format(momentum['rtol'],momentum['atol']))
        self.record('mars')
        self.record('openaccel')

    def compare(self):
        details = self.root / 'history-details'
        details.mkdir()
        public = {}
        profile.compare(self.pair, public, details)
        return public

    def test_both_profiles_verified_at_every_step_and_rank(self):
        self.fixture.enable_momentum(); self.finish()
        public=self.compare()
        self.assertTrue(public['recorded_gpu_momentum_settings_verified'])
        self.assertTrue(public['original_momentum_targets_preserved'])
        self.assertFalse(public['recovery_enabled'])
        self.assertTrue(public['reference_solver_launched'])

    def test_missing_momentum_history_part_fails(self):
        self.fixture.enable_momentum(); self.finish()
        (self.pair/'mars/flow-momentum-step-7-rank-1.settings').unlink()
        self.record('mars')
        with self.assertRaisesRegex(startup.EvidenceError,'runtime_profile_mismatch'): self.compare()

    def test_original_controls_preserved_and_no_recovery_or_audit(self):
        f = self.fixture
        before = {str(p): startup.digest(p) for p in f.baseline.rglob('*') if p.is_file()}
        self.prepare()
        record = startup.pair_inputs(self.pair)
        self.assertEqual(record['steps'], 20)
        self.assertFalse(record['first_step_audit'])
        args = options(startup.solver_arguments(self.pair, 'mars'))
        self.assertEqual((args['--iterations'], args['--snapshot-iterations']), ('20', '20'))
        self.assertEqual((args['--pressure-expansion'], args['--pressure-refinement']), ('0', '0'))
        self.assertNotIn('--first-step-audit', args)
        self.assertEqual((args['--pressure-linear-rtol'], args['--pressure-linear-atol']), ('0.0001', '1e-10'))
        old = startup.load_deck((f.baseline / 'reference/input.i').read_bytes())
        new = startup.load_deck((self.pair / 'reference/input.i').read_bytes())
        for equation in ('coupled_navier_stokes', 'pressure_correction'):
            a, b = [profile.pressure.resolved_solver(deck, equation) for deck in (old, new)]
            a['write_system'] = False
            self.assertEqual(a, b)
        self.assertEqual(before, {str(p): startup.digest(p) for p in f.baseline.rglob('*') if p.is_file()})
        self.assertFalse((self.pair / 'reference/launch.json').exists())

    def test_complete_history_checks_all_steps_and_fresh_launches(self):
        self.finish()
        public = self.compare()
        self.assertEqual(public['outcome'], 'twenty_step_fields_match')
        self.assertEqual(len(public['snapshot_agreement']), 21)
        self.assertTrue(public['comparison']['fresh_launch_records_verified'])
        self.assertTrue(public['comparison']['initial_field_parity_verified'])
        self.assertTrue(public['reference_solver_launched'])
        self.assertFalse(public['recovery_enabled'])
        self.assertFalse(public['nonlinear_convergence_verified'])
        self.assertFalse(public['identical_linear_solvers_verified'])
        self.assertFalse(public['pressure_residuals_recomputed_at_every_step'])
        self.assertNotIn(str(self.root), json.dumps(public))
        self.assertNotIn('first_step_stage_matches', public['comparison'])

    def pressure_log(self, end='1.23e-08'):
        path = self.pair / 'reference/run.log'
        path.write_text(''.join(report(step=i, end=end) for i in range(1, 21)) + 'Simulation is complete\n')
        self.record('openaccel')
        return path

    def test_saved_pressure_inspection_does_not_launch_or_load_runtime(self):
        self.finish()
        self.pressure_log(end='1.23e-04')
        public = {}
        with patch.object(startup, 'launch', side_effect=AssertionError('no launch')), \
             patch.object(startup, 'runtime_libraries', side_effect=AssertionError('no runtime query')):
            profile.inspect_linear(self.pair, public)
        self.assertEqual(public['comparison_status'], 'completed')
        self.assertTrue(public['reference_pressure_log_checks']['any_printed_target_exceeded'])
        self.assertFalse(public['solver_launched'])
        self.assertFalse(public['nonlinear_convergence_verified'])
        self.assertFalse(public['openaccel_field_parity_verified'])
        self.assertNotIn(str(self.root), json.dumps(public))

    def test_inspection_checks_log_hash_and_all_pressure_settings(self):
        self.finish()
        path = self.pressure_log()
        original = path.read_text()
        path.write_text(original + 'changed\n')
        with self.assertRaises(ValueError): profile.inspect_linear(self.pair, {})
        path.write_text(original)
        (self.pair / 'mars/flow-pressure-step-20-rank-1.settings').unlink()
        self.record('mars')
        with self.assertRaisesRegex(startup.EvidenceError, 'runtime_profile_mismatch'):
            profile.inspect_linear(self.pair, {})

    def test_inspection_cli_preserves_evidence_and_existing_output(self):
        self.finish()
        self.pressure_log()
        before = {str(p): startup.digest(p) for p in self.pair.rglob('*') if p.is_file()}
        output = self.root / 'linear-public.json'
        args = ['inspect-linear', '--pair', str(self.pair), '--output', str(output)]
        self.assertEqual(profile.main(args), 0)
        public = startup.read_json(output)
        self.assertTrue(public['reference_pressure_log_checks']['all_printed_relative_targets_met'])
        saved = output.read_bytes()
        self.assertEqual(profile.main(args), 1)
        self.assertEqual(saved, output.read_bytes())
        self.assertEqual(before, {str(p): startup.digest(p) for p in self.pair.rglob('*') if p.is_file()})

    def test_inspection_cli_unknown_format_fails_without_private_text(self):
        self.finish()
        self.pressure_log(end='PRIVATE_ERROR_TEXT')
        output = self.root / 'linear-invalid.json'
        with patch('sys.stdout', new_callable=io.StringIO) as out, patch('sys.stderr', new_callable=io.StringIO) as err:
            code = profile.main(['inspect-linear', '--pair', str(self.pair), '--output', str(output)])
        public = startup.read_json(output)
        self.assertEqual(code, 1)
        self.assertEqual(public['comparison_status'], 'invalid_evidence')
        self.assertEqual(public['failed_check'], 'reference_pressure_log')
        self.assertNotIn('PRIVATE_ERROR_TEXT', out.getvalue() + err.getvalue() + output.read_text())

    def test_inspection_requires_original_twenty_step_history(self):
        self.prepare()
        values = profile.check_inputs(self.pair)
        identity = dict(values[-1])
        del identity['kind']
        with patch.object(profile, 'check_inputs', return_value=values[:-1] + (identity,)), \
             self.assertRaisesRegex(startup.EvidenceError, 'history_required'):
            profile.inspect_linear(self.pair, {})

    def test_late_field_difference_is_reported_without_exact_field_values(self):
        self.finish()
        path = self.pair / 'mars/flow-step-13-fields.csv'
        with path.open() as stream: rows = list(csv.DictReader(stream))
        rows[0]['p'] = str(float(rows[0]['p']) + 123.456789)
        with path.open('w') as stream:
            writer = csv.DictWriter(stream, list(rows[0])); writer.writeheader(); writer.writerows(rows)
        self.record('mars')
        public = self.compare()
        self.assertEqual(public['outcome'], 'twenty_step_fields_differ')
        self.assertEqual(public['comparison']['first_differing_iteration'], 13)
        self.assertEqual(public['snapshot_agreement'][13]['pressure_max_scaled_band'], 'over_5_percent')
        self.assertEqual(public['snapshot_agreement'][14]['pressure_max_scaled_band'], 'within_1e_minus_5')
        self.assertNotIn('123.456789', json.dumps(public))

    def test_missing_or_changed_later_pressure_settings_fail(self):
        self.finish()
        path = self.pair / 'mars/flow-pressure-step-19-rank-1.settings'
        original = path.read_text()
        for text in (None, original.replace('kdim 5\n', 'kdim 100\n')):
            if text is None: path.unlink()
            else: path.write_text(text)
            self.record('mars')
            with self.subTest(missing=text is None), self.assertRaisesRegex(startup.EvidenceError, 'runtime_profile_mismatch'):
                profile.compare(self.pair, {}, self.root)

    def test_missing_last_snapshot_and_reused_one_step_reference_fail(self):
        self.finish()
        path = self.pair / 'mars/flow-step-20-fields.csv'
        saved = path.read_bytes(); path.unlink(); self.record('mars')
        with self.assertRaises((ValueError, OSError)): profile.compare(self.pair, {}, self.root)
        path.write_bytes(saved); self.record('mars')
        for name in ('launch.json', 'launch-start.json'):
            shutil.copyfile(self.fixture.baseline / 'reference' / name, self.pair / 'reference' / name)
        with self.assertRaises(ValueError): profile.compare(self.pair, {}, self.root)

    def test_modified_deck_rejected_even_with_updated_local_hashes(self):
        self.prepare()
        path = self.pair / 'reference/input.i'
        deck = startup.load_deck(path.read_bytes())
        deck['simulation']['solver']['solver_control']['basic_settings']['convergence_controls']['max_iterations'] = 21
        path.write_text(yaml.safe_dump(deck))
        case = startup.read_json(self.pair / 'case.json'); case['deck_sha256'] = startup.digest(path)
        (self.pair / 'case.json').write_text(json.dumps(case))
        record = startup.read_json(self.pair / 'pair.json')
        record.update(deck_sha256=case['deck_sha256'], case_sha256=startup.digest(self.pair / 'case.json'))
        (self.pair / 'pair.json').write_text(json.dumps(record))
        with self.assertRaisesRegex(startup.EvidenceError, 'profile_identity'): profile.check_inputs(self.pair)

    def test_reference_run_uses_its_own_saved_runtime_and_launcher(self):
        self.prepare()
        baseline, launches, values, identity = profile.check_inputs(self.pair)
        launches = copy.deepcopy(launches)
        os.chmod(launches['openaccel']['executable'], 0o700)
        env = dict(OMP_NUM_THREADS='1', OMP_PROC_BIND='close', OMP_PLACES='cores')
        launches['openaccel']['environment'] = env
        with patch.dict(os.environ, {}, clear=True), \
             patch.object(profile, 'check_inputs', return_value=(baseline, launches, values, identity)), \
             patch.object(startup, 'runtime_libraries', return_value=self.fixture.libs), \
             patch.object(startup, 'launch') as launch, patch.object(profile, 'check_launch'):
            profile.run(self.pair, {}, 'openaccel')
            self.assertEqual(launch.call_args[0][1], 'openaccel')
            self.assertEqual(str(launch.call_args[0][2]), launches['openaccel']['executable'])
            self.assertEqual(launch.call_args[0][4], ['launcher'])
            self.assertEqual(launch.call_args[0][-1], env)

    def test_recovery_activation_cannot_pass_as_original_history(self):
        self.finish()
        path = self.pair / 'mars/run.log'
        path.write_text(path.read_text().replace('pressure_expansion=0', 'pressure_expansion=1'))
        self.record('mars')
        with self.assertRaisesRegex(startup.EvidenceError, 'runtime_recovery_mismatch'):
            profile.compare(self.pair, {}, self.root)

    def test_initial_field_mismatch_is_invalid_not_a_completed_comparison(self):
        self.finish()
        path = self.pair / 'mars/flow-step-0-fields.csv'
        with path.open() as stream: rows = list(csv.DictReader(stream))
        rows[0]['p'] = '1'
        with path.open('w') as stream:
            writer = csv.DictWriter(stream, list(rows[0])); writer.writeheader(); writer.writerows(rows)
        self.record('mars')
        output = self.root / 'initial-public.json'
        code = profile.main(['compare', '--pair', str(self.pair), '--detail-dir', str(self.root / 'details'),
                             '--output', str(output)])
        public = startup.read_json(output)
        self.assertEqual(code, 1)
        self.assertEqual(public['comparison_status'], 'invalid_evidence')
        self.assertEqual(public['failed_check'], 'initial_fields_mismatch')

    def test_unknown_history_kind_rejected(self):
        self.prepare()
        path = self.pair / 'pair.json'
        record = startup.read_json(path)
        record['pressure_profile']['kind'] = 'different_experiment'
        path.write_text(json.dumps(record))
        with self.assertRaises(ValueError): profile.check_inputs(self.pair)

    def test_failed_launch_is_not_capture_complete_and_message_stays_private(self):
        self.prepare()
        output = self.root / 'failed-public.json'
        with patch.dict(os.environ, {}, clear=True), \
             patch.object(startup, 'runtime_libraries', return_value=self.fixture.libs), \
             patch.object(startup, 'launch', side_effect=startup.EvidenceError('launcher_exit')), \
             patch.object(profile.pressure, 'failure_summary', return_value={'linear_stage': 'pressure'}):
            code = profile.main(['run', '--pair', str(self.pair), '--output', str(output)])
        public = startup.read_json(output)
        self.assertEqual(code, 1)
        self.assertEqual(public['comparison_status'], 'invalid_evidence')
        self.assertEqual(public['failed_check'], 'launcher_exit')
        self.assertEqual(public['launch_diagnostics']['linear_stage'], 'pressure')
        self.assertNotIn(str(self.root), output.read_text())


if __name__ == '__main__': unittest.main()
