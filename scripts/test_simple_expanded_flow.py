import contextlib
import io
import json
import os
from pathlib import Path
import stat
import sys
import tempfile
import unittest
from unittest import mock

import simple_expanded_flow as flow
from test_simple_convergence_summary import fixture, serialize


class FlowTests(unittest.TestCase):
    def saved(self):
        return {'--mesh': '/private/SECRET.exo', '--inlet-ss': 'SECRET',
                '--iterations': '20', '--pressure-linear-rtol': '1e-10',
                '--pressure-linear-atol': '0', '--snapshot-iterations': '20',
                '--field-output': 'distributed', '--output-prefix': '/private/old',
                '--residual-tol': '1e-6', '--mass-tol': '1e-6', '--change-tol': '1e-6',
                '--rho': '987', '--alpha-u': '.3', '--report-every': '1'}

    def test_only_authorized_arguments_change(self):
        old = self.saved(); original = dict(old)
        new = flow.flow_arguments(old, {'rtol': 1e-10, 'atol': 0}, Path('/private/new'))
        allowed = {'--pressure-expansion', '--snapshot-iterations', '--field-output',
                   '--output-prefix', '--pressure-solver-profile'}
        self.assertEqual({k: v for k, v in new.items() if k not in allowed},
                         {k: v for k, v in old.items() if k not in allowed})
        self.assertEqual(old, original)
        self.assertEqual(new['--pressure-expansion'], '1')
        self.assertEqual(new['--field-output'], 'none')
        self.assertEqual(new['--snapshot-iterations'], '0')
        self.assertNotIn('--pressure-failure-capture', new)

    def test_different_target_rejected(self):
        for values in ({'rtol': 1e-9, 'atol': 0}, {'rtol': 1e-10, 'atol': 1e-8}):
            with self.subTest(values=values), self.assertRaises(Exception):
                flow.flow_arguments(self.saved(), values, Path('/unused'))

    def test_longer_budget_changes_no_numerical_control(self):
        saved = self.saved()
        baseline = flow.flow_arguments(saved, {'rtol': 1e-10, 'atol': 0}, Path('/new'))
        longer = flow.flow_arguments(saved, {'rtol': 1e-10, 'atol': 0}, Path('/new'), 2000)
        self.assertEqual(longer, dict(baseline, **{'--iterations': '2000'}))
        self.assertEqual(saved['--iterations'], '20')
        for invalid in (0, 19, -1, 2**31, 20.5, True):
            with self.subTest(invalid=invalid), self.assertRaises(Exception):
                flow.flow_arguments(saved, {'rtol': 1e-10, 'atol': 0}, Path('/new'), invalid)

    def test_walltime_preserves_rank_binding_and_other_resources(self):
        for time in (['--time=00:03:00'], ['--time', '00:03:00'], ['-t00:03:00'], ['-t', '00:03:00']):
            old = ['srun', '--account', 'test'] + time + ['--nodes=1', '--ntasks-per-node=4',
                   '--export=ALL,MPICH_GPU_SUPPORT_ENABLED=1', '--kill-on-bad-exit=1', '/binding']
            result = flow.extended_launcher(old, '01:00:00')
            self.assertEqual(result[:3], old[:3])
            self.assertEqual(result[-5:], old[-5:])
            self.assertNotIn('00:03:00', ' '.join(result))
            self.assertIn('01:00:00', ' '.join(result))
            self.assertEqual(flow.extended_launcher(old, None), old)
        old = ['srun', '--nodes', '1', '/binding', '--time=do-not-change']
        self.assertEqual(flow.extended_launcher(old, '01:00:00'),
                         ['srun', '--time=01:00:00'] + old[1:])
        for old, time in ((['mpirun', '-np', '4'], '01:00:00'),
                          (['srun', '--time=00:01:00', '-t', '00:02:00'], '01:00:00'),
                          (['srun', '--unknown', 'value'], '01:00:00'),
                          (['srun'], '1:99:00'), (['srun'], '00:00:00'),
                          (['srun'], '01:00:00; SECRET')):
            with self.subTest(old=old, time=time), self.assertRaises(Exception):
                flow.extended_launcher(old, time)

    def test_incompatible_saved_modes_rejected(self):
        for option in ('--pressure-refinement', '--pressure-expansion', '--first-step-audit', '--setup-only'):
            saved = self.saved(); saved[option] = '1'
            with self.subTest(option=option), self.assertRaises(Exception):
                flow.flow_arguments(saved, {'rtol': 1e-10, 'atol': 0}, Path('/unused'))

    def launch(self, root, source, steps=20, ranks=4):
        executable = root / 'source'; executable.write_text('unchanged')
        record = dict(command=[sys.executable, '-c', source], steps=steps, ranks=ranks,
                      inputs=flow.replay.hashes([executable]))
        flow.replay.startup.write_json(root / 'launch-start.json', record)
        public = dict(schema=flow.SCHEMA, comparison_status='invalid_evidence')
        environment = dict(os.environ, EXPANDED_FLOW_TEST='child-only')
        with mock.patch.object(flow, 'prepare', return_value=(record, environment)):
            code = flow.run(Path('/unused'), Path('/unused'), executable, root, public)
        return code, public

    def test_completed_budget_and_early_convergence(self):
        for label, steps, code, converged in (
                ('NOT CONVERGED: iteration limit', 20, 2, False), ('CONVERGED', 7, 0, True)):
            with self.subTest(label=label), tempfile.TemporaryDirectory() as name:
                source = ('import os,sys; assert os.environ["EXPANDED_FLOW_TEST"]=="child-only"; '
                          'print("SECRET"); print({!r}); sys.exit({})').format(
                              '{} iterations={} ranks=4 exchange_rounds=81'.format(label, steps), code)
                status, result = self.launch(Path(name), source)
                self.assertEqual(status, 0)
                self.assertTrue(result['short_run_completed'])
                self.assertEqual(result['nonlinear_convergence_reported'], converged)
                self.assertNotIn('SECRET', json.dumps(result))
                self.assertNotIn('EXPANDED_FLOW_TEST', os.environ)

    def test_false_completions_rejected(self):
        for text, code in (
                ('NOT CONVERGED: iteration limit iterations=19 ranks=4', 2),
                ('NOT CONVERGED: iteration limit iterations=20 ranks=2', 2),
                ('CONVERGED iterations=21 ranks=4', 0),
                ('CONVERGED iterations=20 ranks=4', 2),
                ('CONVERGED iterations=20 ranks=4\nCONVERGED iterations=20 ranks=4', 0),
                ('ERROR: pressure correction failed SECRET', 255), ('SECRET', 0)):
            with self.subTest(text=text), tempfile.TemporaryDirectory() as name:
                source = 'import sys; print({!r}); sys.exit({})'.format(text, code)
                status, public = self.launch(Path(name), source)
                self.assertEqual(status, 1)
                self.assertFalse(public['short_run_completed'])
                self.assertNotIn('SECRET', json.dumps(public))

    def test_changed_input_after_run_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            source = 'from pathlib import Path; Path({!r}).write_text("changed")'.format(str(root / 'source'))
            with self.assertRaises(Exception):
                self.launch(root, source)

    def test_preflight_failure_private_and_no_launch(self):
        with tempfile.TemporaryDirectory() as name:
            output = Path(name) / 'new'
            with mock.patch.object(flow.replay, 'captured_record', side_effect=ValueError('SECRET')), mock.patch.object(flow.diagnostics, 'capture') as launch:
                stream = io.StringIO()
                with contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream):
                    code = flow.main(['--capture-run', '/private/SECRET', '--gpu-profile-pair', '/private/SECRET',
                                      '--executable', '/private/SECRET', '--output-dir', str(output)])
                self.assertEqual(code, 1); launch.assert_not_called()
                public = json.loads((output / 'public.json').read_text())
                self.assertEqual(public['failed_check'], 'saved_capture')
                self.assertNotIn('SECRET', stream.getvalue() + json.dumps(public))
                self.assertEqual(stat.S_IMODE(output.stat().st_mode), 0o700)
                self.assertEqual(stat.S_IMODE((output / 'public.json').stat().st_mode), 0o600)

    def test_inspect_existing_recovery_without_launch(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name); saved = root / 'saved'; saved.mkdir()
            log = ('SIMPLE Tet4, 4 ranks (ElementDomain/cstone), upwind, laminar\n'
                   '[HypreGMRES] rejected: backend=GMRES\n'
                   '[simple-pressure-expansion] rounds=2 correction_iterations=10 hypre_passed=1 mars_passed=1\n'
                   'NOT CONVERGED: iteration limit iterations=20 ranks=4 exchange_rounds=81\nSECRET')
            (saved / 'run.log').write_text(log)
            (saved / 'run.exit').write_text('2\n')
            (saved / 'diagnostics.json').write_text('{"run_status":"failed"}')
            dependency = root / 'dependency'; dependency.write_text('unchanged')
            start = dict(schema=flow.SCHEMA, steps=20, ranks=4, pressure_target_preserved=True,
                         saved_gpu_profile_verified=True, inputs=flow.replay.hashes([dependency]))
            flow.replay.startup.write_json(saved / 'launch-start.json', start)
            record = dict(start, exit_code=2, outputs=flow.replay.hashes(list(saved.iterdir())))
            flow.replay.startup.write_json(saved / 'launch.json', record)
            before = flow.replay.hashes(list(saved.iterdir()))
            with mock.patch.object(flow.diagnostics, 'capture') as launch:
                with contextlib.redirect_stdout(io.StringIO()):
                    code = flow.main(['--inspect-run', str(saved), '--output-dir', str(root / 'inspection')])
                launch.assert_not_called()
            result = json.loads((root / 'inspection/public.json').read_text())
            self.assertEqual(code, 0)
            self.assertTrue(result['short_run_completed'])
            self.assertFalse(result['nonlinear_convergence_reported'])
            self.assertFalse(result['solver_launched'])
            self.assertNotIn('SECRET', json.dumps(result))
            self.assertEqual(before, flow.replay.hashes(list(saved.iterdir())))
            for changed in (saved / 'run.log', dependency):
                old = changed.read_text(); changed.write_text(old + 'changed')
                with self.assertRaises(Exception): flow.inspect_run(saved, {})
                changed.write_text(old)

    def long_record(self, root, count=100, status='iteration_limit', value=lambda i: .97**i):
        rows, _ = fixture(count=count, value=value)
        log, metrics, code = serialize(rows, status)
        log = log.replace('PRIVATE boundary',
            '[HypreGMRES] rejected: backend=GMRES\n'
            '[simple-pressure-expansion] rounds=2 correction_iterations=10 hypre_passed=1 mars_passed=1\n'
            'PRIVATE boundary')
        (root / 'run.log').write_text(log)
        (root / 'run.exit').write_text(code)
        (root / 'flow-metrics.csv').write_text(metrics)
        (root / 'diagnostics.json').write_text('{}')
        targets = {key: 1e-6 for key in flow.TARGETS}
        args = dict(self.saved(), **{'--iterations': str(count), '--pressure-expansion': '1',
                                    '--output-prefix': str(root / 'flow')})
        record = dict(schema=flow.SCHEMA, record_version=2, source_steps=20, steps=count, ranks=4,
            nonlinear_targets=targets, inputs={}, pressure_target_preserved=True,
            saved_gpu_profile_verified=True, launcher=['srun'], source_launcher=['srun'],
            time_limit_override=None, command=['srun', '/synthetic-exe'] + [w for kv in args.items() for w in kv])
        flow.replay.startup.write_json(root / 'launch-start.json', record)
        record.update(exit_code=int(code), outputs=flow.replay.hashes(list(root.iterdir())))
        flow.replay.startup.write_json(root / 'launch.json', record)
        public = dict(diagnostics=flow.summarize(log.splitlines(), code))
        return record, public, int(code)

    def test_extended_completion_reports_tail_without_claiming_convergence(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            record, public, code = self.long_record(root)
            self.assertEqual(flow.completion(record, root, code, public), 0)
            self.assertTrue(public['run_completed'])
            self.assertFalse(public['short_run_completed'])
            self.assertFalse(public['nonlinear_convergence_verified'])
            self.assertEqual(public['convergence']['tail_momentum_trend'], 'decreasing')
            self.assertEqual(public['convergence']['targets_source'], 'verified_launch_arguments')
            self.assertNotIn('PRIVATE', json.dumps(public))
            self.assertNotIn('987.654321', json.dumps(public))
            with mock.patch.object(flow.diagnostics, 'capture') as launch:
                inspected = {}
                self.assertEqual(flow.inspect_run(root, inspected), 0)
                launch.assert_not_called()
            self.assertEqual(inspected['convergence'], public['convergence'])

    def test_extended_child_run_records_metrics_and_can_be_reinspected(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            record, _, _ = self.long_record(root)
            log, metrics = (root / 'run.log').read_text(), (root / 'flow-metrics.csv').read_text()
            for filename in ('run.log', 'run.exit', 'flow-metrics.csv', 'diagnostics.json', 'launch.json'):
                (root / filename).unlink()
            start = {key: value for key, value in record.items() if key not in ('exit_code', 'outputs')}
            source = ('import sys; from pathlib import Path; '
                      'Path({!r}).write_text({!r}); print({!r}); sys.exit(2)').format(
                          str(root / 'flow-metrics.csv'), metrics, log)
            start['launcher'] = [sys.executable, '-c', source]
            start['source_launcher'] = list(start['launcher'])
            start['command'] = start['launcher'] + record['command'][1:]
            (root / 'launch-start.json').write_text(json.dumps(start))
            public = {}
            with mock.patch.object(flow, 'prepare', return_value=(start, os.environ.copy())):
                self.assertEqual(flow.run(root, root, root / 'synthetic-exe', root, public, 100), 0)
            self.assertTrue(public['run_completed'])
            self.assertFalse(public['nonlinear_convergence_verified'])
            saved = json.loads((root / 'launch.json').read_text())
            self.assertIn(str((root / 'flow-metrics.csv').resolve()), saved['outputs'])
            inspected = {}
            self.assertEqual(flow.inspect_run(root, inspected), 0)
            self.assertEqual(inspected['convergence'], public['convergence'])

    def test_early_convergence_within_extended_budget(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            record, public, code = self.long_record(root, count=7, status='converged', value=lambda i: 1e-8)
            record['steps'] = 2000
            self.assertEqual(flow.completion(record, root, code, public), 0)
            self.assertTrue(public['nonlinear_convergence_verified'])
            self.assertEqual(public['convergence']['tail_momentum_trend'], 'insufficient_data')

    def test_converged_label_requires_independent_final_checks(self):
        for value, valid in ((lambda i: 1e-8, True), (lambda i: .01, False)):
            with tempfile.TemporaryDirectory() as name:
                root = Path(name)
                record, public, code = self.long_record(root, status='converged', value=value)
                if valid:
                    self.assertEqual(flow.completion(record, root, code, public), 0)
                    self.assertTrue(public['nonlinear_convergence_verified'])
                else:
                    with self.assertRaises(Exception):
                        flow.completion(record, root, code, public)
                    self.assertFalse(public['nonlinear_convergence_verified'])
                    self.assertEqual(public['comparison_status'], 'invalid_evidence')

    def test_missing_or_changed_metrics_and_mismatched_budget_fail(self):
        for change in ('missing', 'changed', 'unbound', 'budget', 'target'):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as name:
                root = Path(name)
                record, public, code = self.long_record(root)
                metrics = root / 'flow-metrics.csv'
                if change == 'missing': metrics.unlink()
                if change == 'changed': metrics.write_text('SECRET')
                if change == 'unbound': del record['outputs'][str(metrics.resolve())]
                if change == 'budget': record['steps'] = 101
                if change == 'target': record['nonlinear_targets']['--residual-tol'] = 1
                if change in ('budget', 'target', 'unbound'):
                    (root / 'launch.json').write_text(json.dumps(record))
                with self.assertRaises(Exception): flow.inspect_run(root, {})


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.pair = self.root / 'pair'; self.pair.mkdir()
        (self.pair / 'reference').mkdir()
        self.capture = self.root / 'capture'; self.capture.mkdir()
        self.output = self.root / 'new'; self.output.mkdir()
        self.exe = self.root / 'new-exe'; self.exe.write_text('new'); self.exe.chmod(0o700)
        self.old_exe = self.root / 'old-exe'; self.old_exe.write_text('old')
        self.mesh = self.root / 'SECRET.mesh'; self.mesh.write_text('synthetic')
        self.library = self.root / 'lib.so'; self.library.write_text('synthetic library')
        for file in (self.pair / 'pair.json', self.pair / 'case.json', self.pair / 'reference/input.i',
                     self.capture / 'capture.json'):
            file.write_text('synthetic record')
        self.saved = FlowTests().saved()
        self.saved['--mesh'] = str(self.mesh)
        args = dict(self.saved, **{'--output-prefix': str(self.capture / 'flow'),
                    '--field-output': 'none', '--snapshot-iterations': '0',
                    '--pressure-failure-capture': str(self.capture / 'system')})
        self.launcher = ['srun', '--ntasks-per-node=4', '/binding']
        self.libraries = flow.replay.hashes([self.library])
        self.record = dict(pair=str(self.pair), ranks=4, launcher=self.launcher,
            command=self.launcher + [str(self.old_exe)] + [w for kv in args.items() for w in kv],
            libraries=self.libraries, environment={'MARS_HYPRE_SPMV_VENDOR': '0'},
            inputs=flow.replay.hashes([self.old_exe, self.mesh, self.pair / 'case.json', self.library]))
        self.values = {'rtol': 1e-10, 'atol': 0}

    def prepare(self, iterations=None, time_limit=None):
        public = {}
        with contextlib.ExitStack() as stack:
            stack.enter_context(mock.patch.dict(os.environ, {}, clear=True))
            stack.enter_context(mock.patch.object(flow.replay, 'captured_record', return_value=self.record))
            stack.enter_context(mock.patch.object(flow.replay.startup, 'pair_inputs',
                                return_value={'steps': 20, 'mesh': str(self.mesh)}))
            stack.enter_context(mock.patch.object(flow.replay.startup, 'solver_arguments',
                                return_value=[w for kv in self.saved.items() for w in kv]))
            stack.enter_context(mock.patch.object(flow.replay, 'gpu_configuration', return_value=(self.values, {})))
            stack.enter_context(mock.patch.object(flow.replay.startup, 'runtime_libraries', return_value=self.libraries))
            result, env = flow.prepare(self.capture, self.root / 'profile', self.exe, self.output, public,
                                       iterations, time_limit)
        return result, env, public

    def test_rebuilt_executable_allowed_and_launcher_preserved(self):
        self.old_exe.write_text('rebuilt in place')
        result, env, public = self.prepare()
        self.assertEqual(result['command'][:len(self.launcher)], self.launcher)
        self.assertEqual(result['command'][len(self.launcher)], str(self.exe))
        args = flow.options(result['command'][len(self.launcher)+1:])
        self.assertEqual(args['--iterations'], '20')
        self.assertEqual(args['--pressure-linear-rtol'], '1e-10')
        self.assertEqual(args['--pressure-linear-atol'], '0')
        self.assertEqual(env['MARS_HYPRE_SPMV_VENDOR'], '0')
        self.assertTrue(public['saved_input_checks']['matched'])
        self.assertTrue(public['pressure_target_preserved'])

    def test_extended_preparation_preserves_source_case_and_records_direct_launch(self):
        before = flow.replay.hashes([self.pair / 'case.json', self.capture / 'capture.json', self.mesh])
        result, env, public = self.prepare(2000, '01:00:00')
        self.assertEqual(result['steps'], 2000)
        self.assertEqual(result['source_steps'], 20)
        self.assertEqual(result['source_launcher'], self.launcher)
        self.assertEqual(result['launcher'], ['srun', '--time=01:00:00'] + self.launcher[1:])
        self.assertTrue(public['iteration_budget_extended'])
        self.assertFalse(public['launcher_preserved'])
        self.assertTrue(public['launcher_resources_preserved_except_time'])
        self.assertEqual(flow.shlex.split((self.output / 'launch-command.txt').read_text()), result['command'])
        self.assertEqual(before, flow.replay.hashes([self.pair / 'case.json', self.capture / 'capture.json', self.mesh]))

    def test_changed_capture_case_rejected(self):
        (self.pair / 'case.json').write_text('changed')
        with self.assertRaises(Exception): self.prepare()

    def test_changed_library_rejected(self):
        self.library.write_text('changed')
        with self.assertRaises(Exception): self.prepare()

    def test_wrong_captured_command_rejected(self):
        self.record['command'][-1] = '/different-capture'
        with self.assertRaises(Exception): self.prepare()

    def test_profile_target_change_rejected(self):
        self.values['atol'] = 1e-3
        with self.assertRaises(Exception): self.prepare()


if __name__ == '__main__':
    unittest.main()
