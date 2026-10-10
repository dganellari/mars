"""Synthetic systems only: capture identity, solver settings and residual verdicts."""
import contextlib
from fractions import Fraction
import io
import json
import math
import os
from pathlib import Path
import random
import shutil
import struct
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import simple_pressure_replay as replay


def write_system(directory, ranks=1, matrix=None, solution=None, rhs=None, candidate=None, atol=0., rtol=1e-10):
    if matrix is None:
        matrix = [[4. if i == j else (-1. if abs(i-j) == 1 else 0.) for j in range(12)] for i in range(12)]
    n = len(matrix)
    solution = solution or [1.+i/8. for i in range(n)]
    rhs = rhs if rhs is not None else [math.fsum(a*x for a, x in zip(row, solution)) for row in matrix]
    candidate = candidate if candidate is not None else [0.]*n
    directory.mkdir()
    for rank in range(ranks):
        first, last = n*rank//ranks, n*(rank+1)//ranks
        # Solver IDs deliberately differ from local indices, including owned rows.
        mapping = list(reversed(range(n))) + [-1]
        columns, values, offsets = [], [], [0]
        for row in matrix[first:last]:
            for j, a in enumerate(row):
                if a != 0.:
                    columns.append(mapping.index(j)); values.append(a)
            offsets.append(len(columns))
        header = [0x4d53505245535331, 1, rank, ranks, first, last, n, n+1, len(values), 7, 0, 0, 1]
        data = struct.pack('<13Q2d', *(header+[atol, rtol]))
        for fmt, array in (('i', offsets), ('i', columns), ('q', mapping), ('d', values),
                           ('d', rhs[first:last]), ('d', list(reversed(candidate))+[float('nan')])):
            data += struct.pack('<'+str(len(array))+fmt, *array)
        replay.rank_file(directory, rank, '.bin').write_bytes(data)
        replay.rank_file(directory, rank, '.settings').write_text(
            'method 0\nrtol {:.17g}\natol {:.17g}\nmaxiter 200\nrelaxtype 18\ncoarserelax 18\n'.format(rtol, atol))
    (directory/'complete').write_text('mars-pressure-capture-v1\n{}\n'.format(ranks))
    return solution


def write_solution(directory, solution, ranks):
    directory.mkdir()
    for rank in range(ranks):
        owned = solution[len(solution)*rank//ranks:len(solution)*(rank+1)//ranks]
        replay.rank_file(directory, rank, '.solution').write_bytes(struct.pack('<{}d'.format(len(owned)), *owned))
    (directory/'complete').write_text('mars-pressure-replay-v1\n{}\n'.format(ranks))


class ReplayTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        scratch = replay.ROOT/'test-scratch'; scratch.mkdir(exist_ok=True)
        cls.workspace = tempfile.TemporaryDirectory(prefix='pressure-replay-test-', dir=str(scratch))
        if os.environ.get('MARS_TEST_PRESSURE_CHECKER'):
            cls.checker = Path(os.environ['MARS_TEST_PRESSURE_CHECKER']).resolve(strict=True)
        else:
            cls.checker = Path(cls.workspace.name)/'checker'
            subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++17', '-O2', '-ffp-contract=off',
                str(replay.CPP.with_name('pressure_residual_check.cpp')), '-o', str(cls.checker)], check=True)

    @classmethod
    def tearDownClass(cls):
        cls.workspace.cleanup()

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=self.workspace.name)
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.capture = self.root/'capture'; self.result = self.root/'result'

    def check(self, result='-', succeeds=True):
        p = subprocess.run([str(self.checker), str(self.capture), str(result)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(p.returncode == 0, succeeds, p.stderr.decode())
        if succeeds:
            return json.loads(p.stdout)
        self.assertNotIn(str(self.root), p.stderr.decode())

    def test_owned_rows_maps_and_rank_partitions(self):
        for ranks in (1, 2, 4):
            with self.subTest(ranks=ranks):
                solution = write_system(self.capture, ranks)
                self.assertEqual(replay.capture_inputs(self.capture)[0], ranks)
                self.assertTrue(self.check()['residual_failed'])
                write_solution(self.result, solution, ranks)
                result = self.check(self.result)
                self.assertTrue(result['residual_passed'])
                self.assertTrue(result['original_referenced_copies_equal_owners'])
                shutil.rmtree(self.capture); shutil.rmtree(self.result)

    def test_stale_referenced_ghost_detected(self):
        write_system(self.capture, 2)
        p = replay.rank_file(self.capture, 1, '.bin'); data = bytearray(p.read_bytes())
        # Reversed map: local n-1 is solver ID zero, which belongs to rank zero.
        struct.pack_into('<d', data, len(data)-16, 1.)
        p.write_bytes(data)
        # Force that ghost to be referenced by rank one's first row.
        header = struct.unpack('<13Q', data[:104]); start = 120+4*(header[5]-header[4]+1)
        struct.pack_into('<i', data, start, 11); p.write_bytes(data)
        self.assertFalse(self.check()['original_referenced_copies_equal_owners'])

    def test_nonfinite_replay_never_passes(self):
        solution = write_system(self.capture); solution[0] = float('nan')
        write_solution(self.result, solution, 1)
        result = self.check(self.result)
        self.assertFalse(result['finite']); self.assertFalse(result['residual_passed'])

    def test_missing_truncated_and_wrong_rank_parts(self):
        write_system(self.capture, 2)
        path = replay.rank_file(self.capture, 1, '.bin'); original = path.read_bytes()
        for bad in (b'', original[:-1], original[:16]+struct.pack('<Q', 0)+original[24:]):
            with self.subTest(length=len(bad)):
                path.write_bytes(bad)
                with self.assertRaises(Exception): replay.capture_inputs(self.capture)
                self.check(succeeds=False)
        path.unlink()
        with self.assertRaises(Exception): replay.capture_inputs(self.capture)
        self.check(succeeds=False)

    def test_target_mismatch_and_extra_file_rejected(self):
        write_system(self.capture)
        path = replay.rank_file(self.capture, 0, '.settings')
        path.write_text(path.read_text().replace('atol 0', 'atol 1'))
        with self.assertRaises(Exception): replay.capture_inputs(self.capture)
        path.write_text(path.read_text().replace('atol 1', 'atol 0'))
        (self.capture/'stray').write_text('unbound evidence')
        with self.assertRaises(Exception): replay.capture_inputs(self.capture)

    def test_cancellation_uses_compensated_evaluation(self):
        matrix = [[2.**54, 1., -2.**54], [0., 1., 0.], [0., 0., 1.]]
        write_system(self.capture, matrix=matrix, solution=[1.]*3, rhs=[1.]*3, candidate=[1.]*3)
        self.assertTrue(self.check()['residual_passed'])

    def test_boundary_verdict_is_inconclusive(self):
        write_system(self.capture, matrix=[[1.]], rhs=[1.], candidate=[0.], atol=1.)
        result = self.check()
        self.assertFalse(result['residual_passed']); self.assertFalse(result['residual_failed'])
        self.assertTrue(result['residual_inconclusive'])

    def test_exact_rational_residual_verdicts(self):
        random.seed(319)
        for case in range(24):
            n = 4
            matrix = [[math.ldexp(random.uniform(-1, 1), random.randint(-10, 10)) for _ in range(n)] for _ in range(n)]
            candidate = [random.uniform(-1, 1) for _ in range(n)]
            rhs = [math.fsum(a*x for a,x in zip(row, candidate)) for row in matrix]
            residual = [Fraction(b)-sum((Fraction(a)*Fraction(x) for a,x in zip(row,candidate)), Fraction(0))
                        for b,row in zip(rhs,matrix)]
            norm2 = sum(x*x for x in residual); rhs2 = sum(Fraction(b)**2 for b in rhs)
            atol, rtol = 1e-13 if case%2 else 1e-17, 1e-20
            limit2 = max(Fraction(atol)**2, Fraction(rtol)**2*rhs2)
            write_system(self.capture, matrix=matrix, rhs=rhs, candidate=candidate, atol=atol, rtol=rtol)
            result = self.check()
            if result['residual_passed']: self.assertLessEqual(norm2, limit2)
            if result['residual_failed']: self.assertGreater(norm2, limit2)
            shutil.rmtree(self.capture)

    def test_reference_defaults_are_not_replaced_by_mars_defaults(self):
        config = dict(family='Hypre', rtol=1e-8, atol=0., max_iterations=200, min_iterations=15,
                      options=dict(solver=dict(type='GMRES'), precond=dict(type='BoomerAMG', maxiter=3, tol=.1)))
        controls = replay.reference_controls(config, (0., 1e-8))
        self.assertEqual(controls, dict(method=0, maxiter=200., rtol=1e-8, atol=0.))
        config['options']['solver'].update(type='FlexGMRES', kdim=37)
        config['options']['precond'].update(numsweeps=2, strongthreshold=.6)
        controls = replay.reference_controls(config, (0., 1e-8))
        self.assertEqual((controls['method'], controls['kdim'], controls['numsweeps']), (1, 37., 2.))

    def test_unsupported_reference_configuration_rejected(self):
        base = dict(family='Hypre', rtol=1e-8, atol=0., options=dict(solver=dict(type='GMRES'), precond=dict(type='BoomerAMG')))
        for changes in ({'normalize_matrix':True}, {'diagonal_scaling':True}, {'family':'PETSc'}, {'rtol':True}, {'rtol':1e-7}):
            with self.subTest(changes=changes), self.assertRaises(Exception):
                replay.reference_controls(dict(base, **changes), (0., 1e-8))
        base['options']['precond']['unknown'] = 4
        with self.assertRaises(Exception): replay.reference_controls(base, (0., 1e-8))

    def test_existing_output_directory_untouched(self):
        self.result.mkdir(); (self.result/'keep').write_text('keep')
        with contextlib.redirect_stderr(io.StringIO()):
            result = replay.main(['capture', '--pair', str(self.root), '--executable', 'unused', '--output-dir', str(self.result)])
        self.assertEqual(result, 1); self.assertEqual(list(self.result.iterdir()), [self.result/'keep'])

    def test_fixed_error_summary_never_copies_exception(self):
        with patch.object(replay, 'capture', side_effect=ValueError('SECRET geometry')), contextlib.redirect_stdout(io.StringIO()):
            status = replay.main(['capture', '--pair', str(self.root), '--executable', 'unused', '--output-dir', str(self.result)])
        self.assertEqual(status, 1)
        self.assertNotIn('SECRET', (self.result/'public.json').read_text())
        self.assertEqual(json.loads((self.result/'public.json').read_text())['comparison_status'], 'invalid_evidence')

    def test_real_hypre_replay_when_requested(self):
        binaries = os.environ.get('MARS_TEST_PRESSURE_REPLAYS', '').split(os.pathsep)
        if not binaries[0]:
            self.skipTest('set MARS_TEST_PRESSURE_REPLAYS to locally built replay executables')
        solution = write_system(self.capture)
        for binary in binaries:
            for method in (0, 1):
                with self.subTest(binary=binary, method=method):
                    cfg = self.root/'controls'
                    cfg.write_text('method {}\nrtol 1e-10\natol 0\nmaxiter 200\n'.format(method))
                    p = subprocess.run([binary, str(self.capture), str(cfg), str(self.result)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    self.assertEqual(p.returncode, 0, p.stderr.decode())
                    self.assertTrue(self.check(self.result)['residual_passed'])
                    # Reapply the recorded, post-setup control snapshot on a fresh object.
                    report = replay.numeric_file(self.result/'rank-000000.report')
                    self.assertEqual(report['result_converged'], 1)
                    stop = self.stopping_result([report], 'passed')
                    self.assertEqual(stop['assessment'], 'independent_residual_passed')
                    self.assertEqual(stop['backend'], 'FlexGMRES' if method else 'GMRES')
                    self.assertFalse(stop['solve_return_nonzero'])
                    cfg.write_text(''.join('{} {:.17g}\n'.format(k,v) for k,v in report.items() if not k.startswith('result_')))
                    shutil.rmtree(self.result)
                    p = subprocess.run([binary, str(self.capture), str(cfg), str(self.result)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    self.assertEqual(p.returncode, 0, p.stderr.decode())
                    self.assertTrue(self.check(self.result)['residual_passed'])
                    self.assertEqual(replay.numeric_file(self.result/'rank-000000.report'), report)
                    shutil.rmtree(self.result)

    def captured_record(self, ranks=1):
        write_system(self.capture, ranks)
        system = self.root/'captured-run'; system.mkdir()
        self.capture.rename(system/'system')
        self.capture = system/'system'
        files = replay.capture_inputs(self.capture)[2]
        (system/'reference.settings').write_text('method 0\nrtol 1e-10\natol 0\nmaxiter 200\n')
        files.update(replay.hashes([system/'reference.settings']))
        record = dict(schema=replay.SCHEMA, exit_code=255, ranks=ranks, files=files)
        replay.startup.write_json(system/'capture.json', record)
        return system

    def replay_record(self, capture, name, candidate):
        path = self.root/name; path.mkdir()
        write_solution(path/'result', candidate, 1)
        (path/'result/rank-000000.report').write_text(
            'method 0\nrtol 1e-10\natol 0\nmaxiter 200\nrelaxtype 18\ncoarserelax 18\nresult_converged 1\nresult_fatal_error 0\n'
            'result_iterations 12\nresult_solve_error 0\nresult_global_error 0\nresult_reported 1e-12\n')
        record = dict(schema=replay.SCHEMA, backend=name, profile='captured' if name=='mars' else 'reference',
                      capture_sha256=replay.digest(capture/'capture.json'), exit_code=0,
                      inputs=replay.hashes([capture/'capture.json']), files=replay.hashes(list((path/'result').iterdir())))
        replay.startup.write_json(path/'replay.json', record)
        return path

    def compare_arguments(self):
        capture = self.captured_record()
        mars = self.replay_record(capture, 'mars', [1.+i/8. for i in range(12)])
        reference = self.replay_record(capture, 'reference', [0.]*12)
        return ['compare', '--capture-run', str(capture), '--mars-run', str(mars), '--reference-run', str(reference),
                '--checker', str(self.checker), '--output', str(self.root/'public.json')]

    def test_false_backend_success_cannot_override_common_residual(self):
        args = self.compare_arguments()
        with contextlib.redirect_stdout(io.StringIO()): self.assertEqual(replay.main(args), 0)
        public = json.loads((self.root/'public.json').read_text())
        self.assertTrue(public['residual_checks']['mars']['residual_passed'])
        self.assertTrue(public['residual_checks']['reference']['backend_converged_flag'])
        self.assertTrue(public['residual_checks']['reference']['residual_failed'])
        self.assertEqual(public['stopping_checks']['reference']['assessment'], 'failed_before_iteration_limit')
        self.assertTrue(public['stopping_checks']['reference']['convergence_claim_contradicted'])
        self.assertFalse(public['actual_exit_branch_verified'])
        self.assertFalse(public['original_amg_hierarchy_reused'])
        self.assertNotIn(str(self.root), json.dumps(public))

    def stopping_report(self, **changes):
        report = dict(method=0., maxiter=200., rtol=1e-10, atol=0., result_converged=1.,
                      result_fatal_error=0., result_iterations=12., result_solve_error=0.,
                      result_global_error=0., result_reported=1e-12)
        report.update(changes)
        return report

    def stopping_result(self, reports, outcome='failed'):
        residual = dict(finite=outcome != 'nonfinite', residual_passed=outcome == 'passed',
                        residual_failed=outcome == 'failed', residual_inconclusive=outcome == 'inconclusive')
        return replay.stopping_checks(reports, residual)

    def test_stopping_summary_distinguishes_early_exit_cap_and_zero_iterations(self):
        for ranks in (1, 2, 4):
            for iterations, relation, assessment in ((12, 'below', 'failed_before_iteration_limit'),
                    (200, 'at', 'failed_at_or_above_iteration_limit'),
                    (201, 'above', 'failed_at_or_above_iteration_limit'),
                    (0, 'below', 'failed_without_iterations')):
                with self.subTest(ranks=ranks, iterations=iterations):
                    report = self.stopping_report(result_iterations=iterations)
                    public = self.stopping_result([report]*ranks)
                    self.assertEqual(public['iteration_limit_relation'], relation)
                    self.assertEqual(public['assessment'], assessment)
                    self.assertEqual(public['iterations_zero_on_all_ranks'], iterations == 0)
                    self.assertTrue(public['convergence_claim_contradicted'])
                    self.assertTrue(public['exit_metadata_agrees_across_ranks'])
                    self.assertNotIn('stagnation', json.dumps(public))
                    self.assertTrue(all(type(v) in (str, bool, type(None)) for v in public.values()))

    def test_stopping_summary_does_not_confuse_convergence_error_and_fatal_error(self):
        capped = self.stopping_report(result_iterations=200, result_converged=0,
                                     result_solve_error=256, result_global_error=256, result_reported=1e-6)
        public = self.stopping_result([capped])
        self.assertEqual(public['assessment'], 'failed_at_or_above_iteration_limit')
        self.assertTrue(public['solve_return_nonzero'])
        self.assertTrue(public['global_error_nonzero'])
        self.assertFalse(public['fatal_backend_error_seen'])
        self.assertFalse(public['convergence_claim_contradicted'])
        self.assertFalse(public['reported_relative_residual_below_rtol'])
        fatal = self.stopping_report(result_solve_error=1, result_fatal_error=1)
        self.assertEqual(self.stopping_result([fatal])['assessment'], 'fatal_backend_error')

    def test_stopping_summary_keeps_absolute_tolerance_and_true_residual_separate(self):
        report = self.stopping_report(method=1, atol=1., result_reported=1e-4)
        public = self.stopping_result([report], 'passed')
        self.assertEqual(public['backend'], 'FlexGMRES')
        self.assertTrue(public['absolute_tolerance_enabled'])
        self.assertFalse(public['reported_relative_residual_below_rtol'])
        self.assertEqual(public['assessment'], 'independent_residual_passed')
        self.assertFalse(public['convergence_claim_contradicted'])
        for outcome in ('inconclusive', 'nonfinite'):
            with self.subTest(outcome=outcome):
                public = self.stopping_result([report], outcome)
                self.assertEqual(public['assessment'], 'independent_residual_inconclusive' if outcome == 'inconclusive'
                                 else 'nonfinite_candidate_residual')
                self.assertFalse(public['convergence_claim_contradicted'])

    def test_stopping_summary_reports_nonfinite_norm_and_rank_disagreement(self):
        for value in (float('nan'), float('inf')):
            public = self.stopping_result([self.stopping_report(result_reported=value)])
            self.assertEqual(public['assessment'], 'nonfinite_reported_residual')
            self.assertFalse(public['reported_relative_residual_finite'])
            self.assertIsNone(public['reported_relative_residual_below_rtol'])
            json.dumps(public, allow_nan=False)
        for changed in ({'result_iterations': 200}, {'result_converged': 0}, {'maxiter': 300}, {'method': 1}):
            public = self.stopping_result([self.stopping_report(), self.stopping_report(**changed)])
            self.assertFalse(public['exit_metadata_agrees_across_ranks'])
            self.assertEqual(public['assessment'], 'rank_reports_disagree')
        public = self.stopping_result([self.stopping_report(), self.stopping_report(result_iterations=200)])
        self.assertEqual(public['iteration_limit_relation'], 'mixed')

    def test_stopping_summary_rejects_missing_and_invalid_metadata(self):
        for key, value in (('result_iterations', -1.), ('result_iterations', 1.5), ('maxiter', float('inf')),
                           ('result_converged', 2.), ('method', 3.), ('result_solve_error', float('nan')),
                           ('rtol', 0.), ('atol', -1.), ('result_reported', -1.)):
            with self.subTest(key=key, value=value), self.assertRaises(Exception):
                self.stopping_result([self.stopping_report(**{key: value})])
        report = self.stopping_report(); report.pop('result_iterations')
        with self.assertRaises(KeyError): self.stopping_result([report])

    def test_compare_missing_stopping_field_fails_with_private_values_withheld(self):
        args = self.compare_arguments()
        path = self.root/'reference/result/rank-000000.report'
        path.write_text(path.read_text().replace('result_iterations 12\n', ''))
        manifest = self.root/'reference/replay.json'
        record = json.loads(manifest.read_text()); record['files'].update(replay.hashes([path]))
        manifest.write_text(json.dumps(record))
        public = self.compare_saved(args)
        self.assertEqual(public['failed_check'], 'reference_replay_stopping_report')
        self.assertNotIn('stopping_checks', public)

    def test_modified_replay_output_rejected(self):
        args = self.compare_arguments()
        (self.root/'mars/result/rank-000000.solution').write_bytes(b'changed')
        with contextlib.redirect_stdout(io.StringIO()): self.assertEqual(replay.main(args), 1)
        public = json.loads((self.root/'public.json').read_text())
        self.assertEqual(public['failed_check'], 'mars_replay_outputs')
        self.assertEqual(public['replay_evidence_checks']['mars']['outputs']['changed'], ['replay_output'])

    def test_modified_capture_rejected(self):
        args = self.compare_arguments()
        (self.capture/'rank-000000.settings').write_text('rtol 1e-3\natol 0\n')
        with contextlib.redirect_stdout(io.StringIO()): self.assertEqual(replay.main(args), 1)
        self.assertEqual(json.loads((self.root/'public.json').read_text())['failed_check'], 'capture_identity')

    def test_changed_evidence_during_check_rejected(self):
        args = self.compare_arguments(); actual_run = subprocess.run
        def change(command, **kwargs):
            result = actual_run(command, **kwargs)
            if command[-1] != '-': (self.root/'mars/result/rank-000000.report').write_text('changed')
            return result
        with patch.object(replay.subprocess, 'run', side_effect=change), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(replay.main(args), 1)

    def compare_saved(self, args, code=1):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(replay.main(args), code)
        text = Path(args[-1]).read_text()
        self.assertNotIn('SECRET', text)
        self.assertNotIn(str(self.root), text)
        return json.loads(text)

    def add_replay_inputs(self, name, paths, executable=None):
        path = self.root/name/'replay.json'
        record = json.loads(path.read_text())
        record['inputs'].update(replay.hashes(paths))
        if executable is not None:
            record['command'] = ['launcher', str(executable), 'system', 'configuration', 'result']
        path.write_text(json.dumps(record))

    def archive_fixture(self):
        args = self.compare_arguments()
        environment = self.root/'SECRET-uenv'; environment.mkdir()
        deps = [environment/name for name in ('libHYPRE.so', 'libmpi.so', 'libstdc++.so.6',
                                              'SECRET-build.hpp', 'SECRET-replay')]
        for path in deps:
            path.write_text('SECRET bytes for ' + path.name)
        capture_file = self.capture.parent/'capture.json'
        capture = json.loads(capture_file.read_text())
        capture['reference'] = {'libraries': replay.hashes(deps[:3])}
        capture_file.write_text(json.dumps(capture))
        for name in ('mars', 'reference'):
            path = self.root/name/'replay.json'
            record = json.loads(path.read_text())
            record['capture_sha256'] = replay.digest(capture_file)
            record['inputs'].update(replay.hashes([capture_file]))
            record['inputs'].update(capture['files'])
            if name == 'reference':
                libraries = self.root/name/'result/rank-000000.libraries'
                libraries.write_text('{}\n{}\n'.format(*deps[:2]))
                record['files'].update(replay.hashes([libraries]))
                record['inputs'].update(replay.hashes(deps))
                record['command'] = ['launcher', str(deps[-1]), str(self.capture),
                                     str(capture_file.parent/'reference.settings'), str(self.root/name/'result')]
            path.write_text(json.dumps(record))
        archive = self.root/'archive'
        archive_args = ['archive-inputs', '--capture-run', str(capture_file.parent),
                        '--replay-run', str(self.root/'reference'), '--output-dir', str(archive)]
        return args, archive_args, archive, environment, deps

    def archive_saved(self, args, code=0):
        before = replay.hashes(p for name in ('captured-run', 'mars', 'reference')
                               for p in (self.root/name).rglob('*') if p.is_file())
        with patch.object(replay.subprocess, 'run') as run, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(replay.main(args), code)
            run.assert_not_called()
        self.assertEqual(replay.hashes(Path(p) for p in before), before)
        text = (Path(args[-1])/'public.json').read_text()
        self.assertNotIn('SECRET', text)
        self.assertNotIn(str(self.root), text)
        return json.loads(text)

    def test_archived_dependencies_cross_environments_without_changing_residual_verdicts(self):
        args, archive_args, archive, environment, deps = self.archive_fixture()
        public = self.archive_saved(archive_args)
        self.assertEqual(public['comparison_status'], 'input_archive_complete')
        self.assertTrue(public['loaded_library_identity_verified'])
        self.assertFalse(public['solver_launched'])
        self.assertFalse(public['convergence_verified'])
        manifest = json.loads((archive/'archive.json').read_text())
        self.assertEqual(manifest['inputs'], replay.hashes(deps))
        environment.rename(self.root/'unmounted')
        with patch.object(replay.subprocess, 'run') as run:
            failed = self.compare_saved(args)
            run.assert_not_called()
        self.assertEqual(failed['failed_check'], 'reference_replay_inputs')
        args = args[:-2] + ['--reference-input-archive', str(archive), '--output', str(self.root/'archived.json')]
        passed = self.compare_saved(args, 0)
        self.assertEqual(passed['comparison_status'], 'completed')
        self.assertEqual(passed['replay_evidence_checks']['reference']['inputs']['scope'],
                         'archived_dependencies_and_live_capture_inputs')
        self.assertTrue(passed['residual_checks']['mars']['residual_passed'])
        self.assertTrue(passed['residual_checks']['reference']['residual_failed'])
        self.assertFalse(passed['nonlinear_convergence_verified'])

    def test_input_archive_cannot_be_created_from_invalid_evidence(self):
        args, archive_args, archive, environment, deps = self.archive_fixture()
        library_output = self.root/'reference/result/rank-000000.libraries'
        for i, case in enumerate(('missing', 'changed', 'wrong_library', 'changed_output')):
            with self.subTest(case=case):
                originals = {p: p.read_bytes() for p in (deps[1], library_output)}
                if case == 'missing': deps[1].unlink()
                elif case == 'changed': deps[1].write_text('changed')
                else:
                    library_output.write_text('{}\n{}\n'.format(deps[0], deps[0]))
                    if case == 'wrong_library':
                        path = self.root/'reference/replay.json'
                        record = json.loads(path.read_text())
                        record['files'].update(replay.hashes([library_output])); path.write_text(json.dumps(record))
                destination = self.root/('archive-bad-' + str(i)); archive_args[-1] = str(destination)
                public = self.archive_saved(archive_args, 1)
                self.assertEqual(public['failed_check'], 'loaded_library_identity' if case == 'wrong_library'
                                 else 'replay_outputs' if case == 'changed_output' else 'replay_inputs')
                self.assertFalse((destination/'archive.json').exists())
                for p, data in originals.items(): p.write_bytes(data)
                record = json.loads((self.root/'reference/replay.json').read_text())
                record['files'].update(replay.hashes([library_output]))
                (self.root/'reference/replay.json').write_text(json.dumps(record))

    def test_archive_binding_and_bytes_are_required(self):
        args, archive_args, archive, environment, deps = self.archive_fixture()
        self.archive_saved(archive_args)
        manifest = archive/'archive.json'; saved = manifest.read_text()
        obj = archive/'objects'/replay.digest(deps[1]); contents = obj.read_bytes()
        environment.rename(self.root/'unmounted')
        for case in ('binding', 'capture', 'backend', 'omitted', 'invalid_hash', 'public_only', 'missing', 'changed'):
            with self.subTest(case=case):
                record = json.loads(saved)
                if case == 'binding': record['replay_sha256'] = '0'*64
                if case == 'capture': record['capture_sha256'] = '0'*64
                if case == 'backend': record['backend'] = 'mars'
                if case == 'omitted': record['inputs'].pop(str(deps[1]))
                if case == 'invalid_hash': record['inputs'][str(deps[1])] = '../SECRET'
                if case == 'public_only': record = json.loads((archive/'public.json').read_text())
                manifest.write_text(json.dumps(record))
                if case == 'missing': obj.unlink()
                if case == 'changed': obj.write_text('changed')
                compare_args = args[:-2] + ['--reference-input-archive', str(archive), '--output', str(self.root/(case+'.json'))]
                with patch.object(replay.subprocess, 'run') as run:
                    public = self.compare_saved(compare_args)
                    run.assert_not_called()
                self.assertEqual(public['failed_check'], 'reference_replay_inputs' if case in ('missing', 'changed')
                                 else 'reference_input_archive')
                obj.write_bytes(contents)
        manifest.write_text(saved)

    def test_archived_inputs_do_not_hide_changed_live_capture_or_outputs(self):
        args, archive_args, archive, environment, deps = self.archive_fixture()
        self.archive_saved(archive_args)
        environment.rename(self.root/'unmounted')
        for i, path in enumerate((self.root/'reference/result/rank-000000.solution',
                                  self.capture/'rank-000000.bin', self.root/'reference/replay.json')):
            with self.subTest(path=path.name):
                saved = path.read_bytes(); path.write_bytes(saved + b'\n')
                compare_args = args[:-2] + ['--reference-input-archive', str(archive), '--output', str(self.root/('live-'+str(i)+'.json'))]
                with patch.object(replay.subprocess, 'run') as run:
                    public = self.compare_saved(compare_args)
                    run.assert_not_called()
                self.assertEqual(public['failed_check'], ('reference_replay_outputs', 'capture_identity', 'reference_input_archive')[i])
                path.write_bytes(saved)

    def test_input_archive_rechecks_source_after_copy(self):
        args, archive_args, archive, environment, deps = self.archive_fixture()
        copy = replay.shutil.copyfileobj
        def change(src, dst):
            copy(src, dst)
            Path(src.name).write_bytes(b'changed during copy')
        with patch.object(replay.shutil, 'copyfileobj', side_effect=change):
            public = self.archive_saved(archive_args, 1)
        self.assertEqual(public['failed_check'], 'evidence_changed')
        self.assertFalse((archive/'archive.json').exists())

    def test_archived_object_is_rechecked_after_residual_evaluation(self):
        args, archive_args, archive, environment, deps = self.archive_fixture()
        self.archive_saved(archive_args)
        obj = archive/'objects'/replay.digest(deps[1])
        actual_run = subprocess.run
        def change(command, **kwargs):
            result = actual_run(command, **kwargs)
            if command[-1] == str(self.root/'reference/result'): obj.write_bytes(b'changed')
            return result
        args = args[:-2] + ['--reference-input-archive', str(archive), '--output', str(self.root/'changed-during-check.json')]
        with patch.object(replay.subprocess, 'run', side_effect=change):
            public = self.compare_saved(args)
        self.assertEqual(public['failed_check'], 'reference_replay_inputs_changed')

    def test_incomplete_copy_cannot_create_a_valid_archive(self):
        args, archive_args, archive, environment, deps = self.archive_fixture()
        def corrupt(src, dst):
            dst.write(b'incomplete copy')
        with patch.object(replay.shutil, 'copyfileobj', side_effect=corrupt):
            public = self.archive_saved(archive_args, 1)
        self.assertEqual(public['failed_check'], 'input_archive_copy')
        self.assertFalse((archive/'archive.json').exists())

    def test_mars_archive_uses_mars_library_identity(self):
        args, archive_args, archive, environment, deps = self.archive_fixture()
        capture_file = self.capture.parent/'capture.json'
        capture = json.loads(capture_file.read_text())
        capture['libraries'] = capture['reference']['libraries']
        capture['reference']['libraries'] = {str(deps[0]): '0'*64}
        capture_file.write_text(json.dumps(capture))
        record = json.loads((self.root/'reference/replay.json').read_text())
        record.update(backend='mars', profile='captured', capture_sha256=replay.digest(capture_file))
        record['inputs'].update(replay.hashes([capture_file]))
        (self.root/'reference/replay.json').write_text(json.dumps(record))
        public = self.archive_saved(archive_args)
        self.assertEqual(public['backend'], 'mars')
        self.assertTrue(public['loaded_library_identity_verified'])

    def test_compare_checks_both_environments_before_running_checker(self):
        args = self.compare_arguments()
        hypre = self.root/'libHYPRE.so'; hypre.write_text('SECRET hypre')
        mpi = self.root/'libmpi.so'; mpi.write_text('SECRET mpi')
        self.add_replay_inputs('mars', [hypre])
        self.add_replay_inputs('reference', [mpi])
        hypre.unlink(); mpi.write_text('changed')
        with patch.object(replay.subprocess, 'run') as run:
            public = self.compare_saved(args)
            run.assert_not_called()
        self.assertEqual(public['failed_candidate'], 'mars')
        self.assertEqual(public['failed_check'], 'mars_replay_inputs')
        evidence = public['replay_evidence_checks']
        self.assertEqual(evidence['mars']['inputs']['missing'], ['hypre_library'])
        self.assertEqual(evidence['reference']['inputs']['changed'], ['mpi_library'])
        self.assertTrue(evidence['mars']['outputs']['matched'])
        self.assertTrue(evidence['reference']['outputs']['matched'])
        self.assertNotIn('residual_checks', public)

    def test_compare_separates_changed_binary_missing_runtime_and_unreadable_build_input(self):
        args = self.compare_arguments()
        binary, library, source = (self.root/name for name in ('SECRET-exe', 'libstdc++.so.6', 'SECRET-source.hpp'))
        for path in (binary, library, source): path.write_text(path.name)
        self.add_replay_inputs('reference', [binary, library, source], binary)
        binary.write_text('changed'); library.unlink()
        actual_digest = replay.digest
        def digest(path):
            if str(path) == str(source): raise PermissionError('SECRET permission failure')
            return actual_digest(path)
        with patch.object(replay, 'digest', side_effect=digest), patch.object(replay.subprocess, 'run') as run:
            public = self.compare_saved(args)
            run.assert_not_called()
        self.assertEqual(public['failed_check'], 'reference_replay_inputs')
        check = public['replay_evidence_checks']['reference']['inputs']
        self.assertEqual(check['changed'], ['executable'])
        self.assertEqual(check['missing'], ['other_runtime_library'])
        self.assertEqual(check['unreadable'], ['source_or_build_input'])
        self.assertTrue(public['replay_evidence_checks']['mars']['verified'])

    def test_compare_rejects_missing_record_wrong_binding_and_malformed_manifest(self):
        args = self.compare_arguments()
        path = self.root/'reference/replay.json'; saved = path.read_text()
        for case in ('missing', 'binding', 'manifest'):
            with self.subTest(case=case):
                record = json.loads(saved)
                if case == 'binding': record['capture_sha256'] = '0'*64
                if case == 'manifest': record['inputs']['SECRET-relative'] = 'SECRET invalid hash'
                if case == 'missing': path.unlink()
                else: path.write_text(json.dumps(record))
                args[-1] = str(self.root/(case + '.json'))
                with patch.object(replay.subprocess, 'run') as run:
                    public = self.compare_saved(args)
                    run.assert_not_called()
                expected = {'missing': 'replay_record', 'binding': 'replay_binding', 'manifest': 'replay_inputs'}[case]
                self.assertEqual(public['failed_check'], 'reference_' + expected)
                if case == 'manifest':
                    self.assertFalse(public['replay_evidence_checks']['reference']['inputs']['record_valid'])

    def test_compare_checker_failures_are_not_reported_as_replay_identity(self):
        args = self.compare_arguments()
        cases = [
            ('launch', FileNotFoundError('SECRET loader'), 'original_checker_launch'),
            ('permissions', PermissionError('SECRET mode'), 'original_checker_launch'),
            ('loader', subprocess.CompletedProcess([], 127, b'', b'error while loading shared libraries: SECRET'), 'original_checker_exit'),
            ('version', subprocess.CompletedProcess([], 1, b'', b"version `SECRET' not found (required by SECRET)"), 'original_checker_exit'),
            ('format', subprocess.CompletedProcess([], 1, b'', b'ERROR: private pressure residual check failed\n'), 'original_checker_exit'),
            ('signal', subprocess.CompletedProcess([], -9, b'', b'SECRET'), 'original_checker_exit'),
            ('json', subprocess.CompletedProcess([], 0, b'SECRET malformed JSON', b''), 'original_checker_output'),
            ('schema', subprocess.CompletedProcess([], 0, b'{"schema":"SECRET"}', b''), 'original_checker_output')]
        for case, response, expected in cases:
            with self.subTest(case=case):
                args[-1] = str(self.root/(case + '.json'))
                def run(*unused, **kwargs):
                    if isinstance(response, Exception): raise response
                    return response
                with patch.object(replay.subprocess, 'run', side_effect=run):
                    public = self.compare_saved(args)
                self.assertEqual(public['failed_check'], expected)
                self.assertTrue(all(check['verified'] for check in public['replay_evidence_checks'].values()))
                if case in ('loader', 'version'):
                    self.assertTrue(public['checker_diagnostics']['library_load_error_seen'])
                if case == 'format':
                    self.assertTrue(public['checker_diagnostics']['residual_check_error_seen'])
                if case == 'signal': self.assertEqual(public['checker_diagnostics']['process_exit_code'], -9)
                if case == 'permissions': self.assertEqual(public['checker_launch_error'], 'permission_denied')

    def test_compare_distinguishes_invalid_report_from_identity(self):
        args = self.compare_arguments()
        report = self.root/'mars/result/rank-000000.report'
        report.write_text('SECRET malformed\n')
        path = self.root/'mars/replay.json'; record = json.loads(path.read_text())
        record['files'][str(report)] = replay.digest(report)
        path.write_text(json.dumps(record))
        public = self.compare_saved(args)
        self.assertEqual(public['failed_check'], 'mars_replay_report')

    def test_capture_retains_launch_controls_and_rejection(self):
        pair = self.root/'pair'; pair.mkdir()
        for name in ('pair.json', 'case.json'): (pair/name).write_text('{}')
        for name in ('reference', 'mars'):
            (pair/name).mkdir()
            for file in ('launch-start.json', 'launch.json', 'run.log', 'run.exit'): (pair/name/file).write_text('record')
        (pair/'reference/input.i').write_text('SECRET_DECK')
        exe = self.root/'solver'; exe.write_text('binary')
        library = self.root/'libHYPRE.so'; library.write_text('library')
        mpi = self.root/'libmpi.so'; mpi.write_text('mpi')
        arguments = ['--iterations', '20', '--pressure-linear-rtol', '1e-10', '--pressure-linear-atol', '0',
                     '--snapshot-iterations', '20', '--field-output', 'distributed', '--output-prefix', 'old']
        config = dict(family='Hypre', rtol=1e-10, atol=0., options=dict(solver=dict(type='GMRES'), precond=dict(type='BoomerAMG')))
        record = dict(status='finished', environment={}, ranks=1)
        def run(command, stdout, **kwargs):
            actual = replay.options(command[2:])
            self.assertEqual(actual['--iterations'], '20')
            self.assertEqual(actual['--field-output'], 'none')
            self.assertEqual(actual['--snapshot-iterations'], '0')
            write_system(Path(actual['--pressure-failure-capture']))
            stdout.write(('SECRET log\n'+replay.CAPTURE_MARKER+'\n').encode())
            return subprocess.CompletedProcess(command, 255)
        with patch.object(replay.settings, 'saved_configuration', return_value=(config, arguments)), \
             patch.object(replay.settings, 'saved_launch', return_value=record), \
             patch.object(replay.startup, 'solver_arguments', return_value=arguments), \
             patch.object(replay.probe, 'launcher', return_value=['launcher']), \
             patch.object(replay.startup, 'runtime_libraries', return_value=replay.hashes([library, mpi])), \
             patch.object(replay.subprocess, 'run', side_effect=run), patch.dict(os.environ, {}, clear=True), \
             contextlib.redirect_stdout(io.StringIO()):
            status = replay.main(['capture', '--pair', str(pair), '--executable', str(exe), '--output-dir', str(self.result)])
        self.assertEqual(status, 0)
        self.assertEqual(replay.captured_record(self.result)['exit_code'], 255)
        public = (self.result/'public.json').read_text()
        self.assertNotIn('SECRET', public); self.assertIn('capture_complete', public)

    def test_replay_checks_loaded_libraries_and_original_target(self):
        capture = self.captured_record()
        library, mpi, exe = (self.root/name for name in ('libHYPRE.so', 'libmpi.so', 'gpu-replay'))
        for path in (library, mpi, exe): path.write_text(path.name)
        libraries = replay.hashes([library, mpi])
        record = json.loads((capture/'capture.json').read_text())
        record.update(libraries=libraries, reference={'libraries':libraries}, launcher=['launcher'])
        (capture/'capture.json').write_text(json.dumps(record))
        def launch(command, stdout, **kwargs):
            self.assertEqual(command[-2], str(capture/'system'))
            write_solution(Path(command[-1]), [1.+i/8. for i in range(12)], 1)
            (Path(command[-1])/'rank-000000.libraries').write_text('{}\n{}\n'.format(library, mpi))
            stdout.write(b'SECRET replay output\n')
            return subprocess.CompletedProcess(command, 0)
        args = ['replay', '--capture-run', str(capture), '--backend', 'mars', '--executable', str(exe), '--output-dir', str(self.result)]
        with patch.object(replay.startup, 'runtime_libraries', return_value=libraries), \
             patch.object(replay.subprocess, 'run', side_effect=launch), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(replay.main(args), 0)
        public = json.loads((self.result/'public.json').read_text())
        self.assertTrue(public['loaded_library_identity_verified'])
        self.assertFalse(public['convergence_verified'])
        for failure in ('launcher_exit', 'completion_marker', 'loaded_library_identity'):
            with self.subTest(failure=failure):
                args[-1] = str(self.root/failure)
                def failed_launch(command, stdout, **kwargs):
                    result = launch(command, stdout, **kwargs)
                    if failure == 'launcher_exit':
                        return subprocess.CompletedProcess(command, 127)
                    if failure == 'completion_marker':
                        (Path(command[-1])/'complete').unlink()
                    if failure == 'loaded_library_identity':
                        (Path(command[-1])/'rank-000000.libraries').write_text('SECRET malformed paths\n')
                    return result
                with patch.object(replay.startup, 'runtime_libraries', return_value=libraries), \
                     patch.object(replay.subprocess, 'run', side_effect=failed_launch), contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(replay.main(args), 1)
                public = json.loads((Path(args[-1])/'public.json').read_text())
                self.assertEqual(public['failed_check'], failure)
                self.assertNotIn('SECRET', json.dumps(public))
        # Changing a captured dependency must prevent a second launch.
        library.write_text('changed')
        args[-1] = str(self.root/'retry')
        with patch.object(replay.subprocess, 'run') as command, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(replay.main(args), 1)
            command.assert_not_called()

    def test_gpu_profile_replay_verifies_applied_controls_before_completion(self):
        import simple_pressure_profile as profile
        import test_simple_pressure_profile as fixtures
        config = dict(family='hypre', rtol=1e-10, atol=0.,
                      options=dict(solver=dict(type='gmres'), precond=dict(type='boomeramg')), max_iterations=200)
        values, _ = profile.resolve(config, (0., 1e-10), fixtures.measured_defaults())
        capture = self.captured_record()
        pair = self.root / 'profile'; pair.mkdir()
        source = pair / 'pressure.profile'; source.write_text(profile.text_profile(values))
        library, mpi, exe = (self.root/name for name in ('libHYPRE.so', 'libmpi.so', 'gpu-replay'))
        for path in (library, mpi, exe): path.write_text(path.name)
        libraries = replay.hashes([library, mpi])
        record = json.loads((capture/'capture.json').read_text())
        record.update(libraries=libraries, reference={'libraries':libraries}, launcher=['launcher'])
        (capture/'capture.json').write_text(json.dumps(record))
        for fault in ('none', 'cycle', 'target'):
            output = self.root / ('gpu-' + fault)
            def launch(command, stdout, **kwargs):
                settings = replay.numeric_file(Path(command[-2]))
                self.assertEqual(settings, values)
                result = Path(command[-1])
                write_solution(result, [1.+i/8. for i in range(12)], 1)
                (result/'rank-000000.libraries').write_text('{}\n{}\n'.format(library, mpi))
                report = dict(values, effective_levels=3, effective_relax_3=18)
                report['effective_relax_1'] = report.pop('relax_down')
                report['effective_relax_2'] = report.pop('relax_up')
                if fault == 'cycle': report['effective_relax_2'] = 18
                if fault == 'target': report['rtol'] = 1e-4
                (result/'rank-000000.report').write_text(profile.text_profile(report))
                return subprocess.CompletedProcess(command, 0)
            args = ['replay', '--capture-run', str(capture), '--backend', 'mars', '--profile', 'gpu-reference',
                    '--gpu-profile-pair', str(pair), '--executable', str(exe), '--output-dir', str(output)]
            with patch.object(replay, 'gpu_configuration', return_value=(values, replay.hashes([source]))), \
                 patch.object(replay.startup, 'runtime_libraries', return_value=libraries), \
                 patch.object(replay.subprocess, 'run', side_effect=launch), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(replay.main(args), 0 if fault == 'none' else 1)
            public = json.loads((output/'public.json').read_text())
            self.assertEqual(public['failed_check'], 'none' if fault == 'none' else 'gpu_profile_runtime_settings')
            self.assertEqual((output/'replay.json').exists(), fault == 'none')
            if fault == 'none':
                self.assertTrue(public['captured_pressure_target_preserved'])
                self.assertFalse(public['convergence_verified'])

    def test_gpu_profile_comparison_keeps_independent_residual_verdict(self):
        import simple_pressure_profile as profile
        import test_simple_pressure_profile as fixtures
        args = self.compare_arguments()
        config = dict(family='hypre', rtol=1e-10, atol=0., max_iterations=200,
                      options=dict(solver=dict(type='gmres'), precond=dict(type='boomeramg')))
        values, _ = profile.resolve(config, (0., 1e-10), fixtures.measured_defaults())
        path = self.root / 'mars'
        source = self.root / 'pressure.profile'; source.write_text(profile.text_profile(values))
        config_path = path / 'gpu-reference.settings'; config_path.write_text(profile.text_profile(values))
        report_path = path / 'result/rank-000000.report'
        report = replay.numeric_file(report_path)
        report.update(values, effective_levels=3, effective_relax_3=18)
        report['effective_relax_1'] = report.pop('relax_down')
        report['effective_relax_2'] = report.pop('relax_up')
        report_path.write_text(profile.text_profile(report))
        record = json.loads((path/'replay.json').read_text())
        record.update(profile='gpu-reference', gpu_profile_source=str(source),
                      command=['exe', str(self.capture), str(config_path), str(path/'result')])
        record['inputs'].update(replay.hashes([source, config_path]))
        record['files'].update(replay.hashes([report_path]))
        (path/'replay.json').write_text(json.dumps(record))
        public = self.compare_saved(args, code=0)
        self.assertEqual(public['profiles']['mars'], 'gpu-reference')
        self.assertTrue(public['residual_checks']['mars']['residual_passed'])
        self.assertTrue(public['residual_checks']['mars']['captured_pressure_target_preserved'])
        # A backend flag cannot rescue a wrong solution on the frozen system.
        solution = path / 'result/rank-000000.solution'
        solution.write_bytes(bytes(solution.stat().st_size))
        record['files'].update(replay.hashes([solution]))
        (path/'replay.json').write_text(json.dumps(record))
        args[-1] = str(self.root/'bad-candidate.json')
        public = self.compare_saved(args, code=0)
        self.assertTrue(public['residual_checks']['mars']['residual_failed'])
        report['effective_relax_2'] = 18
        report_path.write_text(profile.text_profile(report))
        record['files'].update(replay.hashes([report_path]))
        (path/'replay.json').write_text(json.dumps(record))
        args[-1] = str(self.root/'wrong-settings.json')
        public = self.compare_saved(args, code=1)
        self.assertEqual(public['failed_check'], 'mars_gpu_profile_runtime_settings')

    def test_gpu_profile_invalid_selection_stops_before_launch(self):
        capture = self.captured_record()
        for index, selection in enumerate((['--backend', 'reference', '--profile', 'gpu-reference', '--gpu-profile-pair', 'unused'],
                                          ['--backend', 'mars', '--profile', 'gpu-reference'],
                                          ['--backend', 'mars', '--profile', 'captured', '--gpu-profile-pair', 'unused'])):
            output = self.root / ('selection-' + str(index))
            with patch.object(replay.subprocess, 'run') as launch, contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(replay.main(['replay', '--capture-run', str(capture), '--output-dir', str(output)] + selection), 1)
                launch.assert_not_called()
            self.assertEqual(json.loads((output/'public.json').read_text())['failed_check'], 'profile_selection')

    def inspection_fixture(self, ranks=1):
        capture = self.captured_record(ranks)
        library, mpi, exe = (self.root/name for name in ('libHYPRE.so', 'libmpi.so', 'reference-replay'))
        for path in (library, mpi, exe): path.write_text('SECRET ' + path.name)
        record = json.loads((capture/'capture.json').read_text())
        record['reference'] = {'libraries': replay.hashes([library, mpi])}
        (capture/'capture.json').write_text(json.dumps(record))
        self.result.mkdir()
        write_solution(self.result/'result', [1.+i/8. for i in range(12)], ranks)
        for rank in range(ranks):
            replay.rank_file(self.result/'result', rank, '.libraries').write_text('{}\n{}\n'.format(library, mpi))
            replay.rank_file(self.result/'result', rank, '.report').write_text('SECRET private report\n')
        launch = dict(schema=replay.SCHEMA, backend='reference', profile='reference',
                      command=['launcher', str(exe), str(capture/'system'), str(capture/'reference.settings'), str(self.result/'result')],
                      capture_sha256=replay.digest(capture/'capture.json'), inputs=replay.hashes([exe, library, mpi]))
        replay.startup.write_json(self.result/'launch-start.json', launch)
        (self.result/'run.exit').write_text('0\n')
        (self.result/'run.log').write_text('SECRET private runtime output\n')
        return ['inspect', '--replay-run', str(self.result), '--output', str(self.root/'inspection.json')]

    def test_reference_build_uses_pic_and_still_rejects_executable_identity(self):
        capture = self.captured_record()
        library, mpi, cache = (self.root/name for name in ('libHYPRE.so', 'libmpi.so', 'CMakeCache.txt'))
        for path in (library, mpi, cache): path.write_text('SECRET ' + path.name)
        include = self.root/'include'; include.mkdir()
        libraries = replay.hashes([library, mpi])
        record = json.loads((capture/'capture.json').read_text())
        record.update(reference={'libraries': libraries}, launcher=['launcher'])
        (capture/'capture.json').write_text(json.dumps(record))
        compiler = str(self.root/'compiler')
        for wrong_identity in (False, True):
            with self.subTest(wrong_identity=wrong_identity):
                output = self.root/('reference-wrong' if wrong_identity else 'reference-ok')
                def run(command, **kwargs):
                    if command[0] == compiler:
                        self.assertGreater(command.index('-fPIC'), command.index('-fno-pic'))
                        Path(command[-1]).write_text('SECRET executable')
                    else:
                        self.assertEqual(command[-2], str(capture/'reference.settings'))
                        write_solution(Path(command[-1]), [1.+i/8. for i in range(12)], 1)
                        reported = command[-4] if wrong_identity else str(library)
                        replay.rank_file(Path(command[-1]), 0, '.libraries').write_text('{}\n{}\n'.format(reported, mpi))
                    return subprocess.CompletedProcess(command, 0)
                with patch.object(replay.defaults, 'matching_install', return_value=(library, include)), \
                     patch.object(replay.defaults, 'cached_toolchain', return_value=(compiler, ['-fno-pic'], [], replay.hashes([cache]))), \
                     patch.object(replay.defaults, 'find_compiler', return_value=compiler), \
                     patch.object(replay.startup, 'runtime_libraries', return_value=libraries), \
                     patch.object(replay.subprocess, 'run', side_effect=run), contextlib.redirect_stdout(io.StringIO()):
                    code = replay.main(['replay', '--capture-run', str(capture), '--backend', 'reference',
                                        '--build-cache', str(cache), '--output-dir', str(output)])
                self.assertEqual(code, int(wrong_identity))
                public = json.loads((output/'public.json').read_text())
                self.assertEqual(public['loaded_library_identity_verified'], not wrong_identity)
                self.assertEqual((output/'replay.json').exists(), not wrong_identity)
                self.assertEqual(public['library_identity_checks']['hypre']['failures'],
                                 ['executable_instead_of_library'] if wrong_identity else [])
                self.assertFalse(public.get('convergence_verified', False))
                self.assertNotIn('SECRET', json.dumps(public))
                self.assertNotIn(str(self.root), json.dumps(public))

    def inspect_saved(self, args, expected_status=0):
        before = replay.hashes(p for p in self.result.rglob('*') if p.is_file())
        with patch.object(replay.subprocess, 'run') as command, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(replay.main(args), expected_status)
            command.assert_not_called()
        self.assertEqual(replay.hashes(p for p in self.result.rglob('*') if p.is_file()), before)
        text = Path(args[-1]).read_text()
        self.assertNotIn('SECRET', text)
        self.assertNotIn(str(self.root), text)
        return json.loads(text)

    def test_saved_inspection_does_not_launch_or_claim_convergence(self):
        public = self.inspect_saved(self.inspection_fixture(2))
        self.assertEqual(public['comparison_status'], 'inspection_complete')
        self.assertEqual(public['failed_check'], 'none')
        self.assertTrue(public['loaded_library_identity_verified'])
        self.assertTrue(public['all_report_parts_present'])
        self.assertFalse(public['convergence_verified'])
        self.assertFalse(public['solver_launched_by_inspection'])

    def test_saved_inspection_distinguishes_failed_launch_and_loader(self):
        args = self.inspection_fixture()
        (self.result/'run.exit').write_text('127\n')
        (self.result/'run.log').write_text('/SECRET/exe: error while loading shared libraries: SECRET.so: missing\n')
        shutil.rmtree(self.result/'result')
        public = self.inspect_saved(args)
        self.assertEqual(public['failed_check'], 'launcher_exit')
        self.assertEqual(public['process_exit_code'], 127)
        self.assertTrue(public['library_load_error_seen'])
        self.assertFalse(public['result_directory_present'])

    def test_saved_inspection_identifies_application_abort_and_partial_output(self):
        args = self.inspection_fixture(2)
        (self.result/'run.exit').write_text('1\n')
        (self.result/'run.log').write_text('ERROR: private pressure replay failed\n'
                                         'application called MPI_Abort(MPI_COMM_WORLD, 1) - process 0\n')
        replay.rank_file(self.result/'result', 1, '.solution').unlink()
        public = self.inspect_saved(args)
        self.assertTrue(public['replay_error_seen'])
        self.assertTrue(public['mpi_abort_seen'])
        self.assertTrue(public['any_solution_parts_present'])
        self.assertFalse(public['all_solution_parts_present'])
        self.assertEqual(public['failed_check'], 'launcher_exit')

    def test_saved_inspection_distinguishes_missing_marker(self):
        args = self.inspection_fixture()
        (self.result/'result/complete').unlink()
        public = self.inspect_saved(args)
        self.assertEqual(public['process_exit_code'], 0)
        self.assertEqual(public['failed_check'], 'completion_marker')

    def test_saved_inspection_checks_every_rank_library_identity(self):
        args = self.inspection_fixture(2)
        wrong = self.root/'SECRET-library'; wrong.write_text('different')
        replay.rank_file(self.result/'result', 1, '.libraries').write_text('{}\n{}\n'.format(wrong, self.root/'libmpi.so'))
        public = self.inspect_saved(args)
        self.assertTrue(public['completion_marker_valid'])
        self.assertEqual(public['failed_check'], 'loaded_library_identity')
        self.assertEqual(public['library_identity_checks'], {
            'hypre': {'matched': False, 'failures': ['hash_mismatch']},
            'mpi': {'matched': True, 'failures': []}})

    def test_saved_library_failure_reasons_are_private_and_never_accepted(self):
        args = self.inspection_fixture(2)
        exe = self.root/'reference-replay'
        alias = self.root/'SECRET-alias'; alias.symlink_to(exe)
        identity = replay.rank_file(self.result/'result', 1, '.libraries')
        for i, (path, reason) in enumerate(((str(exe), 'executable_instead_of_library'),
                                           (str(alias), 'executable_instead_of_library'),
                                           ('SECRET-relative.so', 'non_absolute_path'),
                                           (str(self.root/'SECRET-missing.so'), 'file_unreadable'),
                                           (str(self.root/'libmpi.so'), 'hash_mismatch'))):
            with self.subTest(reason=reason):
                identity.write_text('{}\n{}\n'.format(path, self.root/'libmpi.so'))
                args[-1] = str(self.root/'inspect-{}.json'.format(i))
                public = self.inspect_saved(args)
                self.assertEqual(public['failed_check'], 'loaded_library_identity')
                self.assertFalse(public['loaded_library_identity_verified'])
                self.assertEqual(public['library_identity_checks']['hypre']['failures'], [reason])
                self.assertTrue(public['library_identity_checks']['mpi']['matched'])

    def test_saved_library_records_missing_or_malformed_fail_both_identities(self):
        args = self.inspection_fixture(2)
        identity = replay.rank_file(self.result/'result', 1, '.libraries')
        for text, reason in ((None, 'record_unreadable'), ('SECRET malformed\n', 'record_malformed')):
            with self.subTest(reason=reason):
                if text is None: identity.unlink()
                else: identity.write_text(text)
                args[-1] = str(self.root/('inspect-' + reason + '.json'))
                public = self.inspect_saved(args)
                self.assertEqual(public['failed_check'], 'loaded_library_identity')
                self.assertEqual(public['library_identity_checks'], {
                    kind: {'matched': False, 'failures': [reason]} for kind in ('hypre', 'mpi')})

    def test_saved_library_alias_is_accepted_only_by_matching_contents(self):
        args = self.inspection_fixture(2)
        alias = self.root/'SECRET-alias'; alias.symlink_to(self.root/'libHYPRE.so')
        replay.rank_file(self.result/'result', 1, '.libraries').write_text('{}\n{}\n'.format(alias, self.root/'libmpi.so'))
        public = self.inspect_saved(args)
        self.assertEqual(public['failed_check'], 'none')
        self.assertTrue(public['loaded_library_identity_verified'])

    def test_saved_inspection_rejects_changed_inputs_and_record_binding(self):
        args = self.inspection_fixture()
        (self.root/'reference-replay').write_text('changed')
        public = self.inspect_saved(args)
        self.assertEqual(public['failed_check'], 'inputs_changed')
        launch = json.loads((self.result/'launch-start.json').read_text())
        launch['capture_sha256'] = '0'*64
        (self.result/'launch-start.json').write_text(json.dumps(launch))
        args[-1] = str(self.root/'unbound.json')
        public = self.inspect_saved(args, expected_status=1)
        self.assertEqual(public['failed_check'], 'saved_launch_record')

    def test_saved_inspection_exit_status_missing_malformed_and_signal(self):
        args = self.inspection_fixture()
        for text, code in ((None, None), ('SECRET', None), ('999', None), ('-15', -15), ('143', 143)):
            with self.subTest(text=text):
                status = self.result/'run.exit'
                if text is None: status.unlink()
                else: status.write_text(text)
                args[-1] = str(self.root/'inspect-{}.json'.format(text))
                public = self.inspect_saved(args)
                self.assertEqual(public['process_exit_code'], code)
                self.assertEqual(public['failed_check'], 'launcher_exit' if code is not None else 'launcher_exit_missing_or_invalid')

    def test_saved_inspection_classifies_scheduler_and_never_overwrites(self):
        args = self.inspection_fixture()
        (self.result/'run.exit').write_text('143\n')
        (self.result/'run.log').write_text('srun: error: SECRET node: Terminated\n'
                                         'slurmstepd: error: *** STEP SECRET DUE TO TIME LIMIT ***\n')
        public = self.inspect_saved(args)
        self.assertTrue(public['scheduler_signal_seen'])
        self.assertTrue(public['scheduler_time_limit_seen'])
        before = Path(args[-1]).read_bytes()
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(replay.main(args), 1)
        self.assertEqual(Path(args[-1]).read_bytes(), before)


class GpuProfileReplayTests(unittest.TestCase):
    def setUp(self):
        import test_simple_pressure_profile as fixtures
        self.fixture = fixtures.ProfileTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.fixture.finish()
        self.root = self.fixture.root
        self.pair = self.fixture.pair
        self.capture = self.root / 'frozen'
        self.capture.mkdir()
        write_system(self.capture / 'system', 2, rtol=1e-10)
        import simple_pressure_profile as profile
        baseline, launches, values, _ = profile.check_inputs(self.pair)
        self.source_values = values
        self.profile = profile
        explicit = replay.reference_controls(self.fixture.config, (values['atol'], values['rtol']))
        explicit.update(atol=0., rtol=1e-10)
        (self.capture / 'reference.settings').write_text(profile.text_profile(explicit))
        self.record = dict(ranks=2, libraries=self.fixture.libs, reference=launches['openaccel'])

    def test_verified_source_preserves_frozen_target_and_source(self):
        before = replay.hashes(list((self.capture / 'system').iterdir()) + [self.pair / 'pressure.profile'])
        values, inputs = replay.gpu_configuration(self.pair, self.capture, self.record)
        self.assertEqual(values, dict(self.source_values, atol=0., rtol=1e-10))
        replay.verify(before)
        replay.verify(inputs)
        self.assertIn(str(self.pair / 'mars/flow-pressure-rank-1.settings'), inputs)
        self.assertIn(str(self.fixture.baseline / 'reference/launch.json'), inputs)

    def test_wrong_rank_library_reference_or_algorithm_rejected(self):
        import copy
        for fault in ('ranks', 'libraries', 'reference', 'controls'):
            with self.subTest(fault=fault):
                record = copy.deepcopy(self.record)
                path = self.capture / 'reference.settings'
                saved = path.read_text()
                if fault == 'ranks': record['ranks'] = 4
                if fault == 'libraries': record['libraries'] = {}
                if fault == 'reference': record['reference']['executable_sha256'] = '0'*64
                if fault == 'controls': path.write_text(saved.replace('method 0', 'method 1'))
                with self.assertRaises(ValueError):
                    replay.gpu_configuration(self.pair, self.capture, record)
                path.write_text(saved)

    def binding(self):
        values, inputs = replay.gpu_configuration(self.pair, self.capture, self.record)
        output = self.root / 'new-replay'
        output.mkdir()
        config = output / 'gpu-reference.settings'
        config.write_text(self.profile.text_profile(values))
        inputs.update(replay.hashes([config]))
        record = dict(backend='mars', gpu_profile_source=str(self.pair / 'pressure.profile'),
                      inputs=inputs, command=['exe', str(self.capture / 'system'), str(config), str(output / 'result')])
        return values, output, record

    def test_configuration_bound_to_profile_and_original_target(self):
        values, output, record = self.binding()
        self.assertEqual(replay.checked_gpu_configuration(output, record, self.capture), values)
        config = output / 'gpu-reference.settings'
        for key, value in (('rtol', 1e-4), ('coarsentype', 10), ('kdim', 19)):
            with self.subTest(key=key):
                config.write_text(self.profile.text_profile(dict(values, **{key: value})))
                record['inputs'].update(replay.hashes([config]))
                with self.assertRaises(ValueError):
                    replay.checked_gpu_configuration(output, record, self.capture)
        config.write_text(self.profile.text_profile(values))
        record['inputs'].update(replay.hashes([config]))
        del record['inputs'][record['gpu_profile_source']]
        with self.assertRaises(ValueError):
            replay.checked_gpu_configuration(output, record, self.capture)

    def test_actual_cycle_settings_required(self):
        values = dict(self.source_values)
        report = dict(values, effective_levels=3, effective_relax_3=18)
        report['effective_relax_1'] = report.pop('relax_down')
        report['effective_relax_2'] = report.pop('relax_up')
        self.assertTrue(replay.controls_match(report, values))
        for key in ('effective_relax_1', 'effective_relax_2', 'effective_relax_3', 'effective_levels', 'coarserelax', 'kdim', 'rtol'):
            with self.subTest(key=key):
                self.assertFalse(replay.controls_match(dict(report, **{key: -1}), values))
                absent = dict(report); del absent[key]
                self.assertFalse(replay.controls_match(absent, values))

    def test_real_profile_setter_order_when_requested(self):
        binaries = os.environ.get('MARS_TEST_PRESSURE_REPLAYS', '').split(os.pathsep)
        if not binaries[0]:
            self.skipTest('set MARS_TEST_PRESSURE_REPLAYS to locally built replay executables')
        checker = os.environ.get('MARS_TEST_PRESSURE_CHECKER')
        if not checker:
            checker = str(self.root / 'checker')
            subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++17', '-O2', '-ffp-contract=off',
                            str(replay.CPP.with_name('pressure_residual_check.cpp')), '-o', checker], check=True)
        system = self.root / 'one-rank'
        n = 128
        matrix = [[2.01 if i == j else -1. if abs(i-j) == 1 else 0. for j in range(n)] for i in range(n)]
        write_system(system, matrix=matrix, rtol=1e-10)
        config = self.root / 'gpu.settings'
        for binary in binaries:
            for method in (0, 1):
                for levels in (1, 25):
                    with self.subTest(binary=binary, method=method, levels=levels):
                        values = dict(self.source_values, method=method, maxlevels=levels, rtol=1e-10, atol=0., maxiter=200)
                        config.write_text(self.profile.text_profile(values))
                        result = self.root / 'actual-result'
                        process = subprocess.run([binary, str(system), str(config), str(result)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        self.assertEqual(process.returncode, 0, process.stderr.decode())
                        report = replay.numeric_file(result / 'rank-000000.report')
                        self.assertTrue(replay.controls_match(report, values))
                        self.assertEqual(report['effective_levels'] == 1, levels == 1)
                        check = subprocess.run([checker, str(system), str(result)], check=True, stdout=subprocess.PIPE)
                        self.assertTrue(json.loads(check.stdout)['residual_passed'])
                        shutil.rmtree(result)


if __name__ == '__main__':
    unittest.main()
