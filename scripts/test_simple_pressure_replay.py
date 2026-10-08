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
                    cfg.write_text(''.join('{} {:.17g}\n'.format(k,v) for k,v in report.items() if not k.startswith('result_')))
                    shutil.rmtree(self.result)
                    p = subprocess.run([binary, str(self.capture), str(cfg), str(self.result)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    self.assertEqual(p.returncode, 0, p.stderr.decode())
                    self.assertTrue(self.check(self.result)['residual_passed'])
                    self.assertEqual(replay.numeric_file(self.result/'rank-000000.report'), report)
                    shutil.rmtree(self.result)

    def captured_record(self):
        write_system(self.capture)
        system = self.root/'captured-run'; system.mkdir()
        self.capture.rename(system/'system')
        self.capture = system/'system'
        files = replay.capture_inputs(self.capture)[2]
        (system/'reference.settings').write_text('method 0\nrtol 1e-10\natol 0\nmaxiter 200\n')
        files.update(replay.hashes([system/'reference.settings']))
        record = dict(schema=replay.SCHEMA, exit_code=255, ranks=1, files=files)
        replay.startup.write_json(system/'capture.json', record)
        return system

    def replay_record(self, capture, name, candidate):
        path = self.root/name; path.mkdir()
        write_solution(path/'result', candidate, 1)
        (path/'result/rank-000000.report').write_text(
            'method 0\nrtol 1e-10\natol 0\nmaxiter 200\nrelaxtype 18\ncoarserelax 18\nresult_converged 1\nresult_fatal_error 0\n')
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
        self.assertFalse(public['actual_exit_branch_verified'])
        self.assertFalse(public['original_amg_hierarchy_reused'])
        self.assertNotIn(str(self.root), json.dumps(public))

    def test_modified_replay_output_rejected(self):
        args = self.compare_arguments()
        (self.root/'mars/result/rank-000000.solution').write_bytes(b'changed')
        with contextlib.redirect_stdout(io.StringIO()): self.assertEqual(replay.main(args), 1)
        self.assertEqual(json.loads((self.root/'public.json').read_text())['failed_check'], 'replay_identity')

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
        # Changing a captured dependency must prevent a second launch.
        library.write_text('changed')
        args[-1] = str(self.root/'retry')
        with patch.object(replay.subprocess, 'run') as command, contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(replay.main(args), 1)
            command.assert_not_called()


if __name__ == '__main__':
    unittest.main()
