"""Check workload fairness, failure retention and the exact submission wiring."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

import driver
from plan import VARIANTS, cases, depth, dimensions, mars_args, mfem_args, validate_case
from common import save, digest
import analyze

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


class PaperChecks(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=REPO / '.review/bp3-paper-20261107')
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_workload_and_alignment(self):
        self.assertEqual(len(cases('standard')), 14)
        self.assertEqual(len(cases('large')), 2)
        for case in cases('standard') + cases('large'):
            kind, ranks, p, ne = case
            row = validate_case(ranks, p, ne)
            self.assertEqual(row['element_local_dofs_per_gpu'] * (ranks if kind == 'strong' else 1),
                             134217728 if kind == 'strong' else 16777216)
        self.assertEqual(dimensions(256), [8, 8, 4])
        with self.assertRaises(ValueError):
            validate_case(256, 3, 32)
        for p in (3, 7):
            validate_case(256, p, 64)

    def test_matching_and_distinct_preconditioners(self):
        for variant in VARIANTS[:-1]:
            args = mars_args(self.root, variant, 3, 128)
            opts = lambda key: args[args.index(key) + 1]
            self.assertEqual(opts('--h-levels'), '2' if variant == 'best_h' else '0')
            self.assertEqual(opts('--amg-consolidate'), '0' if variant == 'matched' else '-1')
            self.assertEqual(opts('--cheb-interval'), 'mfem' if variant == 'matched' else 'mars')
            self.assertEqual(opts('--solves'), '5')
            self.assertEqual(opts('--true-tol'), '1e-10')
        self.assertIn('--gpu-aware-mpi', mfem_args(3, 128))
        self.assertEqual(depth('best_h', 7), 1)

    def test_mfem_receipt_and_failures(self):
        text = """MFEM 4.8 BP3, p=3, q=5
  4 ranks (2 x 2 x 1), serial mesh 32^3 refined 0 times, GPU-aware MPI on
  Hypre: MARS's settings (SpMV and SpGEMM hypre, GPU-aware MPI, pool 128 MiB)
  CG + V-cycle: 16 iterations, sqrt(<r,z>/<r0,z0>) = 8.0e-14, true residual 9.0e-11
  max |x - u| / max |u| = 1.0e-12 PASS
  solves 500.0 401.0 400.0 399.0 398.0 ms (the last is reported)
  CG iteration 24.875 ms (solve 398.0 ms, + setup 700 ms)
"""
        log = self.root / 'solver.log'
        log.write_text(text)
        row = driver.mfem_parse(self.root)
        self.assertEqual(row['wall_ms'], 398)
        self.assertEqual(row['iterations'], 16)
        for bad in (text.replace('9.0e-11', '1.1e-10'), text.replace('8.0e-14', 'nan'),
                    text.replace('PASS', 'FAIL'), text.replace(' 398.0 ms (the last', ' ms (the last'),
                    text.replace('GPU-aware MPI on', 'GPU-aware MPI off'), text.replace('16 iterations', '2000 iterations')):
            log.write_text(bad)
            with self.assertRaises(ValueError):
                driver.mfem_parse(self.root)

    def test_rotating_rounds_and_retained_failure(self):
        calls = []
        def solve(bundle, prepared, root, name, variant, ranks, p, ne, **kw):
            calls.append(variant)
            if variant == 'mfem':
                raise RuntimeError('exit 139 retained')
            return dict(hosts=['node'] * 4, wall_ms=100., iterations=12,
                        setup_ms=20., wall_solves_ms=[110., 100., 100., 100., 100.], true_residual=1e-12)
        with mock.patch.object(driver, 'solve', solve):
            row = driver.sample_set(self.root, self.root, self.root, 'strong', 4, 3, 128)
        self.assertEqual(calls[:4], list(VARIANTS))
        self.assertEqual(calls.count('mfem'), 1)
        self.assertEqual(len(row['failures']), 1)
        self.assertNotIn('mfem', row['summary'])
        self.assertEqual(len(row['samples']['best_h']), 3)

    def test_large_pilot_blocks_timing_on_failure(self):
        def solve(*args, **kw):
            if args[3] == 'pilot_check_p3':
                raise ValueError('direct check failed')
            return {}
        with mock.patch.object(driver, 'verify_prepared'), mock.patch.object(driver, 'solve', solve), \
                mock.patch.object(driver, 'sample_set') as samples:
            with self.assertRaisesRegex(ValueError, 'direct check'):
                driver.benchmark(self.root, self.root, self.root, 'large')
            samples.assert_not_called()

    def test_large_complete_dispatch(self):
        calls = []
        def solve(bundle, prepared, root, name, variant, ranks, p, ne, checked=False):
            calls.append((name, ranks, p, ne, checked))
            row = dict(iterations=14, hierarchy=[1, 100, 1000, 150 + (12 if checked else 0)])
            if checked:
                row['h_check_worst'] = 1e-16
            out = root / name
            out.mkdir()
            save(out / 'summary.json', row)
            return row
        with mock.patch.object(driver, 'verify_prepared'), mock.patch.object(driver, 'verify_bundle'), \
                mock.patch.object(driver, 'solve', solve), mock.patch.object(driver, 'seal'), \
                mock.patch.object(driver, 'digest', return_value='test'), \
                mock.patch.object(driver, 'sample_set', return_value={'failures': []}) as sample:
            driver.benchmark(self.root, self.root, self.root, 'large')
        self.assertEqual(len(calls), 10)
        self.assertTrue(all(c[1] == 256 and c[3] == 64 for c in calls))
        self.assertEqual(sum(c[4] for c in calls), 2)
        self.assertEqual(sample.call_count, 2)

    def test_summary_seals_and_scaling(self):
        case = dict(workload='weak', ranks=1, p=3, summary={
            'matched': dict(median_ms=100, min_ms=99, max_ms=101, iterations=[14]*3),
            'mfem': dict(median_ms=125, min_ms=124, max_ms=126, iterations=[14]*3)})
        path = self.root / 'results.json'
        save(path, dict(manifest_sha256='same', cases=[case], failures=[]))
        save(self.root / 'evidence.json', dict(manifest_sha256='same', evidence={'results.json': digest(path)}))
        text = analyze.summarize([self.root])
        self.assertIn('100.0%', text)
        self.assertIn('1.250x', text)
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            analyze.summarize([self.root, self.root])
        path.write_text('{}')
        with self.assertRaises(ValueError):
            analyze.summarize([self.root])

    def test_spooled_jobs_and_dependencies(self):
        bundle = self.root / 'bundle'
        bundle.mkdir()
        bindir = self.root / 'bin'
        bindir.mkdir()
        capture = self.root / 'calls'
        for name in ('sbatch', 'python'):
            p = bindir / name
            p.write_text('#!/bin/bash\nprintf "%s\\n" "$*" >> "$CAPTURE"\n' +
                         ('printf "12345\\n"\n' if name == 'sbatch' else ''))
            p.chmod(0o755)
        env = dict(os.environ, PATH=str(bindir) + ':' + os.environ['PATH'], CAPTURE=str(capture), SLURM_JOB_ID='12345')
        for name in ('submit.sh', 'job.sh'):
            text = (HERE / name).read_text().replace('/capstor/scratch/cscs/gandanie/bp3-paper-kit-20261107', str(bundle)).replace('/usr/bin/python3', str(bindir / 'python'))
            p = self.root / name
            p.write_text(text)
            subprocess.run(['bash', '-n', str(p)], check=True)
        subprocess.run(['bash', str(self.root / 'submit.sh')], env=env, check=True)
        spool = self.root / 'slurm_script'
        spool.write_bytes((self.root / 'job.sh').read_bytes())
        subprocess.run(['bash', str(spool), 'standard', 'prepared'], env=env, check=True)
        lines = capture.read_text().splitlines()
        self.assertEqual(len(lines), 5)
        self.assertIn('-N2 -n8', lines[1])
        self.assertIn('-N16 -n64', lines[2])
        self.assertIn('--dependency=afterok:12345,afterany:12345', lines[3])
        self.assertIn('-N64 -n256', lines[3])
        self.assertIn(str(bundle / 'paper/driver.py') + ' standard', lines[4])


if __name__ == '__main__':
    unittest.main()
