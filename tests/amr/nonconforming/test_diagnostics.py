import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

import compare_cells
import run
import verify


def passing_log(element='tet', mode='full'):
    rows = []
    for level in range(4):
        prefix = 'amr-mesh-only' if mode == 'mesh' else 'amr-galerkin'
        rows.append(f'[{prefix}] level={level} elements=100')
        if level and element == 'tet':
            rows += ['[amr-mesh-audit] stage=connectivity-rebuilt count=100', ' expected=100']
        if mode == 'full':
            for field in ('adjoint_relative=', 'symmetry_relative=', 'true_residual_assembled='):
                rows.append(f'[amr-algebra-audit] level={level} {field}0')
            for trial in range(3):
                rows.append(f'[amr-algebra-audit] level={level} trial={trial} operator_relative=0')
            for basis in range(4):
                rows.append(f'[amr-algebra-audit] level={level} basis={basis} max=0')
    return '\n'.join(rows + ['[amr-galerkin] PASS'])


class Diagnostics(unittest.TestCase):
    def test_case_inventory(self):
        cases = run.cases()
        self.assertEqual(len(cases), 9)
        self.assertEqual(len(set(c[0] for c in cases)), 9)
        self.assertEqual(sum(c[1] == 8 and c[2] == 'tet' for c in cases), 6)

    def test_complete_controls(self):
        for element, mode in [('tet', 'mesh'), ('tet', 'full'), ('hex', 'full')]:
            self.assertEqual(run.assess(passing_log(element, mode), 0, element, mode), 'PASS')

    def test_failure_and_truncation(self):
        text = passing_log()
        for code in (1, 77, 78, 143):
            self.assertEqual(run.assess(text, code, 'tet', 'full'), 'FAIL')
        for marker in ('[amr-mesh-audit] FAIL', 'error=nan', 'error=inf'):
            self.assertEqual(run.assess(text + '\n' + marker, 0, 'tet', 'full'), 'FAIL')
        for removed in ('level=2 elements=', 'trial=1 operator_relative=', 'basis=3 max=',
                        'symmetry_relative=', 'stage=connectivity-rebuilt ', '[amr-galerkin] PASS'):
            self.assertEqual(run.assess(text.replace(removed, 'missing'), 0, 'tet', 'full'), 'INCOMPLETE')

    def test_exact_cells_and_incomplete_rank(self):
        with tempfile.TemporaryDirectory(dir=os.environ['MARS_TEST_TMP']) as tmp:
            root = Path(tmp)
            for side in ('before', 'after'):
                for rank in range(2):
                    key = rank if side == 'before' else 1
                    (root / f'{side}-cstone-sync-rank{rank}.csv').write_text(
                        f'key0,key1,key2,key3,volume,invalid\n{key},2,3,4,0.125,0\n')
            result = compare_cells.compare(root, 2)
            self.assertEqual(result['before_count'], result['after_count'])
            self.assertEqual(result['removed'], [{'keys': [0, 2, 3, 4], 'count': 1}])
            self.assertEqual(result['added'], [{'keys': [1, 2, 3, 4], 'count': 1}])
            (root / 'after-cstone-sync-rank1.csv').unlink()
            with self.assertRaises(ValueError):
                compare_cells.compare(root, 2)

    def test_seal_refusals(self):
        with tempfile.TemporaryDirectory(dir=os.environ['MARS_TEST_TMP']) as tmp:
            root = Path(tmp)
            (root / 'source').write_text('original')
            (root / 'manifest.json').write_text(json.dumps({'files': {
                'source': hashlib.sha256(b'original').hexdigest()}}))
            verify.verify(root)
            (root / 'source').write_text('changed')
            with self.assertRaises(ValueError): verify.verify(root)
            (root / 'source').unlink()
            with self.assertRaises(ValueError): verify.verify(root)
            (root / 'source').write_text('original')
            (root / 'extra').write_text('extra')
            with self.assertRaises(ValueError): verify.verify(root)
            (root / 'extra').unlink()
            (root / 'link').symlink_to('source')
            with self.assertRaises(ValueError): verify.verify(root)


    def test_mesh_directory_hash(self):
        with tempfile.TemporaryDirectory(dir=os.environ['MARS_TEST_TMP']) as tmp:
            root = Path(tmp)
            with self.assertRaises(ValueError): run.mesh_hash(root)
            (root / 'x.raw').write_bytes(b'coordinates')
            first = run.mesh_hash(root)
            (root / 'i0.raw').write_bytes(b'connectivity')
            self.assertNotEqual(first, run.mesh_hash(root))
            (root / 'link').symlink_to('x.raw')
            with self.assertRaises(ValueError): run.mesh_hash(root)

    def test_batch_preserves_build_failure(self):
        self.check_batch(build_status=13, run_status=0, expected=13)

    def test_batch_preserves_case_failure(self):
        self.check_batch(build_status=0, run_status=1, expected=1)

    def test_batch_success(self):
        self.check_batch(build_status=0, run_status=0, expected=0)

    def check_batch(self, build_status, run_status, expected):
        with tempfile.TemporaryDirectory(dir=os.environ['MARS_TEST_TMP']) as tmp:
            root = Path(tmp); kit = root / 'kit'; kit.mkdir(); mock = root / 'bin'; mock.mkdir()
            (kit / 'manifest.json').write_text('{}')
            (mock / 'srun').write_text(f'#!/bin/sh\nexit {build_status}\n')
            (mock / 'python3').write_text(f'#!/bin/sh\ncase "$2" in */run.py) exit {run_status};; esac\nexit 0\n')
            for p in mock.iterdir(): p.chmod(0o755)
            env = dict(os.environ, PATH=str(mock) + ':' + os.environ['PATH'], SLURM_JOB_ID='42',
                       SLURM_JOB_NODELIST='mock[1-2]')
            completed = subprocess.run(['bash', str(Path(__file__).with_name('batch.sh')), str(kit)], env=env,
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            self.assertEqual(completed.returncode, expected, completed.stdout)
            self.assertEqual((root / 'kit-results-42/batch-exit.txt').read_text().strip(), str(expected))


if __name__ == '__main__':
    unittest.main()
