import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import probe


def fixture(task):
    passed, buffers = [], []
    for rank in range(task['ranks']):
        receipt = {k: task[k] for k in ('ranks', 'mode', 'count', 'offset', 'words')}
        receipt.update(rank=rank, hostname='node' + str(rank // (task['ranks'] // 2) if task['layout'] == 'multi' else 0),
                       device_id='gpu' + str(rank), kokkos_async=task['build'] == 'original', exchanges=9, errors=0)
        passed.append(receipt)
        for cycle in range(3):
            for label in ('send', 'recv'):
                buffers.append(dict(receipt, cycle=cycle, label=label, pool_query=0, memory_type=2,
                    pooled=task['mode'] == 'async' or (task['mode'] == 'view' and task['build'] == 'original' and task['words'] * 8 + 128 >= 40000)))
    return passed, buffers


def output(passed, buffers):
    return '\n'.join(prefix + json.dumps(row) for prefix, rows in
        (('MARS_ALLOCATOR_PASS ', passed), ('MARS_ALLOCATOR_BUFFER ', buffers)) for row in rows)


class ProbeTests(unittest.TestCase):
    def test_runner_attempts_all_allocators_after_one_fails(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in ('original', 'compatible'):
                (root / name).mkdir()
                (root / name / 'input').write_text('dependency fixture')
            work = root / 'results'
            calls = []
            def execute(argv, log, env, timeout):
                if argv[0] == 'cmake':
                    if argv[1] == '--build':
                        build = Path(argv[2]); build.mkdir()
                        (build / 'mars_mpi_allocator_probe').write_text('binary fixture')
                    log.write_text('build fixture')
                    return 0
                self.assertEqual([a for a in argv if a.startswith('--time=')], ['--time=00:02:00'])
                task = dict(mode=argv[-4], words=int(argv[-3]), count=int(argv[-2]), offset=int(argv[-1]),
                    build=Path(argv[-5]).parent.name, ranks=int(argv[argv.index('-n')+1]),
                    layout='multi' if '--nodes=2' in argv else 'single')
                calls.append(task)
                if task['mode']=='async' and task['layout']=='multi':
                    log.write_text('retained transport failure'); return 97
                log.write_text(output(*fixture(task))); return 0
            argv = ['probe.py', '--original', str(root/'original'), '--compatible', str(root/'compatible'),
                    '--work', str(work), '--cxx', '/cxx', '--gpu-compiler', '/nvcc', '--launcher', '/srun',
                    '--launcher-arg=--time=00:05:00']
            with patch.object(sys, 'argv', argv), patch.object(probe, 'execute', side_effect=execute):
                self.assertEqual(probe.main(), 0)
            data = json.loads((work/'result.json').read_text())
            self.assertEqual(calls, list(probe.schedule()))
            self.assertEqual(len(data['records']), 60)
            self.assertEqual(sum(r['status']=='FAILED' for r in data['records']), 6)
            self.assertTrue(data['compatible_controls_pass'])

    def test_full_matrix_and_bounds(self):
        tasks = list(probe.schedule())
        self.assertEqual(len(tasks), 60)
        for task in tasks:
            self.assertLessEqual(task['offset'] + task['ranks'] * task['count'], task['words'])
            probe.audit(output(*fixture(task)), task)

    def test_wrong_memory_pool_shape_or_completion_rejected(self):
        task = dict(layout='multi', ranks=4, words=8192, count=2048, offset=0, build='original', mode='view')
        original = fixture(task)
        for key, value in [('pool_query', 1), ('pooled', False), ('memory_type', 1), ('words', 4096),
                           ('kokkos_async', False), ('offset', 1), ('count', 1024)]:
            passed, buffers = copy.deepcopy(original)
            buffers[0][key] = value
            with self.assertRaises(ValueError): probe.audit(output(passed, buffers), task)
        for key, value in [('errors', 1), ('exchanges', 8), ('device_id', 'gpu1'), ('rank', 1)]:
            passed, buffers = copy.deepcopy(original)
            passed[0][key] = value
            with self.assertRaises(ValueError): probe.audit(output(passed, buffers), task)
        passed, buffers = original
        for p, b in [(passed[:-1], buffers), (passed, buffers[:-1]), (passed, buffers + buffers[:1])]:
            with self.assertRaises(ValueError): probe.audit(output(p, b), task)


if __name__ == '__main__': unittest.main()
