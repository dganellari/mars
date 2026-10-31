#!/usr/bin/env python3
"""Parser gates, level attribution and trace double-counting controls."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('profile_analysis', Path(__file__).with_name('analyze.py'))
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


class Checks(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self): self.tmp.cleanup()

    def csv(self, report, headers, rows):
        with (self.root / ('rank0.'+report+'.csv')).open('w') as f:
            f.write('Processing report...\n')
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(rows)

    def test_level_accounting(self):
        self.csv('nvtx_pushpop_trace', ['Name', 'Start (ns)', 'End (ns)', 'PID', 'TID'], [
            ['bp3.coarse:bp3.amg/cycle', 0, 1000, 8, 9],
            ['bp3.coarse:bp3.amg/L00', 0, 400, 8, 9],
            ['bp3.coarse:bp3.amg/L00/Relaxation', 10, 200, 8, 9],
            ['bp3.coarse:bp3.amg/L01', 400, 800, 8, 9],
            ['bp3.coarse:bp3.amg/L00', 800, 1000, 8, 9]])
        self.csv('nvtx_kern_sum', ['NVTX Range', 'Total Time (ns)', 'Kern Inst'], [
            ['bp3.amg/L00', 200, 4], ['bp3.amg/L00/Relaxation', 100, 2],
            ['bp3.amg/L01', 80, 2], ['bp3.amg/cycle', 280, 6]])
        self.csv('nvtx_gpu_proj_sum', ['Range', 'Total Proj Time (ns)'], [
            ['bp3.amg/L00', 350], ['bp3.amg/L01', 200]])
        self.csv('mpi_event_trace', ['Start (ns)', 'End (ns)', 'Event', 'Pid', 'Tid'], [
            [410, 460, 'MPI_Waitall', 8, 9], [500, 550, 'MPI_Sendrecv', 8, 9],
            [500, 550, 'MPI_Sendrecv', 8, 9], [900, 920, 'MPI_Allreduce', 8, 9],
            [1001, 1010, 'MPI_Allreduce', 8, 9]])
        data = analysis.trace_summary(self.root, 0)
        self.assertEqual(data['levels'][0]['kernels'], 4)
        self.assertAlmostEqual(data['levels'][0]['kernel_ms_per_amg_cycle'], 200e-6)
        self.assertAlmostEqual(data['levels'][0]['host_ms_per_amg_cycle'], 600e-6)
        self.assertAlmostEqual(data['levels'][1]['mpi_ms_per_amg_cycle'], 100e-6)
        self.assertAlmostEqual(data['levels'][1]['mpi_wait_ms_per_amg_cycle'], 50e-6)
        self.assertEqual(data['mpi_events_outside_levels'], 1)
        (self.root/'rank0.mpi_event_trace.csv').write_text('')
        (self.root/'rank0.mpi_event_trace.stderr').write_text('SKIPPED: report does not contain MPI events')
        single = analysis.trace_summary(self.root, 0, require_mpi=False)
        self.assertEqual(single['levels'][1]['mpi_ms_per_amg_cycle'], 0)
        with self.assertRaises(ValueError): analysis.trace_summary(self.root, 0, require_mpi=True)

    def test_missing_markers_fail(self):
        self.csv('nvtx_pushpop_trace', ['Name', 'Start', 'End', 'PID', 'TID'], [['unrelated', 0, 1, 8, 9]])
        with self.assertRaises(ValueError): analysis.trace_summary(self.root, 0)

    def test_audits(self):
        audit = dict(pushes=100, pops=100, cycles=3, errors=0, levels=[6, 3])
        path = self.root/'rank0.audit.json'
        path.write_text(json.dumps(audit))
        analysis.check_audits(self.root, 1)
        for field, value in [('cycles', 0), ('errors', 1), ('pops', 99)]:
            bad = {**audit, field: value}
            path.write_text(json.dumps(bad))
            with self.assertRaises(ValueError): analysis.check_audits(self.root, 1)
        path.unlink()
        with self.assertRaises(ValueError): analysis.check_audits(self.root, 1)

    def test_accuracy(self):
        good = ('FCG + V-cycle: 15 iterations, sqrt(<r,z>/<r0,z0>) = 5.000e-14, '
                'true assembled residual 5.000e-11\n'
                '  it 10 1.234567e-08\n  max error 1e-12 PASS\n'
                '  FCG iteration 12.3 ms (solve 184.5 ms, + setup 400 ms)\n')
        path = self.root/'solver.log'
        path.write_text(good)
        self.assertEqual(analysis.result(path)['iterations'], 15)
        for bad in (good.replace('5.000e-11', '1.05e-10'), good.replace('PASS', 'FAIL'),
                    good.replace('5.000e-14', 'nan'), good.replace('15 iterations', '2000 iterations')):
            path.write_text(bad)
            with self.assertRaises(ValueError): analysis.result(path)

    def test_profiler_preload_is_preserved(self):
        nsys = self.root/'nsys'
        nsys.write_text('''#!/usr/bin/env python3
import os, subprocess, sys
if sys.argv[1] == 'profile':
    start = sys.argv.index('bash')
    sys.exit(subprocess.call(sys.argv[start:], env={**os.environ, 'LD_PRELOAD':'nsys_fake.so'}))
print('mock report')
''')
        binary = self.root/'solver'
        binary.write_text('''#!/usr/bin/env python3
import os, sys
assert os.environ['LD_PRELOAD'] == 'probe.so:nsys_fake.so'
assert sys.argv[1:] == ['--solves', '3', '--profile']
assert os.environ['BP3_RANGE_AUDIT'].endswith('/rank0.audit.json')
print('preload preserved')
''')
        nsys.chmod(0o755); binary.chmod(0o755)
        run = subprocess.run(['bash', str(Path(__file__).with_name('rank.sh')), 'trace', str(self.root),
                              'probe.so', str(binary), '--solves', '3'], universal_newlines=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             env={**os.environ, 'PATH':str(self.root)+os.pathsep+os.environ['PATH'],
                                  'SLURM_PROCID':'0', 'SLURM_NTASKS':'1'})
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertIn('preload preserved', run.stdout)


if __name__ == '__main__': unittest.main()
