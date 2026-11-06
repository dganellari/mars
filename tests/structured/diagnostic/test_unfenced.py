import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cluster'))
import diagnose
import performance
import setup_trace
import unfenced


class UnfencedTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_timer_audit_requires_all_pairs_and_matches_reported_maximum(self):
        events = [{'kind': 'timer', 'name': 'MPI_Wtime', 'value': i / 100, 'start': i / 100, 'seconds': 0}
                  for i in range(528)]
        data = {'events': events, 'global_fencing': False}
        prefix = self.root / 'sample'
        path = self.root / 'sample.rank0.json'
        path.write_text(json.dumps(data))
        measured = {'cases': {c: {m: .01 for m in performance.METRICS} for c in performance.CASES}}
        self.assertEqual(len(setup_trace.audit(prefix, 1, measured)[0]['cases']), 6)
        measured['cases']['morton-quad-rect']['mesh_seconds'] += .0001
        with self.assertRaisesRegex(ValueError, 'differs from benchmark'):
            setup_trace.audit(prefix, 1, measured)
        with self.assertRaisesRegex(ValueError, '528'):
            setup_trace.intervals(events[:-1])
        with self.assertRaisesRegex(ValueError, 'interception was not active'):
            setup_trace.audit(prefix, 1, measured, True)
        data['global_fencing'] = True
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, 'added global fences'):
            setup_trace.audit(prefix, 1, measured)

    def test_host_mode_rejects_hidden_fencing_and_missing_callbacks(self):
        prefix = self.root / 'sample'
        path = self.root / 'sample.rank0.json'
        data = {'rank': 0, 'status': 'complete', 'dropped': 0, 'tools_loaded': True,
                'tools_finalized': True, 'global_fencing': False, 'hostname': 'fixture',
                'events': [{'kind': 'mpi', 'name': 'MPI_Allreduce', 'seconds': .001},
                           {'kind': 'parallel_for', 'name': 'fixture', 'seconds': .001}]}
        path.write_text(json.dumps(data))
        diagnose.observer(prefix, 1, 'host')
        for change in ({'global_fencing': True}, {'tools_loaded': False}, {'tools_finalized': False}):
            path.write_text(json.dumps(dict(data, **change)))
            with self.assertRaises(ValueError):
                diagnose.observer(prefix, 1, 'host')
        path.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, 'MPI-only'):
            diagnose.observer(prefix, 1, 'mpi-quiet')

    def test_profiles_cannot_pass_controls_or_hide_missing_rounds(self):
        rows = []
        for task in unfenced.schedule():
            value = .002 if task['variant'] == 'reference' else .003
            if task['observer'] != 'control':
                value = .0001
            rows.append(dict(task, status='PASS', audit={'cases': {
                c: {m: value for m in performance.METRICS} for c in performance.CASES}}))
        comparisons = unfenced.comparisons(rows)
        self.assertEqual(len(comparisons), 126)
        self.assertTrue(all(not row['passed'] for row in comparisons))
        self.assertEqual(unfenced.comparisons(rows[1:])[0]['status'], 'INCOMPLETE')

    def test_real_runner_retains_failures_and_finishes_all_modes(self):
        kit = self.root / 'kit'
        kit.mkdir()
        (kit / 'manifest.json').write_text(json.dumps({'files': {}}))
        prepared = self.root / 'prepared'
        binaries = {}
        for variant in ('reference', 'candidate'):
            path = self.root / variant
            path.write_text(variant)
            binaries[variant] = str(path)
        for transport in ('direct', 'staged'):
            config = prepared / ('cuda-' + transport)
            config.mkdir(parents=True)
            (config / 'result.json').write_text(json.dumps({'status': 'PASS', 'environment': {'SLURM_JOB_ID': 'old'}}))
            (config / 'consumer-list.log').write_text(json.dumps({'tests': [
                {'name': 'structured_install_mpi' + str(n), 'command': ['launcher', '/mars_structured_consumer']}
                for n in (1, 2, 4)]}))
            build = prepared / ('cuda-' + transport + '-performance-build')
            build.mkdir()
            (build / 'result.json').write_text(json.dumps({'status': 'BUILT_NOT_RUN', 'binary_paths': binaries}))
        (kit / 'pins.json').write_text(json.dumps({'prepared': str(prepared), 'cuda_include': '/cuda', 'cuda': {'files': {}}}))
        work = self.root / 'work'
        calls = []

        def execute(command, log, env, timeout=330):
            calls.append(command)
            self.assertEqual(env['SLURM_JOB_ID'], 'active')
            log.write_text('fixture')
            if command[:2] == ['cmake', '--build']:
                build = Path(command[2])
                build.mkdir()
                for name in ('libmars_structured_trace.so', 'mars_mpi_trace_probe'):
                    (build / name).write_text('fixture')
            return 9 if log.stem in ('preflight-host', 'direct-multi-np2-control-0-reference') else 0

        audit = {'cases': {c: {m: .001 for m in performance.METRICS} for c in performance.CASES}}
        with patch.object(unfenced, '__file__', str(kit / 'unfenced.py')), \
             patch.object(sys, 'argv', ['unfenced.py', '--work', str(work)]), \
             patch.dict(os.environ, {'SLURM_JOB_ID': 'active'}), \
             patch.object(diagnose, 'execute', side_effect=execute), \
             patch.object(diagnose, 'observer', return_value=[]), \
             patch.object(performance, 'audit', return_value=audit), \
             patch.object(setup_trace, 'audit', return_value=[]):
            self.assertEqual(unfenced.main(), 1)
        result = json.loads((work / 'result.json').read_text())
        self.assertEqual(len(calls), 131)
        self.assertEqual(len(result['records']), 126)
        self.assertEqual(sum(r['status'] == 'FAILED' for r in result['records']), 1)
        self.assertEqual(result['status'], 'COLLECTED_WITH_FAILURES')
        self.assertEqual(result['records'][-1]['status'], 'PASS')
        self.assertEqual(result['preflights'][-1]['status'], 'FAILED')


if __name__ == '__main__':
    unittest.main()
