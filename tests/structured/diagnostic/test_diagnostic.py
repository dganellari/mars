import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cluster'))
import diagnose


def trace(rank, profile=False):
    return {'status': 'complete', 'rank': rank, 'hostname': 'fixture', 'dropped': 0,
            'tools_loaded': profile, 'tools_finalized': profile, 'global_fencing': profile,
            'events': [{'kind': 'mpi', 'name': 'MPI_Allreduce', 'seconds': 0.001},
                       {'kind': 'mpi_error', 'name': 'MPI_Isend: invalid rank', 'seconds': 0},
                       {'kind': 'parallel_for', 'name': 'fixture', 'seconds': 0.002}]}


class DiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_schedule_is_focused_and_keeps_both_variants(self):
        self.assertEqual(len(diagnose.schedule('cpu')), 12)
        tasks = diagnose.schedule('cuda')
        self.assertEqual(len(tasks), 42)
        self.assertEqual([(t['variant'], t['ranks']) for t in tasks[:4]],
                         [('reference', 2), ('candidate', 2), ('reference', 4), ('candidate', 4)])
        for backend in ('cpu', 'cuda'):
            controls = [t for t in diagnose.schedule(backend) if t['mode'] == 'control']
            for index in range(0, len(controls), 2):
                pair = controls[index:index+2]
                self.assertEqual({t['variant'] for t in pair}, {'reference', 'candidate'})
                self.assertEqual(pair[0]['variant'], 'candidate' if pair[0]['repeat'] % 2 else 'reference')

    def test_missing_inactive_truncated_and_incomplete_observers_reject(self):
        prefix = self.root / 'trace'
        path = Path(str(prefix) + '.rank0.json')
        path.write_text(json.dumps(trace(0, True)))
        self.assertEqual(len(diagnose.observer(prefix, 1, 'profile', True)), 1)
        for changes in ({'events': []}, {'tools_loaded': False}, {'tools_finalized': False},
                        {'global_fencing': False}, {'dropped': 1}, {'rank': 1}, {'status': 'aborting'}):
            path.write_text(json.dumps(dict(trace(0, True), **changes)))
            with self.assertRaises(ValueError):
                diagnose.observer(prefix, 1, 'profile', True)
        path.unlink()
        with self.assertRaises(FileNotFoundError):
            diagnose.observer(prefix, 1, 'profile', True)

    def test_original_comparison_limit_and_round_inventory(self):
        records = []
        for task in diagnose.schedule('cpu'):
            if task['mode'] != 'control':
                continue
            record = dict(task, status='PASS', audit={'cases': {}})
            for case in diagnose.performance.CASES:
                record['audit']['cases'][case] = {m: 0.002084 if task['variant'] == 'candidate' else 0.00181
                                                for m in diagnose.performance.METRICS}
            records.append(record)
        compared = diagnose.compare_controls(records)
        self.assertEqual(len(compared), 18)
        self.assertTrue(all(not row['pass'] for row in compared))
        records[0]['repeat'] = 1
        self.assertEqual(diagnose.compare_controls(records)[0]['status'], 'INCOMPLETE')

    def test_real_launcher_and_rank_wrapper_keep_the_active_environment(self):
        stub = self.root / 'srun'
        destination = self.root / 'launcher.json'
        stub.write_text('#!' + sys.executable + '\nimport json,os,sys\n'
                        'open(os.environ["LAUNCH_RECORD"],"w").write(json.dumps(sys.argv))\n'
                        'i=sys.argv.index("/bin/bash"); os.execv(sys.argv[i],sys.argv[i:])\n')
        stub.chmod(0o755)
        binary = self.root / 'benchmark'
        output = self.root / 'environment.json'
        binary.write_text('#!' + sys.executable + '\nimport json,os\n'
                          'open(os.environ["ENV_RECORD"],"w").write(json.dumps(dict(os.environ)))\n')
        binary.chmod(0o755)
        launcher = Path(diagnose.run.__file__).with_name('mpi_launch.py')
        template = [sys.executable, str(launcher), '-n', '4', '--launcher=' + str(stub),
                    '--launcher-arg=--cpu-bind=cores', '/prepared/mars_structured_consumer']
        cmd = diagnose.command(template, Path(__file__).with_name('rank.sh').resolve(),
                               self.root / 'unused.so', self.root / 'control', binary, 'control')
        env = dict(os.environ, LAUNCH_RECORD=str(destination), ENV_RECORD=str(output), SLURM_JOB_ID='live',
                   CUDA_VISIBLE_DEVICES='2,3', MARS_GATE_LAYOUT='multi', KOKKOS_TOOLS_LIBS='wrong')
        subprocess.run(cmd, env=env, check=True, capture_output=True)
        argv = json.loads(destination.read_text())
        for flag in ('--nodes=2', '--ntasks-per-node=2', '--export=ALL', '--time=00:05:00'):
            self.assertIn(flag, argv)
        received = json.loads(output.read_text())
        self.assertEqual(received['SLURM_JOB_ID'], 'live')
        self.assertEqual(received['CUDA_VISIBLE_DEVICES'], '2,3')
        self.assertNotIn('KOKKOS_TOOLS_LIBS', received)

    def test_input_mutation_is_rejected(self):
        path = self.root / 'binary'
        path.write_text('original')
        files = {str(path): diagnose.run.digest(path)}
        diagnose.verify_files(files)
        path.write_text('changed')
        with self.assertRaises(ValueError):
            diagnose.verify_files(files)

    def test_failed_reference_does_not_hide_candidate_or_other_tasks(self):
        kit = self.root / 'kit'
        kit.mkdir()
        (kit / 'manifest.json').write_text(json.dumps({'files': {}}))
        configurations = {}
        for config in ('cuda-direct', 'cuda-staged'):
            configurations[config] = {'threads': 1, 'environment': {},
                'binaries': {'reference': '/reference', 'candidate': '/candidate'},
                'templates': {str(n): ['launcher', '-n', str(n), '/prepared/mars_structured_consumer'] for n in (1, 2, 4)}}
        pins = {'source': '/source', 'source_manifest_sha256': 'source', 'reference': '/reference',
                'reference_manifest_sha256': 'reference', 'kokkos_tools_include': '/include',
                'cuda': {'files': {}, 'compiler': '/compiler', 'configurations': configurations}}
        (kit / 'pins.json').write_text(json.dumps(pins))
        work = self.root / 'work'
        executed = []
        def execute(argv, log, env, timeout=330):
            log.write_text('fixture')
            if log.name == 'build.log':
                build = work / 'observer-build'
                build.mkdir()
                for name in ('libmars_structured_trace.so', 'mars_mpi_trace_probe'):
                    (build / name).write_text('fixture')
            if log.name == 'probe.log':
                for rank in range(2):
                    (work / ('probe.rank' + str(rank) + '.json')).write_text(json.dumps(trace(rank)))
            if log.name.startswith('cuda-'):
                executed.append(log.stem)
                return 15 if '-reference' in log.name and '-mpi-' in log.name else 0
            return 0
        audit = {'cases': {case: {metric: 1.0 for metric in diagnose.performance.METRICS}
                            for case in diagnose.performance.CASES}}
        real_observer = diagnose.observer
        def observed(prefix, ranks, mode, probe=False):
            return real_observer(prefix, ranks, mode, probe) if probe else []
        with patch.object(diagnose, '__file__', str(kit / 'diagnose.py')), \
             patch.object(sys, 'argv', ['diagnose.py', '--backend', 'cuda', '--work', str(work)]), \
             patch.object(diagnose.run, 'verify_snapshot', side_effect=lambda path: path.name), \
             patch.object(diagnose, 'execute', side_effect=execute), \
             patch.object(diagnose, 'observer', side_effect=observed), \
             patch.object(diagnose.performance, 'audit', return_value=audit):
            self.assertEqual(diagnose.main(), 0)
        result = json.loads((work / 'result.json').read_text())
        self.assertEqual(len(executed), 42)
        self.assertEqual(len(result['records']), 42)
        self.assertEqual(result['records'][0]['status'], 'FAILED')
        self.assertEqual(result['records'][1]['variant'], 'candidate')
        self.assertEqual(result['records'][1]['status'], 'PASS')
        self.assertEqual(result['status'], 'DIAGNOSTIC_COLLECTED')


if __name__ == '__main__':
    unittest.main()
