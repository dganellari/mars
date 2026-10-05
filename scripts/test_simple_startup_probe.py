"""Startup evidence checks use synthetic fields only."""
import contextlib
import copy
import csv
import io
import json
import os
from pathlib import Path
import shutil
import sys
import unittest
from unittest.mock import patch

import numpy as np
from netCDF4 import Dataset

import simple_startup_probe as probe
import test_simple_snapshot_compare as fixtures


class StartupTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.SnapshotTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.tearDown)
        f = self.fixture
        self.pair = f.root / 'startup'
        probe.prepare(f.case, f.reference, self.pair)
        self.mars = self.pair / 'mars'
        self.mars.mkdir()
        f.log.write_text(f.log.read_text().replace('iterations=100', 'iterations=20'))
        shutil.copy(f.log, self.mars / 'run.log')
        shutil.copy(f.exit, self.mars / 'run.exit')
        # Preserve a complete 0..20 metrics history and both gathered/partitioned readers.
        lines = f.metrics.read_text().splitlines()
        (self.mars / 'flow-metrics.csv').write_text('\n'.join(lines[:22]) + '\n')
        fields_stream = io.StringIO()
        writer = csv.writer(fields_stream)
        writer.writerow(['node', 'x', 'y', 'z', 'u', 'v', 'w', 'p'])
        for n in range(len(f.ids)):
            writer.writerow([n] + list(f.xyz[n]) + list(f.fields[n]))
        fields = fields_stream.getvalue().encode()
        for iteration in range(21):
            (self.mars / ('flow-step-' + str(iteration) + '-fields.csv')).write_bytes(fields)
        (self.mars / 'flow-fields.csv').write_bytes(fields)
        reference = self.pair / 'reference'
        (reference / 'run.log').write_text('\n'.join('Iter = ' + str(i) for i in range(1, 21)) + '\nSimulation is complete\n')
        (reference / 'run.exit').write_text('0\n')
        f.write_exodus(reference / 'results.e.2.0', [4, 0, 2, 1], times=list(range(21)))
        f.write_exodus(reference / 'results.e.2.1', [3, 1, 5, 2], times=list(range(21)), combined=True)
        for rank, nodes in enumerate(([4, 0, 2, 1], [3, 1, 5, 2])):
            with Dataset(str(reference / ('results.e.2.' + str(rank))), 'a') as ds:
                if rank:
                    ds.variables['vals_nod_var'][:] = np.tile(f.fields[nodes].T, (21, 1, 1))
                else:
                    for k in range(4):
                        ds.variables['vals_nod_var' + str(k+1)][:] = np.tile(f.fields[nodes, k], (21, 1))
        for solver in ('mars', 'openaccel'):
            self.record(solver)
        self.public = f.root / 'startup-public.json'

    def record(self, solver):
        directory = self.pair / ('reference' if solver == 'openaccel' else 'mars')
        executable = self.fixture.root / 'fake-executable'
        executable.write_text('synthetic executable identity')
        record = dict(schema=probe.SCHEMA, solver=solver, ranks=2,
                      command=['launcher', str(executable)] + probe.solver_arguments(self.pair, solver),
                      executable=str(executable), executable_sha256=probe.digest(executable),
                      pair_sha256=probe.digest(self.pair / 'pair.json'), libraries={'synthetic-library': 'a'*64},
                      environment={}, status='started')
        (directory / 'launch-start.json').write_text(json.dumps(record))
        record.update(status='finished', exit_code=int((directory / 'run.exit').read_text()),
                      files=probe.run_files(directory, solver))
        (directory / 'launch.json').write_text(json.dumps(record))

    def compare(self):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            code = probe.main(['compare', '--pair', str(self.pair), '--output', str(self.public)])
        result = json.loads(self.public.read_text())
        self.assertNotIn(str(self.fixture.root), self.public.read_text())
        return code, result

    def test_equal_histories(self):
        code, result = self.compare()
        self.assertEqual(code, 0)
        self.assertTrue(result['all_twenty_snapshots_match'])
        self.assertIsNone(result['first_differing_iteration'])
        self.assertFalse(result['identical_linear_solvers_verified'])
        self.assertTrue(result['initial_field_parity_verified'])

    def test_first_difference_not_last(self):
        for step in (3, 7):
            path = self.mars / ('flow-step-' + str(step) + '-fields.csv')
            lines = path.read_text().splitlines()
            row = lines[1].split(',')
            row[-1] = str(float(row[-1]) + .01)
            lines[1] = ','.join(row)
            path.write_text('\n'.join(lines) + '\n')
        self.record('mars')
        code, result = self.compare()
        self.assertEqual(code, 0)
        self.assertEqual(result['first_differing_iteration'], 3)
        self.assertEqual(result['first_differing_fields'], ['pressure'])

    def test_initial_mismatch_is_not_attributed_to_first_iteration(self):
        path = self.mars / 'flow-step-0-fields.csv'
        lines = path.read_text().splitlines()
        row = lines[1].split(',')
        row[-1] = '1'
        lines[1] = ','.join(row)
        path.write_text('\n'.join(lines) + '\n')
        self.record('mars')
        code, result = self.compare()
        self.assertEqual(code, 0)
        self.assertEqual(result['first_differing_iteration'], 0)
        self.assertFalse(result['initial_field_parity_verified'])

    def test_reference_initial_snapshot_required(self):
        for path in probe.result_files(self.pair / 'reference'):
            with Dataset(str(path), 'a') as ds:
                ds.variables['time_whole'][0] = -.5
        self.record('openaccel')
        self.assertEqual(self.compare()[0], 1)

    def test_missing_or_tampered_evidence(self):
        files = ['pair.json', 'case.json', 'reference/input.i', 'reference/results.e.2.1',
                 'reference/launch-start.json', 'mars/launch.json', 'mars/run.log',
                 'mars/run.exit', 'mars/flow-step-8-fields.csv']
        for name in files:
            with self.subTest(name=name):
                path = self.pair / name
                original = path.read_bytes()
                path.unlink()
                try:
                    code, result = self.compare()
                    self.assertEqual(code, 1)
                    self.assertEqual(result['comparison_status'], 'invalid_evidence')
                finally:
                    self.public.unlink()
                    path.write_bytes(original)

    def test_launch_identity_and_exit(self):
        record = self.mars / 'launch.json'
        original = json.loads(record.read_text())
        for key, value in [('exit_code', 137), ('ranks', 4), ('command', ['wrong-executable']),
                           ('pair_sha256', '0'*64), ('status', 'started')]:
            with self.subTest(key=key):
                record.write_text(json.dumps(dict(original, **{key: value})))
                self.assertEqual(self.compare()[0], 1)
                self.public.unlink()
        record.write_text(json.dumps(original))

    def test_incomplete_history_does_not_become_first_difference(self):
        path = self.mars / 'flow-step-9-fields.csv'
        path.unlink()
        self.record('mars')
        code, result = self.compare()
        self.assertEqual(code, 1)
        self.assertNotIn('first_differing_iteration', result)

    def test_nan_and_missing_node_are_rejected(self):
        path = self.mars / 'flow-step-9-fields.csv'
        original = path.read_text()
        for text in (original.replace('0.03', 'nan'), '\n'.join(original.splitlines()[:-1]) + '\n'):
            path.write_text(text)
            self.record('mars')
            self.assertEqual(self.compare()[0], 1)
            self.public.unlink()

    def test_preparation_preserves_equations_and_linear_settings(self):
        modified = probe.load_deck((self.pair / 'reference/input.i').read_bytes())
        original = copy.deepcopy(self.fixture.doc)
        solver = original['simulation']['solver']
        solver['output_control'] = modified['simulation']['solver']['output_control']
        conv = solver['solver_control']['basic_settings']['convergence_controls']
        conv.update(min_iterations=20, max_iterations=20)
        original['mesh']['file_path'] = str(self.fixture.mesh.resolve())
        self.assertEqual(modified, original)

    def test_changed_mesh_is_rejected(self):
        with self.fixture.mesh.open('ab') as stream:
            stream.write(b'changed')
        self.assertEqual(self.compare()[0], 1)

    def test_no_overwrite(self):
        self.assertEqual(self.compare()[0], 0)
        before = self.public.read_bytes()
        self.assertEqual(self.compare()[0], 1)
        self.assertEqual(self.public.read_bytes(), before)

    def test_launch_records_real_subprocess(self):
        pair = self.fixture.root / 'launched'
        probe.prepare(self.fixture.case, self.fixture.reference, pair)
        exe = self.fixture.root / 'synthetic-solver.py'
        exe.write_text("from pathlib import Path\nPath('results.e').write_text('synthetic')\nprint('synthetic completion')\n")
        exe.chmod(0o700)
        with patch.object(probe, 'runtime_libraries', return_value={'test': 'a'*64}), patch.dict(
                os.environ, {'MARS_OPENACCEL_EXPORT_DIR': 'must-not-reach-child'}):
            probe.launch(pair, 'openaccel', exe, 1, [sys.executable])
        record = probe.verified_launch(pair, 'openaccel')
        self.assertEqual(record['exit_code'], 0)
        self.assertNotIn('MARS_OPENACCEL_EXPORT_DIR', record['environment'])
        self.assertIn('synthetic completion', (pair / 'reference/run.log').read_text())


if __name__ == '__main__':
    unittest.main()
