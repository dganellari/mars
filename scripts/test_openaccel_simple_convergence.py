"""Transport and comparison tests with synthetic fields, not CFD validation."""
import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from netCDF4 import Dataset
import openaccel_simple_convergence as gate

DECK = ('max_iterations: 2\nresidual_target: 1.0e-6\n'
        'rtol: 1.0e-8\natol: 1.0e-12\n')


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.reference = self.root / 'reference'
        self.mars = self.root / 'mars'
        self.reference.mkdir(); self.mars.mkdir()
        self.output = self.root / 'comparison.json'
        self.ids = np.arange(101, 526)
        self.xyz = np.column_stack([np.linspace(0, 4, 425), np.zeros((425, 2))])
        self.fields = np.column_stack([np.linspace(0, .18, 425), np.zeros((425, 2)), np.linspace(1, 0, 425)])
        (self.reference / 'input.i').write_text(gate.convergence_deck(DECK))
        (self.reference / 'channel.exo').write_text('synthetic public test input')
        (self.reference / 'run.log').write_text('Git hash: ' + gate.PIN + '\nIter = 1\nIter = 2\nConverged.\n')
        self.mesh_hash = gate.digest(self.reference / 'channel.exo')
        self.write_exodus()
        self.manifest()
        (self.mars / 'channel.txt').write_text('MARS_PUBLIC_SIMPLE_MESH_V1 425 1536 576\n')
        (self.mars / 'channel.json').write_text(json.dumps(dict(node_global_ids=self.ids.tolist(),
             packed_sha256=gate.digest(self.mars / 'channel.txt'))))
        with (self.mars / 'channel-fields.csv').open('w') as f:
            f.write('node,x,y,z,u,v,w,p\n')
            for i in range(425):
                f.write(','.join(map(str, [i] + self.xyz[i].tolist() + self.fields[i].tolist())) + '\n')
        (self.mars / 'channel-metrics.csv').write_text(
            'iteration,momentum,continuity,mass_balance,du,dp,dflux,cancellation,changed_faces\n'
            '1277,9.9e-7,3e-11,2e-12,2e-10,1e-9,4e-12,1e-16,0\n')
        (self.mars / 'run.log').write_text('CONVERGED iterations=1277\n')

    def write_exodus(self, combined=False, offset=0, drift=0, masked=False, wrong_id=False, stale=False):
        # Reverse file order: matching by row instead of global ID must fail this fixture.
        order = np.arange(424, -1, -1)
        with Dataset(str(self.reference / 'results.e'), 'w') as ds:
            for k, n in [('num_nodes', 425), ('num_dim', 3), ('time_step', 2), ('num_nod_var', 4), ('len_name', 33)]:
                ds.createDimension(k, n)
            ids = self.ids[order].copy()
            if wrong_id: ids[0] = ids[1]
            ds.createVariable('node_num_map', 'i8', ('num_nodes',))[:] = ids
            if combined:
                ds.createVariable('coord', 'f8', ('num_dim', 'num_nodes'))[:] = self.xyz[order].T
            else:
                for j, a in enumerate('xyz'):
                    ds.createVariable('coord' + a, 'f8', ('num_nodes',))[:] = self.xyz[order, j]
            ds.createVariable('time_whole', 'f8', ('time_step',))[:] = [0, 1] if stale else [1, 2]
            names = np.zeros((4, 33), dtype='S1')
            for i, s in enumerate(('velocity_x', 'velocity_y', 'velocity_z', 'pressure')):
                names[i, :len(s)] = np.frombuffer(s.encode(), dtype='S1')
            ds.createVariable('name_nod_var', 'S1', ('num_nod_var', 'len_name'))[:] = names
            values = np.stack([self.fields[order], self.fields[order]])
            values[:, :, 3] += offset
            values[0, :, 0] += drift
            if masked: values[1, 0, 0] = np.nan
            if combined:
                ds.createVariable('vals_nod_var', 'f8', ('time_step', 'num_nod_var', 'num_nodes'))[:] = values.transpose(0, 2, 1)
            else:
                for k in range(4):
                    ds.createVariable('vals_nod_var' + str(k+1), 'f8', ('time_step', 'num_nodes'))[:] = values[:, :, k]

    def manifest(self):
        record = dict(fixture='public_simple_convergence_v1', status='native_convergence_reported',
                      returncode=0, iteration=2, result_file='results.e', source_deck_sha256=gate.DECK_SHA256)
        for name, key in [('input.i', 'deck_sha256'), ('channel.exo', 'mesh_sha256'),
                          ('run.log', 'log_sha256'), ('results.e', 'result_sha256')]:
            record[key] = gate.digest(self.reference / name)
        (self.reference / 'comparison-run.json').write_text(json.dumps(record))

    def compare(self, native_mesh=None):
        with patch.object(gate, 'MESH_SHA256', self.mesh_hash), contextlib.redirect_stdout(io.StringIO()):
            gate.compare(self.reference, self.mars, self.output, native_mesh)

    def native_input(self, explicit_ids=True):
        path = self.reference / 'channel.exo'
        if not explicit_ids:
            self.ids = np.arange(1, 426)
            self.write_exodus()
        with Dataset(str(path), 'w') as mesh:
            mesh.createDimension('num_nodes', 425)
            mesh.createDimension('num_elem', 1536)
            if explicit_ids:
                mesh.createVariable('node_num_map', 'i8', ('num_nodes',))[:] = self.ids
        self.mesh_hash = gate.digest(path)
        self.manifest()
        (self.mars / 'run.log').write_text('Native SIMPLE Exodus input: ' + str(path.resolve())
                                         + '\nCONVERGED iterations=1277\n')
        # Native runtime order differs from the Exodus storage row written in column node.
        fields = self.mars / 'channel-fields.csv'
        lines = fields.read_text().splitlines()
        fields.write_text('\n'.join(lines[:1] + lines[:0:-1]) + '\n')
        (self.mars / 'channel.json').unlink()
        (self.mars / 'channel.txt').unlink()
        return path

    def test_native_permuted_rows_and_explicit_ids(self):
        self.compare(self.native_input())
        report = json.loads(self.output.read_text())
        self.assertTrue(report['passed']); self.assertEqual(report['input_format'], 'exodus')

    def test_native_default_node_ids(self):
        self.compare(self.native_input(explicit_ids=False))

    def test_native_requires_run_path(self):
        path = self.native_input()
        (self.mars / 'run.log').write_text('CONVERGED iterations=1277\n')
        with self.assertRaisesRegex(ValueError, 'native input path'):
            self.compare(path)

    def test_native_rejects_other_mesh(self):
        path = self.native_input()
        other = self.root / 'wrong.exo'; other.write_bytes(b'not the pinned mesh')
        with self.assertRaisesRegex(ValueError, 'pinned public reference'):
            self.compare(other)

    def test_permuted_ids_pass(self):
        self.compare()
        r = json.loads(self.output.read_text())
        self.assertTrue(r['passed']); self.assertEqual(r['pressure_max_scaled'], 0)

    def test_combined_layout_pass(self):
        self.write_exodus(combined=True); self.manifest(); self.compare()

    def test_pressure_offset_is_not_hidden(self):
        self.write_exodus(offset=.001); self.manifest()
        with self.assertRaisesRegex(ValueError, 'parity'): self.compare()
        self.assertAlmostEqual(json.loads(self.output.read_text())['pressure_mean_shift_pa'], -.001)

    def test_final_state_drift_rejected(self):
        self.write_exodus(drift=1e-5); self.manifest()
        with self.assertRaisesRegex(ValueError, 'reference change'): self.compare()

    def test_duplicate_ids(self):
        self.write_exodus(wrong_id=True); self.manifest()
        with self.assertRaisesRegex(ValueError, 'global IDs'): self.compare()

    def test_stale_output(self):
        self.write_exodus(stale=True); self.manifest()
        with self.assertRaisesRegex(ValueError, 'final two'): self.compare()

    def test_nonfinite(self):
        self.write_exodus(masked=True); self.manifest()
        with self.assertRaisesRegex(ValueError, 'nonfinite'): self.compare()

    def test_changed_result(self):
        self.write_exodus(offset=1)
        with self.assertRaisesRegex(ValueError, 'changed'): self.compare()

    def test_missing_mars_convergence(self):
        (self.mars / 'run.log').write_text('NOT CONVERGED\n')
        with self.assertRaisesRegex(ValueError, 'completion'): self.compare()

    def test_mars_residual_gate(self):
        path = self.mars / 'channel-metrics.csv'
        path.write_text(path.read_text().replace('9.9e-7', '9.9e-3'))
        with self.assertRaisesRegex(ValueError, 'MARS gate failed'): self.compare()

    def test_missing_exodus_values(self):
        with Dataset(str(self.reference / 'results.e'), 'r+') as ds:
            ds.variables['vals_nod_var1'][1, 0] = np.ma.masked
        self.manifest()
        with self.assertRaisesRegex(ValueError, 'missing/unwritten'): self.compare()

    def test_velocity_difference(self):
        self.fields[0, 1] += .001
        self.write_exodus(); self.manifest()
        with self.assertRaisesRegex(ValueError, 'parity'): self.compare()

    def test_mars_mesh_change(self):
        with (self.mars / 'channel.txt').open('a') as f: f.write('changed')
        with self.assertRaisesRegex(ValueError, 'mesh checksum'): self.compare()

    def test_coordinate_mismatch(self):
        self.xyz[0, 0] += .01; self.write_exodus(); self.manifest()
        with self.assertRaisesRegex(ValueError, 'coordinates'): self.compare()

    def test_preserve_existing_report(self):
        self.output.write_text('old report')
        with self.assertRaisesRegex(ValueError, 'exists'): self.compare()
        self.assertEqual(self.output.read_text(), 'old report')

    def test_runner_disables_capture_and_preserves_exit(self):
        capture = self.root / 'capture'; capture.mkdir()
        (capture / 'input.i').write_text(DECK)
        (capture / 'channel.exo').write_text('fake mesh; no CFD')
        exe = self.root / 'fake-solver'
        exe.write_text('#!/bin/sh\n'
                       'test -z "$MARS_OPENACCEL_EXPORT_DIR" || exit 91\n'
                       'test "$OMP_NUM_THREADS" = 1 || exit 92\n'
                       "printf 'Git hash: " + gate.PIN + "\\nIter = 1\\nIter = 2\\nConverged.\\n'\n"
                       'touch results.e\nexit 0\n')
        exe.chmod(0o700)
        dh, mh = gate.digest(capture / 'input.i'), gate.digest(capture / 'channel.exo')
        (capture / 'run.json').write_text(json.dumps(dict(fixture='public_channel', returncode=0,
            status='update_capture_completed', binary_sha256=gate.digest(exe),
            bundle=dict(sha256={'input.i':dh, 'channel.exo':mh}))))
        out = self.root / 'run'
        with patch.object(gate, 'DECK_SHA256', dh), patch.object(gate, 'MESH_SHA256', mh), \
             patch.dict(os.environ, {'MARS_OPENACCEL_EXPORT_DIR':'must-not-use'}, clear=True), \
             contextlib.redirect_stdout(io.StringIO()):
            gate.run(capture, exe, out)
            self.assertEqual(json.loads((out / 'comparison-run.json').read_text())['status'], 'native_convergence_reported')
            exe.write_text('#!/bin/sh\nexit 7\n'); exe.chmod(0o700)
            record = json.loads((capture / 'run.json').read_text()); record['binary_sha256'] = gate.digest(exe)
            (capture / 'run.json').write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError, 'failed'): gate.run(capture, exe, self.root / 'failed')
            self.assertEqual(json.loads((self.root / 'failed/comparison-run.json').read_text())['returncode'], 7)

    def test_iteration_limit_not_convergence(self):
        with self.assertRaisesRegex(ValueError, 'did not report'):
            gate.last_iteration('Iter = 1\nIter = 2\nSimulation is complete\n')


if __name__ == '__main__':
    unittest.main()
