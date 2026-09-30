"""Public probe generation and failure checks; these do not establish solver parity."""
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from netCDF4 import Dataset

import openaccel_simple_geometry_probe as probe
from prepare_simple_deck import load_deck, translate

ROOT = Path(__file__).resolve().parents[1]
DECK = (ROOT / 'tests/reference/openaccel/distributed_simple/preparation_fixture.i').read_text()


class ProbeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def test_only_outlet_and_iteration_controls_change(self):
        for case in probe.CASES:
            doc = load_deck(probe.probe_deck(DECK, case))
            args, mesh = translate(doc)
            options = dict(zip(args[::2], args[1::2]))
            self.assertEqual(mesh, 'channel.exo')
            self.assertEqual(options['--advection'], 'high-resolution')
            self.assertEqual(options['--velocity-interpolation'], 'linear-linear')
            self.assertEqual(float(options['--outlet-beta']), 1 if case.endswith('static') else .05)
            self.assertEqual(float(options['--inlet-velocity']), .1)
            self.assertEqual(float(options['--rho']), 1)
            self.assertEqual(float(options['--mu']), .1)
            c = doc['simulation']['solver']['solver_control']['basic_settings']['convergence_controls']
            self.assertEqual((c['min_iterations'], c['max_iterations']), (20, 20))

    def test_straight_and_warped_decks_are_identical(self):
        for outlet in ('average', 'static'):
            self.assertEqual(probe.probe_deck(DECK, 'straight-' + outlet),
                             probe.probe_deck(DECK, 'warped-' + outlet))

    def test_unknown_case_and_changed_controls_fail(self):
        with self.assertRaises(ValueError):
            probe.probe_deck(DECK, 'private')
        with self.assertRaises(ValueError):
            probe.probe_deck(DECK.replace('min_iterations: 2', 'min_iterations: 7'), probe.CASES[0])

    def test_warp_bends_inlet_and_preserves_input(self):
        xyz = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 1], [0, .5, .5]], dtype=float)
        old = xyz.copy()
        actual = probe.warp(xyz)
        np.testing.assert_array_equal(xyz, old)
        self.assertGreater(abs(np.linalg.det(actual[1:] - actual[0])), .01)

    def test_nonpublic_mesh_rejected_without_output(self):
        source = self.root / 'mesh.exo'; source.write_bytes(b'not a public mesh')
        out = self.root / 'out.exo'
        with self.assertRaises(ValueError):
            probe.prepare_mesh(source, out, probe.CASES[0])
        self.assertFalse(out.exists())

    def identity(self):
        (self.root / 'input.i').write_text(DECK)
        (self.root / 'channel.exo').write_bytes(b'synthetic test identity')
        exe = self.root / 'reference'; exe.write_text('synthetic executable identity'); exe.chmod(0o700)
        record = dict(fixture='public_channel', returncode=0, status='update_capture_completed',
                      binary_sha256=probe.digest(exe), bundle=dict(sha256={
                          'input.i': probe.digest(self.root / 'input.i'),
                          'channel.exo': probe.digest(self.root / 'channel.exo')}))
        (self.root / 'run.json').write_text(json.dumps(record))
        return exe, record

    def test_binary_and_input_identity_are_required(self):
        exe, record = self.identity()
        with patch.object(probe, 'DECK_SHA256', record['bundle']['sha256']['input.i']), \
             patch.object(probe, 'MESH_SHA256', record['bundle']['sha256']['channel.exo']):
            self.assertEqual(probe.capture_identity(self.root, exe), record['binary_sha256'])
            exe.write_text('another binary')
            with self.assertRaisesRegex(ValueError, 'executable changed'):
                probe.capture_identity(self.root, exe)
            (self.root / 'input.i').write_text('another deck')
            with self.assertRaisesRegex(ValueError, 'inputs changed'):
                probe.capture_identity(self.root, exe)

    def test_reference_failure_is_recorded_and_stops_suite(self):
        exe, _ = self.identity()
        output = self.root / 'probe'
        def fake_mesh(source, destination, case):
            destination.write_bytes(b'synthetic test mesh')
        with patch.object(probe, 'capture_identity', return_value='binary'), \
             patch.object(probe, 'prepare_mesh', side_effect=fake_mesh), \
             patch.dict(os.environ, dict(SLURM_NTASKS='1', OMPI_COMM_WORLD_SIZE='1', PMI_SIZE='1', PMIX_SIZE='1')), \
             patch.object(probe.subprocess, 'run') as launch:
            launch.return_value.returncode = 9
            with self.assertRaisesRegex(ValueError, 'OpenAccel failed'):
                probe.run(self.root, exe, output)
            record = json.loads((output / probe.CASES[0] / 'probe.json').read_text())
            self.assertEqual((record['returncode'], record['status']), (9, 'failed'))
            self.assertEqual(launch.call_count, 1)
            self.assertEqual((output / probe.CASES[0] / 'run.exit').read_text(), '9\n')

    def test_complete_suite_checks_and_hashes_outputs(self):
        exe, _ = self.identity()
        output = self.root / 'probe'
        def launch(command, cwd, env, stdout, stderr):
            self.assertNotIn('MARS_OPENACCEL_EXPORT_DIR', env)
            self.assertNotIn('MARS_OPENACCEL_PUBLIC_FIXTURE', env)
            stdout.write('Git hash: ' + probe.PIN + '\n' +
                         ''.join('Iter = {}\n'.format(i) for i in range(1, 21)))
            with Dataset(str(Path(cwd) / 'results.e'), 'w') as ds:
                for name, size in [('num_nodes', 425), ('time_step', 20), ('num_nod_var', 4), ('len_name', 33)]:
                    ds.createDimension(name, size)
                ds.createVariable('node_num_map', 'i4', ('num_nodes',))[:] = np.arange(1, 426)
                ds.createVariable('time_whole', 'f8', ('time_step',))[:] = np.arange(1, 21)
                for a in 'xyz':
                    ds.createVariable('coord' + a, 'f8', ('num_nodes',))[:] = 0
                names = np.zeros((4, 33), dtype='S1')
                for j, name in enumerate(('velocity_x', 'velocity_y', 'velocity_z', 'pressure')):
                    names[j, :len(name)] = np.frombuffer(name.encode(), dtype='S1')
                    ds.createVariable('vals_nod_var' + str(j + 1), 'f8', ('time_step', 'num_nodes'))[:] = 0
                ds.createVariable('name_nod_var', 'S1', ('num_nod_var', 'len_name'))[:] = names
            return probe.subprocess.CompletedProcess(command, 0)
        with patch.object(probe, 'capture_identity', return_value='binary'), \
             patch.object(probe, 'prepare_mesh', side_effect=lambda source, dest, case: dest.write_bytes(b'mesh')), \
             patch.dict(os.environ, dict(SLURM_NTASKS='1', OMPI_COMM_WORLD_SIZE='1', PMI_SIZE='1', PMIX_SIZE='1',
                                         MARS_OPENACCEL_EXPORT_DIR='stale', MARS_OPENACCEL_PUBLIC_FIXTURE='stale')), \
             patch.object(probe.subprocess, 'run', side_effect=launch):
            probe.run(self.root, exe, output)
        for case in probe.CASES:
            record = json.loads((output / case / 'probe.json').read_text())
            self.assertEqual(record['status'], 'snapshots_completed')
            for name, key in [('input.i', 'deck_sha256'), ('channel.exo', 'mesh_sha256'),
                              ('results.e', 'result_sha256'), ('run.log', 'log_sha256')]:
                self.assertEqual(probe.digest(output / case / name), record[key])

    def test_incomplete_iteration_sequence_fails(self):
        (self.root / 'run.log').write_text('Git hash: ' + probe.PIN + '\nIter = 1\nIter = 20\n')
        with self.assertRaisesRegex(ValueError, 'requested iterations'):
            probe.check_result(self.root)

    def test_missing_output_state_fails(self):
        (self.root / 'run.log').write_text('Git hash: ' + probe.PIN + '\n' +
                                         ''.join('Iter = {}\n'.format(i) for i in range(1, 21)))
        with Dataset(str(self.root / 'results.e'), 'w') as ds:
            ds.createDimension('time_step', 2)
            ds.createVariable('time_whole', 'f8', ('time_step',))[:] = [1, 20]
        with self.assertRaisesRegex(ValueError, 'missing or repeated output'):
            probe.check_result(self.root)


if __name__ == '__main__':
    unittest.main()
