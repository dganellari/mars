"""Startup evidence checks use synthetic fields only."""
import contextlib
import copy
import csv
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import unittest
from unittest.mock import patch

import numpy as np
from netCDF4 import Dataset

import simple_startup_probe as probe
import test_simple_snapshot_compare as fixtures


class ReferenceLogTests(unittest.TestCase):
    def scan(self, text):
        state = probe.ReferenceLogState()
        for line in text.splitlines():
            state.feed(line)
        result = state.result()
        self.assertNotIn('PRIVATE', json.dumps(result))
        return result

    def test_uncaught_yaml_error_during_controls(self):
        result = self.scan("""Reading controls ..
terminate called after throwing an instance of 'YAML::TypedBadConversion<double>'
  what(): yaml-cpp: error at line 45, column 7: PRIVATE
srun: error: PRIVATE: task 0: Aborted (core dumped)
""")
        self.assertEqual(result['reference_stages_seen'], ['controls_read'])
        self.assertEqual(result['reference_error_categories'], ['yaml'])
        self.assertTrue(result['reference_cpp_termination_seen'])
        self.assertTrue(result['reference_exception_message_seen'])
        self.assertTrue(result['reference_abort_seen'])

    def test_multiline_io_error_hides_mesh_name(self):
        result = self.scan("""Finished reading controls ..
Reading mesh ..
terminate called after throwing an instance of 'std::runtime_error'
  what(): Ioss::DatabaseIO PRIVATE
Could not open PRIVATE/results.e.4.0: No such file or directory
""")
        self.assertEqual(result['reference_error_categories'], ['file_access', 'mesh_io'])
        self.assertEqual(result['reference_stages_seen'], ['controls_ready', 'mesh_read'])

    def test_known_failure_categories(self):
        messages = (
            ('std::bad_alloc', 'allocation'),
            ('MPI_Init_thread PRIVATE', 'mpi_initialization'),
            ('Provided MPI thread-level support is not sufficient', 'mpi_thread_support'),
            ('Kokkos::Cuda::initialize PRIVATE', 'kokkos'),
            ('unsupported decomposition method PRIVATE', 'decomposition'),
            ('Disk quota exceeded PRIVATE', 'disk_space'),
            ('invalid boundary part PRIVATE', 'input_validation'),
            ('Assertion PRIVATE failed', 'assertion'),
            ('Belos:: PRIVATE', 'linear_solver'),
        )
        for message, category in messages:
            with self.subTest(category=category):
                result = self.scan('terminate called\n  what(): ' + message)
                self.assertEqual(result['reference_error_categories'], [category])

    def test_normal_banners_and_private_values_are_not_errors(self):
        result = self.scan("""Command line: PRIVATE/Kokkos/yaml-cpp
Automatic domain decomposition: input Exodus file must be a serial file
Validating YAML input against Exodus file
Finished validating YAML input
[3] Initializing equation `PRIVATE` on realm `PRIVATE`
Iter = 12
""")
        self.assertEqual(result['reference_error_categories'], [])
        self.assertFalse(result['reference_cpp_termination_seen'])
        self.assertEqual(result['reference_stages_seen'],
                         ['mesh_validation', 'mesh_validated', 'equation_initialization', 'iteration'])

    def test_unknown_exception_is_reported_without_exporting_text(self):
        result = self.scan("terminate called after throwing an instance of 'PRIVATE'\nwhat(): PRIVATE 456.78")
        self.assertTrue(result['reference_cpp_termination_seen'])
        self.assertEqual(result['reference_error_categories'], [])
        self.assertNotIn('456.78', json.dumps(result))
        self.assertEqual(result['reference_exception_classes'], ['other'])

    def test_setup_exception_classes_and_source_signatures_are_allowlisted(self):
        result = self.scan("""Finished reading mesh ..
terminate called after throwing an instance of 'std::runtime_error'
what(): stk::mesh::impl::FieldRepository PRIVATE field restriction incompatible
  at /PRIVATE/FieldRepository.cpp:456
  at /PRIVATE/meshGeometry.cpp:123
  at /PRIVATE/PRIVATE.cpp:789
""")
        self.assertEqual(result['reference_exception_classes'], ['std::runtime_error'])
        self.assertEqual(result['reference_error_categories'], ['field_registration', 'stk'])
        self.assertEqual(result['reference_source_signatures'], ['FieldRepository.cpp', 'meshGeometry.cpp'])
        self.assertNotIn('456', json.dumps(result))

    def test_out_of_range_is_not_misreported_as_field_registration(self):
        result = self.scan("terminate called after throwing an instance of 'std::out_of_range'\nwhat(): map::at")
        self.assertEqual(result['reference_exception_classes'], ['std::out_of_range'])
        self.assertEqual(result['reference_error_categories'], ['container_lookup'])

    def test_master_element_failure_hides_topology_and_paths(self):
        result = self.scan("""terminate called after throwing an instance of 'std::logic_error'
what(): Expr 'theElem != nullptr' eval'd to false
location /PRIVATE/MasterElementFactory.C:179
PRIVATE topology 456
""")
        self.assertEqual(result['reference_error_categories'], ['master_element'])
        self.assertEqual(result['reference_source_signatures'], ['MasterElementFactory.C'])
        self.assertNotIn('456', json.dumps(result))

    def test_signature_like_private_names_do_not_escape(self):
        result = self.scan("""terminate called after throwing an instance of 'PRIVATE::runtime_error'
what(): PRIVATE_FieldRepository.cpp.tmp
""")
        self.assertEqual(result['reference_exception_classes'], ['other'])
        self.assertEqual(result['reference_source_signatures'], [])

    def test_libcpp_exception_class(self):
        result = self.scan('libc++abi: terminating due to uncaught exception of type std::length_error: PRIVATE')
        self.assertEqual(result['reference_exception_classes'], ['std::length_error'])

    def test_cpp_runtime_assertion_and_mpi_interleaving(self):
        result = self.scan("""[2] Reading mesh ..
[1] Finished reading controls ..
libc++abi: terminating due to uncaught exception of type PRIVATE
what(): Assertion PRIVATE failed
""")
        self.assertEqual(result['reference_stages_seen'], ['controls_ready', 'mesh_read'])
        self.assertEqual(result['reference_progress_scope'], 'any_logged_rank')
        self.assertEqual(result['reference_error_categories'], ['assertion'])


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

    def test_failed_process_without_results_preserves_actual_exit(self):
        pair = self.fixture.root / 'failed-launch'
        probe.prepare(self.fixture.case, self.fixture.reference, pair)
        exe = self.fixture.root / 'failed-solver.py'
        exe.write_text("import sys\nprint('private-path: error while loading shared libraries: private-lib')\nsys.exit(127)\n")
        exe.chmod(0o700)
        with patch.object(probe, 'runtime_libraries', return_value={'test': 'a'*64}):
            with self.assertRaisesRegex(probe.EvidenceError, '^launcher_exit$'):
                probe.launch(pair, 'openaccel', exe, 1, [sys.executable])
            record = probe.read_json(pair / 'reference/launch.json')
            self.assertEqual(record['status'], 'failed')
            self.assertEqual(record['exit_code'], 127)
            result = probe.inspect_launch(pair, 'openaccel', exe)
        self.assertEqual(result['process_exit_code'], 127)
        self.assertTrue(result['library_load_error_seen'])
        self.assertEqual(result['outputs_check'], 'reference_outputs')
        self.assertNotIn('private-path', json.dumps(result))
        self.assertNotIn('private-lib', json.dumps(result))

    def test_inspection_before_launch_is_read_only(self):
        pair = self.fixture.root / 'before-launch'
        probe.prepare(self.fixture.case, self.fixture.reference, pair)
        before = {str(p): probe.digest(p) for p in pair.rglob('*') if p.is_file()}
        exe = self.fixture.root / 'fake-executable'
        exe.chmod(0o700)
        with patch.object(probe, 'runtime_libraries', side_effect=probe.EvidenceError('runtime_libraries_unresolved')):
            result = probe.inspect_launch(pair, 'openaccel', exe)
        self.assertFalse(result['launch_start_present'])
        self.assertFalse(result['log_present'])
        self.assertIsNone(result['process_exit_code'])
        self.assertEqual(result['input_check'], 'passed')
        self.assertEqual(result['runtime_check'], 'runtime_libraries_unresolved')
        self.assertEqual(before, {str(p): probe.digest(p) for p in pair.rglob('*') if p.is_file()})

    def test_inspection_of_reference_abort_without_final_record(self):
        reference = self.pair / 'reference'
        (reference / 'launch.json').unlink()
        (reference / 'run.exit').write_text('134\n')
        (reference / 'run.log').write_text("Reading mesh ..\nterminate called after throwing an instance of 'std::runtime_error'\nwhat(): Ioss:: PRIVATE\n")
        for path in probe.result_files(reference):
            path.unlink()
        before = {str(p): probe.digest(p) for p in self.pair.rglob('*') if p.is_file()}
        with patch.object(probe, 'runtime_libraries', return_value={'synthetic': 'a'*64}), contextlib.redirect_stdout(io.StringIO()):
            code = probe.main(['inspect', '--pair', str(self.pair), '--solver', 'openaccel',
                               '--executable', str(self.fixture.root / 'fake-executable'), '--output', str(self.public)])
        result = probe.read_json(self.public)
        self.assertEqual(code, 0)
        self.assertEqual(result['process_exit_code'], 134)
        self.assertFalse(result['launch_record_present'])
        self.assertTrue(result['reference_cpp_termination_seen'])
        self.assertEqual(result['reference_error_categories'], ['mesh_io'])
        self.assertEqual(result['outputs_check'], 'reference_outputs')
        self.assertNotIn('PRIVATE', self.public.read_text())
        self.assertEqual(before, {str(p): probe.digest(p) for p in self.pair.rglob('*') if p.is_file()})

    def test_inspection_matches_public_source_without_sharing_error_text(self):
        reference = self.pair / 'reference'
        (reference / 'run.log').write_text("terminate called after throwing an instance of 'std::runtime_error'\nwhat(): Error in the expression provided at boundary PRIVATE\n")
        catalog = {tuple('error in the expression provided at boundary'.split()): {'src/model/model.cpp:62'}}
        with patch.object(probe, 'read_catalog', return_value=catalog), contextlib.redirect_stdout(io.StringIO()):
            code = probe.main(['inspect', '--pair', str(self.pair), '--solver', 'openaccel',
                               '--executable', str(self.fixture.root / 'fake-executable'),
                               '--reference-source', '/synthetic/source', '--output', str(self.public)])
        result = probe.read_json(self.public)
        self.assertEqual(code, 0)
        self.assertEqual(result['reference_catalog_check'], 'passed')
        self.assertEqual(result['reference_catalog_revision'], probe.REFERENCE_REVISION)
        self.assertEqual(result['reference_message_candidates'], ['src/model/model.cpp:62'])
        self.assertNotIn('PRIVATE', self.public.read_text())
        self.assertNotIn('/synthetic/source', self.public.read_text())

    def test_unavailable_source_does_not_export_process_error(self):
        with patch.object(probe, 'read_catalog', side_effect=OSError('PRIVATE')):
            result = probe.inspect_launch(self.pair, 'openaccel', self.fixture.root / 'fake-executable', Path('/synthetic/source'))
        self.assertEqual(result['reference_catalog_check'], 'rejected')
        self.assertNotIn('reference_message_candidates', result)
        self.assertNotIn('PRIVATE', json.dumps(result))

    def test_library_probe_labels_hide_paths(self):
        for output, code, label in (
                (b'private-lib => not found\n', 0, 'runtime_libraries_unresolved'),
                (b'private-executable: failure\n', 1, 'runtime_library_probe'),
                (b'linux-vdso.so.1 (0x0)\n', 0, 'runtime_library_paths'),
                (b'libfoo => /nonexistent/private-library (0x0)\n', 0, 'runtime_library_unreadable')):
            with self.subTest(label=label), patch.object(probe.subprocess, 'run', return_value=subprocess.CompletedProcess([], code, output)):
                with self.assertRaisesRegex(probe.EvidenceError, '^' + label + '$'):
                    probe.runtime_libraries(self.fixture.root / 'fake-executable')

    def test_library_probe_hashes_resolved_targets(self):
        library = self.fixture.root / 'library.so'
        library.write_bytes(b'synthetic library')
        output = ('libfoo => ' + str(library) + ' (0x0)\n').encode()
        with patch.object(probe.subprocess, 'run', return_value=subprocess.CompletedProcess([], 0, output)):
            self.assertEqual(probe.runtime_libraries(library), {str(library): probe.digest(library)})


if __name__ == '__main__':
    unittest.main()
