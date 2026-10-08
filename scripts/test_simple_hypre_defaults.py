"""Synthetic capture records exercise library identity and the public output boundary."""

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import simple_hypre_defaults as probe
import simple_startup_probe as startup


def raw_defaults(library, mpi):
    values = {key: 1 for key in probe.INT_KEYS}
    values.update({key: .25 for key in probe.REAL_KEYS})
    return dict(schema=probe.RAW_SCHEMA, version='2.33.0', build_cuda=False,
                library_path=str(library), mpi_library_path=str(mpi),
                memory_location=1, execution_policy=0, header_version_matches=True,
                getter_layout_checks_passed=True, defaults=values)


class DefaultsTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.pair = self.root / 'PRIVATE_PAIR'
        reference = self.pair / 'reference'
        reference.mkdir(parents=True)
        (self.pair / 'pair.json').write_text('{}')
        (reference / 'input.i').write_text('PRIVATE_DECK_WITHOUT_MESH_OR_FIELDS')
        (reference / 'run.log').write_text('PRIVATE_LOG')
        (reference / 'run.exit').write_text('0\n')
        self.exe = self.root / 'PRIVATE_EXECUTABLE'
        self.exe.write_text('synthetic executable')
        self.library = self.root / 'PRIVATE_INSTALL/lib/libHYPRE.so.2.33.0'
        self.library.parent.mkdir(parents=True)
        self.library.write_text('synthetic Hypre library')
        self.mpi = self.library.with_name('libmpi_gnu.so.12')
        self.mpi.write_text('synthetic MPI library')
        self.include = self.library.parent.parent / 'include'
        self.include.mkdir()
        for name in ('HYPRE.h', 'HYPRE_config.h', '_hypre_parcsr_ls.h'):
            (self.include / name).write_text('// synthetic installed header\n')
        self.libraries = {str(path): startup.digest(path) for path in (self.library, self.mpi)}
        self.record = dict(schema=startup.SCHEMA, solver='openaccel', status='started', ranks=4,
            command=['PRIVATE_LAUNCHER', str(self.exe), '-i', 'input.i'],
            executable=str(self.exe), executable_sha256=startup.digest(self.exe),
            libraries=self.libraries, environment={}, pair_sha256=startup.digest(self.pair / 'pair.json'))
        (reference / 'launch-start.json').write_text(json.dumps(self.record))
        self.record.update(status='finished', exit_code=0,
            files={name: startup.digest(reference / name) for name in ('input.i', 'run.log', 'run.exit')})
        (reference / 'launch.json').write_text(json.dumps(self.record))
        self.output = self.root / 'result'
        self.raw = raw_defaults(self.library, self.mpi)
        self.compilation_exit = 0
        self.launch_exit = 0
        self.after_launch = lambda: None
        self.calls = []

    def subprocess(self, command, stdout, **kwargs):
        self.calls.append(command)
        if '-o' in command:
            Path(command[-1]).write_text('synthetic probe executable')
            stdout.write(b'PRIVATE compiler output\n')
            code = self.compilation_exit
        else:
            stdout.write(b'PRIVATE binding output\n')
            stdout.write((json.dumps(self.raw, separators=(',', ':')) + '\n').encode())
            self.after_launch()
            code = self.launch_exit
        return type('ProcessResult', (), dict(returncode=code))()

    def invoke(self, linked=None):
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream), contextlib.redirect_stderr(stream), \
                patch.object(probe.shutil, 'which', side_effect=lambda name: '/PRIVATE_BIN/' + name), \
                patch.object(probe.subprocess, 'run', side_effect=self.subprocess), \
                patch.object(startup, 'runtime_libraries', return_value=linked or self.libraries):
            code = probe.main(['--pair', str(self.pair), '--output-dir', str(self.output),
                               '--launcher', 'srun', '--nodes=1', '-n', '1'])
        text = (self.output / 'public.json').read_text()
        self.assertNotIn('PRIVATE', text + stream.getvalue())
        self.assertNotIn(str(self.root), text + stream.getvalue())
        return code, json.loads(text)

    def test_queries_captured_install_and_exports_only_fresh_defaults(self):
        before = {p: p.read_bytes() for p in self.pair.rglob('*') if p.is_file()}
        code, result = self.invoke()
        self.assertEqual(code, 0)
        self.assertEqual(result['defaults'], self.raw['defaults'])
        self.assertTrue(result['loaded_hypre_library_matches_capture'])
        self.assertTrue(result['mpi_library_identity_verified'])
        self.assertFalse(result['application_effective_settings_verified'])
        self.assertFalse(result['matrix_dependent_amg_setup_verified'])
        self.assertEqual(len(self.calls), 2)
        self.assertIn(str(self.library), self.calls[0])
        self.assertIn('-I' + str(self.include), self.calls[0])
        self.assertEqual(self.calls[1], ['srun', '--nodes=1', '-n', '1', str(self.output / 'probe')])
        self.assertNotIn(str(self.exe), json.dumps(self.calls))
        self.assertEqual(before, {p: p.read_bytes() for p in self.pair.rglob('*') if p.is_file()})
        self.assertEqual(self.output.stat().st_mode & 0o777, 0o700)

    def test_missing_or_stale_capture_stops_before_compile(self):
        (self.pair / 'reference/run.log').write_text('PRIVATE_STALE')
        _, result = self.invoke()
        self.assertEqual(result['failed_check'], 'saved_capture')
        self.assertEqual(self.calls, [])

    def test_changed_recorded_library_stops_before_compile(self):
        self.library.write_text('PRIVATE_CHANGED')
        _, result = self.invoke()
        self.assertEqual(result['failed_check'], 'capture_libraries')
        self.assertEqual(self.calls, [])

    def test_changed_recorded_executable_stops_before_compile(self):
        self.exe.write_text('PRIVATE_CHANGED')
        _, result = self.invoke()
        self.assertEqual(result['failed_check'], 'capture_libraries')
        self.assertEqual(self.calls, [])

    def test_missing_installed_header_stops_before_compile(self):
        (self.include / '_hypre_parcsr_ls.h').unlink()
        _, result = self.invoke()
        self.assertEqual(result['failed_check'], 'installed_headers')
        self.assertEqual(self.calls, [])

    def test_compile_failure_does_not_launch(self):
        self.compilation_exit = 1
        code, result = self.invoke()
        self.assertEqual(code, 1)
        self.assertEqual(result['failed_check'], 'probe_compile')
        self.assertEqual(len(self.calls), 1)

    def test_mpi_link_mismatch_does_not_launch(self):
        _, result = self.invoke(dict(self.libraries, **{str(self.mpi): 'a'*64}))
        self.assertEqual(result['failed_check'], 'probe_libraries')
        self.assertEqual(len(self.calls), 1)

    def test_sequential_library_does_not_claim_mpi_verification(self):
        (self.include / 'HYPRE_config.h').write_text('#define HYPRE_SEQUENTIAL 1\n')
        self.raw['mpi_library_path'] = ''
        code, result = self.invoke({str(self.library): startup.digest(self.library)})
        self.assertEqual(code, 0)
        self.assertFalse(result['mpi_enabled'])
        self.assertFalse(result['mpi_library_identity_verified'])

    def test_wrong_runtime_hypre_is_rejected(self):
        self.raw['library_path'] = str(self.exe)
        _, result = self.invoke()
        self.assertEqual(result['failed_check'], 'loaded_library_identity')
        self.assertNotIn('defaults', result)

    def test_wrong_runtime_mpi_is_rejected(self):
        self.raw['mpi_library_path'] = str(self.exe)
        _, result = self.invoke()
        self.assertEqual(result['failed_check'], 'loaded_library_identity')

    def test_nonzero_launch_exit_rejects_even_complete_output(self):
        self.launch_exit = 137
        _, result = self.invoke()
        self.assertEqual(result['failed_check'], 'probe_launch')
        self.assertNotIn('defaults', result)

    def test_header_change_during_probe_is_rejected(self):
        self.after_launch = lambda: (self.include / 'HYPRE.h').write_text('PRIVATE_CHANGED')
        _, result = self.invoke()
        self.assertEqual(result['failed_check'], 'inputs_changed')
        self.assertNotIn('defaults', result)

    def test_library_change_during_probe_is_rejected(self):
        self.after_launch = lambda: self.mpi.write_text('PRIVATE_CHANGED')
        _, result = self.invoke()
        self.assertEqual(result['failed_check'], 'loaded_library_identity')

    def test_existing_output_directory_is_not_reused(self):
        self.output.mkdir()
        marker = self.output / 'keep'
        marker.write_text('original')
        with contextlib.redirect_stderr(io.StringIO()):
            code = probe.main(['--pair', str(self.pair), '--output-dir', str(self.output),
                               '--launcher', 'srun'])
        self.assertEqual(code, 1)
        self.assertEqual(marker.read_text(), 'original')

    def test_public_projection_rejects_malformed_values_and_omits_unknown_text(self):
        self.raw['PRIVATE_EXTRA'] = 'PRIVATE_SECRET'
        self.assertNotIn('PRIVATE', json.dumps(probe.public_defaults(self.raw)))
        for key, value in (('version', 'PRIVATE_VERSION'), ('header_version_matches', False),
                           ('getter_layout_checks_passed', False), ('build_cuda', 'PRIVATE_BOOL'),
                           ('execution_policy', True), ('memory_location', 99)):
            with self.subTest(key=key), self.assertRaises(ValueError):
                probe.public_defaults(dict(self.raw, **{key: value}))
        for value in (True, 'PRIVATE_VALUE', float('nan'), float('inf'), .5, 2**50):
            bad = dict(self.raw, defaults=dict(self.raw['defaults'], gmres_restart_dimension=value))
            with self.subTest(value=value), self.assertRaises(ValueError):
                probe.public_defaults(bad)
        del self.raw['defaults']['amg_relax_coarse']
        with self.assertRaises(ValueError):
            probe.public_defaults(self.raw)

    def test_ambiguous_or_missing_raw_record_fails(self):
        log = self.root / 'raw.log'
        line = json.dumps(self.raw, separators=(',', ':')) + '\n'
        for content in ('PRIVATE_NO_RECORD', line + line):
            log.write_text(content)
            with self.assertRaises(ValueError):
                probe.raw_output(log)


if __name__ == '__main__':
    unittest.main()
