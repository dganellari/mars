"""Local checks for the actual CLI, build guards and user-run job recipe."""
import ast
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import common
import export
import rank
import run
import reuse

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
EVIDENCE = REPO / '.review/bp3-integration/local'
HISTORY = REPO.parent / 'bp3-coarse-architect-20261031/finish-results'


class Integration(unittest.TestCase):
    def setUp(self):
        EVIDENCE.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=str(EVIDENCE))
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_cli_from_real_driver(self):
        text = (REPO / 'examples/distributed/unstructured/mars_marsir_bp.cu').read_text()
        options = text[text.index('struct Options {'):text.index('\ntemplate <int P>\nint run(const Options& o)')]
        start = text.index('    Options o;', text.index('int main(int argc, char** argv)'))
        body = text[start:text.index('    const int status = o.p == 1 ? run<1>(o)', start)]
        source = self.root / 'cli.cpp'
        source.write_text('#include <string>\n#include <cstdio>\n#include <cstring>\n#include <cstdlib>\n'
            '#define MARS_ENABLE_HYPRE 1\n#define MPI_MAX_LIBRARY_VERSION_STRING 128\n'
            'void MPI_Get_library_version(char*s,int*n){std::strcpy(s,"MPICH");*n=5;}\nvoid MPI_Finalize(){}\n' +
            (REPO / 'tests/bp3_coarse/integration/cli_stubs.hpp').read_text() + options + '\nint main(int argc,char**argv){int rank=0,ranks=64;\n' + body +
            'std::printf("%d %d %d\\n",o.consolidate_rows,o.p,o.ne);return 0;}\n')
        records = []
        for enabled in (False, True):
            binary = self.root / ('on' if enabled else 'off')
            args = ['c++', '-std=c++20', '-Wall', '-Wextra', '-Werror', str(source), '-o', str(binary)]
            if enabled:
                args.append('-DMARS_BP3_COARSE_CONSOLIDATION')
            subprocess.run(args, check=True)
            cases = [(p, ne, c, enabled or c == 0) for _, p, _, ne in run.timing_cases() for c in (0, -1)]
            cases += [(p, 8, c, enabled or c == 0) for p in (1, 3, 7) for c in (0, -1)]
            cases += [(3, 8, c, False) for c in (-3, 1, 9, 50001, '1x', '', '99999999999999999999999')]
            cases += [(3, 8, c, enabled) for c in (-2, "root", 10, 50000)]
            for p, ne, cutoff, expected in cases:
                args = [str(binary)] + common.solver_args(Path('/generated'), p, ne, cutoff)
                result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        env=dict(os.environ, MPICH_GPU_SUPPORT_ENABLED='1'), universal_newlines=True)
                self.assertEqual(result.returncode == 0, expected, result.stderr)
                records.append(dict(enabled=enabled, args=args[1:], passed=result.returncode == 0))
            base = common.solver_args(Path('/generated'), 3, 8, 0)
            index = base.index('--amg-consolidate')
            del base[index:index + 2]
            for extra, values, selected in (
                ([], {}, -1),
                (['--amg-consolidate', '0'], {}, 0),
                (['--amg-agg', '1'], {}, 0),
                (['--amg-theta', '0.5'], {}, 0),
                (['--amg-coarse-size', '500'], {}, 0),
                ([], {'MARS_HYPRE_GPU_AWARE': '0'}, 0),
                ([], {'MARS_HYPRE_SPMV_VENDOR': '1'}, 0),
                ([], {'MARS_HYPRE_SPGEMM_VENDOR': '1'}, 0)):
                env = dict(os.environ, MPICH_GPU_SUPPORT_ENABLED='1', MARS_HYPRE_GPU_AWARE='1',
                           MARS_HYPRE_SPMV_VENDOR='0', MARS_HYPRE_SPGEMM_VENDOR='0')
                env.update(values)
                result = subprocess.run([str(binary)] + base + extra, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE, env=env, universal_newlines=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(int(result.stdout.splitlines()[-1].split()[0]), selected if enabled else 0)
                records.append(dict(enabled=enabled, automatic=True, extra=extra, environment=values,
                                    selected=selected if enabled else 0))
        common.save(EVIDENCE / 'cli-checks.json', records)

    def test_config_guards(self):
        module = REPO / 'cmake/MarsOptions.cmake'
        required = ('MARSIR', 'FEM_EXAMPLES', 'UNSTRUCTURED', 'CUDA', 'HYPRE', 'MPI')
        for missing in (None,) + required + ('Linux',):
            lines = ['set(CMAKE_SYSTEM_NAME {})'.format('Darwin' if missing == 'Linux' else 'Linux')]
            lines += ['set(MARS_ENABLE_BP3_COARSE_CONSOLIDATION ON CACHE BOOL "")']
            lines += ['set(MARS_ENABLE_{} {} CACHE BOOL "")'.format(x, 'OFF' if missing == x else 'ON') for x in required]
            lines += ['set(MARS_ENABLE_HIP OFF CACHE BOOL "")', 'include("{}")'.format(module)]
            script = self.root / 'options.cmake'
            script.write_text('\n'.join(lines) + '\n')
            result = subprocess.run(['cmake', '-P', str(script)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(result.returncode == 0, missing is None, result.stderr.decode())
            script.write_text(script.read_text().replace('CONSOLIDATION ON CACHE BOOL', 'CONSOLIDATION AUTO CACHE STRING') + '\nif(NOT MARS_BP3_COARSE_REQUESTED STREQUAL \"' + ('TRUE' if missing is None else 'FALSE') + '\")\n message(FATAL_ERROR \"wrong automatic selection\")\nendif()\n')
            if missing == 'CUDA':
                script.write_text(script.read_text().replace('MARS_ENABLE_UNSTRUCTURED ON', 'MARS_ENABLE_UNSTRUCTURED OFF'))
            result = subprocess.run(['cmake', '-P', str(script)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(result.returncode, 0, result.stderr.decode())

    def test_callback_abi_guards(self):
        # These synthetic inputs exercise CMake rejection, not GPU ABI compatibility.
        module = (REPO / 'cmake/MarsBP3Coarse.cmake').read_text().split('add_library(mars_bp3_coarse OBJECT')[0]
        module = module.replace('${CMAKE_CURRENT_LIST_DIR}', str(REPO / 'cmake'))
        include = self.root / 'include'
        include.mkdir()
        library = self.root / 'libHYPRE.so'
        library.touch()
        reader = self.root / 'readelf'
        reader.write_text('#!/bin/sh\ncat "' + str(self.root / 'relocs') + '"\n')
        reader.chmod(0o755)
        header = '#define HYPRE_RELEASE_NUMBER 23300\n#define HYPRE_USING_CUDA\n'
        relocations = '0 R_AARCH64_JUMP_SLOT 0 hypre_seqAMGSetup + 0\n0 R_AARCH64_JUMP_SLOT 0 hypre_seqAMGCycle + 0\n'
        for mutation in ('valid', 'version', 'backend', 'static', 'unresolved', 'setup', 'cycle', 'single', 'complex', 'bigint', 'long_double'):
            config = header.replace('23300', '23400') if mutation == 'version' else header
            if mutation == 'backend':
                config = config.replace('#define HYPRE_USING_CUDA', '/* #undef HYPRE_USING_CUDA */')
            if mutation in ('single', 'complex', 'bigint', 'long_double'):
                config += '#define HYPRE_' + mutation.upper() + '\n'
            (include / 'HYPRE_config.h').write_text(config)
            text = relocations
            if mutation in ('setup', 'cycle'):
                text = '\n'.join(l for l in text.splitlines() if ('Setup' if mutation == 'setup' else 'Cycle') not in l)
            (self.root / 'relocs').write_text(text)
            path = str(library).replace('.so', '.a') if mutation == 'static' else 'HYPRE' if mutation == 'unresolved' else str(library)
            script = self.root / 'abi.cmake'
            script.write_text('set(_mars_bp3_mode ON)\nset(HYPRE_INCLUDE_DIRS "{}")\nset(HYPRE_LIBRARIES "{}")\n'
                              'set(_mars_bp3_readelf "{}")\n'.format(include, path, reader) + module)
            result = subprocess.run(['cmake', '-P', str(script)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(result.returncode == 0, mutation == 'valid', result.stderr.decode())
            script.write_text(script.read_text().replace('set(_mars_bp3_mode ON)', 'set(_mars_bp3_mode AUTO)'))
            result = subprocess.run(['cmake', '-P', str(script)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.assertEqual(result.returncode, 0, result.stderr.decode())
            self.assertIn(b'enabled' if mutation == 'valid' else b'disabled', result.stdout)

    def test_native_log_reader(self):
        if not (HISTORY / '5253440').is_dir():
            self.skipTest('Historical native result files are not available in this checkout')
        for job, name, ranks, cutoff in ((5253439, 'smoke_p3_r4_stock1', 4, 0),
                                         (5253440, 'weak_p3_r64_adaptive1', 64, -1)):
            old = HISTORY / str(job) / name
            out = self.root / name
            out.mkdir()
            for source in [old / 'solver.log'] + list(old.glob('rank*.json')):
                shutil.copy2(str(source), str(out / source.name))
            result = common.parse(out, ranks, cutoff)
            original = json.loads((old / 'summary.json').read_text())
            for key in ('iterations', 'true_residual', 'wall_ms'):
                self.assertEqual(result[key], original[key])
            path = out / 'rank0.json'
            audit = path.read_text()
            bad = json.loads(audit)
            bad['destroyed'] += 1
            path.write_text(json.dumps(bad))
            with self.assertRaises(ValueError):
                common.parse(out, ranks, cutoff)
            path.write_text(audit)
            log = (out / 'solver.log').read_text()
            (out / 'solver.log').write_text(log.replace('PASS', 'FAIL'))
            with self.assertRaises(ValueError):
                common.parse(out, ranks, cutoff)

    def test_complete_schedule_and_preparation_seal(self):
        bundle = self.root / 'bundle'
        bundle.mkdir()
        common.save(bundle / 'manifest.json', dict(files={}))
        prepared = self.root / 'prepared'
        prepared.mkdir()
        records = []

        def commands(args, log, env=None, expected=0, **kwargs):
            args = [str(a) for a in args]
            records.append(args)
            text = ''
            if str(run.HERE / 'build.py') in args:
                folder = prepared / 'build'
                folder.mkdir()
                outputs = {}
                for mode in ('off', 'on'):
                    target = folder / mode / 'examples/distributed/unstructured/mars_marsir_bp'
                    target.parent.mkdir(parents=True)
                    target.write_text('synthetic build receipt for workflow check only')
                    outputs[str(target.relative_to(folder))] = common.digest(target)
                common.save(folder / 'build.json', dict(binaries=outputs))
            elif '--show-only=json-v1' in args:
                text = json.dumps(dict(tests=[dict(name='marsBP3CoarseSingle')]))
            elif 'ctest' in args:
                text = '100% tests passed, 0 tests failed out of 1'
            elif args[-2] in ('balanced', 'empty-root', 'empty-tail'):
                ranks = next(int(a[2:]) for a in args if a.startswith('-n') and a[2:].isdigit())
                text = 'PASS {} ranks={}'.format(args[-2], ranks)
            elif expected == 1:
                text = 'MARS_ENABLE_BP3_COARSE_CONSOLIDATION=ON invalid --amg-consolidate expects -1'
            else:
                raise AssertionError('Unexpected launch: ' + repr(args))
            Path(log).write_text(text)
            return text

        def solve(bundle, preparation, root, name, mode, ranks, p, ne, cutoff):
            out = root / name
            out.mkdir()
            result = dict(iterations=15, hierarchy=[3, 200, 600, 144], wall_ms=100.0)
            common.save(out / 'summary.json', result)
            records.append(dict(solve=name, mode=mode, ranks=ranks, p=p, ne=ne, cutoff=cutoff))
            return result

        with mock.patch.object(run, 'command', commands), mock.patch.object(run, 'solve', solve), mock.patch.object(run, 'verify_accepted'):
            run.prepare(bundle, prepared)
            run.verify_prepared(prepared, bundle)
            trial = self.root / 'trial'
            trial.mkdir()
            run.trial(bundle, prepared, trial)
            self.assertTrue(json.loads((trial / 'verdict.json').read_text())['passed'])
            self.assertEqual(len([r for r in records if isinstance(r, dict)]), 78)
            target = prepared / 'smoke_p3_r4_on/summary.json'
            target.write_text('{}')
            with self.assertRaises(ValueError):
                run.verify_prepared(prepared, bundle)
        common.save(EVIDENCE / 'workflow-launches.json', records)

    def test_rank_export_and_failure(self):
        program = self.root / 'child.py'
        program.write_text('import os, pathlib, sys\npathlib.Path(os.environ["BP3_CONSOLIDATION_AUDIT"]).write_text("{}")\nsys.exit(int(sys.argv[1]))\n')
        original = tempfile.TemporaryDirectory
        def scratch(**kwargs):
            return original(dir=str(self.root))
        with mock.patch.dict(os.environ, SLURM_PROCID='0'), mock.patch.object(rank.tempfile, 'TemporaryDirectory', scratch):
            for code in (0, 7):
                out = self.root / str(code)
                out.mkdir()
                self.assertEqual(rank.run(out, True, [sys.executable, str(program), str(code)]), code)
                receipt = json.loads((out / 'process0.json').read_text())
                self.assertEqual(receipt['returncode'], code)
                self.assertEqual(receipt['audit_sha256'], common.digest(out / 'rank0.json'))

    def test_inherited_mpi_descriptor_and_signal(self):
        program = self.root / 'descriptor.py'
        program.write_text('import os,signal,sys\n'
            'assert os.read(int(sys.argv[1]),4)==b"PMI!"\n'
            'assert signal.SIGCHLD not in signal.pthread_sigmask(signal.SIG_BLOCK,[])\n')
        read_fd, write_fd = os.pipe()
        os.set_inheritable(read_fd, True)
        os.write(write_fd, b'PMI!')
        os.close(write_fd)
        previous = signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGCHLD})
        try:
            with mock.patch.dict(os.environ, SLURM_PROCID='0'):
                original = subprocess.run
                def old_descriptor_policy(*args, **kwargs):
                    kwargs['close_fds'] = True
                    kwargs['stderr'] = subprocess.PIPE
                    return original(*args, **kwargs)
                with mock.patch.object(rank.subprocess, 'run', old_descriptor_policy):
                    self.assertNotEqual(rank.run(self.root, False, [sys.executable, str(program), str(read_fd)], self.root), 0)
                self.assertEqual(rank.run(self.root, False, [sys.executable, str(program), str(read_fd)], self.root), 0)
        finally:
            signal.pthread_sigmask(signal.SIG_SETMASK, previous)
            os.close(read_fd)
        self.assertIn(int(signal.SIGCHLD), json.loads((self.root / 'process0.json').read_text())['parent_blocked_signals'])

    def test_ctest_uses_production_mpi_registration(self):
        # Exercise real CTest and its launcher ordering, environment and exit status.
        bundle = self.root / 'bundle'
        registration = bundle / reuse.REGISTRATION
        registration.parent.mkdir(parents=True)
        shutil.copyfile(str(REPO / 'tests/bp3_coarse/register.cmake'), str(registration))
        fixture = self.root / 'fixture'
        fixture.write_text('#!' + sys.executable + '\nimport os,sys\n'
            'assert sys.argv[1:]==["balanced","511"]\n'
            'assert os.environ["BP3_CONSOLIDATE_ROWS"]=="100"\n'
            'assert os.environ["BP3_FIXTURE_PROFILE"]=="mars"\n'
            'assert "SLURM_PROCID" not in os.environ\n'
            'print("TEST_FIXTURE_ENTERED",flush=True)\n'
            'sys.exit(int(os.environ.get("TEST_EXIT","0")))\n')
        fixture.chmod(0o755)
        launcher = self.root / 'mpi'
        launcher.write_text('#!' + sys.executable + '\nimport os,sys\n'
            'assert sys.argv[1:4]==["-n","1","--site-flag"]\n'
            'os.execv(sys.argv[4],sys.argv[4:])\n')
        launcher.chmod(0o755)
        directory = reuse.configure_test(bundle, self.root, fixture,
            ['-DMPIEXEC_EXECUTABLE=' + str(launcher), '-DMPIEXEC_NUMPROC_FLAG=-n', '-DMPIEXEC_PREFLAGS=--site-flag'])
        inventory = json.loads(common.command(['ctest', '--test-dir', directory, '--show-only=json-v1'], self.root / 'inventory.json'))
        properties = {p['name']: p['value'] for p in inventory['tests'][0]['properties']}
        self.assertEqual(properties['TIMEOUT'], 180)
        for code in (0, 7):
            text = common.command(['ctest', '--test-dir', directory, '--verbose'], self.root / ('ctest%d.log' % code),
                                  env=dict(os.environ, TEST_EXIT=str(code)), expected=0 if code == 0 else 8, timeout=15)
            self.assertIn('TEST_FIXTURE_ENTERED', text)

    def test_resume_copies_verified_native_binaries(self):
        evidence = REPO / '.review/bp3-integration'
        old = evidence / 'bundle-final'
        bundle = evidence / 'bundle-retry-20261105'
        previous = evidence / 'results-20261105-2/runs/5266939'
        if not bundle.is_dir() or not previous.is_dir():
            self.skipTest('Native build evidence is not available in this checkout')
        # Fake only the installed Hypre path; copy and verify the five real downloaded binaries.
        hypre = self.root / 'hypre/lib/libHYPRE.so'
        hypre.parent.mkdir(parents=True)
        hypre.write_bytes(b'local runtime guard input')
        header = self.root / 'hypre/include/HYPRE_config.h'
        header.parent.mkdir()
        header.write_bytes(b'local header guard input')
        root = self.root / 'resumed'
        root.mkdir()
        with mock.patch.object(reuse, 'OLD_BUNDLE', old), mock.patch.object(reuse, 'HYPRE', hypre), \
                mock.patch.object(reuse, 'HYPRE_SHA', common.digest(hypre)), \
                mock.patch.object(reuse, 'CONFIG_SHA', common.digest(header)):
            directory = reuse.resume(bundle, previous, root)
            self.assertTrue((directory / 'CTestTestfile.cmake').is_file())
            receipt = json.loads((root / 'reused-build.json').read_text())
            self.assertFalse(receipt['recompiled'])
            self.assertEqual(len(receipt['binaries']), 5)
            common.verify_files(root / 'build', receipt['binaries'])
            with mock.patch.object(reuse, 'BUILD_SHA', 'wrong'):
                with self.assertRaisesRegex(ValueError, 'independently checked'):
                    reuse.resume(bundle, previous, self.root / 'reject')
        # A source edit with a correspondingly updated manifest still invalidates reuse.
        altered = self.root / 'altered'
        shutil.copytree(str(bundle), str(altered))
        name = 'source/backend/distributed/unstructured/solvers/mars_cellwise_amg_consolidate.cu'
        (altered / name).write_text((altered / name).read_text() + '\n// changed\n')
        manifest = json.loads((altered / 'manifest.json').read_text())
        manifest['files'][name] = common.digest(altered / name)
        common.save(altered / 'manifest.json', manifest)
        with self.assertRaisesRegex(ValueError, 'identical numerical'):
            reuse.check_sources(old, altered)

    def test_submit_and_slurm_spool_paths(self):
        # Replace only site paths/interpreter; execute the real wrappers from a spool directory.
        bundle = self.root / 'bundle'
        (bundle / 'kit').mkdir(parents=True)
        fake = self.root / 'bin'
        fake.mkdir()
        calls = self.root / 'calls.jsonl'
        py = fake / 'python3'
        py.write_text('#!' + sys.executable + '\nimport json,os,sys\n'
                      'with open(os.environ["CALLS"],"a") as f:f.write(json.dumps(["python"]+sys.argv[1:])+"\\n")\n')
        batch = fake / 'sbatch'
        batch.write_text('#!' + sys.executable + '\nimport json,os,sys\n'
                         'with open(os.environ["CALLS"],"a") as f:f.write(json.dumps(["sbatch"]+sys.argv[1:])+"\\n")\n'
                         'print("1234;cluster")\n')
        py.chmod(0o755)
        batch.chmod(0o755)
        for name in ('job.sh', 'submit.sh'):
            text = (HERE / name).read_text().replace('/capstor/scratch/cscs/gandanie/bp3-coarse-integration-retry-20261105', str(bundle))
            text = text.replace('/usr/bin/python3', str(py))
            (bundle / 'kit' / name).write_text(text)
        spool = self.root / 'spool'
        spool.mkdir()
        shutil.copyfile(str(bundle / 'kit/job.sh'), str(spool / 'slurm_script'))
        env = dict(os.environ, PATH=str(fake) + ':' + os.environ['PATH'], CALLS=str(calls), SLURM_JOB_ID='999')
        for phase in ('prepare', 'resume', 'trial'):
            subprocess.run(['bash', str(spool / 'slurm_script'), phase, '/prepared'], env=env, cwd=str(spool), check=True)
        subprocess.run(['bash', str(bundle / 'kit/submit.sh')], env=env, cwd=str(spool), check=True, stdout=subprocess.PIPE)
        records = [json.loads(s) for s in calls.read_text().splitlines()]
        for record in records[:3]:
            self.assertEqual(record[1], str(bundle / 'kit/run.py'))
        sbatch = [r for r in records if r[0] == 'sbatch']
        self.assertEqual(len(sbatch), 2)
        self.assertIn('--dependency=afterok:1234', sbatch[1])
        for record in sbatch:
            self.assertIn('--uenv-passthrough=ignore', record)
            self.assertIn('--uenv=' + common.IMAGE + ':/user-environment', record)
            self.assertIn(str(bundle / 'kit/job.sh'), record)
        calls.write_text('')
        subprocess.run(['bash', str(bundle / 'kit/submit.sh'), '/previous'], env=env, cwd=str(spool), check=True, stdout=subprocess.PIPE)
        resumed = [json.loads(s) for s in calls.read_text().splitlines() if json.loads(s)[0] == 'sbatch']
        self.assertEqual(resumed[0][-2:], ['resume', '/previous'])
        self.assertIn('--time=00:45:00', resumed[0])
        common.save(EVIDENCE / 'wrapper-checks.json', records)

    def test_integrity_and_export_boundary(self):
        self.assertFalse(export.eligible('tests/data/private_mesh.e'))
        self.assertFalse(export.eligible('internal-notes/pump/report.md'))
        self.assertFalse(export.eligible('examples/distributed/unstructured/mars_ns_pump_solver.cu'))
        self.assertTrue(export.eligible('examples/distributed/unstructured/mars_marsir_bp.cu'))
        (self.root / 'value').write_text('original')
        common.save(self.root / 'manifest.json', dict(files={'value': common.digest(self.root / 'value')}))
        common.verify_bundle(self.root)
        (self.root / 'value').write_text('changed')
        with self.assertRaises(ValueError):
            common.verify_bundle(self.root)

    def test_python_and_shell_compatibility(self):
        for path in HERE.glob('*.py'):
            ast.parse(path.read_text(), filename=str(path), feature_version=(3, 6))
        for path in HERE.glob('*.sh'):
            subprocess.run(['bash', '-n', str(path)], check=True)
        self.assertIn('--time=01:55:00', common.step(1, 115))


if __name__ == '__main__':
    unittest.main(verbosity=2)
