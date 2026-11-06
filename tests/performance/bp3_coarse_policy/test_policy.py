"""Local gates for candidate scope, automatic selection and retained trial failures."""
import ast
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import candidate
import policy
from common import IMAGE, digest, save

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
EVIDENCE = REPO / '.review/bp3-auto'


class Policy(unittest.TestCase):
    def setUp(self):
        EVIDENCE.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=str(EVIDENCE))
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_candidate_preserves_adaptive_policy(self):
        out = self.root / 'private'
        candidate.generate(REPO, out)
        old = (REPO / 'backend/distributed/unstructured/solvers/mars_cellwise_amg_consolidate.cu').read_text()
        new = (out / 'solvers/mars_cellwise_amg_consolidate.cu').read_text()
        adaptive = 'return ranks>=50 ? 50000 : (ranks<5 ? 5000 : ranks*1000);'
        self.assertIn(adaptive, new)
        self.assertIn('hypre_ParCSRMatrixGlobalNumRows(a)<=limit', new)
        self.assertIn('BP3 root memory:', new)
        self.assertNotIn('after setup\n"', new)
        # Reversing the bounded edits must recover the production adapter exactly.
        before, after = new.split('#ifdef BP3_CONSOLIDATE_CUDA\n        if (rank==0) {\n            size_t free_bytes', 1)
        after = after.split('#endif\n', 1)[1]
        self.assertEqual((before + after).replace('500000', '50000'), old)
        broken = self.root / 'broken'
        adapter = broken / 'backend/distributed/unstructured/solvers/mars_cellwise_amg_consolidate.cu'
        adapter.parent.mkdir(parents=True)
        adapter.write_text(old.replace('requested<=50000', 'requested<=60000'))
        with self.assertRaisesRegex(ValueError, 'cutoff guards'):
            candidate.generate(broken, self.root / 'reject')

    def test_real_private_cli(self):
        private = self.root / 'private'
        candidate.generate(REPO, private)
        text = (private / 'driver/mars_marsir_bp.cu').read_text()
        options = text[text.index('struct Options {'):text.index('\n\ntemplate <int P>\nint run(const Options& o)')]
        start = text.index('    Options o;', text.index('int main(int argc, char** argv)'))
        body = text[start:text.index('    const int status = o.p == 1 ? run<1>(o)', start)]
        source = self.root / 'cli.cpp'
        source.write_text('#include <string>\n#include <cstdio>\n#include <cstring>\n#include <cstdlib>\n'
            '#define MARS_ENABLE_HYPRE 1\n#define MARS_BP3_COARSE_CONSOLIDATION 1\n'
            '#define MPI_MAX_LIBRARY_VERSION_STRING 128\n'
            'void MPI_Get_library_version(char*s,int*n){std::strcpy(s,"MPICH");*n=5;}\nvoid MPI_Finalize(){}\n' +
            options + '\nint main(int argc,char**argv){int rank=0,ranks=64;\n' + body +
            'std::printf("selected %d\\n",o.consolidate_rows);return 0;}\n')
        binary = self.root / 'cli'
        subprocess.run(['c++', '-std=c++20', '-Wall', '-Wextra', '-Werror', str(source), '-o', str(binary)], check=True)
        env = dict(os.environ, MPICH_GPU_SUPPORT_ENABLED='1', MARS_HYPRE_GPU_AWARE='1',
                   MARS_HYPRE_SPMV_VENDOR='0', MARS_HYPRE_SPGEMM_VENDOR='0')
        for cutoff, expected in ((500000, 0), (500001, 1), (-2, 1), (0, 0), (None, 0)):
            args = [str(binary), '--ptx-dir', '/generated', '--mg', '--amg']
            if cutoff is not None:
                args += ['--amg-consolidate', str(cutoff)]
            result = subprocess.run(args, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            self.assertEqual(result.returncode, expected, result.stderr)
            if not expected:
                self.assertIn('selected %d' % (-1 if cutoff is None else cutoff), result.stdout)

    def test_real_deferred_cmake_targets(self):
        folder = self.root / 'policy'
        folder.mkdir()
        shutil.copyfile(str(HERE / 'extra.cmake'), str(folder / 'extra.cmake'))
        files = [('solvers/mars_cellwise_amg_consolidate.cu', 'BP3_CONSOLIDATE_CUDA'),
                 ('driver/mars_marsir_bp.cu', 'MARS_BP3_COARSE_CONSOLIDATION'),
                 ('fixture.cu', 'BP3_CONSOLIDATE_CUDA')]
        for name, macro in files:
            path = folder / 'private' / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('#ifndef %s\n#error missing definition\n#endif\nint main(){return 0;}\n' % macro)
        source = self.root / 'source'
        source.mkdir()
        adapter = folder / 'private' / files[0][0]
        lines = ['cmake_minimum_required(VERSION 3.24)', 'project(registration LANGUAGES CXX)',
                 'set(CMAKE_EXPORT_COMPILE_COMMANDS ON)', 'add_library(mars INTERFACE)',
                 'add_library(mars_bp3_coarse INTERFACE)',
                 'target_include_directories(mars_bp3_coarse INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}")',
                 'set_property(TARGET mars_bp3_coarse PROPERTY INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}")',
                 'set(HYPRE_LIBRARY mars)']
        lines += ['add_library(%s INTERFACE IMPORTED)' % name for name in ('MPI::MPI_CXX', 'CUDA::cudart', 'CUDA::cuda_driver')]
        lines += ['include("%s")' % (folder / 'extra.cmake')]
        lines += ['set_source_files_properties("%s" PROPERTIES LANGUAGE CXX)' % (folder / 'private' / name) for name, _ in files]
        (source / 'CMakeLists.txt').write_text('\n'.join(lines) + '\n')
        subprocess.run(['cmake', '-S', str(source), '-B', str(self.root / 'build')], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        subprocess.run(['cmake', '--build', str(self.root / 'build'), '--target', 'bp3_policy_adapter'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        db = json.loads((self.root / 'build/compile_commands.json').read_text())
        self.assertEqual(len(db), 3)
        for entry in db:
            is_driver = Path(entry['file']).name == 'mars_marsir_bp.cu'
            self.assertEqual('MARS_BP3_COARSE_CONSOLIDATION' in entry['command'], is_driver)
        # This is CMake registration and host compilation, not CUDA validation.

    @staticmethod
    def sample(ms=100, earlier=False):
        return dict(iterations=15, hierarchy=[2, 267766, 26500002, 144] if earlier else [3, 28660, 2000000, 144],
                    wall_ms=ms, setup_ms=300.0, wall_solves_ms=[ms + 10, ms, ms])

    def test_performance_gates_reject_regressions_and_overlap(self):
        rows = dict(auto=[self.sample() for _ in range(3)], private=[self.sample() for _ in range(3)],
                    large=[self.sample(90, True) for _ in range(3)])
        self.assertTrue(policy.root_verdict(rows)['warm_candidate_passed'])
        rows['large'][0]['wall_ms'] = 101
        self.assertFalse(policy.root_verdict(rows)['warm_candidate_passed'])
        rows['large'] = [self.sample(96, True) for _ in range(3)]
        self.assertFalse(policy.root_verdict(rows)['warm_candidate_passed'])
        rows['large'] = [self.sample(90, True) for _ in range(3)]
        rows['private'][0]['iterations'] = 16
        self.assertFalse(policy.root_verdict(rows)['warm_candidate_passed'])
        rows = dict(accepted=[self.sample() for _ in range(3)], auto=[self.sample(104) for _ in range(3)])
        self.assertFalse(policy.auto_verdict(rows)['passed'])
        rows['auto'] = [self.sample(100) for _ in range(3)]
        self.assertTrue(policy.auto_verdict(rows)['passed'])
        self.assertFalse(policy.root_memory('BP3 root memory: 50 used of 100 bytes after setup')['peak_measured'])
        with self.assertRaises(ValueError):
            policy.root_memory('BP3 root memory: 110 used of 100 bytes after setup')

    def test_complete_trial_retains_failed_candidate(self):
        for fail in (False, True):
            out = self.root / str(fail)
            out.mkdir()
            calls = []
            def solve(bundle, prepared, root, name, mode, ranks, p, ne, cutoff, executable=None, automatic=False):
                calls.append(dict(name=name, mode=mode, ranks=ranks, cutoff=cutoff, automatic=automatic))
                self.assertEqual(automatic, '_auto' in name)
                if fail and name == 'weak_p3_r64_large1':
                    raise ValueError('injected numerical rejection')
                folder = root / name
                folder.mkdir()
                (folder / 'solver.log').write_text('BP3 root memory: 50 used of 100 bytes after setup\n')
                return self.sample(90 if cutoff == 500000 else 100, earlier=cutoff == 500000)
            with mock.patch.object(policy.integration, 'solve', solve), \
                    mock.patch.object(policy.integration, 'verify_prepared'), \
                    mock.patch.object(policy.integration, 'verify_accepted'), mock.patch.object(policy, 'verify_bundle'):
                if fail:
                    with self.assertRaises(RuntimeError):
                        policy.trial(self.root, self.root, out)
                else:
                    policy.trial(self.root, self.root, out)
            self.assertEqual(len(calls), 72)
            report = json.loads((out / 'verdict.json').read_text())
            self.assertTrue(report['auto_passed'])
            self.assertFalse(report['larger_root_promoted'])
            self.assertEqual(len(report['failures']), int(fail))
            if fail:
                self.assertTrue(report['cases'][4]['larger_root']['incomplete'])
            save(EVIDENCE / ('policy-dispatch-%s.json' % fail), calls)

    def test_preparation_seals_all_extra_checks(self):
        bundle, root = self.root / 'bundle', self.root / 'prepared'
        bundle.mkdir()
        root.mkdir()
        save(bundle / 'manifest.json', dict(files={}))
        launches = []
        def command(args, log, env=None, **kwargs):
            args = [str(a) for a in args]
            launches.append(args)
            if 'build' in args:
                folder = root / 'build'
                folder.mkdir()
                binary = folder / 'policy-binary'
                binary.write_text('workflow input only')
                save(folder / 'build.json', dict(binaries={'policy-binary': digest(binary)}))
                text = 'synthetic build receipt'
            else:
                ranks = next(int(a[2:]) for a in args if a.startswith('-n') and a[2:].isdigit())
                text = 'BP3 consolidation: level 1, 113235 rows\nPASS %s ranks=%d' % (args[-2], ranks)
            Path(log).write_text(text)
            return text
        def validate(bundle, root, build):
            save(root / 'prepared.json', dict(old_seal=True))
            for ranks in (1, 4):
                for p in (1, 3, 7):
                    folder = root / ('smoke_p%d_r%d_on' % (p, ranks))
                    folder.mkdir()
                    save(folder / 'summary.json', self.sample())
        def solve(bundle, prepared, root, name, mode, ranks, p, ne, cutoff, **kwargs):
            self.assertFalse((root / 'prepared.json').exists())
            launches.append([name, cutoff, kwargs])
            folder = root / name
            folder.mkdir()
            save(folder / 'summary.json', self.sample())
            return self.sample()
        with mock.patch.object(policy, 'command', command), mock.patch.object(policy.integration, 'validate', validate), \
                mock.patch.object(policy.integration, 'solve', solve):
            policy.prepare(bundle, root)
        policy.integration.verify_prepared(root, bundle)
        self.assertEqual(len(launches), 24)  # build, five fixtures, six defaults, twelve private controls
        (root / 'private_p7_r4_c500000/summary.json').write_text('{}')
        with self.assertRaisesRegex(ValueError, 'Changed or missing'):
            policy.integration.verify_prepared(root, bundle)

    def test_spooled_scripts_and_dependency(self):
        bundle, fake, spool = (self.root / name for name in ('bundle', 'bin', 'spool'))
        (bundle / 'policy').mkdir(parents=True)
        fake.mkdir()
        spool.mkdir()
        calls = self.root / 'calls.jsonl'
        for name in ('python3', 'sbatch'):
            path = fake / name
            path.write_text('#!' + sys.executable + '\nimport json,os,sys\n'
                'with open(os.environ["CALLS"],"a") as f:f.write(json.dumps(sys.argv)+"\\n")\n' +
                ('print("1234;cluster")\n' if name == 'sbatch' else ''))
            path.chmod(0o755)
        for name in ('job.sh', 'submit.sh'):
            text = (HERE / name).read_text().replace('/capstor/scratch/cscs/gandanie/bp3-coarse-policy-kit', str(bundle))
            text = text.replace('/usr/bin/python3', str(fake / 'python3'))
            (bundle / 'policy' / name).write_text(text)
        shutil.copyfile(str(bundle / 'policy/job.sh'), str(spool / 'slurm_script'))
        env = dict(os.environ, PATH=str(fake) + ':' + os.environ['PATH'], CALLS=str(calls), SLURM_JOB_ID='999')
        for phase in ('prepare', 'trial'):
            subprocess.run(['bash', str(spool / 'slurm_script'), phase, '/prepared'], env=env, cwd=str(spool), check=True)
        subprocess.run(['bash', str(bundle / 'policy/submit.sh')], env=env, cwd=str(spool), check=True, stdout=subprocess.PIPE)
        records = [json.loads(s) for s in calls.read_text().splitlines()]
        self.assertTrue(all(r[1] == str(bundle / 'policy/policy.py') for r in records[:2]))
        submissions = [r for r in records if Path(r[0]).name == 'sbatch']
        self.assertEqual(len(submissions), 2)
        self.assertIn('--dependency=afterok:1234', submissions[1])
        self.assertIn('--kill-on-invalid-dep=yes', submissions[1])
        for record in submissions:
            self.assertIn('--uenv-passthrough=ignore', record)
            self.assertIn('--uenv=' + IMAGE + ':/user-environment', record)
            self.assertIn(str(bundle / 'policy/job.sh'), record)
        save(EVIDENCE / 'policy-wrapper-checks.json', records)

    def test_site_syntax(self):
        for path in HERE.glob('*.py'):
            ast.parse(path.read_text(), filename=str(path), feature_version=(3, 6))
        for path in HERE.glob('*.sh'):
            subprocess.run(['bash', '-n', str(path)], check=True)


if __name__ == '__main__':
    unittest.main(verbosity=2)
