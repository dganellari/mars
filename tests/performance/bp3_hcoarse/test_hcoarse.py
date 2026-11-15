"""Exercise the actual CLI, native workflow, receipt guards and submission wiring locally."""
import ast
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock

import integrate
from common import digest, save

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
EVIDENCE = REPO / '.review/bp3-h-integration'


class Integration(unittest.TestCase):
    def setUp(self):
        EVIDENCE.mkdir(parents=True, exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=str(EVIDENCE))
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def cli(self, macro=True):
        text = (REPO / 'examples/distributed/unstructured/mars_marsir_bp.cu').read_text()
        options = text[text.index('struct Options {'):text.index('\ntemplate <int P>\nint run(const Options& o)')]
        start = text.index('    Options o;', text.index('int main('))
        body = text[start:text.index('    const int status = o.p == 1 ? run<1>(o)', start)]
        source = self.root / ('cli_%s.cpp' % macro)
        source.write_text('''#include <string>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#define MARS_ENABLE_HYPRE 1
#define MPI_MAX_LIBRARY_VERSION_STRING 128
#define MARS_CELLWISE_MPI(x) x
#define BP_CK(x) x
#define MPI_IN_PLACE nullptr
#define MPI_INT 0
#define MPI_MIN 0
#define MPI_COMM_WORLD 0
struct cudaDeviceProp { int major=9, minor=0; };
[[maybe_unused]] void cudaGetDevice(int*d){*d=0;}
[[maybe_unused]] void cudaGetDeviceProperties(cudaDeviceProp*p,int){if(std::getenv("OLD_GPU"))p->major=8;}
[[maybe_unused]] void MPI_Allreduce(void*,int*s,int,int,int,int){if(std::getenv("MIXED_GPU"))*s=0;}
void MPI_Get_library_version(char*s,int*n){std::strcpy(s,"MPICH");*n=5;}
void MPI_Finalize(){}
void MPI_Dims_create(int ranks,int,int*d){
    d[0]=d[1]=d[2]=1;
    for(int f=2;f<=ranks;){if(ranks%f){++f;continue;}
        int k=d[0]<=d[1]?(d[0]<=d[2]?0:2):(d[1]<=d[2]?1:2);d[k]*=f;ranks/=f;}}
''' + ('#define MARS_BP3_COARSE_CONSOLIDATION 1\n' if macro else '') + options +
            '\nint main(int argc,char**argv){int rank=0,ranks=std::atoi(std::getenv("TEST_RANKS"));\n' + body + 'return 0;}\n')
        binary = source.with_suffix('')
        subprocess.run(['c++', '-std=c++20', '-Wall', '-Wextra', '-Werror', str(source), '-o', str(binary)], check=True)
        return binary

    def invoke(self, binary, ranks, p, ne, flags=(), env=None):
        base = dict(os.environ, MPICH_GPU_SUPPORT_ENABLED='1', MARS_HYPRE_GPU_AWARE='1',
                    MARS_HYPRE_SPMV_VENDOR='0', MARS_HYPRE_SPGEMM_VENDOR='0', TEST_RANKS=str(ranks))
        base.update(env or {})
        return subprocess.run([str(binary), '--ptx-dir', '/generated', '--p', str(p), '--ne', str(ne), '--mg', '--amg'] + list(flags),
                              env=base, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)

    def test_real_cli_selection_and_fallbacks(self):
        binary = self.cli()
        for ranks, p, ne, solves, depth in ((1,3,64,3,0), (4,7,64,5,0), (8,3,128,1,2),
            (8,7,64,1,0), (8,7,64,4,0), (8,7,64,5,1), (64,3,128,1,2), (64,3,256,1,2),
            (64,7,128,1,1), (64,7,64,5,0), (16,3,128,5,0), (64,3,512,5,0), (64,1,128,5,0)):
            result = self.invoke(binary,ranks,p,ne,['--solves',str(solves)])
            self.assertEqual(result.returncode,0,result.stdout)
            self.assertIn('requested=auto levels=%d' % depth,result.stdout)
        for flags, env in ((['--amg-consolidate','0'],{}), (['--amg-consolidate','10000'],{}),
            (['--amg-cycles','2'],{}), (['--deform','0'],{}), (['--stretch','2'],{}),
            (['--pre','1'],{}), (['--post','3'],{}), (['--power-steps','30'],{}),
            (['--cheb-interval','mfem'],{}), (['--spec','bp3_tri'],{}),
            ([],{'OLD_GPU':'1'}), ([],{'MIXED_GPU':'1'}), ([],{'MARS_HYPRE_SPMV_VENDOR':'1'}),
            ([],{'MARS_HYPRE_GPU_AWARE':'0'}), ([],{'MARS_HYPRE_SPGEMM_VENDOR':'1'})):
            result=self.invoke(binary,64,3,128,flags,env)
            self.assertEqual(result.returncode,0,result.stdout)
            self.assertIn('requested=auto levels=0',result.stdout)
        for depth in (0,1,2,3):
            result=self.invoke(binary,64,3,128,['--h-levels',str(depth)])
            self.assertEqual(result.returncode,0,result.stdout)
            self.assertIn('requested=explicit levels=%d' % depth,result.stdout)
        for flags, needle in ((['--h-levels','4'],'expects -1'), (['--h-levels','-2'],'expects -1'),
            (['--h-levels','1x'],'expects -1'), (['--h-levels',''],'expects -1'), (['--h-check'],'an option needs a value'),
            (['--h-levels','1','--p','1'],'require unblocked'), (['--h-levels','2','--mesh','unread.e'],'require unblocked'),
            (['--h-levels','3','--ne','32'],'need aligned'), (['--h-levels','1','--post','3'],'require unblocked')):
            result=self.invoke(binary,64,3,128,flags)
            self.assertEqual(result.returncode,1,result.stdout)
            self.assertIn(needle,result.stdout)
        result=self.invoke(self.cli(False),64,3,128)
        self.assertEqual(result.returncode,0,result.stdout)
        self.assertIn('requested=auto levels=0',result.stdout)

    @staticmethod
    def sample(ms=100):
        return dict(iterations=15,hierarchy=[2,100,1600,144],wall_ms=ms,wall_solves_ms=[110,ms,ms],setup_ms=300)

    def test_timing_rejections(self):
        rows={m:[self.sample() for _ in range(3)] for m in ('on','reference')}
        self.assertTrue(integrate.verdict(rows)['passed'])
        rows['on']=[self.sample(104) for _ in range(3)]
        self.assertFalse(integrate.verdict(rows)['passed'])
        rows['on']=[self.sample(96) for _ in range(3)]
        self.assertFalse(integrate.verdict(rows)['passed'])
        rows['on']=[self.sample() for _ in range(3)]
        rows['on'][0]['iterations']=14
        self.assertFalse(integrate.verdict(rows)['passed'])
        rows['on'].pop()
        self.assertFalse(integrate.verdict(rows)['passed'])

    def test_direct_checks_have_exact_extra_cycles(self):
        import copy
        normal=self.sample()
        checked=copy.deepcopy(normal)
        checked['h_check_worst']=2.65e-16
        checked['hierarchy'][-1]+=12
        self.assertTrue(integrate.checked_equivalent(normal,checked))
        self.assertFalse(integrate.equivalent([normal,checked]))
        for delta in (-12,-1,1,12):
            wrong=copy.deepcopy(checked);wrong['hierarchy'][-1]+=delta
            self.assertFalse(integrate.checked_equivalent(normal,wrong))
        for field in ('iterations','h_check_worst'):
            wrong=copy.deepcopy(checked);wrong[field]=None
            self.assertFalse(integrate.checked_equivalent(normal,wrong))
        wrong=copy.deepcopy(checked);wrong['hierarchy'][1]+=1
        self.assertFalse(integrate.checked_equivalent(normal,wrong))
        normal['hierarchy']=[-1,0,0,0];checked['hierarchy']=[-1,0,0,0]
        self.assertTrue(integrate.checked_equivalent(normal,checked))
        checked['hierarchy'][-1]=12
        self.assertFalse(integrate.checked_equivalent(normal,checked))
        self.assertEqual((HERE/'check.hpp').read_text().count('hierarchy.solve('),4)

    def test_trial_dispatch_and_retained_failure(self):
        for fail in (False,True):
            root=self.root/str(fail);root.mkdir();calls=[]
            def solve(*args,**kw):
                calls.append((args[3:9],kw))
                if fail and args[3]=='weak_p3_64_on_1':raise ValueError('injected numerical failure')
                return self.sample()
            with mock.patch.object(integrate,'solve',solve), mock.patch.object(integrate.standard_run,'verify_prepared'), \
                 mock.patch.object(integrate,'verify_reference'), mock.patch.object(integrate,'verify_bundle'):
                if fail:
                    with self.assertRaises(RuntimeError):integrate.trial(root,root,root)
                else:integrate.trial(root,root,root)
            self.assertEqual(len(calls),54)
            report=json.loads((root/'verdict.json').read_text())
            self.assertEqual(len(report['failures']),int(fail))
            self.assertEqual(report['passed'],not fail)
            self.assertEqual(sum(c[1]['solves']==5 for c in calls),6)

    def test_preparation_seal(self):
        bundle=self.root/'bundle';bundle.mkdir();save(bundle/'manifest.json',dict(files={}))
        root=self.root/'prepared';root.mkdir();calls=[]
        def command(args,log,**kw):
            args=list(map(str,args));calls.append(args)
            if 'build' in args:
                (root/'build').mkdir();(root/'build/binary').write_text('test binary receipt')
                save(root/'build/build.json',dict(binaries={'binary':digest(root/'build/binary')}));text='mock build'
            elif '--show-only=json-v1' in args:text=json.dumps(dict(tests=[dict(name='marsBP3CoarseSingle')]))
            elif args[0]=='ctest':
                self.assertNotIn('srun',args);text='100% tests passed, 0 tests failed out of 1'
            else:
                self.assertEqual(kw['expected'],1)
                text={'range':'expects -1','negative':'expects -1','p1':'require unblocked',
                      'rotated':'require unblocked','alignment':'need aligned','check':'an option needs a value'}[Path(log).stem[7:]]
            Path(log).write_text(text);return text
        def solve(*args,**kw):
            calls.append(args[3:]);folder=root/args[3];folder.mkdir();row=self.sample()
            if kw.get('checked'):
                row['hierarchy'][-1]+=12
                row['h_check_worst']=1e-16
            save(folder/'summary.json',row)
            return row
        with mock.patch.object(integrate,'command',command),mock.patch.object(integrate,'solve',solve), \
             mock.patch.object(integrate,'verify_bundle'),mock.patch.object(integrate,'verify_reference'):
            integrate.prepare(bundle,root)
        self.assertEqual(len(calls),45)
        integrate.standard_run.verify_prepared(root,bundle)
        (root/'explicit_r4_h3_check/summary.json').write_text('{}')
        with self.assertRaises(ValueError):integrate.standard_run.verify_prepared(root,bundle)

    def test_submission(self):
        # Run a spooled copy against recording commands, without any scheduler or remote access.
        bundle=self.root/'bundle';(bundle/'kit').mkdir(parents=True);(bundle/'hcoarse').mkdir()
        bindir=self.root/'bin';bindir.mkdir()
        capture=self.root/'calls'
        for name in ('sbatch','python'):
            f=bindir/name
            f.write_text('#!/bin/bash\nprintf "%s\\n" "$*" >> "$CAPTURE"\n' + ('echo "12345;daint"\n' if name=='sbatch' else ''))
            f.chmod(0o755)
        replacements={'/capstor/scratch/cscs/gandanie/bp3-hcoarse-resume-kit-20261106':str(bundle),
                      '/usr/bin/python3':str(bindir/'python')}
        for name in ('submit.sh','job.sh'):
            text=(HERE/name).read_text()
            for a,b in replacements.items():text=text.replace(a,b)
            (bundle/'hcoarse'/name).write_text(text)
            subprocess.run(['bash','-n',str(bundle/'hcoarse'/name)],check=True)
        env=dict(os.environ,PATH=str(bindir)+':'+os.environ['PATH'],CAPTURE=str(capture),SLURM_JOB_ID='111')
        subprocess.run(['bash',str(bundle/'hcoarse/submit.sh')],check=True,env=env,stdout=subprocess.PIPE)
        spooled=self.root/'slurm_script';shutil.copyfile(bundle/'hcoarse/job.sh',spooled)
        subprocess.run(['bash',str(spooled),'prepare'],check=True,env=env)
        subprocess.run(['bash',str(spooled),'resume','/previous/prepared'],check=True,env=env)
        subprocess.run(['bash',str(spooled),'trial','/sealed/prepared'],check=True,env=env)
        lines=capture.read_text().splitlines()
        self.assertEqual(len(lines),6)
        self.assertIn('--uenv-passthrough=ignore',lines[1])
        self.assertIn('--dependency=afterok:12345',lines[2])
        self.assertIn('--kill-on-invalid-dep=yes',lines[2])
        self.assertIn(str(bundle/'hcoarse/integrate.py')+' prepare',lines[3])
        self.assertIn('resume /capstor/scratch/cscs/gandanie/bp3-hcoarse-integrated/5279173',lines[1])
        self.assertIn(str(bundle/'hcoarse/integrate.py')+' resume',lines[4])
        self.assertIn('--prepared /previous/prepared',lines[4])
        self.assertIn('--prepared /sealed/prepared',lines[5])

    def test_registered_check_target(self):
        folder=self.root/'hcoarse';(folder/'private').mkdir(parents=True)
        shutil.copyfile(HERE/'extra.cmake',folder/'extra.cmake')
        private=folder/'private/mars_marsir_bp.cu'
        private.write_text('int main(){return 0;}\n')
        adapter=self.root/'adapter.cpp';adapter.write_text('int fixture(){return 0;}\n')
        (self.root/'CMakeLists.txt').write_text('\n'.join([
            'cmake_minimum_required(VERSION 3.24)','project(registration LANGUAGES CXX)',
            'set(CMAKE_EXPORT_COMPILE_COMMANDS ON)','add_library(mars INTERFACE)',
            'add_library(CUDA::cuda_driver INTERFACE IMPORTED)',
            'add_library(mars_bp3_coarse OBJECT "'+str(adapter)+'")',
            'include("'+str(folder/'extra.cmake')+'")',
            'set_source_files_properties("'+str(private)+'" PROPERTIES LANGUAGE CXX)'])+'\n')
        subprocess.run(['cmake','-S',str(self.root),'-B',str(self.root/'build')],check=True,
                       stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        db=json.loads((self.root/'build/compile_commands.json').read_text())
        self.assertEqual(len(db),2)
        for entry in db:
            for macro in ('MARS_BP3_COARSE_CONSOLIDATION','MARS_BP3_HCOARSE_CHECK'):
                self.assertEqual(macro in entry['command'],Path(entry['file'])==private)

    def test_receipts_and_five_solves(self):
        # Synthetic native-format input exercises the parser, not numerical correctness.
        out=self.root/'receipt';out.mkdir()
        log="""Hypre: SpMV hypre, SpGEMM hypre, GPU-aware MPI on
  FCG + V-cycle: 15 iterations, sqrt(<r,z>/<r0,z0>) = 7.398e-14, true assembled residual 1.544e-12
  max |x - u| / max |u| = 2.951e-14  PASS
  host wall solves 44.805 41.613 41.613 41.613 41.613 ms (max over ranks; last is reported)
  FCG iteration 2.774 ms (solve 41.6 ms, + setup 166 ms)
  coarse solve 1.667 ms per V-cycle
BP3 geometric levels: 2, global cells per axis 32 16 8; bottom rows 343; smoother bounds 1.9 1.9
PASS h validation: levels=2, adjoint/Galerkin/symmetry/linearity/continuity/repeat worst 1e-15
"""
        for rank in range(4):
            save(out/('rank%d.json'%rank),dict(version=1,setups=1,destroyed=1,cycles=240,empty_setups=0,
                skipped_setups=0,threshold_rows=5000,level=1,rows=66,nnz=2172))
        (out/'solver.log').write_text(log)
        result=integrate.parse_h(out,4,32,2,True,solves=5)
        self.assertEqual(len(result['wall_solves_ms']),5)
        with self.assertRaises(ValueError):integrate.parse_h(out,4,32,2,True)
        for a,b in (('bottom rows 343','bottom rows 344'),('worst 1e-15','worst nan'),
                    ('true assembled residual 1.544e-12','true assembled residual 1e-8'),
                    ('smoother bounds 1.9 1.9','smoother bounds 1.9')):
            (out/'solver.log').write_text(log.replace(a,b))
            with self.assertRaises(ValueError):integrate.parse_h(out,4,32,2,True,solves=5)

    def test_python_compatibility(self):
        for path in HERE.glob('*.py'):
            ast.parse(path.read_text(),filename=str(path),feature_version=(3,6))


if __name__=='__main__':unittest.main()
