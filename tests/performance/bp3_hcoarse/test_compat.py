"""Local tests for combined-source dispatch, operator ABI pairing and job submission."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest import mock
import compat
import export_compat
import export_hcoarse
from common import digest,save,verify_bundle

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[2]


class Compatibility(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=str(REPO/'.review/bp3-h-compat'))
        self.root=Path(self.temp.name)

    def tearDown(self):self.temp.cleanup()

    def test_all_dispatches_and_failure(self):
        for fail in (False,True):
            root=self.root/str(fail);root.mkdir();calls=[]
            def command(args,log,**kw):
                self.assertIn('build',args);self.assertIn('srun',args)
                (root/'build').mkdir();save(root/'build/build.json',dict(binaries={}))
            def solve(*args,**kw):
                name=args[3];calls.append((name,args[4:],kw))
                if fail and name=='check_p7_h3':raise ValueError('injected GPU validation failure')
                row=dict(iterations=15,hierarchy=[1,100,1000,144],true_residual=1e-12,relative=1e-14)
                if kw.get('checked'):row.update(h_check_worst=1e-16,hierarchy=[1,100,1000,156])
                return row
            with mock.patch.object(compat,'command',command),mock.patch.object(compat.i,'solve',solve), \
                 mock.patch.object(compat,'verify_bundle'),mock.patch.object(compat,'digest',return_value='test'):
                if fail:
                    with self.assertRaises(ValueError):compat.check(root,root)
                    self.assertFalse((root/'passed.json').exists())
                else:
                    compat.check(root,root)
                    report=json.loads((root/'passed.json').read_text())
                    self.assertEqual(report['solves'],24)
                    self.assertEqual(sum(bool(c[2].get('checked')) for c in calls),8)
                    self.assertEqual(sum(c[1][1]==8 for c in calls),4)

    def test_current_operator_export_and_legacy_rejection(self):
        with self.assertRaisesRegex(ValueError,'historical reference'):
            export_hcoarse.freeze(REPO,self.root,self.root,self.root/'reject')
        def freeze(repo,deps,accepted,target,**kw):
            self.assertTrue(kw['current_operators'])
            (target/'generated').mkdir(parents=True)
            save(target/'manifest.json',dict(files={}))
        target=self.root/'bundle'
        with mock.patch.object(export_compat,'freeze',freeze):
            export_compat.export_compat(REPO,self.root,self.root,target)
        record=verify_bundle(target)
        self.assertEqual(len(record['operator_ptx']),12)
        for name,sha in record['operator_ptx'].items():
            self.assertEqual(digest(REPO/'marsir-mlir/generated'/name),sha)
            self.assertEqual(digest(target/'generated'/name),sha)

    def test_spooled_compat_job(self):
        bundle=self.root/'bundle';bundle.mkdir();bindir=self.root/'bin';bindir.mkdir()
        calls=self.root/'calls'
        for name in ('sbatch','python'):
            p=bindir/name;p.write_text('#!/bin/bash\nprintf "%s\\n" "$*" >> "$CAPTURE"\n');p.chmod(0o755)
        for name in ('submit-compat.sh','compat-job.sh'):
            text=(HERE/name).read_text().replace('/capstor/scratch/cscs/gandanie/bp3-hcoarse-compat-kit-20261106',str(bundle)).replace('/usr/bin/python3',str(bindir/'python'))
            p=self.root/name;p.write_text(text);subprocess.run(['bash','-n',str(p)],check=True)
        env=dict(os.environ,PATH=str(bindir)+':'+os.environ['PATH'],CAPTURE=str(calls),SLURM_JOB_ID='123')
        subprocess.run(['bash',str(self.root/'submit-compat.sh')],env=env,check=True)
        spool=self.root/'slurm_script';shutil.copyfile(self.root/'compat-job.sh',spool)
        subprocess.run(['bash',str(spool)],env=env,check=True)
        lines=calls.read_text().splitlines();self.assertEqual(len(lines),3)
        self.assertIn('--uenv-passthrough=ignore',lines[1]);self.assertIn('-N2 -n8',lines[1])
        self.assertIn(str(bundle/'hcoarse/compat.py')+' check',lines[2])


if __name__=='__main__':unittest.main()
