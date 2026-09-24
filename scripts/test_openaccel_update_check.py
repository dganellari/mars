"""Transport and state-order failures using hand-worked, non-CFD fixtures."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from openaccel_update_check import CONTRACT, PHASES, load_updates, pack_updates


def fixture_rows(iteration):
    head=dict(kind='header',schema=1,producer='openaccel',fixture='public_channel',ranks=1,iteration=iteration,
              reference_revision=CONTRACT['reference']['revision'],solver_revision=CONTRACT['reference']['solver_gitlink']['revision'])
    fields={
        0:([2,1,.3],[2.3]), 1:([1],[1]), 2:([3,0,0,2,3,4,.1,.2,.3,0],[2.8,-.6,-1.2]),
        3:([1,3,0,0]+[0]*15+[1/3,0,0,1,.75],[1]),
        4:([1,-3,0,0,1/3,0,0,-1,.75],[-1]),
        5:([1,3,0,0]+[0]*15+[1/3,0,0,1,.75,0],[1]),
        6:([1]*3+[0]*3+[3,0,0]*3+[2]*6+[1/3,0,0]*3+[0],[1]*3+[0]*3),
        7:([2]*3+[1/3]*3+[1/3,0,0,0],[2/3,1/3]),
        8:([2,2,2,.05,2,0],[2]), 9:([2,1],[2]),
    }
    rows=[head]
    for phase in range(8):
        rows.append(dict(kind='phase',value=phase))
        for stage in (0,1,2,3,4,5,6,7,9,8):
            if PHASES[stage]!=phase:continue
            entity=1 if stage==9 else 50 if stage==3 else 20 if stage==4 else 30 if stage>=5 else 10
            samples=6 if stage==3 else 3 if stage in (4,5,7,8) else 1
            for sample in range(samples):
                x,y=fields[stage]
                rows.append(dict(kind='update',stage=stage,entity=entity,sample=sample,phase=phase,
                                 inputs=x[:],outputs=y[:]))
    rows.append(dict(kind='end',records=sum(r['kind']=='update' for r in rows)))
    return rows


def write_fixture(directory):
    for i in (1,2):
        (directory/('iteration'+str(i)+'.jsonl')).write_text(''.join(json.dumps(r)+'\n' for r in fixture_rows(i)))


class UpdateCaptureTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.data=self.root/'updates';self.data.mkdir();write_fixture(self.data)

    def change(self,action):
        path=self.data/'iteration2.jsonl';rows=[json.loads(l) for l in path.read_text().splitlines()]
        action(rows);path.write_text(''.join(json.dumps(r)+'\n' for r in rows))

    def test_pack_and_no_overwrite(self):
        pack_updates(self.data,self.root/'inputs.txt')
        self.assertEqual(len(load_updates(self.data)[0]),46)
        self.assertTrue((self.root/'inputs.txt').read_text().startswith('MARS_PUBLIC_UPDATE_REPLAY_V1 46\n'))
        with self.assertRaises(ValueError):pack_updates(self.data,self.root/'inputs.txt')

    def test_bad_records(self):
        cases=[(0,'inputs',2,0),(2,'inputs',9,2),(3,'inputs',0,0),(3,'inputs',22,2),
               (4,'inputs',7,2),(5,'outputs',0,2),(5,'inputs',24,1),(6,'inputs',3,1),
               (6,'inputs',18,3),(6,'inputs',21,float('nan')),(7,'outputs',0,0),
               (8,'inputs',2,3),(8,'inputs',3,2),(9,'inputs',1,0)]
        for stage,group,index,value in cases:
            with self.subTest(stage=stage,group=group,index=index):
                write_fixture(self.data)
                self.change(lambda rows:next(r for r in rows if r.get('stage')==stage)[group].__setitem__(index,value))
                with self.assertRaises(ValueError):load_updates(self.data)

    def test_phase_order_and_trace_order(self):
        self.change(lambda rows:next(r for r in rows if r.get('stage')==2).__setitem__('phase',3))
        with self.assertRaises(ValueError):load_updates(self.data)
        write_fixture(self.data)
        def swap(rows):
            a=next(i for i,r in enumerate(rows) if r.get('stage')==9)
            b=next(i for i,r in enumerate(rows) if r.get('stage')==8)
            rows[a],rows[b]=rows[b],rows[a]
        self.change(swap)
        with self.assertRaises(ValueError):load_updates(self.data)

    def test_duplicate_and_truncation(self):
        def duplicate(rows):
            rows.insert(4,copy.deepcopy(rows[3]));rows[-1]['records']+=1
        self.change(duplicate)
        with self.assertRaises(ValueError):load_updates(self.data)
        write_fixture(self.data);self.change(lambda rows:rows.pop())
        with self.assertRaises(ValueError):load_updates(self.data)

    def test_pin_fixture_and_rank(self):
        for key,value in [('reference_revision','wrong'),('fixture','private'),('ranks',4)]:
            write_fixture(self.data);self.change(lambda rows:rows[0].__setitem__(key,value))
            with self.assertRaises(ValueError):load_updates(self.data)

    def test_missing_sample(self):
        def remove(rows):
            rows.remove(next(r for r in rows if r.get('stage')==3 and r['sample']==5));rows[-1]['records']-=1
        self.change(remove)
        with self.assertRaises(ValueError):load_updates(self.data)

    def test_empty_directory(self):
        for p in self.data.iterdir():p.unlink()
        with self.assertRaises(ValueError):load_updates(self.data)

if __name__=='__main__':unittest.main()
