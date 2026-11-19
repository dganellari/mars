import copy,json,tempfile,unittest
from pathlib import Path
from checks import parse,compare_fields,phases,overhead

LOG='''CG + V-cycle: 17 iterations, sqrt(<r,z>/<r0,z0>) = 2e-14, true residual 3e-12
  max error 1e-14 PASS
  solves 12.0 10.0 ms (the last is reported)
  CG iteration 0.588 ms (solve 10.0 ms, + setup 20 ms)
'''
NAMES=('solve','coarse_entry_wait','coarse_complete','hypre_complete','operator_complete','outer_dss_complete','smoother_complete','cycle_complete')

def receipt(rank,total,coarse):
    regions={n:dict(inclusive_ns=0,exclusive_ns=0,calls=0) for n in NAMES}
    for name,ns in [('solve',total),('coarse_entry_wait',2),('coarse_complete',coarse),('hypre_complete',coarse-1)]:
        regions[name]=dict(inclusive_ns=ns,exclusive_ns=0,calls=1)
    regions['solve']['exclusive_ns']=total-coarse-2
    regions['coarse_complete']['exclusive_ns']=1
    regions['hypre_complete']['exclusive_ns']=coarse-1
    regions['coarse_entry_wait']['exclusive_ns']=2
    return dict(version=2,rank=rank,detail=False,regions=regions)

class Checks(unittest.TestCase):
    def test_parse_both_drivers(self):
        a=parse(LOG,'mfem');self.assertEqual(a['setup_plus_first_ms'],32.)
        b=LOG.replace('CG +','FCG +').replace('true residual','true assembled residual').replace('  solves','  host wall solves')
        self.assertEqual(parse(b,'mars'),a)
    def test_reject_failed_incomplete_or_nonfinite(self):
        for text in [LOG.replace('PASS','FAIL'),LOG.replace('3e-12','nan'),LOG.replace('3e-12','2e-10'),LOG.replace('17 iterations','2000 iterations'),LOG.replace('  solves','  absent'),LOG+'\nFAIL']:
            with self.subTest(text=text),self.assertRaises(ValueError):parse(text,'mfem')
    def test_same_rank_subtraction(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            for r in [receipt(0,100,60),receipt(1,120,40)]:
                (d/('phase-rank%d.json'%r['rank'])).write_text(json.dumps(r))
            out=phases(d,2)
            self.assertEqual(out['critical_rank'],1)
            self.assertEqual(out['remainder_ms'],80/1e6)
            self.assertNotEqual(out['remainder_ms'],60/1e6)
    def test_phase_schema_and_call_counts(self):
        with tempfile.TemporaryDirectory() as d:
            f=Path(d)/'phase-rank0.json';a=receipt(0,100,50)
            for change in ['rank','version','calls','negative','accounting']:
                b=copy.deepcopy(a)
                if change=='rank':b['rank']=1
                elif change=='version':b['version']=1
                elif change=='calls':b['regions']['hypre_complete']['calls']=2
                elif change=='negative':b['regions']['solve']['exclusive_ns']=-1
                else:b['regions']['hypre_complete']['inclusive_ns']=80
                f.write_text(json.dumps(b))
                with self.subTest(change=change),self.assertRaises(ValueError):phases(Path(d),1)
    def test_fine_system_and_refusal_controls(self):
        with tempfile.TemporaryDirectory() as d:
            a=Path(d)/'a';b=Path(d)/'b';a.mkdir();b.mkdir()
            row='0.25 0.5 0.75 0 1 2 3 4 1\n'
            (a/'fine-rank0.txt').write_text(row)
            (b/'fine-rank0.txt').write_text(row)
            self.assertEqual(compare_fields(a,b)['nodes'],1)
            for column in range(9):
                v=list(map(float,row.split()));v[column]+=.1
                (b/'fine-rank0.txt').write_text(' '.join(map(str,v))+'\n')
                with self.subTest(column=column),self.assertRaises(ValueError):compare_fields(a,b)
            (b/'fine-rank0.txt').write_text(row*2)
            with self.assertRaises(ValueError):compare_fields(a,b)
    def test_probe_gate_and_iteration_control(self):
        a=[dict(warm_ms=100.,iterations=10)]*3
        for time,it,expected in [(102.,10,True),(104.,10,False),(96.,10,False),(100.,11,False)]:
            b=[dict(warm_ms=time,iterations=it)]*3
            self.assertEqual(overhead(a,b)['accepted'],expected)
        with self.assertRaises(ValueError):overhead(a,[])

if __name__=='__main__':unittest.main()
