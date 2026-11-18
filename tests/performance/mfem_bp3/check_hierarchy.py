"""CPU/MPI fine-system gate; the native GPU kit applies the same comparisons."""
from pathlib import Path
import argparse,json,os,subprocess,sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'bp3_phases'))
from checks import compare_fields,parse,phases


def run(executable,output,launcher):
    output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for p,ne,depths in [(3,8,(0,1,2)),(7,4,(0,1))]:
        reference=None
        for ranks in (1,4,8):
            for h in depths:
                case=output/('p%d_r%d_h%d'%(p,ranks,h));case.mkdir()
                cmd=launcher+['-n',str(ranks),str(executable),'--device','cpu','--amg-settings','mfem','--p',str(p),'--ne',str(ne),'--h-levels',str(h),'--tol','1e-13','--true-tol','1e-10','--solves','2','--reps','1','--audit-prefix',str(case/'fine')]
                # Exercise detailed nesting at every degree/rank/depth in this small gate.
                env=dict(os.environ,BP3_PHASE_PREFIX=str(case/'phase'),BP3_PHASE_DETAIL='1')
                (case/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
                with (case/'solver.log').open('w') as f:
                    subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,check=True,timeout=180)
                result=parse((case/'solver.log').read_text(),'mfem')
                result['phases']=phases(case,ranks)
                if reference is None:reference=case
                result['fine_system']=compare_fields(reference,case)
                result.update(p=p,ranks=ranks,h=h);rows.append(result)
                print('PASS',case.name,flush=True)
    (output/'result.json').write_text(json.dumps(dict(status='PASS',cases=rows),indent=2)+'\n')

if __name__=='__main__':
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('executable',type=Path);a.add_argument('output',type=Path)
    a.add_argument('--launcher',nargs='+',default=['mpiexec','--oversubscribe'])
    opt=a.parse_args();run(opt.executable.resolve(),opt.output.resolve(),opt.launcher)
