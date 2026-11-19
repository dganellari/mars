from pathlib import Path
import sys,json,subprocess,importlib.util,os,tempfile
import numpy as np
from netCDF4 import Dataset
import argparse
args=argparse.ArgumentParser(description="CPU/MPI import parity gate (requires NumPy and netCDF4)")
args.add_argument("--build",type=Path,required=True)
args.add_argument("--mpiexec",default="mpiexec")
args=args.parse_args()
B=args.build.resolve();T=Path(tempfile.mkdtemp(prefix="import-check-",dir=B))
def fixture(path,packed=False,bad=None):
 ne=4;nodes=np.array([(x/ne,y/ne,z/ne) for z in range(ne+1) for y in range(ne+1) for x in range(ne+1)],dtype=float);nodes+=.05*np.prod(np.sin(np.pi*nodes),axis=1)[:,None]
 idx=lambda x,y,z:1+x+(ne+1)*(y+(ne+1)*z)
 cells=np.array([[idx(x+dx,y+dy,z+dz) for dx,dy,dz in [(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]] for z in range(ne) for y in range(ne) for x in range(ne)],dtype=int)
 if bad=='id':cells[0,0]=0
 if bad=='nan':nodes[0,0]=np.nan
 with Dataset(str(path),'w',format='NETCDF3_64BIT_OFFSET') as f:
  for name,size in [('num_dim',3),('num_nodes',len(nodes)),('num_elem',len(cells)),('num_el_blk',3),('num_el_in_blk1',20),('num_el_in_blk3',44),('num_nod_per_el1',8),('num_nod_per_el3',8)]:f.createDimension(name,size)
  f.createVariable('eb_status','i4',('num_el_blk',))[:]=[1,0,1]
  for b,rows in [(1,cells[:20]),(3,cells[20:])]:
   v=f.createVariable('connect%d'%b,'i4',('num_el_in_blk%d'%b,'num_nod_per_el%d'%b));v.elem_type='TETRA' if bad=='type' else 'HEX8';v[:]=rows
  if packed:f.createVariable('coord','f8',('num_dim','num_nodes'))[:]=nodes.T
  else:
   for d,k in enumerate(('coordx','coordy','coordz')):f.createVariable(k,'f8',('num_nodes',))[:]=nodes[:,d]
 return nodes,cells

def run(cmd,name,expect=0):
 with (T/(name+'.log')).open('w') as f:r=subprocess.run(list(map(str,cmd)),stdout=f,stderr=subprocess.STDOUT,timeout=120)
 assert r.returncode==expect,(name,r.returncode,(T/(name+'.log')).read_text()[-1500:])
 return (T/(name+'.log')).read_text()
for packed in (False,True):
 fixture(T/('packed.e' if packed else 'split.e'),packed)
for bad in ('id','nan','type'):
 fixture(T/(bad+'.e'),bad=bad)
 run([B/'mfem_partition',T/(bad+'.e'),T/('bad-'+bad),1],'reject-'+bad,1)
# MeshPartitioner preserves shared topology, including vertices shared by >2 ranks.
base=[B/'mfem_bp3','--device','cpu','--p','3','--ne','4','--tol','1e-13','--true-tol','1e-10','--solves','2','--reps','1','--maxit','2000']
report=[]
for n in (1,2,4,8):
 prefix=T/('part%d.'%n);run([B/'mfem_partition',T/'split.e',prefix,n],'partition%d'%n)
 imported=run([args.mpiexec,'--oversubscribe','-n',n]+base+['--mesh-prefix',prefix,'--deform','0','--audit-prefix',T/('import%d'%n)],'import%d'%n)
 cube=run([args.mpiexec,'--oversubscribe','-n',n]+base+['--deform','.05','--audit-prefix',T/('cube%d'%n)],'cube%d'%n)
 assert 'PASS' in imported and 'PASS' in cube
 def read(prefix):
  a=np.vstack([np.loadtxt(p) for p in T.glob(prefix+'-rank*.txt')]);keys=[tuple(np.round(v[:3],11)) for v in a];assert len(keys)==len(set(keys))==2197;return dict(zip(keys,a))
 a,b=read('import%d'%n),read('cube%d'%n);assert a.keys()==b.keys()
 errors=np.max(np.abs(np.array([a[k]-b[k] for k in a])),axis=0);scales=np.maximum(1,np.max(np.abs(np.array(list(a.values()))),axis=0));relative=errors/scales
 assert np.max(relative[:8])<1e-10 and relative[8]<1e-9,(n,relative)
 report.append(dict(ranks=n,relative_errors=relative.tolist()))
run([B/'mfem_partition',T/'packed.e',T/'packedpart.',1],'partition-packed')
assert (T/'packedpart.000000').read_bytes()==(T/'part1.000000').read_bytes()
run([B/'mfem_partition',T/'split.e',T/'part1.',1],'reject-existing-output',1)
# Missing part-count metadata must fail collectively instead of hanging.
run([args.mpiexec,'--oversubscribe','-n','2']+base+['--mesh-prefix',T/'part1.','--deform','0'],'wrong-parts',1)
(T/'cpu-audit.json').write_text(json.dumps(dict(status='PASS',cases=report,conversion_controls=5),indent=2)+'\n');print('Evidence:',T);print('PASS imported/cube discrete systems on 1/2/4/8 ranks, conversion and rejection gates')
