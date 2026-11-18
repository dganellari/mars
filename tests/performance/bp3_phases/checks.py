"""Numerical, provenance and timing gates; no failures become successful samples."""
from pathlib import Path
import json,math,re,statistics

def parse(text,kind):
    if re.search(r'\bFAIL\b',text) or not re.search(r'\bPASS\b',text):raise ValueError('Numerical failure')
    m=re.findall(r'true (?:assembled )?residual (\S+)',text)
    if len(m)!=1 or not(0<=float(m[0])<=1e-10):raise ValueError('True residual gate')
    it=re.search(r'(?:FCG|CG) \+ V-cycle: (\d+) iterations',text)
    if not it or not(0<int(it[1])<2000):raise ValueError('Iteration gate')
    line=re.findall(r'^\s+(?:run time|solves)\s+([\d. eE+\-]+)ms',text,re.M)
    if kind=='mars':
        line=re.findall(r'^\s+host wall solves\s+([\d. eE+\-]+)ms',text,re.M)
    if len(line)!=1:raise ValueError('Missing timing series')
    values=[float(x) for x in line[0].split()]
    if not values or not all(math.isfinite(x) and x>0 for x in values):raise ValueError('Invalid times')
    setup=re.findall(r'\+ setup ([\d.]+) ms\)',text)
    if len(setup)!=1 or not math.isfinite(float(setup[0])) or float(setup[0])<0:raise ValueError('Missing setup time')
    return dict(true_residual=float(m[0]),iterations=int(it[1]),times_ms=values,warm_ms=values[-1],reported_setup_ms=float(setup[0]),setup_plus_first_ms=float(setup[0])+values[0])

def fields(case):
    rows={}
    for file in sorted(Path(case).glob('fine-rank*.txt')):
        for line in file.read_text().splitlines():
            x=list(map(float,line.split()))
            if len(x)!=9 or not all(map(math.isfinite,x)):raise ValueError('Bad audit row')
            key=tuple(round(v,10) for v in x[:3])
            if key in rows:raise ValueError('Duplicate true-node ownership')
            rows[key]=x
    if not rows:raise ValueError('No fine-system audit')
    return rows

def compare_fields(reference,candidate):
    a,b=fields(reference),fields(candidate)
    if a.keys()!=b.keys():raise ValueError('Fine geometry changed')
    errors=[0.]*9;scales=[1.]*9
    for key in a:
        for j in range(9):
            errors[j]=max(errors[j],abs(a[key][j]-b[key][j]));scales[j]=max(scales[j],abs(a[key][j]))
    rel=[x/y for x,y in zip(errors,scales)]
    if max(rel[:8])>1e-10 or rel[8]>1e-9:raise ValueError('Fine operator/RHS/boundary/solution changed')
    return dict(nodes=len(a),relative_errors=rel)

def phases(case,ranks):
    files=sorted(Path(case).glob('phase-rank*.json'))
    if len(files)!=ranks:raise ValueError('Missing phase ranks')
    records=[json.loads(p.read_text()) for p in files]
    if sorted(x['rank'] for x in records)!=list(range(ranks)):raise ValueError('Wrong phase rank set')
    for x in records:
        if x.get('version')!=2 or not isinstance(x.get('detail'),bool):raise ValueError('Wrong phase schema')
        s=x['regions']
        if set(s)!={'solve','coarse_entry_wait','coarse_complete','hypre_complete','operator_complete','outer_dss_complete','smoother_complete','cycle_complete'}:raise ValueError('Missing phase regions')
        for v in s.values():
            if any(not isinstance(v[k],int) or v[k]<0 for k in ['inclusive_ns','exclusive_ns','calls']) or v['exclusive_ns']>v['inclusive_ns']:raise ValueError('Bad phase accounting')
        if x['detail'] and any(s[n]['calls']==0 for n in ['operator_complete','smoother_complete','cycle_complete']):raise ValueError('Missing detailed intervals')
        if sum(v['exclusive_ns'] for v in s.values())!=s['solve']['inclusive_ns']:raise ValueError('Exclusive intervals do not partition the solve')
        total=s['solve']['inclusive_ns'];coarse=s['coarse_complete']['inclusive_ns'];hypre=s['hypre_complete']['inclusive_ns'];entry=s['coarse_entry_wait']['inclusive_ns']
        if not(s['solve']['calls']==1 and s['coarse_complete']['calls']>0 and s['hypre_complete']['calls']==s['coarse_complete']['calls']==s['coarse_entry_wait']['calls'] and 0<=hypre<=coarse and coarse+entry<=total):raise ValueError('Invalid completed regions')
    # All differences below use the same rank, selected by its complete solve interval.
    x=max(records,key=lambda r:r['regions']['solve']['inclusive_ns']);s=x['regions']
    return dict(detail=x['detail'],regions=s,critical_rank=x['rank'],solve_ms=s['solve']['inclusive_ns']/1e6,coarse_ms=s['coarse_complete']['inclusive_ns']/1e6,hypre_ms=s['hypre_complete']['inclusive_ns']/1e6,entry_ms=s['coarse_entry_wait']['inclusive_ns']/1e6,remainder_ms=(s['solve']['inclusive_ns']-s['coarse_complete']['inclusive_ns'])/1e6)

def overhead(control,probe):
    if not control or len(control)!=len(probe):raise ValueError('Missing controls')
    same_iterations=len({x['iterations'] for x in control+probe})==1
    a=statistics.median(x['warm_ms'] for x in control);b=statistics.median(x['warm_ms'] for x in probe)
    ratio=b/a
    return dict(control_ms=a,probe_ms=b,ratio=ratio,same_iterations=same_iterations,accepted=same_iterations and abs(ratio-1)<=.03)
