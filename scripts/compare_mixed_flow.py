#!/usr/bin/env python3
"""Compare original-node CSV fields; physical-boundary masking is explicit."""
import argparse
import csv
import glob
import json
import math
from pathlib import Path

FIELDS=('velocity_x','velocity_y','velocity_z','pressure')

def read(pattern):
    result={}
    files=sorted(glob.glob(pattern))
    if not files:raise ValueError('no field files: '+pattern)
    for name in files:
        with open(name,newline='') as source:
            for row in csv.DictReader(source):
                node=int(row['node_id'])
                if node in result:raise ValueError('duplicate owned node '+str(node))
                values={key:float(row[key]) for key in ('x','y','z')+FIELDS}
                if not all(math.isfinite(x) for x in values.values()):raise ValueError('nonfinite field')
                result[node]=values
    return result

def compare(reference,candidate,interior,velocity_atol,pressure_atol):
    if reference.keys()!=candidate.keys():raise ValueError('owned node IDs differ')
    errors={key:[] for key in FIELDS};included=0
    for node,a in reference.items():
        b=candidate[node]
        if max(abs(a[key]-b[key]) for key in ('x','y','z'))>1e-12:raise ValueError('coordinates differ')
        if interior and not (1e-9<a['x']<2-1e-9 and 1e-9<a['y']<1-1e-9 and 1e-9<a['z']<1-1e-9):continue
        included+=1
        for key in FIELDS:errors[key].append((abs(a[key]-b[key]),node,a[key],b[key]))
    if not included:raise ValueError('comparison has no nodes')
    report={'nodes':len(reference),'compared_nodes':included,'scope':'physical interior including both interface sides' if interior else 'all raw nodes','fields':{},'pass':True}
    for key,values in errors.items():
        worst=max(values);limit=pressure_atol if key=='pressure' else velocity_atol
        rms=math.sqrt(sum(x[0]*x[0] for x in values)/included)
        scale=math.sqrt(sum(x[2]*x[2] for x in values)/included)
        report['fields'][key]={'max_abs':worst[0],'rms_abs':rms,'relative_l2':rms/scale if scale else None,'worst_node':worst[1],'reference':worst[2],'candidate':worst[3],'absolute_tolerance':limit,'pass':worst[0]<=limit}
        report['pass']=report['pass'] and worst[0]<=limit
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('reference');p.add_argument('candidate');p.add_argument('--interior',action='store_true')
    p.add_argument('--velocity-atol',type=float,default=3e-7);p.add_argument('--pressure-atol',type=float,default=4.5e-5)
    p.add_argument('--output',type=Path)
    a=p.parse_args()
    if not all(math.isfinite(x) and x>0 for x in (a.velocity_atol,a.pressure_atol)):p.error('positive finite tolerances required')
    result=compare(read(a.reference),read(a.candidate),a.interior,a.velocity_atol,a.pressure_atol)
    text=json.dumps(result,indent=2)+'\n'
    if a.output:a.output.write_text(text)
    print(text,end='');return 0 if result['pass'] else 2

if __name__=='__main__':
    raise SystemExit(main())
