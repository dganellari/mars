#!/usr/bin/env python3
"""Compare owned channel fields without requiring a scientific Python environment."""
import argparse
import csv
import json
import math
from pathlib import Path

FIELDS=('velocity_x','velocity_y','velocity_z','pressure')

def read(paths):
    rows={}
    for path in paths:
        with Path(path).open() as f:
            for row in csv.DictReader(f):
                node=int(row['node_id'])
                if node<=0 or node in rows: raise ValueError('invalid or duplicate owned node')
                value={k:float(row[k]) for k in ('x','y','z')+FIELDS}
                if not all(map(math.isfinite,value.values())): raise ValueError('non-finite field')
                rows[node]=value
    if not rows: raise ValueError('empty fields')
    return rows

def compare(reference,candidate,atol=1e-9,rtol=1e-7):
    if reference.keys()!=candidate.keys(): raise ValueError('owned identity sets differ')
    maximum=dict.fromkeys(FIELDS,0.)
    scaled=dict.fromkeys(FIELDS,0.)
    for node,a in reference.items():
        b=candidate[node]
        if any(a[k]!=b[k] for k in ('x','y','z')): raise ValueError('node coordinates differ')
        for k in FIELDS:
            error=abs(a[k]-b[k])
            maximum[k]=max(maximum[k],error)
            scaled[k]=max(scaled[k],error/(atol+rtol*max(abs(a[k]),abs(b[k]))))
    return {'passed':max(scaled.values())<=1,'max_absolute':maximum,'max_scaled':scaled,'atol':atol,'rtol':rtol}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',required=True,type=Path)
    p.add_argument('--results',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args()
    try:
        ref=read([a.reference]); results={}
        for ranks in (1,2,4):
            prefix=a.results/f'channel-{ranks}'
            paths=[Path(str(prefix)+f'.rank{rank}.csv') for rank in range(ranks)]
            log=Path(str(prefix)+'.log').read_text()
            if 'PASS: stationary non-conformal channel iteration' not in log:
                raise ValueError('channel did not finish successfully')
            results[str(ranks)]=compare(ref,read(paths))
        report={'schema':'mars-nonconformal-field-comparison-v1',
                'passed':all(r['passed'] for r in results.values()),'rank_comparisons':results}
        with a.output.open('x') as f: json.dump(report,f,indent=2); f.write('\n')
    except (ValueError,OSError,KeyError) as e: p.exit(1,f'Field comparison rejected: {e}\n')
    print(json.dumps(report,indent=2))
    raise SystemExit(0 if report['passed'] else 1)

if __name__=='__main__': main()
