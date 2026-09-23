#!/usr/bin/env python3
"""Compile extracted native arithmetic/hooks and replay a synthetic ordered capture."""
import argparse
import os
from pathlib import Path
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'scripts'))
from openaccel_update_check import pack_updates
from extract_oracle import extract


def run(source,output):
    output.mkdir(parents=True,exist_ok=False)
    extract(source,output/'mars_update_reference.hpp')
    common=['-std=c++20','-Wall','-Wextra','-Werror','-Wno-unused-variable','-Wno-unused-parameter',
            '-I'+str(ROOT/'backend/distributed/unstructured/fem/segregated'),
            '-I'+str(ROOT/'tests/reference/openaccel/update_replay'),
            '-I'+str(ROOT/'tests/reference/openaccel'),'-I'+str(output)]
    environment=os.environ.copy()
    environment.pop('MARS_OPENACCEL_EXPORT_DIR',None)
    environment.pop('MARS_OPENACCEL_PUBLIC_FIXTURE',None)
    environment['OMPI_MCA_btl']='self'
    environment['OMPI_MCA_pml']='ob1'
    for name in ('native_check','capture_fixture'):
        subprocess.run(['mpicxx']+common+[str(Path(__file__).with_name(name+'.cpp')),'-o',str(output/name)],check=True)
        if name=='capture_fixture':
            environment.update(MARS_OPENACCEL_EXPORT_DIR=str(output),MARS_OPENACCEL_PUBLIC_FIXTURE='public_channel')
        subprocess.run([str(output/name)],env=environment,check=True)
    pack_updates(output/'updates',output/'inputs.txt')
    subprocess.run(['c++','-x','c++']+common+[str(ROOT/'examples/distributed/unstructured/mars_segregated_update_replay.cu'),
                   '-o',str(output/'replay')],check=True)
    subprocess.run([str(output/'replay'),str(output/'inputs.txt')],check=True)
    # Fail closed on corrupt transport before a user spends a GPU allocation.
    content=(output/'inputs.txt').read_text()
    for name,data in [('truncated',content[:-30]),('trailing',content+'unexpected\n')]:
        path=output/(name+'.txt');path.write_text(data)
        result=subprocess.run([str(output/'replay'),str(path)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if result.returncode==0:raise ValueError('accepted '+name+' transport')
    print('PASS: native expressions, capture hooks, ordered transport and host replay; STK/CUDA pending')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.source.resolve(),a.output.resolve())
