#!/usr/bin/env python3
"""Extract actual update arithmetic and its instrumentation from the pinned source."""
import argparse
from pathlib import Path
import subprocess
import sys
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'scripts'))
from openaccel_update_instrumentation import FLOW, BASE, PRESSURE, SEQUENCE, add_update_edits


def block(text,start):
    brace=text.index('{',start);end=brace+1;depth=1
    while depth:
        depth+=(text[end]=='{')-(text[end]=='}');end+=1
    return text[start:end]


def function(text,name,occurrence=0):
    start=-1
    for _ in range(occurrence+1):start=text.index(name,start+1)
    return block(text,start)


def capture_end(text,start,stage):
    call=text.index('mars_reference::update_record('+str(stage)+',',start)
    return text.index('\n                }',call)+len('\n                }')

def extract(source,output):
    pin='0d69041ba1afda63e9e4328d9e0d9834bba37756'
    edits={};add_update_edits(source,edits,ROOT/'tests/reference/openaccel')
    for name in (FLOW,BASE,PRESSURE,SEQUENCE):
        pinned=subprocess.check_output(['git','-C',str(source),'show',pin+':'+name])
        if (source/name).read_bytes()!=pinned:raise ValueError('update source differs from pin: '+name)
    text=edits[FLOW][1]
    common='''#pragma once
// Arithmetic extracted from pinned OpenAccel; original source attribution follows.
'''+text[:text.index('#include')]+'''
#include "evaluate.hpp"
#include "update_export.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#define SPATIAL_DIM 3
namespace update_reference {
using scalar=double;using label=int;
constexpr double SMALL=std::numeric_limits<double>::epsilon();
struct Bulk { unsigned long long identifier(int n) const { return n; } struct Bucket { bool owned() const {return true;} }; Bucket bucket(int) const {return {};}};
'''
    bodies=[]
    base=edits[BASE][1];start=base.index('                    const scalar export_old')
    end=capture_end(base,start,0)
    bodies.append(('''double fieldVal[]={x[0]}, effectiveCorrection[]={x[1]};
const double effectiveRelaxValue=x[2],offset=0,lowerBoundValue=0,upperBoundValue=1e99,clipFactor=1;
const int FIELD_DIM=1,BLOCKSIZE=1,STRIDE=0,OFFSET=0,CLIP=0,i=0,k=0,id=0;
''',base[start:end],'out.values[0]=fieldVal[0];'))
    part=edits[PRESSURE][1];start=part.index('                pCorrVal[i] = correction[');end=capture_end(part,start,1)
    bodies.append(('double pCorrVal[1]{}, correction[]={x[0]}; const int i=0,row=0,BLOCKSIZE=1,STRIDE=0; int bucket[]={entity};\n',part[start:end],'out.values[0]=pCorrVal[0];'))
    part=edits[SEQUENCE][1];start=part.index('                        double export_velocity');end=capture_end(part,start,2)
    bodies.append(('''double Ub[]={x[0],x[1],x[2]};const double* dub=x+3;const double* dpCorrdxb=x+6;
const bool consistent=x[9]!=0;const int iNode=0;struct Bucket {int entity;bool owned()const{return true;}int operator[](int)const{return entity;}} nodeBucket{entity};
struct Mesh {Bulk& b; Bulk& bulkDataRef(){return b;}} mesh{bulkData};auto meshRef=[&]() -> Mesh& {return mesh;};
''',part[start:end],'std::copy(Ub,Ub+3,out.values);'))
    for name,stage in (('Interior_',3),('BoundaryFieldInletSpecifiedVelocity_',4),('BoundaryFieldOutletSpecifiedPressure_',5)):
        part=function(text,'void flowModel::updateMassFlowRate'+name+'(',1)
        start=part.index('                scalar tmDot =');end=capture_end(part,start,stage)
        prefix='const int ip=0; double mDot[]={x['+('7' if stage==4 else '22')+']}; const double mDotURF=x['+('8' if stage==4 else '23')+'];\n'
        if stage==3:
            prefix+='''const double rhoHR=x[0]; const double *p_uIp=x+1,*p_duIp=x+4,*p_dpdxIp=x+7,*p_GpdxIp=x+10,*p_FOrigIp=x+13,*p_FIp=x+16,*p_scs_areav=x+19;
int elementBucket[]={entity};const int iElement=0;
'''
        elif stage==4:prefix+='const double rhoBip=x[0];const double *UbcVec=x+1,*areaVec=x+4;\n'
        else:
            prefix+='''const double rhoBip=x[0];const double *p_uBip=x+1,*p_duBip=x+4,*p_dpdxBip=x+7,*p_GpdxBip=x+10,*p_FOrigBip=x+13,*p_FBip=x+16,*areaVec=x+19;
const int rfflag[]={int(x[24])};
'''
        loop=part[start:end]
        if stage==5:
            # Preserve the early continue by wrapping the actual skip branch and final arithmetic in one iteration.
            skip=block(part,part.index('                if (rfflag[ip] == 1)'))
            loop='for(int once=0;once<1;++once) {\n'+skip+'\n'+loop+'\n}'
        bodies.append((prefix,loop,'out.values[0]=mDot[0];'))
    part=function(text,'void flowModel::updateFlowReversalFlag_(')
    part=part[part.index('            case boundaryPhysicalType::outlet:'):]
    start=part.index('                                        double export_old_flux')
    end=capture_end(part,start,6)
    prefix='''double mDot[]={x[0],x[1],x[2]};int revf_val[]={int(x[3]),int(x[4]),int(x[5])};
const double *p_U=x+6,*p_p=x+15,*pbc=x+18,*areaVec=x+21;
const bool ignoreFlagUpdate=x[30]!=0;const int numScsBip=3,nodesPerSide=3;const double f=1.0/3;
double p_uAvg[3]{},p_faceAreaVec[3]{};
'''
    bodies.append((prefix,part[start:end],'for(int j=0;j<3;++j){out.values[j]=mDot[j];out.values[3+j]=revf_val[j];}'))
    part=function(text,'void flowModel::updatePressureBoundarySideFieldAverageStaticPressure_(')
    start=part.index('                            for (label ip = 0; ip < numScsBip; ++ip)')
    loop=block(part,start)
    prefix='''const double *p_p=x,*p_face_shape_function=x+3,*areaVec=x+6;const int rfflag[]={int(x[9])};
const int numScsBip=1,nodesPerSide=3; double p_estimate=0,area=0;
'''
    bodies.append((prefix,loop,'out.values[0]=p_estimate;out.values[1]=area;'))
    start=part.index('                            for (label ip = 0; ip < numScsBip; ++ip)',start+len(loop))
    loop=block(part,start)
    prefix='''const double p_p[]={x[0]},pAvg=x[1],p_estimate=x[2],beta=x[3];double pbc[]={x[4]};
const int rfflag[]={int(x[5])},ipNodeMap[]={0};const int numScsBip=1;
'''
    bodies.append((prefix,loop,'out.values[0]=pbc[0];'))
    start=part.index('                    const double export_moment');end=capture_end(part,start,9)
    bodies.append(('double p_estimate=x[0],area=x[1];\n',part[start:end],'out.values[0]=p_estimate;'))
    for stage,(prefix,body,suffix) in enumerate(bodies):
        common+='inline void stage'+str(stage)+'(const update_replay::Input& in, update_replay::Output& out, int entity=1, int export_sample=0) {\n'
        common+='out=update_replay::Output{};const double* x=in.values;Bulk bulkData; const int side=entity;\n'+prefix+body.replace(', ip,\n', ', export_sample,\n')+'\n'+suffix+'\n}\n'
    common+='inline void evaluate(const update_replay::Input& in,update_replay::Output& out,int entity=1, int export_sample=0) {\n switch(in.stage) {\n'
    for stage in range(10):common+='case '+str(stage)+':stage'+str(stage)+'(in,out,entity,export_sample);break;\n'
    common+='}\n}\n}\n#undef SPATIAL_DIM\n';output.write_text(common)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();extract(a.source,a.output)
