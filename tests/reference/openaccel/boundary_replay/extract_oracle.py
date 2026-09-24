#!/usr/bin/env python3
"""Compile-test the pinned reference's actual arithmetic loops without STK mesh plumbing."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'scripts'))
from openaccel_boundary_instrumentation import function_range, NAMES


def extract(source, output):
    contract=json.loads((ROOT/'tests/data/public_openaccel_reference/contract_v1.json').read_text())
    texts={}
    for family in ('navierStokes','pressureCorrection'):
        name='src/assemble/flow/segregatedFlow/{0}/{0}AssemblerElemBoundaryConditions.cpp'.format(family)
        text=(source/name).read_text()
        if hashlib.sha256((source/name).read_bytes()).hexdigest()!=contract['reference']['inspected_source_sha256'][name]:
            raise ValueError('boundary reference source differs from inspected pin')
        texts[family]=text
    mom=texts['navierStokes']
    macros=mom[mom.index('#define IP_EXPLICIT_ADVECTIVE_FLUX__'):mom.index('void navierStokesAssembler::assembleElemTermsBoundaryInletSpecifiedVelocity_')]
    common='''#pragma once
// Arithmetic extracted from the pinned OpenAccel source; original attribution follows.
'''+mom[:mom.index('#include')]+'''
#include "mars_segregated_boundary.hpp"
#include <vector>
#include <cmath>
#include <limits>
#define SPATIAL_DIM 3
namespace boundary_reference {
using label=int; using scalar=double;
constexpr double SMALL=std::numeric_limits<double>::epsilon();
'''+macros+'\n'
    for stage in range(6):
        family='pressureCorrection' if stage<3 else 'navierStokes'
        text=texts[family];a,b=function_range(text,family,NAMES[stage%3]);section=text[a:b]
        marker='// loop over face nodes' if stage==5 else '// loop over boundary ips'
        a=section.index(marker);a=section.index('            for (label ip',a)
        b=section.index('            this->applyCoeff_',a)
        loop=section[a:b]
        # Exclude debug-only SCL fields: mesh motion is outside this fixed-frame fixture.
        import re
        loop=re.sub(r'#ifndef NDEBUG.*?#endif[^\n]*\n','',loop,flags=re.S)
        prefix='''inline void stageSTAGE(const mars::segregated::BoundaryInput& x, mars::segregated::BoundaryOutput& out) {
    out = mars::segregated::BoundaryOutput{};
    const int nodesPerElement=4, nodesPerSide=3, numScsBip=3, faceOrdinal=0;
    const int* ipNodeMap=x.nearest;
    const int* faceIpNodeMap=x.nearest;
    const int* faceNodeOrdinals=x.face_nodes;
    int localFaceMap[3];
    for(int s=0;s<3;++s) for(int f=0;f<3;++f) if(x.face_nodes[f]==x.nearest[s]) localFaceMap[s]=f;
    struct Master { const int* opp; int opposingNodes(int,int ip) const {return opp[ip];} } master{x.opposing};
    const auto* meSCS=&master;
    const double *areaVec=x.area, *rfflag_dummy=nullptr;
    const int* rfflag=x.reversal;
    const double* p_velocity_face_shape_function=x.shape;
    const double* p_coordinate_face_shape_function=x.shape;
    const double* p_dndx=x.gradient;
    const double* p_bcMultiplier=x.bc_multiplier;
    const double* p_muEff=x.viscosity;
    const double* p_p=x.pressure;
    const double* p_rho=x.density;
    const double* UbcVec=x.boundary_velocity;
    const double* mDot=x.stored_flux;
    const double* uWallCoeffsBip=x.wall_coefficient;
    const double* p_Gpdx_elem=x.pressure_gradient;
    const double *p_du=x.influence_lhs, *p_duRhs=x.influence_rhs;
    double face_u[9], face_g[9];
    for(int f=0;f<3;++f) for(int j=0;j<3;++j) { face_u[3*f+j]=x.velocity[3*x.face_nodes[f]+j]; face_g[3*f+j]=x.pressure_gradient[3*x.face_nodes[f]+j]; }
    const double* p_U=STAGE==1 ? face_u : x.velocity;
    const double* p_Gpdx=face_g;
    double zero[12]{};
    const double *p_psi=zero, *p_Um=zero, *p_mat=zero, *p_ori=zero, *p_coordinates=zero, *p_F=zero, *p_F_elem=zero, *B_el=zero;
    const double comp=0, cvpgHarm=0;
    double nx[3]{}, coordBip[3]{}, uBip[3]{}, umBip[3]{}, GpdxBip[3]{}, dpdxBip[3]{}, duBip[3]{}, duRhsBip[3]{}, FBip[3]{}, FOrigBip[3]{};
    double *p_nx=nx, *p_coordBip=coordBip, *p_uBip=uBip, *p_umBip=umBip, *p_GpdxBip=GpdxBip, *p_dpdxBip=dpdxBip, *p_duBip=duBip, *p_duRhsBip=duRhsBip, *p_FBip=FBip, *p_FOrigBip=FOrigBip;
    double *p_lhs=out.lhs, *p_rhs=out.rhs;
'''.replace('STAGE',str(stage))
        # The reference has separate parent-local and face-local IP maps.
        if stage==1: prefix=prefix.replace('const int* faceIpNodeMap=x.nearest;', 'const int* faceIpNodeMap=nullptr;').replace('    struct Master', '    faceIpNodeMap=localFaceMap;\n    struct Master')
        common+=prefix+loop+'}\n'
    common+='}\n#undef SPATIAL_DIM\n'
    output.write_text(common)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();extract(a.source,a.output)
