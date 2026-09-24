#!/usr/bin/env python3
"""Compile exact capture expressions with mock mesh plumbing, then replay their output.

This executes the pinned arithmetic, not STK. Actual STK capture is a separate gate.
"""
import argparse
from pathlib import Path
import subprocess
import sys

from extract_oracle import ROOT, extract
sys.path.insert(0, str(ROOT/'scripts'))
from openaccel_boundary_instrumentation import function_range, instrument_function, NAMES
from openaccel_boundary_check import pack_boundary


def generate(source, output):
    text = '''#include "mars_boundary_reference.hpp"
#include "boundary_export.hpp"
#include <cstdlib>
struct Bulk {
    MPI_Comm parallel() const { return MPI_COMM_WORLD; }
    int parallel_rank() const { return 0; }
    int parallel_owner_rank(int) const { return 0; }
    unsigned long long identifier(int node) const { return node; }
};
struct Domain {
    bool isMaterialCompressible() const { return false; }
    bool frameRotating() const { return false; }
    bool meshMoving() const { return false; }
    Domain* zonePtr() { return this; }
};
'''
    stages = ('pressure.inlet','pressure.outlet','pressure.wall','momentum.inlet','momentum.outlet','momentum.wall')
    for stage, name in enumerate(stages):
        family = 'pressureCorrection' if stage < 3 else 'navierStokes'
        path = source/'src/assemble/flow/segregatedFlow'/family/(family+'AssemblerElemBoundaryConditions.cpp')
        original = path.read_text()
        a,b = function_range(original,family,NAMES[stage%3])
        instrumented = instrument_function(original[a:b],stage)
        a = instrumented.index('            if (export_boundary.active())')
        b = instrumented.index('            this->applyCoeff_(',a)
        capture = instrumented[a:b]
        text += '''void captureSTAGE(Bulk& bulkData, const mars::segregated::BoundaryInput& x,
                        const mars::segregated::BoundaryOutput& y) {
    mars_reference::BoundaryExport export_boundary(bulkData, "NAME", 3);
    Domain domain_object; auto* domain = &domain_object;
    const int nodesPerSide=3, numScsBip=3, nodesPerElement=4, faceOrdinal=0, side=21+STAGE%3;
    FACE_NODE_DECLARATION
    const int *ipNodeMap=x.nearest, *faceIpNodeMap=x.nearest, *rfflag=x.reversal;
    struct Master {
        const int *opp, *face;
        int opposingNodes(int,int ip) const { return opp[ip]; }
        const int* side_node_ordinals(int) const { return face; }
    } master{x.opposing,x.face_nodes};
    auto* meSCS=&master;
    const double *areaVec=x.area, *UbcVec=x.boundary_velocity, *mDot=x.stored_flux, *uWallCoeffsBip=x.wall_coefficient;
    std::vector<int> connectedNodes = STAGE==5 ? std::vector<int>{20,30,40} : std::vector<int>{10,20,30,40};
    const int sideNodeRels[]={20,30,40};
    const double cvpgHarm=0;
    const int width = STAGE==5 ? 9 : (STAGE<3 ? 4 : 12);
    std::vector<double> lhs(y.lhs,y.lhs+width*width), rhs(y.rhs,y.rhs+width);
    std::vector<double> ws_velocity_face_shape_function(x.shape,x.shape+9), ws_rho(x.density,x.density+3);
    std::vector<double> ws_U(x.velocity,x.velocity+(STAGE==5?9:12));
    if(STAGE==1) { ws_U.resize(9); for(int f=0;f<3;++f)for(int j=0;j<3;++j)ws_U[3*f+j]=x.velocity[3*x.face_nodes[f]+j]; }
    std::vector<double> ws_p(x.pressure,x.pressure+4), ws_dndx(x.gradient,x.gradient+36);
    std::vector<double> ws_Gpdx_elem(x.pressure_gradient,x.pressure_gradient+12);
    std::vector<double> ws_Gpdx(9);
    for(int f=0;f<3;++f)for(int j=0;j<3;++j)ws_Gpdx[3*f+j]=x.pressure_gradient[3*x.face_nodes[f]+j];
    std::vector<double> ws_du(x.influence_lhs,x.influence_lhs+9), ws_duRhs(x.influence_rhs,x.influence_rhs+9);
    std::vector<double> ws_bcMultiplier(x.bc_multiplier,x.bc_multiplier+4), ws_muEff(x.viscosity,x.viscosity+3);
    std::vector<double> ws_F(9,0), ws_FOrig_elem(12,0), ws_F_elem(12,0);
'''.replace('STAGE',str(stage)).replace('NAME',name).replace(
            'FACE_NODE_DECLARATION', '' if stage == 4 else 'const int* faceNodeOrdinals=x.face_nodes;')+capture+'}\n'
    text += '''int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    if(argc!=2) { MPI_Finalize(); return 1; }
    setenv("MARS_OPENACCEL_EXPORT_DIR",argv[1],1);
    setenv("MARS_OPENACCEL_PUBLIC_FIXTURE","public_channel",1);
    Bulk bulk;
    using namespace mars::segregated;
    using Fn=void(*)(const BoundaryInput&,BoundaryOutput&);
    const Fn reference[]={boundary_reference::stage0,boundary_reference::stage1,boundary_reference::stage2,
                          boundary_reference::stage3,boundary_reference::stage4,boundary_reference::stage5};
    using Capture=void(*)(Bulk&,const BoundaryInput&,const BoundaryOutput&);
    const Capture capture[]={capture0,capture1,capture2,capture3,capture4,capture5};
    for(int call=0;call<2;++call)for(int stage=0;stage<6;++stage) {
        BoundaryInput x; x.stage=stage;
        const double grad[]={-1,-1,-1,1,0,0,0,1,0,0,0,1};
        for(int s=0;s<3;++s) {
            x.face_nodes[s]=x.nearest[s]=stage==5?s:s+1;
            x.density[s]=1+.2*s; x.viscosity[s]=.5+s; x.wall_coefficient[s]=1+s;
            x.stored_flux[s]=.1*(s+1);x.reversal[s]=call && s==1;
            for(int j=0;j<3;++j) {
                x.area[3*s+j]=1.0/6; x.shape[3*s+j]=s==j?11.0/18:7.0/36;
                x.influence_lhs[3*s+j]=.5+.1*j+.2*s; x.influence_rhs[3*s+j]=2+x.influence_lhs[3*s+j];
                x.boundary_velocity[3*s+j]=stage==2?0:.1*(3*s+j+1);
            }
            for(int j=0;j<12;++j)x.gradient[12*s+j]=grad[j];
        }
        for(int n=0;n<4;++n) {
            x.pressure[n]=.2*n*n;
            for(int j=0;j<3;++j) {x.velocity[3*n+j]=.1*(n+1)*(j-1); x.pressure_gradient[3*n+j]=.5*n-.2*j;}
        }
        x.bc_multiplier[0]=1;
        BoundaryOutput y; reference[stage](x,y); capture[stage](bulk,x,y);
    }
    MPI_Finalize();
}
'''
    output.write_text(text)


def run(args):
    args.output.mkdir(parents=True,exist_ok=False)
    extract(args.source,args.output/'mars_boundary_reference.hpp')
    generate(args.source,args.output/'capture.cpp')
    common = ['-std=c++20','-Wall','-Wextra','-Werror','-I'+str(ROOT/'backend/distributed/unstructured/fem/segregated')]
    subprocess.check_call([args.cxx]+common+['-Wno-unused-variable','-I'+str(ROOT/'tests/reference/openaccel'),
                          str(args.output/'capture.cpp'),'-o',str(args.output/'capture')])
    (args.output/'exports').mkdir()
    subprocess.check_call([str(args.output/'capture'),str(args.output/'exports')])
    pack_boundary(args.output/'exports/boundary',args.output/'inputs.txt')
    subprocess.check_call([args.cxx,'-x','c++']+common+[
        str(ROOT/'examples/distributed/unstructured/mars_segregated_boundary_replay.cu'),'-o',str(args.output/'replay')])
    subprocess.check_call([str(args.output/'replay'),str(args.output/'inputs.txt')])
    print('PASS: capture expressions/writer/packer/host replay; actual STK and CUDA still pending')


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cxx',default='mpicxx')
    a=p.parse_args(); a.source=a.source.resolve(); a.output=a.output.resolve(); run(a)
