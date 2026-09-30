#!/usr/bin/env python3
"""Replay production Hex8 arithmetic on invented cells; no CUDA/MPI validation.

CUDA qualifiers, thread indices and atomicAdd are replaced for serial host execution.
The numerical helpers and tensor assembly body are extracted unchanged from the repo.
"""
import argparse
from pathlib import Path
import shlex
import subprocess
import tempfile


PREAMBLE = r'''
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#define __device__
#define __constant__
#define __global__
#define __forceinline__ inline
struct Dim { int x; } blockIdx{0}, blockDim{1}, threadIdx{0};
template<class T> void atomicAdd(T* p, T v) { *p += v; }
'''

CHECK = r'''
} }
using namespace mars::fem;
int checks=0, failures=0;
void near(double actual, double expected, double tolerance) {
    ++checks;
    if (!std::isfinite(actual) || std::abs(actual-expected)>tolerance) {
        if (failures++<8) std::printf("FAIL check %d: %.17g expected %.17g tolerance %.3g\n",
                                    checks,actual,expected,tolerance);
    }
}
void cell(const double T[3][3]) {
    const double ref[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                          {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    double c[8][3]={}, x[8],y[8],z[8],one[8],zero[8]={},rhs[8]={};
    double av[3][12],mdot[12]={},values[64]={},cached[3][96]={};
    unsigned nodes[8]; int dofs[8],rp[9],ci[64],diag[8]; uint8_t own[8];
    for(int n=0;n<8;++n) {
        for(int i=0;i<3;++i) for(int j=0;j<3;++j) c[n][i]+=T[i][j]*ref[n][j];
        x[n]=c[n][0]; y[n]=c[n][1]; z[n]=c[n][2]; one[n]=1; own[n]=1;
        nodes[n]=n; dofs[n]=n; rp[n]=n*8; diag[n]=n*8+n;
        for(int j=0;j<8;++j) ci[n*8+j]=j;
    }
    rp[8]=64;
    precomputeShapeDerivativesKernel<unsigned,double>(
        nodes,nodes+1,nodes+2,nodes+3,nodes+4,nodes+5,nodes+6,nodes+7,1,
        x,y,z,cached[0],cached[1],cached[2]);
    for(int ip=0;ip<12;++ip) {
        double grad[8][3]; computeShapeDerivatives(ip,c,grad);
        for(int component=0;component<3;++component) {
            double sum=0,cache_sum=0;
            for(int n=0;n<8;++n){sum+=grad[n][component];cache_sum+=cached[component][ip*8+n];}
            near(sum,0,1e-9); near(cache_sum,0,1e-9);
            for(int field=0;field<3;++field) {
                double g=0,cg=0;
                for(int n=0;n<8;++n){g+=c[n][field]*grad[n][component];cg+=c[n][field]*cached[component][ip*8+n];}
                near(g,field==component?1:0,1e-9);
                near(cg,field==component?1:0,1e-9);
            }
        }
        double a[3]; computeAreaVector(ip,c,a);
        for(int j=0;j<3;++j) av[j][ip]=a[j];
    }
    CSRMatrix<double> mat{rp,ci,values,diag,8,64,8};
    cvfem_hex_assembly_kernel_tensor<unsigned,double>(
        nodes,nodes+1,nodes+2,nodes+3,nodes+4,nodes+5,nodes+6,nodes+7,1,
        x,y,z,one,zero,zero,zero,zero,zero,mdot,av[0],av[1],av[2],dofs,own,&mat,rhs);
    const double volume=T[0][0]*(T[1][1]*T[2][2]-T[1][2]*T[2][1])
                       -T[0][1]*(T[1][0]*T[2][2]-T[1][2]*T[2][0])
                       +T[0][2]*(T[1][0]*T[2][1]-T[1][1]*T[2][0]);
    for(int i=0;i<8;++i) {
        double row=0,scale=0;
        for(int j=0;j<8;++j){row+=values[8*i+j];scale+=std::abs(values[8*i+j]);}
        near(row,0,1e-12*std::max(1.,scale));
    }
    // For an affine physical field, integrated diffusion energy is V*|grad q|^2.
    for(int field=0;field<3;++field) {
        double energy=0;
        for(int i=0;i<8;++i)for(int j=0;j<8;++j) energy+=c[i][field]*values[8*i+j]*c[j][field];
        near(energy,volume,5e-9*volume);
    }
}
int main() {
    const double maps[][3][3]={
        {{2,0,0},{0,3,0},{0,0,4}},
        {{0,10./149.,0},{-.0001,0,0},{0,0,.06}},
        {{1,.25,.1},{-.2,2,.3},{.1,-.15,.7}},
        {{0,0,2},{3,0,0},{0,4,0}}
    };
    for(const auto& map:maps) cell(map);
    std::printf("%s: %d host arithmetic checks, %d failures; CUDA/MPI not executed\n",
                failures?"FAIL":"PASS",checks,failures);
    return failures?1:0;
}
'''


def source(root):
    fem = root / 'backend/distributed/unstructured/fem'
    header = (fem / 'mars_cvfem_hex_kernel.hpp').read_text()
    tensor = (fem / 'mars_cvfem_hex_kernel_tensor.hpp').read_text()
    utils = (fem / 'mars_cvfem_utils.hpp').read_text()
    numerical = header[header.index('namespace mars'):header.index('// CVFEM assembly kernel for hex elements')]
    numerical += tensor[tensor.index('template<typename KeyType'):tensor.rindex('} // namespace fem')]
    start = utils.index('template<typename RealType>\n__device__ inline void invert3x3_generic')
    numerical += utils[start:utils.index('// Host function to pre-compute shape derivatives')]
    return PREAMBLE + numerical + CHECK


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--cxx', default='c++')
    args = parser.parse_args()
    scratch = args.root / '.local-worktrees' / 'poiseuille-host'
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='gradient-', dir=str(scratch)) as tmp:
        cpp, executable = Path(tmp)/'gate.cpp', Path(tmp)/'gate'
        cpp.write_text(source(args.root))
        # Existing CUDA bodies have unused locals, mixed index signs and CUDA pragmas.
        subprocess.run(shlex.split(args.cxx) + ['-std=c++17', '-O2', '-Wall', '-Wextra', '-Werror',
                       '-Wno-unused-variable', '-Wno-sign-compare', '-Wno-unknown-pragmas',
                       str(cpp), '-o', str(executable)], check=True)
        return subprocess.run([str(executable)]).returncode


if __name__ == '__main__':
    raise SystemExit(main())
