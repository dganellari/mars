#!/usr/bin/env python3
"""Compile production channel scatter/mask/lift arithmetic on invented boxes.

Serial host replay only: this does not compile CUDA launch sites or exercise MPI.
"""
from pathlib import Path
import subprocess
import tempfile
import cvfem_gradient_host_check as gradient

ROOT = Path(__file__).resolve().parents[1]
FEM = ROOT / 'backend/distributed/unstructured/fem'


def extract(text, name):
    start = text.rfind('template<', 0, text.index(name))
    begin = text.index('{', text.index(name))
    depth = 1
    end = begin + 1
    while depth:
        depth += (text[end] == '{') - (text[end] == '}')
        end += 1
    return text[start:end] + '\n'


CHECK = r'''
using Vec=std::vector<double>;
int checks=0;
void near(double a,double b,double tol=2e-10) {
    ++checks;
    if(!std::isfinite(a)||std::abs(a-b)>tol) {
        std::fprintf(stderr,"FAIL %d: %.17g != %.17g (tol %.3g)\n",checks,a,b,tol);
        std::exit(1);
    }
}
struct Grid {
    static constexpr int nx=5,ny=5,n=nx*ny*2,ne=(nx-1)*(ny-1);
    std::array<std::vector<unsigned>,8> c;
    std::array<Vec,3> area;
    Vec mass=Vec(n),x=Vec(n),y=Vec(n),z=Vec(n),ain=Vec(n),aout=Vec(n),target=Vec(n);
    std::vector<uint8_t> fixed=std::vector<uint8_t>(n),pressure=std::vector<uint8_t>(n),own=std::vector<uint8_t>(n,1);
    std::vector<int> active;
    static int id(int i,int j,int k){return (i*ny+j)*2+k;}
    Grid(bool stretched) {
        const double ys[ny]={0,.0001,.2,.8,1};
        for(int i=0;i<nx;++i)for(int j=0;j<ny;++j)for(int k=0;k<2;++k){
            int a=id(i,j,k); x[a]=.5*i; y[a]=stretched?ys[j]:double(j)/(ny-1);z[a]=.06*k;
            fixed[a]=(i==0||j==0||j==ny-1);
            pressure[a]=(i==nx-1||(i==0&&(j==0||j==ny-1)));
            target[a]=(i==0&&j!=0&&j!=ny-1)?1:0;
            if(!pressure[a])active.push_back(a);
        }
        for(auto& v:c)v.resize(ne);
        for(auto& v:area)v.resize(ne*12);
        int e=0;
        for(int i=0;i<nx-1;++i)for(int j=0;j<ny-1;++j,++e){
            unsigned ns[8]={unsigned(id(i,j,0)),unsigned(id(i+1,j,0)),unsigned(id(i+1,j+1,0)),unsigned(id(i,j+1,0)),
                            unsigned(id(i,j,1)),unsigned(id(i+1,j,1)),unsigned(id(i+1,j+1,1)),unsigned(id(i,j+1,1))};
            double coords[8][3];
            double vol=.5*(y[ns[3]]-y[ns[0]])*.06;
            for(int a=0;a<8;++a){c[a][e]=ns[a];coords[a][0]=x[ns[a]];coords[a][1]=y[ns[a]];coords[a][2]=z[ns[a]];mass[ns[a]]+=vol/8;}
            near(channel_rectilinear_hex(coords,1e-10),1);
            coords[0][0]+=.01;near(channel_rectilinear_hex(coords,1e-10),0);coords[0][0]-=.01;
            for(int ip=0;ip<12;++ip){double av[3];computeAreaVector(ip,coords,av);for(int d=0;d<3;++d)area[d][e*12+ip]=av[d];}
            if(i==0)for(int a:{0,3,4,7})ain[ns[a]]-=vol/.5/4;
            if(i==nx-2)for(int a:{1,2,5,6})aout[ns[a]]+=vol/.5/4;
        }
    }
    Vec div(const std::array<Vec,3>& u,int begin=0,int end=ne)const{
        Vec r(n);
        for(int e=begin;e<end;++e){blockIdx.x=e-begin;
            computeDivergencePerNodeKernel<unsigned,double,HexTag>(c[0].data(),c[1].data(),c[2].data(),c[3].data(),c[4].data(),c[5].data(),c[6].data(),c[7].data(),
            u[0].data(),u[1].data(),u[2].data(),area[0].data(),area[1].data(),area[2].data(),r.data(),begin,end-begin);}
        return r;
    }
    std::array<Vec,3> grad(const Vec& p,bool project=true,bool legacy=false,int begin=0,int end=ne)const{
        std::array<Vec,3> g{Vec(n),Vec(n),Vec(n)};
        for(int e=begin;e<end;++e){blockIdx.x=e-begin;
            if(legacy)computeGradientPerNodeKernel<unsigned,double,HexTag>(c[0].data(),c[1].data(),c[2].data(),c[3].data(),c[4].data(),c[5].data(),c[6].data(),c[7].data(),
            p.data(),area[0].data(),area[1].data(),area[2].data(),g[0].data(),g[1].data(),g[2].data(),begin,end-begin);
            else applyDivTransposePerNodeKernel<unsigned,double,HexTag>(c[0].data(),c[1].data(),c[2].data(),c[3].data(),c[4].data(),c[5].data(),c[6].data(),c[7].data(),
            p.data(),area[0].data(),area[1].data(),area[2].data(),g[0].data(),g[1].data(),g[2].data(),begin,end-begin);}
        for(int i=0;i<n;++i){for(int d=0;d<3;++d)g[d][i]/=mass[i];
            if(project){blockIdx.x=i;project_channel_gradient_kernel(fixed.data(),g[0].data(),g[1].data(),g[2].data(),n);}}
        return g;
    }
    Vec advect(const std::array<Vec,3>& u,const Vec& q,int mode=3,bool boundary=true)const{
        Vec r(n);
        for(int e=0;e<ne;++e){blockIdx.x=e;
            explicitAdvectionFluxScatterPerNodeKernel<unsigned,double,HexTag>(c[0].data(),c[1].data(),c[2].data(),c[3].data(),c[4].data(),c[5].data(),c[6].data(),c[7].data(),
            u[0].data(),u[1].data(),u[2].data(),q.data(),area[0].data(),area[1].data(),area[2].data(),r.data(),0,ne,mode,
            nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr);}
        if(boundary)for(int i=0;i<n;++i){blockIdx.x=i;add_channel_advection_boundary_kernel(ain.data(),aout.data(),target.data(),u[0].data(),q.data(),own.data(),r.data(),n);}
        return r;
    }
    Vec residual(const std::array<Vec,3>& u)const{
        Vec r=div(u);
        for(int i=0;i<n;++i){blockIdx.x=i;add_channel_opening_flux_kernel(ain.data(),aout.data(),target.data(),u[0].data(),own.data(),r.data(),n);}
        return r;
    }
};
void run(bool stretch){
    Grid g(stretch);const int n=Grid::n,m=int(g.active.size());
    auto constant=g.grad(Vec(n,3),false),old=g.grad(Vec(n,3),false,true);
    double old_error=0;
    for(int i=0;i<n;++i)for(int d=0;d<3;++d){near(constant[d][i],0);old_error=std::max(old_error,std::abs(old[d][i]));}
    near(old_error>1,1);
    // -M^-1 D^T is the physical gradient, including boundary nodes.
    auto affine=g.grad(g.x,false);
    for(int i=0;i<n;++i){near(affine[0][i],-1,1e-9);near(affine[1][i],0);near(affine[2][i],0);}
    Vec p(n);std::array<Vec,3> u{Vec(n),Vec(n),Vec(n)};
    for(int i=0;i<n;++i){p[i]=g.pressure[i]?0:std::sin(.71*i);
        u[0][i]=g.fixed[i]?g.target[i]:1+.04*std::sin(i);
        u[1][i]=g.fixed[i]?0:.03*std::cos(i);}
    auto advection=g.advect(u,p); double adv_energy=0,bound_energy=0;
    for(int i=0;i<n;++i){adv_energy+=p[i]*advection[i];bound_energy-=.5*p[i]*p[i]*channel_opening_flux(g.ain[i],g.aout[i],g.target[i],u[0][i]);}
    near(adv_energy,bound_energy);
    auto old_adv=g.advect(u,p,0);double old_energy=0;for(int i=0;i<n;++i)old_energy+=p[i]*old_adv[i];
    near(std::abs(old_energy-bound_energy)>1e-8,1);
    auto saved_target=g.target;g.target=Vec(n,1);
    std::array<Vec,3> uniform{Vec(n,1),Vec(n),Vec(n)};
    auto uniform_adv=g.advect(uniform,Vec(n,1));for(double a:uniform_adv)near(a,0);
    g.target=saved_target;
    auto gp=g.grad(p),raw=g.grad(p,false);auto du=g.div(u);
    double lhs=0,rhs=0;
    for(int i=0;i<n;++i){lhs+=p[i]*du[i];for(int d=0;d<3;++d)rhs+=raw[d][i]*g.mass[i]*u[d][i];}
    near(lhs,rhs);
    // Emulate additive element partitions before the owner applies M^-1 and Q.
    auto split_a=g.grad(p,true,false,0,Grid::ne/2),split_b=g.grad(p,true,false,Grid::ne/2,Grid::ne);
    for(int i=0;i<n;++i)for(int d=0;d<3;++d)near(split_a[d][i]+split_b[d][i],gp[d][i],1e-8);
    Vec A(m*m),L(m*m);
    for(int j=0;j<m;++j){Vec e(n);e[g.active[j]]=1;auto a=g.div(g.grad(e));for(int i=0;i<m;++i)A[i*m+j]=a[g.active[i]];}
    for(int i=0;i<m;++i)for(int j=0;j<m;++j)near(A[i*m+j],A[j*m+i],1e-9);
    // Cholesky must succeed after only outlet and zero-row pressure elimination.
    for(int i=0;i<m;++i)for(int j=0;j<=i;++j){double a=A[i*m+j];for(int k=0;k<j;++k)a-=L[i*m+k]*L[j*m+k];
        if(i==j){near(a>0,1);L[i*m+j]=std::sqrt(a);}else L[i*m+j]=a/L[j*m+j];}
    const auto before=g.residual(u);
    for(double h:{.01,2*.01/3}){
        Vec b(m),phi(n);
        for(int i=0;i<m;++i){double a=-before[g.active[i]]/h;for(int j=0;j<i;++j)a-=L[i*m+j]*b[j];b[i]=a/L[i*m+i];}
        for(int i=m-1;i>=0;--i){double a=b[i];for(int j=i+1;j<m;++j)a-=L[j*m+i]*phi[g.active[j]];phi[g.active[i]]=a/L[i*m+i];}
        auto q=u;auto corr=g.grad(phi);auto action=g.div(corr);
        for(int i=0;i<n;++i)for(int d=0;d<3;++d)q[d][i]+=h*corr[d][i];
        auto after=g.residual(q);
        for(int i=0;i<n;++i){near(q[2][i],0);if(g.fixed[i]){near(q[0][i],g.target[i]);near(q[1][i],0);}
            if(!g.pressure[i]){near(after[i],0,1e-10);near(after[i]-before[i],h*action[i],1e-10);}
            if(g.x[i]==0&&(g.y[i]==0||g.y[i]==1))near(after[i],0);}
        double total=0,flux=0;for(int i=0;i<n;++i){total+=after[i];flux+=channel_opening_flux(g.ain[i],g.aout[i],g.target[i],q[0][i]);}near(total,flux);
        // Omitting Q in the pressure action must be detectable.
        auto wrong=g.div(g.grad(phi,false));double gap=0;for(int i:g.active)gap=std::max(gap,std::abs(h*(wrong[i]-action[i])));near(gap>1e-8,1);
    }
    // Nonzero Dirichlet elimination must preserve a constant solution.
    const int rp[4]={0,2,5,7},ci[7]={0,1,0,1,2,1,2},dn[3]={2,0,1};
    const uint8_t fixed[3]={0,1,1};
    const double a[7]={3,-1,-1,4,-2,-2,4},t[3]={0,2,2};double lift[3]={};
    for(int i=0;i<3;++i){blockIdx.x=i;build_channel_velocity_lift_kernel(rp,ci,dn,fixed,a,t,lift,3);}
    near(lift[0],0);near(lift[1],6);near(lift[2],0);
    // A_ff*2 = mass_rhs + lift = 2+6; omission would give .5 instead of2.
    near((2+lift[1])/4,2);
}
void validation_checks(){
    near(channel_state_converged(1e-8,1e-9,1e-8,1e-6,1e-6),1);
    for(double bad:{-1.,1e-3,double(INFINITY),double(NAN)}){
        near(channel_state_converged(bad,0,0,1e-6,1e-6),0);
        near(channel_state_converged(0,bad,0,1e-6,1e-6),0);
        near(channel_state_converged(0,0,bad,1e-6,1e-6),0);
    }
    const int dof[3]={2,0,1};const uint8_t own[3]={1,1,0},fixed[3]={1,0,0};
    const double mass[3]={2,4,8},q[3]={1,2,3},old[3]={.9,1.8,2.7},
        adv[3]={.1,.2,.3},adv_old[3]={.01,.02,.03},grad[3]={1.8,2,3},target[3]={7,9,11};
    for(double dt:{.01,.125}){
        double star[3]={-99,-99,-99};
        for(int i=0;i<3;++i){blockIdx.x=i;
            applyPredictorPerNodeKernel(star,q,adv,grad,mass,mass,target,fixed,dof,own,dt,.5,.6,3,3);}
        near(star[0],1+dt*(.1/8-.5*1.8+.5*.6));near(star[1],9);near(star[2],-99);
        for(int i=0;i<3;++i){blockIdx.x=i;
            applyPredictorBdf2PerNodeKernel(star,q,old,adv,adv_old,grad,mass,target,fixed,dof,own,dt,.5,.6,3,3);}
        near(star[0],4./3-.9/3+(2*dt/3)*((.2-.01)/8-.5*1.8+.5*.6));
        near(star[1],9);near(star[2],-99);
        double result[3]={-99,-99,-99};
        for(int i=0;i<3;++i){blockIdx.x=i;
            applyCorrectorPerNodeKernel(result,star,grad,target,fixed,dof,own,2*dt/3,.5,3,3);}
        near(result[0],star[0]-(2*dt/3)*.5*grad[0]);near(result[1],9);near(result[2],-99);
    }
}
int main(){run(false);run(true);validation_checks();std::printf("PASS: %d production host projection checks; CUDA/MPI not executed\n",checks);}
'''


def source():
    solver=(FEM/'mars_ns_channel_solver.hpp').read_text()
    header=(FEM/'mars_cvfem_hex_kernel.hpp').read_text()
    utils=(FEM/'mars_cvfem_utils.hpp').read_text()
    geom=header[header.index('namespace mars'):header.index('// CVFEM assembly kernel for hex elements')]
    text=gradient.PREAMBLE+'\n#define __host__\n#include <array>\n#include <vector>\n#include <cstdlib>\n#include <type_traits>\n'+geom+'\n} }\nusing namespace mars::fem;\n'
    text+='struct HexTag {}; struct TetTag {}; template<class T> struct ElemTraits {static constexpr int NodesPerElem=8, ScsPerElem=12;};\nint d_tetLRSCV[12]={};\n'
    text+=utils[utils.index('__device__ __constant__ int d_hexLRSCV'):utils.index('// GPU kernel: Count NNZ')]
    text+=extract(solver,'void scsLR(')
    text+='#include "backend/distributed/unstructured/fem/mars_channel_projection.hpp"\n'
    for name in ('computeDivergencePerNodeKernel','applyDivTransposePerNodeKernel','computeGradientPerNodeKernel',
                 'project_channel_gradient_kernel','build_channel_velocity_lift_kernel','add_channel_opening_flux_kernel',
                 'applyPredictorPerNodeKernel','applyPredictorBdf2PerNodeKernel','applyCorrectorPerNodeKernel',
                 'explicitAdvectionFluxScatterPerNodeKernel','add_channel_advection_boundary_kernel'):
        text+=extract(solver,'void '+name+'(')
    return text+CHECK


def main():
    scratch = ROOT / '.local-worktrees' / 'poiseuille-host'
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='projection-', dir=str(scratch)) as folder:
        cpp=Path(folder)/'check.cpp';exe=Path(folder)/'check'
        cpp.write_text(source())
        subprocess.run(['c++','-std=c++17','-O2','-Wall','-Wextra','-Werror','-Wno-unknown-pragmas','-I',str(ROOT),str(cpp),'-o',str(exe)],check=True)
        subprocess.run([str(exe)],check=True)


if __name__=='__main__':
    main()
