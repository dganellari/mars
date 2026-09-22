#include "mars_segregated_boundary.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

using mars::segregated::BoundaryInput;
using mars::segregated::BoundaryOutput;
using mars::segregated::boundary_block;
int checks = 0;
void close(double a, double b, double tolerance=1e-12) {
    ++checks;
    if (!std::isfinite(a) || std::abs(a-b)>tolerance*std::max(1.0,std::abs(b)))
        throw std::runtime_error("boundary algebra mismatch: "+std::to_string(a)+" versus "+std::to_string(b));
}
BoundaryInput tetrahedron(int stage) {
    BoundaryInput x; x.stage=stage;
    const double gradient[12]={-1,-1,-1,1,0,0,0,1,0,0,0,1};
    for (int s=0;s<3;++s) {
        x.face_nodes[s]=s+1; x.nearest[s]=s+1; x.opposing[s]=0;
        x.density[s]=1; x.viscosity[s]=2; x.wall_coefficient[s]=4;
        for(int j=0;j<3;++j) {
            x.area[3*s+j]=1.0/6;
            x.shape[3*s+j]=(s==j ? 11.0/18 : 7.0/36);
            x.influence_lhs[3*s+j]=x.influence_rhs[3*s+j]=2;
        }
        for(int j=0;j<12;++j) x.gradient[12*s+j]=gradient[j];
    }
    x.bc_multiplier[0]=1;
    return x;
}
int main() {
    try {
        BoundaryOutput base, trial;
        auto x=tetrahedron(1);
        // p=2x-3y+5z; compact and reconstructed gradients agree exactly.
        x.pressure[1]=2; x.pressure[2]=-3; x.pressure[3]=5;
        for(int n=0;n<4;++n) {
            x.velocity[3*n]=2; x.velocity[3*n+1]=3; x.velocity[3*n+2]=1;
            x.pressure_gradient[3*n]=2; x.pressure_gradient[3*n+1]=-3; x.pressure_gradient[3*n+2]=5;
        }
        boundary_block(x,base);
        for(int s=0;s<3;++s) { close(base.flux[s],1); close(base.lhs[4*(s+1)],1); close(base.rhs[s+1],-1); }
        auto p=x;p.pressure[0]+=.2;boundary_block(p,trial);
        double delta=0;for(int s=0;s<3;++s)delta+=trial.flux[s]-base.flux[s];close(delta,.6);
        p=x;for(int n=1;n<4;++n)p.pressure[n]+=.2;boundary_block(p,trial);
        delta=0;for(int s=0;s<3;++s)delta+=trial.flux[s]-base.flux[s];close(delta,-.6);
        // Constant volume increment has nonzero action because trace columns are masked.
        for(int r=1;r<4;++r) { double sum=0;for(int c=0;c<4;++c)sum+=base.lhs[4*r+c];close(sum,1); }
        p=x;for(double& d:p.influence_lhs)d*=2;boundary_block(p,trial);
        for(int k=0;k<16;++k)close(trial.lhs[k],2*base.lhs[k]);
        for(int s=0;s<3;++s)close(trial.flux[s],base.flux[s]);
        // Nonuniform face coefficients: each sample uses its own shape weights.
        p=tetrahedron(1);const double ds[3]={1,2,4};
        for(int s=0;s<3;++s)for(int j=0;j<3;++j)p.influence_lhs[3*s+j]=ds[s];
        boundary_block(p,trial);
        close(trial.lhs[4],.5*(11.0/18+7.0/36*6));
        close(trial.lhs[8],.5*(22.0/18+7.0/36*5));
        close(trial.lhs[12],.5*(44.0/18+7.0/36*3));
        p=tetrahedron(1);p.pressure_gradient[3]=6;boundary_block(p,trial);
        close(trial.flux[0],1); close(trial.flux[1],0);close(trial.flux[2],0);
        // Flags suppress the whole open sample; no slip-wall replacement diagonal.
        p=x;for(int& f:p.reversal)f=1;boundary_block(p,trial);
        for(double v:trial.lhs)close(v,0);for(double v:trial.rhs)close(v,0);
        p=tetrahedron(0);for(int s=0;s<3;++s)p.boundary_velocity[3*s]=-6;
        boundary_block(p,trial);
        for(int s=0;s<3;++s)close(trial.flux[s],-1);
        for(double v:trial.lhs)close(v,0);
        // Finite differences of the momentum residual -rhs, holding side values fixed.
        for(int stage:{3,4,5}) {
            p=tetrahedron(stage);
            if(stage==5)for(int s=0;s<3;++s)p.face_nodes[s]=p.nearest[s]=s;
            const int nodes=stage==5?3:4, width=3*nodes;
            for(int i=0;i<width;++i)p.velocity[i]=.1*(i+1);
            for(int s=0;s<3;++s)p.stored_flux[s]=.2*(s+1);
            boundary_block(p,base);
            for(int c=0;c<width;++c) {
                if(stage==3 && c>=3)continue; // prescribed inlet face values are frozen
                auto z=p;z.velocity[c]+=1e-5;boundary_block(z,trial);
                for(int r=0;r<width;++r)close((base.rhs[r]-trial.rhs[r])/1e-5,base.lhs[r*width+c],1e-9);
            }
            if(stage==3)for(int r=0;r<width;++r)for(int c=3;c<width;++c)close(base.lhs[r*width+c],0);
        }
        // Wall: normal velocity is unconstrained by this tangential law; tangential work dissipates.
        p=tetrahedron(5);for(int s=0;s<3;++s){p.face_nodes[s]=p.nearest[s]=s;for(int j=0;j<3;++j)p.velocity[3*s+j]=1;}
        boundary_block(p,trial);for(double v:trial.rhs)close(v,0);
        for(int s=0;s<3;++s){p.velocity[3*s]=1;p.velocity[3*s+1]=-1;p.velocity[3*s+2]=0;}
        boundary_block(p,trial);double dissipation=0;for(int i=0;i<9;++i)dissipation-=p.velocity[i]*trial.rhs[i];close(dissipation,24);
        // Outlet viscous traction has no normal component, including cross-component stress.
        p=tetrahedron(4);for(int i=0;i<12;++i)p.velocity[i]=.2*i*i;
        boundary_block(p,trial);for(int s=0;s<3;++s)close(trial.rhs[3*(s+1)]+trial.rhs[3*(s+1)+1]+trial.rhs[3*(s+1)+2],0);
        std::cout<<"PASS: "<<checks<<" independent boundary algebra checks\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
