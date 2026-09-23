#include "mars_segregated_simple.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace mars::segregated;
namespace {
int checks=0;
void check(bool ok) { ++checks; if (!ok) throw std::runtime_error("independent SIMPLE algebra mismatch"); }
void near(double a,double b) { check(std::isfinite(a) && std::abs(a-b)<=2e-13*std::max(1.,std::abs(b))); }
template<int C> void csr() {
    int offsets[]={0,2,5,7},columns[]={0,1,0,1,2,1,2};
    std::vector<double> blocks(7*C*C),rhs(3*C,0),values(7*C*C);
    for (int k=0;k<7;++k) for (int i=0;i<C;++i) for (int j=0;j<C;++j) blocks[C*C*k+C*i+j]=100*k+10*i+j+1;
    BlockCsrView<C> a{3,offsets,columns,blocks.data(),rhs.data()};
    std::vector<int> out_rows(3*C+1,-1),out_columns(7*C*C,-1);
    for (int r=0;r<3*C;++r) simple_scalar_row(a,r,out_rows.data(),out_columns.data(),values.data());
    check(out_rows.front()==0 && out_rows.back()==7*C*C);
    for (int r=0;r<3*C;++r) {
        check(out_rows[r+1]-out_rows[r]==C*(offsets[r/C+1]-offsets[r/C]));
        int previous=-1;
        for (int k=out_rows[r];k<out_rows[r+1];++k) {
            check(out_columns[k]>previous); previous=out_columns[k];
            const int block=block_position(a,r/C,out_columns[k]/C);
            near(values[k],100*block+10*(r%C)+out_columns[k]%C+1);
        }
    }
    double solution[3*C]; for (int i=0;i<3*C;++i) solution[i]=i+1;
    std::vector<double> squares(6*C);
    SimpleResidual residual{out_rows.data(),out_columns.data(),values.data(),rhs.data(),solution,squares.data()};
    for (int r=0;r<3*C;++r) {
        double expected=0;
        for (int k=offsets[r/C];k<offsets[r/C+1];++k) for (int j=0;j<C;++j) expected+=(100*k+10*(r%C)+j+1)*solution[C*columns[k]+j];
        residual(r); near(squares[2*r],expected*expected); near(squares[2*r+1],0);
    }
}
}
int main() {
    try {
        csr<1>(); csr<3>();
        double xyz[]={0,0,0,1,0,0,0,1,0,0,0,1}; TetGeometry<double> g; check(tet_geometry(xyz,g));
        int nodes[]={0,1,2,3}; double u[]={1,2,3,4,5,6,7,8,9,10,11,12},p[]={2,3,4,5},trace[]={7,8,9},flux[]={.1,.2,.3};
        SimpleControls c; SimpleFace face{0,1,0};
        auto inlet=simple_boundary(true,face,nodes,g,u,p,trace,flux,c,false);
        for (int j=0;j<3;++j) { near(inlet.values.velocity[j],u[j]); for (int f=0;f<3;++f) near(inlet.values.boundary_velocity[3*f+j],-.1/std::sqrt(3.)); }
        face.kind=1;
        auto outlet=simple_boundary(false,face,nodes,g,u,p,trace,flux,c,false);
        near(outlet.values.pressure[0],2);
        for (int f=0;f<3;++f) near(outlet.values.pressure[f+1],trace[f]);
        double pg[12]{},d[12]; std::fill(d,d+12,2.);
        check(native_boundary(outlet.values,outlet,g,nodes,pg,d));
        BoundaryOutput out; boundary_block(outlet.values,out);
        for (int f=0;f<3;++f) { near(out.lhs[4*(f+1)],1.); for (int k=1;k<4;++k) near(out.lhs[4*(f+1)+k],0); }
        // K annihilates constants, but the open-boundary pressure partial anchors them.
        double row_sum=0; for (double a:out.lhs) row_sum+=a; near(row_sum,3);
        face.kind=2;
        for (bool initialized:{false,true}) {
            auto wall=simple_boundary(true,face,nodes,g,u,p,trace,flux,c,initialized);
            for (int f=0;f<3;++f) { near(wall.values.wall_coefficient[f],initialized?.2:0); check(wall.nodes[f]==f+1); }
        }
        // Pressure relaxation must not enter velocity correction a second time.
        double increment[]={2,3,5},influence[]={.2,.3,.4},velocity[]={7,11,13},corrected[3];
        velocity_update(velocity,influence,increment,corrected);
        for (int j=0;j<3;++j) near(corrected[j],velocity[j]-influence[j]*increment[j]);
        near(pressure_update(10,4,.3),11.2);
        double area[]={1,0,0},prescribed[]={2,0,0}; near(inlet_flux_update(3,prescribed,area,2,.75),5);
        // A closed outlet must fail until artificial-wall selection is integrated.
        int one[]={0},two[]={1},three[]={2},four[]={3},error=0;
        double px[]={0,1,0,0},py[]={0,0,1,0},pz[]={0,0,0,1};
        double backward[12],zero_pressure[4]{},zero_gradient[36]{},zero_d[12]{},zero_trace[3]{},bf[3]{},div[4]{};
        std::fill(backward,backward+12,-1.); SimpleFace opening{0,1,1};
        SimpleMesh mesh{4,1,1,{one,two,three,four},px,py,pz,&opening,&g};
        SimpleState state{}; state.velocity=backward; state.pressure=zero_pressure;
        state.pressure_gradient=zero_gradient; state.influence=zero_d; state.trace=zero_trace;
        state.boundary_flux=bf; state.mass_divergence=div; state.error=&error;
        SimpleBoundary<1>{mesh,state,c,{},false,true}(0);
        check(error==1); for (double q:bf) near(q,0);
        std::cout<<"PASS: "<<checks<<" independent SIMPLE integration checks\n";
    } catch (const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1; }
}
