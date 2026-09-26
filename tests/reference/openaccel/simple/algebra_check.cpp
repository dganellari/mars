#include "mars_segregated_simple.hpp"
#include "mars_segregated_simple_metrics.hpp"
#include "mars_segregated_native_mapping.hpp"
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
void native_mapping() {
    const int source[]={17,3,29,8}; int permutation[]={0,1,2,3};
    const auto identity=simple_cell_key(source);
    do {
        int cell[4]; for (int j=0;j<4;++j) cell[j]=source[permutation[j]];
        check(simple_cell_key(cell)==identity);
        for (int f=0;f<4;++f) {
            int face[3]; for (int j=0;j<3;++j) face[j]=source[tet_face_node(f,j)];
            std::sort(face,face+3);
            do {
                const int mapped=simple_native_face(face,cell); check(mapped>=0);
                // Compare the excluded vertex, independently of face numbering/orientation.
                int absent=-1;
                for (int value:source) if (std::find(face,face+3,value)==face+3) absent=value;
                check(cell[tet_opposite_node(mapped)]==absent);
            } while (std::next_permutation(face,face+3));
        }
    } while (std::next_permutation(permutation,permutation+4));
    const int absent[]={17,3,99},duplicate[]={17,17,3};
    check(simple_native_face(absent,source)==-1); check(simple_native_face(duplicate,source)==-1);
    const int keys[]={3,8,17,29};
    for (int j=0;j<4;++j) check(simple_find_key(keys,4,keys[j])==j);
    for (int key:{-1,0,4,30}) check(simple_find_key(keys,4,key)==-1);
    check(simple_find_key(keys,0,3)==-1);
}
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
        csr<1>(); csr<3>(); native_mapping();
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
        // A caller without persistent reversal storage must reject closure.
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
        // Persist closure, omit both opening blocks, retain trace, then reopen.
        int flags[3]{}; state.reversal=flags; error=0;
        SimpleBoundary<1>{mesh,state,c,{},false,true}(0);
        check(error==0); for (int flag:flags) check(flag==1);
        int offsets[]={0,4,8,12,16},columns[]={0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
        double pa[16]{},pb[4]{},ma[144]{},mb[12]{};
        SimpleBoundary<1>{mesh,state,c,{4,offsets,columns,pa,pb},true,false}(0);
        SimpleBoundary<3>{mesh,state,c,{4,offsets,columns,ma,mb},true,false}(0);
        for (double v:pa) near(v,0); for (double v:pb) near(v,0);
        for (double v:ma) near(v,0); for (double v:mb) near(v,0);
        double moments[2]{};
        std::fill(zero_trace,zero_trace+3,2.);
        SimpleTraceMoment{mesh,state,moments}(0); near(moments[0],0); near(moments[1],0);
        SimpleTrace{mesh,state,c,moments}(0);
        for (double t:zero_trace) near(t,2); check(error==0);
        std::fill(backward,backward+12,1.);
        // Outward velocity alone cannot reopen against the pressure condition.
        SimpleBoundary<1>{mesh,state,c,{},true,true}(0);
        for (int flag:flags) check(flag==1);
        std::fill(zero_pressure,zero_pressure+4,3.);
        SimpleBoundary<1>{mesh,state,c,{},true,true}(0);
        for (int flag:flags) check(flag==0);
        for (double q:bf) near(q,0);
        SimpleTraceMoment{mesh,state,moments}(0);
        near(moments[1],std::sqrt(3.)/2); near(moments[0]/moments[1],3);
        SimpleTrace{mesh,state,c,moments}(0); for (double t:zero_trace) near(t,0);
        SimpleBoundary<1>{mesh,state,c,{},true,true}(0);
        for (double q:bf) near(q,.375); near(div[1]+div[2]+div[3],1.125);
        // The newly open pressure derivative restores its constant-mode anchor.
        std::fill(zero_d,zero_d+12,2.);
        SimpleBoundary<1>{mesh,state,c,{4,offsets,columns,pa,pb},true,false}(0);
        double anchor=0; for (double v:pa) anchor+=v; near(anchor,3);
        SimpleSums sums; sums.volume=2; sums.momentum2=8; sums.continuity2=2;
        sums.velocity_change2=.02; sums.pressure_change2=.0002;
        sums.inlet=-1; sums.outlet=.9; sums.continuity=-.1; sums.inlet_area=1;
        SimpleControls scales; scales.density=2; scales.inlet_speed=.5;
        auto metrics=simple_metrics(sums,scales);
        near(metrics.momentum,4); near(metrics.continuity,1); near(metrics.flux,.1);
        near(metrics.velocity_change,.2); near(metrics.pressure_change,.02); near(metrics.cancellation,0);
        check(metrics.finite); check(!simple_converged(metrics,100,0,1e-6,1e-6,1e-6));
        SimpleMetrics good{0,0,0,0,0,0,true};
        check(simple_converged(good,2,0,1e-6,1e-6,1e-6));
        check(!simple_converged(good,1,0,1e-6,1e-6,1e-6));
        check(!simple_converged(good,2,1,1e-6,1e-6,1e-6));
        for (double* value:{&good.momentum,&good.continuity,&good.flux,&good.velocity_change,&good.pressure_change,&good.cancellation,&good.flux_change}) {
            *value=1.; check(!simple_converged(good,10,0,1e-6,1e-6,1e-6)); *value=0;
        }
        double fresh_flux[]={1,2},prior_flux[]={.9,2.2};
        near(SimpleFluxChange{fresh_flux,prior_flux}(1).flux_change,.2);
        good.finite=false; check(!simple_converged(good,10,0,1e-6,1e-6,1e-6));
        // Boundary RHS relaxation must not hide a physical momentum defect.
        double volumes[4]={2,2,2,2},factors[4]={1,0,0,0},rhs[12]={1.5,0,0};
        double old_u[12]{},old_p[4]{}; state.volume=volumes; state.boundary_factor=factors;
        state.mass_divergence=div;
        auto node=SimpleNodeSums{state,rhs,old_u,old_p}(0); near(node.momentum2,2);
        rhs[0]=std::numeric_limits<double>::quiet_NaN();
        check(SimpleNodeSums{state,rhs,old_u,old_p}(0).invalid);
        std::cout<<"PASS: "<<checks<<" independent SIMPLE integration checks\n";
    } catch (const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1; }
}
