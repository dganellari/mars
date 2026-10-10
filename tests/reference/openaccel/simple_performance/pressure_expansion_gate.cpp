#include "channel.hpp"
#include "mars_segregated_simple_expansion.hpp"
#include <iostream>
#ifdef __CUDACC__
#define MARS_EXPANSION_TEST_HD __host__ __device__
#else
#define MARS_EXPANSION_TEST_HD
#endif
using namespace mars::segregated;
using namespace mars::segregated::runtime;
using namespace dsimple_gate;

void check(bool ok,const char* what) { if(!ok) throw std::runtime_error(what); }
struct HiddenField {
    SimpleMesh mesh; double *high,*low; bool constant;
    MARS_EXPANSION_TEST_HD void operator()(int n) const {
        high[n]=0x1p60; low[n]=constant?.25:mesh.x[n]*.125+mesh.y[n]*.25-mesh.z[n]*.0625;
    }
};
struct ExpandedActionTest {
    int* failed;
    MARS_EXPANSION_TEST_HD void operator()(int) const {
        const auto p=PressureValue{0x1p60,.25}+PressureValue{0,.125}*.3;
        const auto change=p-PressureValue{0x1p60,.25};
        bool ok=p.high==0x1p60 && fabs(change.rounded()-.0375)<4e-17;
        // The compact and reconstructed high gradients cancel, leaving a low pressure force.
        TetInteriorInput x{}; x.stage=0; x.pressure_expansion=true;
        for(int i=0;i<6;++i) {
            x.edges[2*i]=0; x.edges[2*i+1]=1; x.velocity_shape[4*i]=1;
            x.area[3*i]=1; x.shape_gradient[12*i]=1; x.shape_gradient[12*i+3]=-1;
        }
        for(int k=0;k<4;++k) { x.density[k]=1; x.influence_rhs[3*k]=1; }
        x.pressure[0]=0x1p60; x.pressure[1]=0x1p60; x.pressure_low[0]=.25;
        TetInteriorOutput y; tet_interior(x,y);
        ok=ok && y.flux[0]==-.25;
        x.pressure_expansion=false; tet_interior(x,y); ok=ok && y.flux[0]==0;
        BoundaryInput b{}; b.stage=1; b.pressure_expansion=true;
        for(int i=0;i<3;++i) {
            b.face_nodes[i]=i; b.nearest[i]=i; b.opposing[i]=3; b.density[i]=1;
            b.shape[3*i+i]=1; b.area[3*i]=1; b.influence_rhs[3*i]=1;
            b.gradient[12*i]=1; b.gradient[12*i+3]=-1;
        }
        b.pressure[0]=b.pressure[1]=0x1p60; b.pressure_low[0]=.25;
        BoundaryOutput out; boundary_block(b,out); ok=ok && out.flux[0]==-.25;
        double velocity[3]={1,0,0},influence[3]={1,1,1},gradient[3]={1,0,0},gradient_low[3]={0x1p-54,0,0};
        SimpleState state{}; state.velocity=velocity; state.influence=influence;
        ExpandedVelocityUpdate{state,gradient,gradient_low}(0);
        ok=ok && velocity[0]==-0x1p-54;
        double volume[1]={1},div[1]={0},factor[1]={0},pressure_low[1]={0},matrix_values[9]{},rhs[3]={1,0,0};
        int offsets[2]={0,1},columns[1]={0};
        state.volume=volume; state.mass_divergence=div; state.boundary_factor=factor; state.error=failed;
        state.pressure_gradient=gradient; state.pressure_gradient_low=gradient_low; state.pressure_low=pressure_low;
        SimpleControls ctl; ctl.alpha_u=1;
        SimpleMomentumNode{state,ctl,BlockCsrView<3>{1,offsets,columns,matrix_values,rhs}}(0);
        ok=ok && rhs[0]==-0x1p-54;
        // A closed outlet must not reopen when its positive mean difference lives only in the low part.
        int n0[1]={0},n1[1]={1},n2[1]={2},n3[1]={3},flags[3]={1,1,1};
        double cx[4]={0,1,0,0},cy[4]={0,0,1,0},cz[4]={0,0,0,1};
        const double xyz[12]={0,0,0,1,0,0,0,1,0,0,0,1};
        TetGeometry<double> geometry; ok=ok && tet_geometry(xyz,geometry);
        SimpleFace face{0,0,1}; SimpleMesh mesh{4,1,1,{n0,n1,n2,n3},cx,cy,cz,&face,&geometry};
        double u[12],d[12]{},pg[12]{},pfield[4]={10,10,10,10},plow[4]{},trace[3]={11,9,10},tl[3]={0x1p-54,0,0},flux[3]{},balance[4]{};
        double area[3]; tet_boundary_area(geometry,0,area);
        for(int k=0;k<4;++k) for(int j=0;j<3;++j) u[3*k+j]=area[j];
        SimpleState outlet{}; outlet.velocity=u; outlet.pressure=pfield; outlet.pressure_low=plow;
        outlet.pressure_gradient=pg; outlet.pressure_gradient_low=pg; outlet.influence=d;
        outlet.trace=trace; outlet.trace_low=tl; outlet.boundary_flux=flux; outlet.reversal=flags;
        outlet.mass_divergence=balance; outlet.error=failed;
        SimpleBoundary<1>{mesh,outlet,ctl,{},false,true}(0);
        ok=ok && flags[0]==1 && flags[1]==1 && flags[2]==1;
        if(!ok) *failed=1;
    }
};
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    int rank=0,ranks=0; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try {
#ifdef MARS_REPLAY_CUDA
        int local=0,devices=0; MPI_Comm node;
        MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node);
        MPI_Comm_rank(node,&local); MPI_Comm_free(&node);
        assembly_cuda_check(cudaGetDeviceCount(&devices)); check(devices>0,"no GPU");
        assembly_cuda_check(cudaSetDevice(local%devices));
#endif
        {
        auto input=channel(4,2,2); SimpleRunner run(input);
        Array<double> low(run.n),grad(3*run.n),grad_low(3*run.n),old_low(run.n),trace_low(3*run.b),outlet_low(run.n);
        run.state.pressure_low=low.data(); run.state.pressure_gradient_low=grad_low.data(); run.state.trace_low=trace_low.data();
        run.state.old_pressure_low=old_low.data();
        PressureIncidence cells(run.mesh,false),faces(run.mesh,true);
        launch(1,ExpandedActionTest{run.error.data()}); run.check("pair arithmetic or pressure flux action failed");
        launch(run.n,HiddenField{run.mesh,run.pressure.data(),low.data(),true});
        launch(run.n,ExpandedGradient{run.mesh,cells.offsets.data(),cells.entries.data(),run.pressure.data(),low.data(),run.volume.data(),grad.data(),grad_low.data(),run.error.data()});
        for(double v:grad.host()) check(v==0,"constant pressure gradient nonzero");
        launch(run.n,HiddenField{run.mesh,run.pressure.data(),low.data(),false});
        launch(run.n,ExpandedGradient{run.mesh,cells.offsets.data(),cells.entries.data(),run.pressure.data(),low.data(),run.volume.data(),grad.data(),grad_low.data(),run.error.data()});
        // The ordinary operator on the small field is an independent linearity oracle.
        gradient<1>(run.mesh,run.state,low.data(),run.sum,run.pg.data());
        auto g=grad.host(),gl=grad_low.host(),reference=run.pg.host();
        for(int i=0;i<3*run.n;++i) check(std::abs(g[i]+gl[i]-reference[i])<1e-15,"hidden field gradient lost");
        std::vector<int> owned; for(int i=rank;i<run.b;i+=ranks) owned.push_back(i);
        Array<int> device_owned(owned); PressureMomentReduction reduction;
        const auto mean=reduction.mean(MPI_COMM_WORLD,int(owned.size()),ExpandedMoment{run.mesh,run.state,device_owned.data()});
        // At 2^60, two-component arithmetic has a spacing of order 1e-14.
        check(mean.high==0x1p60 && std::abs(mean.low-(.25+.125-.03125))<1e-12,"paired outlet mean lost");
        launch(run.b,ExpandedTrace{run.mesh,run.state,run.controls,mean});
        launch(run.n,ExpandedOutletPressure{run.mesh,run.state,faces.offsets.data(),faces.entries.data(),run.outlet_area.data(),run.outlet_pressure.data(),outlet_low.data()});
        auto trace=run.trace.host(),tl=trace_low.host(),small=low.host();
        for(int i=0;i<run.b;++i) if(input.faces[i].kind==1) for(int j=0;j<3;++j) {
            const auto f=input.faces[i]; const int n=input.nodes[tet_face_node(f.ordinal,j)][f.element];
            check(std::abs(trace[3*i+j]+tl[3*i+j]-.95*(small[n]-(.25+.125-.03125)))<1e-12,"outlet trace lost low pressure");
        }
        // Repeated updates retain the pressure change used by nonlinear stopping.
        run.old_pressure.copy_from(run.pressure); old_low.copy_from(low);
        Array<double> increment(run.n),increment_low(std::vector<double>(run.n,.125));
        for(int step=0;step<2;++step) launch(run.n,ExpandedPressureUpdate{run.pressure.data(),low.data(),increment.data(),increment_low.data(),.3});
        const auto sums=reduce_sums(run.n,SimpleNodeSums{run.state,run.momentum.rhs.data(),run.old_velocity.data(),run.old_pressure.data()});
        check(std::abs(sums.pressure_change2/sums.volume-.075*.075)<1e-16,"pressure change ignores low history");
        run.check("expanded kernels failed");
        }
        if(!rank) std::cout<<"PASS: retained pressure actions, outlet mean and repeated update; ranks="<<ranks<<'\n';
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Finalize();
}
