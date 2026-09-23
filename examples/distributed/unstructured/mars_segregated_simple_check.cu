#ifdef MARS_REPLAY_CUDA
#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/solvers/mars_hypre_gmres_solver.hpp"
#include "mars_segregated_assembly_device.hpp"
#include <thrust/iterator/permutation_iterator.h>
#endif
#include "fixture.hpp"
#include <algorithm>
#include <numeric>
using namespace mars::segregated;
using namespace simple_check;

#ifdef MARS_REPLAY_CUDA
template<class Function> __global__ void simple_apply(int n,Function f) { int i=blockIdx.x*blockDim.x+threadIdx.x; if (i<n) f(i); }
template<class Function> void launch(int n,Function f) { if (n) { simple_apply<<<(n+63)/64,64>>>(n,f); assembly_cuda_check(cudaGetLastError()); } }
template<class T> struct Array {
    thrust::device_vector<T> values;
    explicit Array(std::size_t n):values(n) {}
    explicit Array(const std::vector<T>& x):values(x) {}
    T* data() { return device_data(values); }
    void zero() { thrust::fill(values.begin(),values.end(),T{}); }
    std::vector<T> host() { std::vector<T> x(values.size()); assembly_cuda_check(cudaMemcpy(x.data(),data(),x.size()*sizeof(T),cudaMemcpyDeviceToHost)); return x; }
};
#else
template<class Function> void launch(int n,Function f) { for (int i=0;i<n;++i) f(i); }
template<class T> struct Array {
    std::vector<T> values;
    explicit Array(std::size_t n):values(n) {}
    explicit Array(const std::vector<T>& x):values(x) {}
    T* data() { return values.data(); }
    void zero() { std::fill(values.begin(),values.end(),T{}); }
    std::vector<T> host() { return values; }
};
#endif
struct Graph {
#ifdef MARS_REPLAY_CUDA
    DeviceBlockGraph graph;
    Graph(SimpleMesh m) { TetMeshView<int,double> v{{m.nodes[0],m.nodes[1],m.nodes[2],m.nodes[3]},m.x,m.y,m.z,m.node_count,m.element_count}; graph.build(v); }
    int blocks() const { return int(graph.columns.size()); }
    template<int C> BlockCsrView<C> view(double* a,double* b) { return graph.view<C>(a,b); }
#else
    std::vector<int> offsets,columns;
    explicit Graph(SimpleMesh m) {
        std::vector<std::set<int>> rows(m.node_count);
        for (int e=0;e<m.element_count;++e) for (int i=0;i<4;++i) for (int j=0;j<4;++j) rows[m.nodes[i][e]].insert(m.nodes[j][e]);
        offsets.push_back(0);
        for (int n=0;n<m.node_count;++n) { rows[n].insert(n); columns.insert(columns.end(),rows[n].begin(),rows[n].end()); offsets.push_back(int(columns.size())); }
    }
    int blocks() const { return int(columns.size()); }
    template<int C> BlockCsrView<C> view(double* a,double* b) { return {int(offsets.size())-1,offsets.data(),columns.data(),a,b}; }
#endif
};
template<int C> struct LinearSystem {
    int rows;
    Array<double> blocks,rhs,increment,residual;
#ifdef MARS_REPLAY_CUDA
    using Solver=mars::fem::HypreGMRESSolver<double,int,cstone::GpuTag>;
    typename Solver::Matrix matrix;
    typename Solver::Vector b,x;
    thrust::device_vector<HYPRE_BigInt> mapping;
    Solver solver;
    LinearSystem(int n,int nnz):rows(C*n),blocks(nnz*C*C),rhs(rows),increment(rows),residual(2*rows),mapping(rows),solver(MPI_COMM_WORLD,2000,1e-12,Solver::BOOMERAMG,100) {
        matrix.allocate(rows,rows,nnz*C*C); b.resize(rows); x.resize(rows);
        thrust::sequence(mapping.begin(),mapping.end(),HYPRE_BigInt(0)); solver.setVerbose(false); solver.setPointBlock(C);
    }
    void solve(BlockCsrView<C> view) {
        launch(rows,SimpleScalarRows<C>{view,matrix.rowOffsetsPtr(),matrix.colIndicesPtr(),matrix.valuesPtr()});
        assembly_cuda_check(cudaMemcpy(b.data(),rhs.data(),rows*sizeof(double),cudaMemcpyDeviceToDevice));
        assembly_cuda_check(cudaMemset(x.data(),0,rows*sizeof(double)));
        require(solver.solve(matrix,b,x,0,rows,0,rows,mapping),"Hypre solve failed");
        assembly_cuda_check(cudaMemcpy(increment.data(),x.data(),rows*sizeof(double),cudaMemcpyDeviceToDevice));
        launch(rows,SimpleResidual{matrix.rowOffsetsPtr(),matrix.colIndicesPtr(),matrix.valuesPtr(),rhs.data(),increment.data(),residual.data()});
        // Only scalar reductions leave the device during the solve.
        auto r=thrust::make_permutation_iterator(residual.values.begin(),thrust::make_transform_iterator(thrust::make_counting_iterator(0),Even{}));
        double r2=thrust::reduce(r,r+rows,0.,thrust::plus<double>());
        auto b2it=thrust::make_permutation_iterator(residual.values.begin()+1,thrust::make_transform_iterator(thrust::make_counting_iterator(0),Even{}));
        double b2=thrust::reduce(b2it,b2it+rows,0.,thrust::plus<double>());
        check_residual(r2,b2); std::cout<<"linear components="<<C<<" iterations="<<solver.getLastIterations()<<'\n';
    }
    struct Even { __host__ __device__ int operator()(int i) const { return 2*i; } };
#else
    Array<int> offsets,columns; Array<double> values;
    LinearSystem(int n,int nnz):rows(C*n),blocks(nnz*C*C),rhs(rows),increment(rows),residual(2*rows),offsets(rows+1),columns(nnz*C*C),values(nnz*C*C) {}
    void solve(BlockCsrView<C> view) {
        launch(rows,SimpleScalarRows<C>{view,offsets.data(),columns.data(),values.data()});
        // Independent host oracle only; production calls Hypre on device CSR.
        std::vector<double> a(std::size_t(rows)*rows,0),b=rhs.values;
        for (int i=0;i<rows;++i) for (int k=offsets.values[i];k<offsets.values[i+1];++k) a[std::size_t(i)*rows+columns.values[k]]=values.values[k];
        for (int k=0;k<rows;++k) {
            int pivot=k;
            for (int i=k+1;i<rows;++i) if (std::abs(a[std::size_t(i)*rows+k])>std::abs(a[std::size_t(pivot)*rows+k])) pivot=i;
            require(std::isfinite(a[std::size_t(pivot)*rows+k]) && std::abs(a[std::size_t(pivot)*rows+k])>1e-20,"singular host system");
            if (pivot!=k) { for (int j=k;j<rows;++j) std::swap(a[std::size_t(k)*rows+j],a[std::size_t(pivot)*rows+j]); std::swap(b[k],b[pivot]); }
            for (int i=k+1;i<rows;++i) {
                const double ratio=a[std::size_t(i)*rows+k]/a[std::size_t(k)*rows+k];
                for (int j=k+1;j<rows;++j) a[std::size_t(i)*rows+j]-=ratio*a[std::size_t(k)*rows+j];
                b[i]-=ratio*b[k];
            }
        }
        for (int i=rows-1;i>=0;--i) { double v=b[i]; for (int j=i+1;j<rows;++j) v-=a[std::size_t(i)*rows+j]*increment.values[j]; increment.values[i]=v/a[std::size_t(i)*rows+i]; }
        launch(rows,SimpleResidual{offsets.data(),columns.data(),values.data(),rhs.data(),increment.data(),residual.data()});
        double r2=0,b2=0; for (int i=0;i<rows;++i) { r2+=residual.values[2*i]; b2+=residual.values[2*i+1]; } check_residual(r2,b2);
    }
#endif
    void check_residual(double r2,double b2) {
        const double absolute=std::sqrt(r2),relative=absolute/std::max(std::sqrt(b2),1e-300);
        std::cout<<"linear components="<<C<<" true_relative="<<relative<<" absolute="<<absolute<<'\n';
        require(std::isfinite(absolute) && std::isfinite(b2) && absolute<=1e-13+1e-10*std::sqrt(b2),"true linear residual failed");
    }
};
template<int C> void gradient(SimpleMesh m,SimpleState s,const double* field,Array<double>& sum,double* output) {
    sum.zero(); launch(m.element_count,SimpleGradientInterior<C>{m,field,sum.data()});
    launch(m.face_count,SimpleGradientBoundary<C>{m,field,sum.data()});
    launch(m.node_count*C*3,SimpleGradientFinish<C>{s.volume,sum.data(),output,s.error});
}
void execute(const Fixture& f) {
    const int n=int(f.x.size()),e=int(f.nodes[0].size()),b=int(f.faces.size());
    Array<double> x(f.x),y(f.y),z(f.z),velocity(3*n),pressure(n),vg(9*n),pg(3*n),d(3*n),volume(n),div(n),eflux(6*e),bflux(3*b),trace(3*b),factor(n),sum(9*n),gp(3*n),moment(2);
    Array<int> n0(f.nodes[0]),n1(f.nodes[1]),n2(f.nodes[2]),n3(f.nodes[3]),error(1);
    Array<SimpleFace> faces(f.faces); Array<TetGeometry<double>> geometry(e);
    SimpleMesh mesh{n,e,b,{n0.data(),n1.data(),n2.data(),n3.data()},x.data(),y.data(),z.data(),faces.data(),geometry.data()};
    SimpleState s{velocity.data(),pressure.data(),vg.data(),pg.data(),d.data(),volume.data(),div.data(),eflux.data(),bflux.data(),trace.data(),factor.data(),error.data()};
    SimpleControls controls;
    launch(e,SimpleGeometry{mesh,s}); require(error.host()[0]==0,"native geometry failed");
    launch(b,SimpleBoundaryFactor{mesh,factor.data()}); Graph graph(mesh);
    LinearSystem<3> momentum(n,graph.blocks()); LinearSystem<1> poisson(n,graph.blocks());
    for (int it=0;it<2;++it) {
        const auto& expected=f.expected[it]; std::string label="iteration="+std::to_string(it+1)+" ";
        gradient<3>(mesh,s,s.velocity,sum,s.velocity_gradient);
        gradient<1>(mesh,s,s.pressure,sum,s.pressure_gradient);
        momentum.blocks.zero(); momentum.rhs.zero(); auto am=graph.view<3>(momentum.blocks.data(),momentum.rhs.data());
        launch(e,SimpleInterior<3>{mesh,s,controls,am});
        launch(b,SimpleBoundary<3>{mesh,s,controls,am,it>0,false});
        launch(n,SimpleMomentumNode{s,controls,am}); require(error.host()[0]==0,"momentum assembly failed");
        momentum.solve(am); launch(3*n,SimpleAddIncrement{s.velocity,momentum.increment.data(),1});
        compare(velocity.host(),expected.momentum,label+"momentum"); compare(d.host(),expected.influence,label+"influence");
        moment.zero(); launch(b,SimpleTraceMoment{mesh,s,moment.data()}); launch(b,SimpleTrace{mesh,s,controls,moment.data()});
        compare(trace.host(),expected.trace,label+"trace");
        poisson.blocks.zero(); poisson.rhs.zero(); auto ap=graph.view<1>(poisson.blocks.data(),poisson.rhs.data());
        launch(e,SimpleInterior<1>{mesh,s,controls,ap}); launch(b,SimpleBoundary<1>{mesh,s,controls,ap});
        require(error.host()[0]==0,"pressure assembly failed"); poisson.solve(ap);
        compare(poisson.increment.host(),expected.increment,label+"raw_pressure_increment");
        launch(n,SimpleAddIncrement{s.pressure,poisson.increment.data(),controls.alpha_p});
        gradient<1>(mesh,s,poisson.increment.data(),sum,gp.data());
        // Fluxes use corrected p, uncorrected u, and the old reconstructed grad(p).
        div.zero(); launch(e,SimpleInterior<1>{mesh,s,controls,ap,true}); launch(b,SimpleBoundary<1>{mesh,s,controls,ap,it>0,true});
        require(error.host()[0]==0,"unsupported reversed outlet or flux update failure");
        launch(n,SimpleCorrectVelocity{s,gp.data()});
        compare(pressure.host(),expected.pressure,label+"pressure"); compare(velocity.host(),expected.velocity,label+"velocity");
        compare(eflux.host(),expected.interior_flux,label+"interior_flux"); compare(bflux.host(),expected.boundary_flux,label+"boundary_flux");
        compare(div.host(),expected.mass_divergence,label+"mass_divergence");
        auto hdiv=div.host(),hb=bflux.host(); double balance=std::accumulate(hb.begin(),hb.end(),0.),total=std::accumulate(hdiv.begin(),hdiv.end(),0.);
        require(std::abs(balance-total)<=1e-12,"global flux cancellation failed");
        std::cout<<label<<"net_boundary_mass="<<balance<<" sum_continuity="<<total<<'\n';
    }
}
int main(int argc,char** argv) {
#ifdef MARS_REPLAY_CUDA
    MPI_Init(&argc,&argv); int ranks=0; MPI_Comm_size(MPI_COMM_WORLD,&ranks);
#endif
    int result=0;
    try {
#ifdef MARS_REPLAY_CUDA
        require(ranks==1,"this SIMPLE gate supports exactly one rank");
#endif
        require(argc==2,"usage: mars_segregated_simple_check PUBLIC_SIMPLE_FILE"); execute(load(argv[1]));
#ifdef MARS_REPLAY_CUDA
        const char* backend="CUDA/Hypre";
#else
        const char* backend="host/direct";
#endif
        std::cout<<"PASS: "<<backend<<" two native SIMPLE iterations scalar_checks="<<checks
                 <<"; no injected stage fields; MPI, reversal selection and nonlinear convergence not tested\n";
    } catch (const std::exception& error) { std::cerr<<"FAIL: "<<error.what()<<'\n'; result=1; }
#ifdef MARS_REPLAY_CUDA
    MPI_Finalize();
#endif
    return result;
}
