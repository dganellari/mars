// Distributed (multi-rank, multi-GPU) gate for the high-order tet Laplacian: Kuhn tet mesh,
// C0 nodal DOF (barycentric-key handler), elements partitioned across ranks, device CG apply per
// rank, MPI_Allreduce assembly. Gates: global A.1=0, A.(x,y,z)=0 at interior DOF, and parity with
// the host-assembled reference. Tet counterpart of mars_ho_dist_apply_test.
//
// Replicated numbering + Allreduce assembly is the correctness stage; the targeted HoHalo
// exchange (ghost-DOF only) is the perf follow-up once this gate is green.
//
//   srun -N2 --ntasks-per-node=4 ./examples/distributed/unstructured/mars_ho_dist_apply_tet_test --ncells=8 --p=3
#include "backend/distributed/unstructured/fem/mars_ho_laplacian_tet.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler_tet.hpp"
#include <mpi.h>
#include <cuda_runtime.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <cmath>
using namespace mars::fem;

#define CK(x) do{ cudaError_t e_=(x); if(e_!=cudaSuccess){ \
    std::printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e_), __FILE__, __LINE__); MPI_Abort(MPI_COMM_WORLD,1);} }while(0)

template <typename T> static T* up(const std::vector<T>& h){ T* d; CK(cudaMalloc(&d,h.size()*sizeof(T)));
    CK(cudaMemcpy(d,h.data(),h.size()*sizeof(T),cudaMemcpyHostToDevice)); return d; }

template <int P>
static bool run_order(int ncells, int rank, int nranks) {
    auto o  = buildHoTetOps<double>(P);
    auto nd = buildHoTetNodal<double>(o);
    std::vector<std::array<double,3>> coords;
    std::vector<int> elemCorners;
    buildKuhnTetMesh(ncells, coords, elemCorners);
    const size_t nElem = elemCorners.size()/4;
    HoTetDofHandler dh; dh.build(elemCorners, coords, nd);
    const int Np = nd.Np; const long nDof = dh.numDof;

    // contiguous element partition
    const size_t e0 = nElem*rank/nranks, e1 = nElem*(rank+1)/nranks, nOwn = e1-e0;

    // owned metric + dof tables
    std::vector<double> G((size_t)nOwn*o.nq*6);
    for (size_t e=0;e<nOwn;++e) {
        double corners[4][3];
        for (int c=0;c<4;++c) for (int d=0;d<3;++d) corners[c][d]=coords[elemCorners[(e0+e)*4+c]][d];
        auto g = computeHoTetMetric<double>(corners, o);
        std::memcpy(&G[(size_t)e*o.nq*6], g.data(), g.size()*sizeof(double));
    }
    std::vector<int> eDof((size_t)nOwn*Np);
    std::memcpy(eDof.data(), &dh.elemDof[e0*Np], eDof.size()*sizeof(int));
    std::vector<int> modes((size_t)o.nModes*3);
    for (int m=0;m<o.nModes;++m){ modes[m*3]=o.modes[m][0]; modes[m*3+1]=o.modes[m][1]; modes[m*3+2]=o.modes[m][2]; }

    double *d_A=up(o.A),*d_Ad=up(o.Ad),*d_B=up(o.B),*d_Bd=up(o.Bd),*d_C=up(o.C),*d_Cd=up(o.Cd),
           *d_G=up(G),*d_Vi=up(nd.Vinv);
    int *d_eDof=up(eDof),*d_modes=up(modes);
    double *d_u,*d_y; CK(cudaMalloc(&d_u,nDof*sizeof(double))); CK(cudaMalloc(&d_y,nDof*sizeof(double)));

    const int bs=64, nb=(int)((nOwn+bs-1)/bs);
    auto applyDist=[&](const std::vector<double>& u, std::vector<double>& y){
        CK(cudaMemcpy(d_u,u.data(),nDof*sizeof(double),cudaMemcpyHostToDevice));
        CK(cudaMemset(d_y,0,nDof*sizeof(double)));
        if (nOwn) ho_laplacian_tet_apply_cg<double,P><<<nb,bs>>>(d_u,d_y,d_eDof,d_Vi,d_modes,
                      d_A,d_Ad,d_B,d_Bd,d_C,d_Cd,d_G,nOwn);
        CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
        CK(cudaMemcpy(y.data(),d_y,nDof*sizeof(double),cudaMemcpyDeviceToHost));
        MPI_Allreduce(MPI_IN_PLACE,y.data(),(int)nDof,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    };

    // gate B: A.1 = 0
    std::vector<double> one(nDof,1.0), y(nDof);
    applyDist(one,y);
    double a1=0; for(long i=0;i<nDof;++i) a1=std::max(a1,std::fabs(y[i]));

    // gate C: A.(x,y,z)=0 at interior DOF
    const double L=ncells, eps=1e-9;
    double cons=0;
    for (int d=0;d<3;++d) {
        std::vector<double> ul(nDof);
        for(long i=0;i<nDof;++i) ul[i]=dh.dofPos[i][d];
        applyDist(ul,y);
        for(long i=0;i<nDof;++i){ const auto&x=dh.dofPos[i];
            if(x[0]>eps&&x[0]<L-eps&&x[1]>eps&&x[1]<L-eps&&x[2]>eps&&x[2]<L-eps)
                cons=std::max(cons,std::fabs(y[i])); }
    }

    // gate A: GPU distributed apply == host-assembled reference on a random field
    std::vector<double> u(nDof);
    for(long i=0;i<nDof;++i) u[i]=std::sin(0.9*i+0.4);
    applyDist(u,y);
    std::vector<double> yref(nDof,0.0), ue(Np), ye(Np);
    for (size_t e=0;e<nElem;++e) {
        double corners[4][3];
        for (int c=0;c<4;++c) for (int d=0;d<3;++d) corners[c][d]=coords[elemCorners[e*4+c]][d];
        auto g = computeHoTetMetric<double>(corners, o);
        for (int i=0;i<Np;++i) ue[i]=u[dh.elemDof[e*Np+i]];
        applyHoLaplacianTetElementNodal<double>(o, nd, g.data(), ue.data(), ye.data());
        for (int i=0;i<Np;++i) yref[dh.elemDof[e*Np+i]] += ye[i];
    }
    double par=0; for(long i=0;i<nDof;++i) par=std::max(par,std::fabs(y[i]-yref[i]));

    cudaFree(d_A);cudaFree(d_Ad);cudaFree(d_B);cudaFree(d_Bd);cudaFree(d_C);cudaFree(d_Cd);
    cudaFree(d_G);cudaFree(d_Vi);cudaFree(d_eDof);cudaFree(d_modes);cudaFree(d_u);cudaFree(d_y);

    const bool ok = (a1<1e-9) && (cons<1e-8) && (par<1e-9);
    if (rank==0)
        std::printf(" p=%d | elems=%zu | nDof=%ld | ranks=%d | A.1=%.2e | A.x_int=%.2e | gpu-host=%.2e | %s\n",
                    P, nElem, nDof, nranks, a1, cons, par, ok?"PASS":"FAIL");
    return ok;
}

int main(int argc, char** argv) {
    MPI_Init(&argc,&argv);
    int rank,nranks; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&nranks);
    int ndev=0; cudaGetDeviceCount(&ndev);
    if (ndev) CK(cudaSetDevice(rank%ndev));   // tet drivers must bind the device before first CUDA call

    int ncells=4, p=3;
    for (int i=1;i<argc;++i){
        if (!std::strncmp(argv[i],"--ncells=",9)) ncells=std::atoi(argv[i]+9);
        if (!std::strncmp(argv[i],"--p=",4))      p=std::atoi(argv[i]+4);
    }
    if (rank==0) std::printf("=== HO tet distributed apply gate: ncells=%d p=%d ranks=%d ===\n",ncells,p,nranks);

    bool ok=false;
    switch(p){
        case 2: ok=run_order<2>(ncells,rank,nranks); break;
        case 3: ok=run_order<3>(ncells,rank,nranks); break;
        case 4: ok=run_order<4>(ncells,rank,nranks); break;
        case 5: ok=run_order<5>(ncells,rank,nranks); break;
        default: if(rank==0) std::printf("unsupported --p=%d (2..5)\n",p);
    }
    if (rank==0) std::printf(ok ? "==== DISTRIBUTED TET HO GATE PASS ====\n" : "==== FAILURE ====\n");
    MPI_Finalize();
    return ok?0:1;
}
