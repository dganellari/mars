// GPU + distributed gate for the high-order tet CVFEM (box-partition, collapsed sum-factorization).
// Sorted Kuhn mesh, merged collapsed DOF (C0), elements partitioned across ranks, device apply,
// MPI_Allreduce assembly. Gates: assembled A.1=0, A.(x,y,z)=0 at domain-interior DOF (EXACT -- the
// scheme is quadrature-exact), and parity with the host sum-fac reference.
//
//   srun -N1 --ntasks-per-node=4 ./examples/distributed/unstructured/mars_cvfem_ho_tet_test --ncells=4 --p=3
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_tet.hpp"
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
    constexpr int M = (P + 3) / 2;
    auto o = buildHoCvfemTetOps<double>(P);
    if (o.m != M) { if (rank==0) std::printf("m mismatch\n"); return false; }
    std::vector<std::array<double,3>> coords;
    std::vector<int> elemCorners;
    buildKuhnTetMesh(ncells, coords, elemCorners);
    const size_t nElem = elemCorners.size()/4;
    std::vector<double> Zd(o.Z.begin(), o.Z.end());
    HoCvfemTetDofHandler dh; dh.build(elemCorners, coords, Zd);
    const int NN = o.NN; const long nDof = dh.numDof;

    const size_t e0 = nElem*rank/nranks, e1 = nElem*(rank+1)/nranks, nOwn = e1-e0;

    std::vector<double> geom(nOwn*10);
    for (size_t e=0;e<nOwn;++e) {
        double corners[4][3];
        for (int c=0;c<4;++c) for (int d=0;d<3;++d) corners[c][d]=coords[dh.sortedCorners[(e0+e)*4+c]][d];
        double MMt[9], detB; computeTetGeom<double>(corners, MMt, detB);
        for (int i=0;i<9;++i) geom[e*10+i]=MMt[i]; geom[e*10+9]=detB;
    }
    std::vector<int> eDof(nOwn*NN);
    std::memcpy(eDof.data(), &dh.elemDof[e0*NN], eDof.size()*sizeof(int));

    std::vector<double> hxs(o.xs.begin(),o.xs.end()), hyq(o.yq.begin(),o.yq.end()), hwq(o.wqv.begin(),o.wqv.end()),
                        hLs(o.Ls.begin(),o.Ls.end()), hDs(o.Ds.begin(),o.Ds.end()),
                        hLq(o.Lq.begin(),o.Lq.end()), hDq(o.Dq.begin(),o.Dq.end());
    double *d_xs=up(hxs),*d_yq=up(hyq),*d_wq=up(hwq),*d_Ls=up(hLs),*d_Ds=up(hDs),*d_Lq=up(hLq),*d_Dq=up(hDq),
           *d_geom=up(geom);
    int *d_eDof=up(eDof);
    double *d_u,*d_y; CK(cudaMalloc(&d_u,nDof*sizeof(double))); CK(cudaMalloc(&d_y,nDof*sizeof(double)));

    const int bs=64, nb=(int)((nOwn+bs-1)/bs);
    auto applyDist=[&](const std::vector<double>& u, std::vector<double>& y){
        CK(cudaMemcpy(d_u,u.data(),nDof*sizeof(double),cudaMemcpyHostToDevice));
        CK(cudaMemset(d_y,0,nDof*sizeof(double)));
        if (nOwn) ho_cvfem_tet_apply_cg<double,P,M><<<nb,bs>>>(d_u,d_y,d_eDof,
                      d_xs,d_yq,d_wq,d_Ls,d_Ds,d_Lq,d_Dq,d_geom,nOwn);
        CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
        CK(cudaMemcpy(y.data(),d_y,nDof*sizeof(double),cudaMemcpyDeviceToHost));
        MPI_Allreduce(MPI_IN_PLACE,y.data(),(int)nDof,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    };

    // A.1 = 0
    std::vector<double> one(nDof,1.0), y(nDof);
    applyDist(one,y);
    double a1=0; for(long i=0;i<nDof;++i) a1=std::max(a1,std::fabs(y[i]));

    // A.(x,y,z)=0 at domain-interior DOF (exact consistency)
    const double L=ncells, eps=1e-9;
    double lin=0;
    for (int d=0;d<3;++d) {
        std::vector<double> ul(nDof);
        for(long i=0;i<nDof;++i) ul[i]=dh.dofPos[i][d];
        applyDist(ul,y);
        for(long i=0;i<nDof;++i){ const auto&x=dh.dofPos[i];
            if(x[0]>eps&&x[0]<L-eps&&x[1]>eps&&x[1]<L-eps&&x[2]>eps&&x[2]<L-eps)
                lin=std::max(lin,std::fabs(y[i])); }
    }

    // parity: GPU distributed apply == host sum-fac assembled reference
    std::vector<double> u(nDof);
    for(long i=0;i<nDof;++i) u[i]=std::sin(0.9*i+0.4);
    applyDist(u,y);
    std::vector<double> yref(nDof,0.0), ue(NN), ye(NN);
    for (size_t e=0;e<nElem;++e) {
        double corners[4][3];
        for (int c=0;c<4;++c) for (int d=0;d<3;++d) corners[c][d]=coords[dh.sortedCorners[e*4+c]][d];
        double MMt[9], detB; computeTetGeom<double>(corners, MMt, detB);
        for (int i=0;i<NN;++i) ue[i]=u[dh.elemDof[e*NN+i]];
        applyHoCvfemTetElement<double>(o, MMt, detB, ue.data(), ye.data());
        for (int i=0;i<NN;++i) yref[dh.elemDof[e*NN+i]] += ye[i];
    }
    double par=0; for(long i=0;i<nDof;++i) par=std::max(par,std::fabs(y[i]-yref[i]));

    cudaFree(d_xs);cudaFree(d_yq);cudaFree(d_wq);cudaFree(d_Ls);cudaFree(d_Ds);cudaFree(d_Lq);cudaFree(d_Dq);
    cudaFree(d_geom);cudaFree(d_eDof);cudaFree(d_u);cudaFree(d_y);

    const bool ok = (a1<1e-9) && (lin<1e-9) && (par<1e-9);
    if (rank==0)
        std::printf(" p=%d | elems=%zu | nDof=%ld | ranks=%d | A.1=%.2e | A.lin_int=%.2e | gpu-host=%.2e | %s\n",
                    P, nElem, nDof, nranks, a1, lin, par, ok?"PASS":"FAIL");
    return ok;
}

int main(int argc, char** argv) {
    MPI_Init(&argc,&argv);
    int rank,nranks; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&nranks);
    int ndev=0; cudaGetDeviceCount(&ndev);
    if (ndev) CK(cudaSetDevice(rank%ndev));   // tet drivers must bind the device before first CUDA call

    int ncells=2, p=3;
    for (int i=1;i<argc;++i){
        if (!std::strncmp(argv[i],"--ncells=",9)) ncells=std::atoi(argv[i]+9);
        if (!std::strncmp(argv[i],"--p=",4))      p=std::atoi(argv[i]+4);
    }
    if (rank==0) std::printf("=== HO tet CVFEM gate (box-partition, quadrature-exact): ncells=%d p=%d ranks=%d ===\n",
                             ncells,p,nranks);
    bool ok=false;
    switch(p){
        case 2: ok=run_order<2>(ncells,rank,nranks); break;
        case 3: ok=run_order<3>(ncells,rank,nranks); break;
        case 4: ok=run_order<4>(ncells,rank,nranks); break;
        case 5: ok=run_order<5>(ncells,rank,nranks); break;
        default: if(rank==0) std::printf("unsupported --p=%d (2..5)\n",p);
    }
    if (rank==0) std::printf(ok ? "==== TET HO CVFEM GATE PASS ====\n" : "==== FAILURE ====\n");
    MPI_Finalize();
    return ok?0:1;
}
