// Stage 1c: topological HO DOF handler + CG apply on a real ElementDomain mesh
// (single rank). Same handler as Stage 1b, now fed the SFC-ordered connectivity
// and real node numbering from ElementDomain instead of a hand-built lattice.
// Validates the handler against the actual domain pipeline.
//
// Geometry assumption (axis-aligned cube elements): per-element diagonal metric
// from edge lengths. General curved/unstructured hexes need the full Jacobian +
// face-orientation table (next stages). Use a cube mesh here.
//
// Gates: nullspace A*1=0, symmetry, continuity A*(physical x)=0 at interior.

#include "mars.hpp"
#include "backend/distributed/unstructured/domain.hpp"
#include "backend/distributed/unstructured/fem/mars_spectral_basis.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_dof_handler.hpp"
#include "backend/distributed/unstructured/fem/mars_ho_laplacian.hpp"

#include <cuda_runtime.h>
#include <mpi.h>
#include <array>
#include <cstdio>
#include <cmath>
#include <random>
#include <string>
#include <vector>

using namespace mars;
using namespace mars::fem;

#define CK(call) do { cudaError_t e=(call); if(e!=cudaSuccess){ \
  printf("CUDA error %s at %s:%d\n", cudaGetErrorString(e), __FILE__, __LINE__); return false; } } while(0)

template<int P>
bool runDomain(const std::string& meshFile)
{
    constexpr int n  = P + 1;
    constexpr int N3 = n * n * n;
    using KeyType = uint64_t;
    using Real    = double;

    ElementDomain<HexTag, Real, KeyType, cstone::GpuTag> domain(meshFile, 0, 1, true, 64, 8);
    (void)domain.getNodeOwnershipMap();
    const auto& d_conn = domain.getElementToNodeConnectivity();
    domain.cacheNodeCoordinates();
    const auto& d_x = domain.getNodeX();
    const auto& d_y = domain.getNodeY();
    const auto& d_z = domain.getNodeZ();

    size_t nE = std::get<0>(d_conn).size();
    size_t nN = d_x.size();

    // Pull connectivity + coords to host.
    std::vector<KeyType> c[8];
    for (int t = 0; t < 8; ++t) {
        c[t].resize(nE);
        const KeyType* src = (t==0?std::get<0>(d_conn):t==1?std::get<1>(d_conn):t==2?std::get<2>(d_conn):
                              t==3?std::get<3>(d_conn):t==4?std::get<4>(d_conn):t==5?std::get<5>(d_conn):
                              t==6?std::get<6>(d_conn):std::get<7>(d_conn)).data();
        CK(cudaMemcpy(c[t].data(), src, nE*sizeof(KeyType), cudaMemcpyDeviceToHost));
    }
    std::vector<Real> X(nN), Y(nN), Z(nN);
    CK(cudaMemcpy(X.data(), d_x.data(), nN*sizeof(Real), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(Y.data(), d_y.data(), nN*sizeof(Real), cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(Z.data(), d_z.data(), nN*sizeof(Real), cudaMemcpyDeviceToHost));

    std::vector<std::array<int,8>> ec(nE);
    for (size_t e = 0; e < nE; ++e)
        for (int t = 0; t < 8; ++t) ec[e][t] = (int)c[t][e];

    HODofHandler dh;
    dh.build(ec, (long)nN, P);
    long nDof = dh.numDof;

    std::vector<double> gx, gw, gD;
    gllBasis(P, gx, gw, gD);
    std::vector<Real> Dr(gD.begin(), gD.end());

    // Domain extent for boundary classification.
    Real xmn=X[0],xmx=X[0],ymn=Y[0],ymx=Y[0],zmn=Z[0],zmx=Z[0];
    for (size_t i=0;i<nN;++i){ xmn=std::min(xmn,X[i]);xmx=std::max(xmx,X[i]); ymn=std::min(ymn,Y[i]);ymx=std::max(ymx,Y[i]); zmn=std::min(zmn,Z[i]);zmx=std::max(zmx,Z[i]); }
    Real ext = std::max({xmx-xmn, ymx-ymn, zmx-zmn});
    Real eps = 1e-9 * ext;

    // Per-element diagonal metric (axis-aligned) + per-DOF physical x + boundary flag.
    std::vector<Real>    G((size_t)nE * N3 * 6, 0.0);
    std::vector<Real>    dofX(nDof, 0.0);
    std::vector<uint8_t> isB(nDof, 0);
    for (size_t e = 0; e < nE; ++e)
    {
        int c0=ec[e][0], c1=ec[e][1], c3=ec[e][3], c4=ec[e][4];
        Real hx = std::abs(X[c1]-X[c0]), hy = std::abs(Y[c3]-Y[c0]), hz = std::abs(Z[c4]-Z[c0]);
        Real x0 = X[c0], y0 = Y[c0], z0 = Z[c0];
        for (int i=0;i<n;++i) for (int j=0;j<n;++j) for (int k=0;k<n;++k)
        {
            int l = i*n*n + j*n + k;
            Real wp = gw[i]*gw[j]*gw[k];
            size_t gp = (e*N3 + l)*6;
            G[gp+0] = (hy*hz)/(2*hx) * wp;
            G[gp+3] = (hx*hz)/(2*hy) * wp;
            G[gp+5] = (hx*hy)/(2*hz) * wp;
            Real px = x0 + 0.5*(gx[i]+1.0)*hx;
            Real py = y0 + 0.5*(gx[j]+1.0)*hy;
            Real pz = z0 + 0.5*(gx[k]+1.0)*hz;
            int dof = dh.elemDof[e*N3 + l];
            dofX[dof] = px;
            bool b = (px<xmn+eps||px>xmx-eps||py<ymn+eps||py>ymx-eps||pz<zmn+eps||pz>zmx-eps);
            isB[dof] = b ? 1 : 0;
        }
    }

    int *dElemDof; Real *dD, *dG, *dU, *dY;
    CK(cudaMalloc(&dElemDof, dh.elemDof.size()*sizeof(int)));
    CK(cudaMalloc(&dD, Dr.size()*sizeof(Real)));
    CK(cudaMalloc(&dG, G.size()*sizeof(Real)));
    CK(cudaMalloc(&dU, nDof*sizeof(Real)));
    CK(cudaMalloc(&dY, nDof*sizeof(Real)));
    CK(cudaMemcpy(dElemDof, dh.elemDof.data(), dh.elemDof.size()*sizeof(int), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dD, Dr.data(), Dr.size()*sizeof(Real), cudaMemcpyHostToDevice));
    CK(cudaMemcpy(dG, G.data(), G.size()*sizeof(Real), cudaMemcpyHostToDevice));

    int block = 128, grid = int((nE + block - 1)/block);
    auto apply = [&](const std::vector<Real>& u, std::vector<Real>& y) -> bool {
        CK(cudaMemcpy(dU, u.data(), nDof*sizeof(Real), cudaMemcpyHostToDevice));
        CK(cudaMemset(dY, 0, nDof*sizeof(Real)));
        ho_laplacian_apply_cg<Real, P><<<grid, block>>>(dU, dY, dElemDof, dD, dG, nE);
        CK(cudaDeviceSynchronize());
        y.resize(nDof);
        CK(cudaMemcpy(y.data(), dY, nDof*sizeof(Real), cudaMemcpyDeviceToHost));
        return true;
    };

    std::vector<Real> ones(nDof, 1.0), y;
    if (!apply(ones, y)) return false;
    Real nullMax = 0; for (Real v : y) nullMax = std::max(nullMax, std::abs(v));

    std::mt19937 rng(11); std::uniform_real_distribution<Real> uni(-1,1);
    std::vector<Real> uu(nDof), vv(nDof), Au, Av;
    for (auto& z : uu) z = uni(rng);
    for (auto& z : vv) z = uni(rng);
    if (!apply(uu, Au)) return false;
    if (!apply(vv, Av)) return false;
    Real uAv=0, vAu=0; for (long p=0;p<nDof;++p){ uAv+=uu[p]*Av[p]; vAu+=vv[p]*Au[p]; }
    Real symErr = std::abs(uAv-vAu)/(std::abs(uAv)+1.0);

    std::vector<Real> Ax;
    if (!apply(dofX, Ax)) return false;
    Real interMax = 0;
    for (long p=0;p<nDof;++p) if (!isB[p]) interMax = std::max(interMax, std::abs(Ax[p]));

    bool ok = (nullMax<1e-9) && (symErr<1e-10) && (interMax<1e-8*ext);
    printf("P=%d | elems=%zu nodes=%zu | nDof=%ld nEdge=%ld nFace=%ld | null=%.2e | sym=%.2e | continuity=%.2e | %s\n",
           P, nE, nN, nDof, dh.nEdge, dh.nFace, nullMax, symErr, interMax, ok?"PASS":"FAIL");

    cudaFree(dElemDof); cudaFree(dD); cudaFree(dG); cudaFree(dU); cudaFree(dY);
    return ok;
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int rank=0, nr=1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nr);
    int dev=0; cudaGetDeviceCount(&dev); if (dev>0) cudaSetDevice(rank % dev);

    std::string mesh; int order = 4;
    for (int i=1;i<argc;++i){ std::string a=argv[i];
        if (a.find("--mesh=")==0) mesh=a.substr(7);
        else if (a.find("--order=")==0) order=std::stoi(a.substr(8)); }
    if (mesh.empty()) { if(rank==0) printf("need --mesh=<cube dir> [--order=2|4|7]\n"); MPI_Finalize(); return 1; }

    bool ok=false;
    if      (order==2) ok=runDomain<2>(mesh);
    else if (order==4) ok=runDomain<4>(mesh);
    else if (order==7) ok=runDomain<7>(mesh);
    else { if(rank==0) printf("order must be 2,4,7\n"); MPI_Finalize(); return 1; }

    if (rank==0) printf("Stage 1c (HO DOF on ElementDomain): %s\n", ok?"PASS":"FAIL");
    MPI_Finalize();
    return ok?0:1;
}
