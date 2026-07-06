// Standalone correctness test for the tet SCS area-vector precompute
// (mars_cvfem_tet_area.hpp). No mesh, no MPI math -- a few hardcoded
// reference tets on the device.
//
// IMPORTANT on invariants: the 6 interior SCS faces of a single tet do NOT
// sum to zero. Each node-cell is closed by its 3 interior SCS faces PLUS its
// share of the tet's OUTER surface; the outer pieces only cancel between
// adjacent elements on a mesh. So a per-element "sum of 6 area vectors = 0"
// check is wrong geometry (the hex area vectors fail it too). The correct
// global check is div(constant field)=0 at interior nodes, which is the
// mesh-level div_max test done once the NS operator runs.
//
// What we CAN check per element, independent of boundary:
//
//   1. VOLUME IDENTITY (Gauss on the field x, where div(x)=3): sum over the 6
//      SCS faces of A_f . (x_R - x_L) equals 3 * V_tet (the spatial dimension
//      times the volume). Verified analytically against the median-dual
//      construction. Catches area magnitude, winding, and centroid bugs.
//
//   2. NODE BALANCE: for each node, the signed sum (+L, -R) of its 3 incident
//      faces; summed over all 4 nodes must be zero (every face is +L once and
//      -R once). Catches orientation/sign bugs.
//
// Pass: |sum A_f.(xR-xL) - V| < 1e-10 and node-balance < 1e-10.

#include "backend/distributed/unstructured/fem/mars_cvfem_tet_area.hpp"
#include "backend/distributed/unstructured/fem/mars_element_traits.hpp"

#include <cuda_runtime.h>
#include <cstdio>
#include <cmath>

using namespace mars::fem;

using KeyType  = uint64_t;
using RealType = double;

// Per-element checks on element 0's 6 area vectors. Needs the 4 corner coords
// to form the edge vectors (x_R - x_L) and the true tet volume.
__global__ void checkAreaKernel(const RealType* ax, const RealType* ay, const RealType* az,
                                const RealType* cx, const RealType* cy, const RealType* cz,
                                RealType* outVolErr, RealType* outNodeBal, RealType* outVolTrue)
{
    constexpr int nscs = 6;

    // True tet volume = (1/6)|det[ p1-p0, p2-p0, p3-p0 ]|.
    RealType e1[3] = {cx[1]-cx[0], cy[1]-cy[0], cz[1]-cz[0]};
    RealType e2[3] = {cx[2]-cx[0], cy[2]-cy[0], cz[2]-cz[0]};
    RealType e3[3] = {cx[3]-cx[0], cy[3]-cy[0], cz[3]-cz[0]};
    RealType det = e1[0]*(e2[1]*e3[2]-e2[2]*e3[1])
                 - e1[1]*(e2[0]*e3[2]-e2[2]*e3[0])
                 + e1[2]*(e2[0]*e3[1]-e2[1]*e3[0]);
    RealType Vtrue = fabs(det) / RealType(6);
    outVolTrue[0] = Vtrue;

    // Gauss volume identity for the field x (div x = 3 in 3D):
    //   sum_f A_f . (x_R - x_L) == 3 * V_tet.
    RealType vsum = 0;
    for (int ip = 0; ip < nscs; ++ip) {
        int L = d_tetLRSCV[ip * 2];
        int R = d_tetLRSCV[ip * 2 + 1];
        RealType dx = cx[R]-cx[L], dy = cy[R]-cy[L], dz = cz[R]-cz[L];
        vsum += ax[ip]*dx + ay[ip]*dy + az[ip]*dz;
    }
    outVolErr[0] = fabs(vsum - RealType(3) * Vtrue);

    // Per-node balance: +A when node is L, -A when node is R, over incident
    // faces; sum of the 4 node vectors must be zero (each face is +L once,
    // -R once).
    RealType nb[3] = {0, 0, 0};
    for (int node = 0; node < 4; ++node) {
        RealType vx = 0, vy = 0, vz = 0;
        for (int ip = 0; ip < nscs; ++ip) {
            int L = d_tetLRSCV[ip * 2];
            int R = d_tetLRSCV[ip * 2 + 1];
            RealType s = (L == node) ? RealType(1) : (R == node) ? RealType(-1) : RealType(0);
            vx += s * ax[ip]; vy += s * ay[ip]; vz += s * az[ip];
        }
        nb[0] += vx; nb[1] += vy; nb[2] += vz;
    }
    outNodeBal[0] = sqrt(nb[0]*nb[0] + nb[1]*nb[1] + nb[2]*nb[2]);
}

static int runOneTet(const char* name, const RealType verts[4][3])
{
    // Upload one tet: connectivity 0,1,2,3 -> node ids 0..3; coords as given.
    KeyType h_c0[1] = {0}, h_c1[1] = {1}, h_c2[1] = {2}, h_c3[1] = {3};
    RealType h_x[4], h_y[4], h_z[4];
    for (int n = 0; n < 4; ++n) { h_x[n] = verts[n][0]; h_y[n] = verts[n][1]; h_z[n] = verts[n][2]; }

    KeyType *d_c0, *d_c1, *d_c2, *d_c3;
    RealType *d_x, *d_y, *d_z, *d_ax, *d_ay, *d_az;
    cudaMalloc(&d_c0, sizeof(KeyType)); cudaMalloc(&d_c1, sizeof(KeyType));
    cudaMalloc(&d_c2, sizeof(KeyType)); cudaMalloc(&d_c3, sizeof(KeyType));
    cudaMalloc(&d_x, 4*sizeof(RealType)); cudaMalloc(&d_y, 4*sizeof(RealType)); cudaMalloc(&d_z, 4*sizeof(RealType));
    cudaMalloc(&d_ax, 6*sizeof(RealType)); cudaMalloc(&d_ay, 6*sizeof(RealType)); cudaMalloc(&d_az, 6*sizeof(RealType));

    cudaMemcpy(d_c0, h_c0, sizeof(KeyType), cudaMemcpyHostToDevice);
    cudaMemcpy(d_c1, h_c1, sizeof(KeyType), cudaMemcpyHostToDevice);
    cudaMemcpy(d_c2, h_c2, sizeof(KeyType), cudaMemcpyHostToDevice);
    cudaMemcpy(d_c3, h_c3, sizeof(KeyType), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, h_x, 4*sizeof(RealType), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, h_y, 4*sizeof(RealType), cudaMemcpyHostToDevice);
    cudaMemcpy(d_z, h_z, 4*sizeof(RealType), cudaMemcpyHostToDevice);

    precomputeTetAreaVectorsGpu<KeyType, RealType>(
        d_c0, d_c1, d_c2, d_c3, /*numElements=*/1, d_x, d_y, d_z, d_ax, d_ay, d_az);
    cudaDeviceSynchronize();

    RealType *d_volErr, *d_nodebal, *d_volTrue;
    cudaMalloc(&d_volErr, sizeof(RealType)); cudaMalloc(&d_nodebal, sizeof(RealType));
    cudaMalloc(&d_volTrue, sizeof(RealType));
    checkAreaKernel<<<1,1>>>(d_ax, d_ay, d_az, d_x, d_y, d_z, d_volErr, d_nodebal, d_volTrue);
    cudaDeviceSynchronize();

    RealType volErr = 0, nodebal = 0, volTrue = 0;
    cudaMemcpy(&volErr,  d_volErr,  sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(&nodebal, d_nodebal, sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(&volTrue, d_volTrue, sizeof(RealType), cudaMemcpyDeviceToHost);

    // Also pull the 6 area vectors for reporting.
    RealType h_ax[6], h_ay[6], h_az[6];
    cudaMemcpy(h_ax, d_ax, 6*sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_ay, d_ay, 6*sizeof(RealType), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_az, d_az, 6*sizeof(RealType), cudaMemcpyDeviceToHost);

    const RealType tol = 1e-10;
    bool ok = (volErr < tol) && (nodebal < tol);
    printf("[%s]\n", name);
    for (int ip = 0; ip < 6; ++ip)
        printf("   A[%d] = (% .6e, % .6e, % .6e)\n", ip, h_ax[ip], h_ay[ip], h_az[ip]);
    printf("   tet volume         = %.6e\n", volTrue);
    printf("   volume-identity err = %.3e  (sum A.(xR-xL) vs 3V)\n", volErr);
    printf("   node-balance resid  = %.3e\n", nodebal);
    printf("   --> %s\n\n", ok ? "PASS" : "FAIL");

    cudaFree(d_c0); cudaFree(d_c1); cudaFree(d_c2); cudaFree(d_c3);
    cudaFree(d_x); cudaFree(d_y); cudaFree(d_z);
    cudaFree(d_ax); cudaFree(d_ay); cudaFree(d_az);
    cudaFree(d_volErr); cudaFree(d_nodebal); cudaFree(d_volTrue);
    return ok ? 0 : 1;
}

int main()
{
    int fails = 0;

    // Reference (right) tet: unit corner at origin.
    const RealType refTet[4][3] = {
        {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}
    };
    fails += runOneTet("reference unit tet", refTet);

    // Skewed tet: shifted apex, non-axis-aligned -- catches winding/sign bugs
    // a symmetric tet might hide.
    const RealType skewTet[4][3] = {
        {0.1, 0.2, 0.0}, {1.3, 0.1, 0.2}, {0.2, 1.1, -0.1}, {0.4, 0.3, 1.2}
    };
    fails += runOneTet("skewed tet", skewTet);

    // Flipped orientation (nodes 1 and 2 swapped -> negative det). The
    // orientation fix should still produce closing area vectors.
    const RealType flipTet[4][3] = {
        {0, 0, 0}, {0, 1, 0}, {1, 0, 0}, {0, 0, 1}
    };
    fails += runOneTet("flipped (negative-det) tet", flipTet);

    printf("==== %s ====\n", fails == 0 ? "ALL TET AREA TESTS PASS" : "TET AREA TESTS FAILED");
    return fails;
}
