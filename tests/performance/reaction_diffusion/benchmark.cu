// Element-local FP64 comparison. Each arm starts with its own prepacked coefficient arrays.
#include "backend/distributed/unstructured/fem/mars_cvfem_ho_basis.hpp"
#include "backend/distributed/unstructured/marsir/mars_marsir_ptx_galerkin.hpp"
#include <thrust/device_vector.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Vec = std::vector<double>;
using Device = thrust::device_vector<double>;
constexpr int p = 3, n = p + 1, q = p + 2, nodes = n*n*n, points = q*q*q;

double* ptr(Device& v) { return thrust::raw_pointer_cast(v.data()); }
void check(cudaError_t code) {
    if (code != cudaSuccess) throw std::runtime_error(cudaGetErrorString(code));
}

__global__ void add(const double* a, const double* b, double* y, int size) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < size) y[i] = a[i] + b[i];
}

double field(double x, double y, double z) { return 1 + x + x*y + 0.2*z*z; }

struct Inputs {
    Vec z, xq, wq, B, Bt, dB, D, Dt, u, metric, reaction, packed;
    explicit Inputs(int elements)
      : B(q*n), Bt(n*q), dB(q*n), D(q*q), Dt(q*q), u(elements*nodes),
        metric(elements*6*points), reaction(elements*points), packed(elements*7*points) {
        Vec weights, nodal_deriv, v, d;
        mars::fem::gllBasis(p, z, weights, nodal_deriv);
        mars::fem::gaussLegendre(q, xq, wq);
        for (int i = 0; i < q; ++i) {
            mars::fem::lagrangeValDeriv(z, xq[i], v, d);
            for (int j = 0; j < n; ++j) {
                B[i*n+j] = Bt[j*q+i] = v[j];
                dB[i*n+j] = d[j];
            }
            mars::fem::lagrangeValDeriv(xq, xq[i], v, d);
            for (int j = 0; j < q; ++j) D[i*q+j] = Dt[j*q+i] = d[j];
        }
        for (int e = 0; e < elements; ++e) {
            const double a = warp(e);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    for (int k = 0; k < n; ++k)
                        u[e*nodes+(i*n+j)*n+k] = field(z[i]+a*z[j]*z[k], z[j], z[k]);
            for (int i = 0; i < q; ++i)
                for (int j = 0; j < q; ++j)
                    for (int k = 0; k < q; ++k) {
                        int t = (i*q+j)*q+k;
                        double y = xq[j], zc = xq[k], w = wq[i]*wq[j]*wq[k];
                        // x = xi + a*eta*zeta, y = eta, z = zeta. det(J) = 1.
                        double g[] = {1+a*a*(y*y+zc*zc), -a*zc, -a*y, 1, 0, 1};
                        for (int c = 0; c < 6; ++c)
                            packed[(e*7+c)*points+t] = metric[(e*6+c)*points+t] = w*g[c];
                        packed[(e*7+6)*points+t] = reaction[e*points+t] = w*sigma(xq[i]+a*y*zc, y);
                    }
        }
    }

    static double warp(int e) { return 0.06 + 0.012*(e % 17); }
    static double sigma(double x, double y) { return 1.2 + 0.3*x + 0.1*y; }

    // Direct weak form from the analytic physical gradient of the prescribed field.
    double expected(int e, int node) const {
        int i = node/(n*n), j = node/n%n, k = node%n;
        double total = 0, a = warp(e);
        for (int l = 0; l < q; ++l)
            for (int m = 0; m < q; ++m)
                for (int h = 0; h < q; ++h) {
                    double y = xq[m], zc = xq[h], x = xq[l]+a*y*zc;
                    double bx = B[l*n+i], by = B[m*n+j], bz = B[h*n+k];
                    double dx = dB[l*n+i]*by*bz;
                    double dy = bx*dB[m*n+j]*bz-a*zc*dx;
                    double dz = bx*by*dB[h*n+k]-a*y*dx;
                    total += wq[l]*wq[m]*wq[h]*(dx*(1+y)+dy*x+dz*0.4*zc
                             + bx*by*bz*sigma(x,y)*field(x,y,zc));
                }
        return total;
    }
};

template<class Apply>
double timed(Apply apply, cudaEvent_t begin, cudaEvent_t end) {
    constexpr int repeats = 20;
    check(cudaEventRecord(begin));
    for (int i = 0; i < repeats; ++i) apply();
    check(cudaEventRecord(end));
    check(cudaEventSynchronize(end));
    float ms;
    check(cudaEventElapsedTime(&ms, begin, end));
    if (!(ms > 0 && std::isfinite(ms))) throw std::runtime_error("invalid event interval");
    return ms / repeats;
}
}

int main() {
    constexpr int E = 32769, size = E*nodes;
    Inputs h(E);
    Device u(h.u), B(h.B), Bt(h.Bt), D(h.D), Dt(h.Dt), G(h.metric), S(h.reaction), GS(h.packed);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Device y(size,nan), split(size,nan), mass(size,nan), diff(size,nan);
    mars::marsir::GalerkinPtx fused("reaction_p3_sm90.ptx","reaction_apply_p3",p,q,true,7),
        bp1("bp1_p3_sm90.ptx","bp1_apply_p3",p,q,false,1), bp3("bp3_p3_sm90.ptx","bp3_apply_p3",p,q,true,6);
    auto apply = [&](mars::marsir::GalerkinPtx& op, Device& g, Device& out) {
        op.apply(ptr(u),ptr(B),ptr(Bt),ptr(D),ptr(Dt),ptr(g),ptr(out),E);
    };
    auto one = [&] { apply(fused,GS,y); };
    auto three = [&] {
        apply(bp3,G,diff); apply(bp1,S,mass);
        add<<<(size+255)/256,256>>>(ptr(diff),ptr(mass),ptr(split),size);
        check(cudaGetLastError());
    };
    one(); three(); check(cudaDeviceSynchronize());
    Vec yf(size), ys(size);
    check(cudaMemcpy(yf.data(),ptr(y),size*sizeof(double),cudaMemcpyDeviceToHost));
    check(cudaMemcpy(ys.data(),ptr(split),size*sizeof(double),cudaMemcpyDeviceToHost));
    double error=0, scale=0, oracle_error=0;
    for (int i=0; i<size; ++i) {
        if (!std::isfinite(yf[i]) || !std::isfinite(ys[i])) throw std::runtime_error("nonfinite output");
        error=std::max(error,std::abs(yf[i]-ys[i])); scale=std::max(scale,std::abs(ys[i]));
    }
    for (int e : {0, 1, 16, E-1})
        for (int i=0; i<nodes; ++i) {
            double ref=h.expected(e,i);
            if (!std::isfinite(ref)) throw std::runtime_error("nonfinite oracle");
            oracle_error=std::max({oracle_error,std::abs(yf[e*nodes+i]-ref),std::abs(ys[e*nodes+i]-ref)});
        }
    if (!(scale>0 && std::max(error,oracle_error)<2e-12*scale)) throw std::runtime_error("field comparison failed");
    std::printf("Correctness PASS: %d entries; split difference %.4e; direct quadrature difference %.4e; scale %.4e\n",
                size,error,oracle_error,scale);
    cudaEvent_t begin,end; check(cudaEventCreate(&begin)); check(cudaEventCreate(&end));
    for (int r=0; r<3; ++r) { one(); three(); }
    check(cudaDeviceSynchronize());
    std::array<double,6> ratios;
    for (int r=0; r<6; ++r) {
        double a,b;
        if (r%2==0) { a=timed(one,begin,end); b=timed(three,begin,end); }
        else { b=timed(three,begin,end); a=timed(one,begin,end); }
        ratios[r]=b/a;
        std::printf("round %d fused_ms %.6f split_plus_add_ms %.6f ratio %.4f\n",r+1,a,b,b/a);
    }
    std::sort(ratios.begin(),ratios.end());
    std::printf("Median paired split/fused %.4f; range %.4f to %.4f\n",(ratios[2]+ratios[3])/2,ratios.front(),ratios.back());
    check(cudaEventDestroy(end)); check(cudaEventDestroy(begin));
}
