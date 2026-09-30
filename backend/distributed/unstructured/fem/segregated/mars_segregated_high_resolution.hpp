#pragma once
#include "mars_segregated_simple.hpp"

#if defined(__CUDACC__)
#define MARS_BLEND_HD __host__ __device__
#else
#define MARS_BLEND_HD
#endif

namespace mars::segregated {

// OpenAccel's ordinary high-resolution scheme, blend_factor_max=1, without
// gradient limiting or shock damping. Both sides of the local range constrain it.
MARS_BLEND_HD inline double velocity_blend_bound(double value,double lower,double upper,double projection) {
    constexpr double small=std::numeric_limits<double>::epsilon();
    const double r=projection+(projection>0?small:-small);
    if (fabs(r)<=100*small) return 1;
    const double y=fmin(upper-value,value-lower)/fabs(r);
    // Equivalent to (y*y+2*y)/(y*y+y+2), capped at one. Avoid overflow for tiny slopes.
    if (y>=2) return 1;
    return (y*y+2*y)/(y*y+y+2);
}
MARS_BLEND_HD inline double relax_velocity_blend(double bound,double previous) {
    return fmin(bound,.75*previous+.25*bound);
}
MARS_BLEND_HD inline void blend_min(double* address,double value) {
#if defined(__CUDA_ARCH__)
    // All candidates lie in [0,1], where positive IEEE double bits are ordered.
    atomicMin(reinterpret_cast<unsigned long long*>(address),static_cast<unsigned long long>(__double_as_longlong(value)));
#else
    *address=fmin(*address,value);
#endif
}

struct SimpleBlendBounds {
    BlockCsrView<3> graph; const double* velocity;
    double *lower,*upper,*candidate; int* error;
    MARS_BLEND_HD void operator()(int node) const {
        for (int c=0;c<3;++c) {
            double lo=velocity[3*node+c],hi=lo;
            for (int k=graph.offsets[node];k<graph.offsets[node+1];++k) {
                const double v=velocity[3*graph.columns[k]+c];
                if (!geometry_finite(v)) simple_error(error);
                lo=fmin(lo,v); hi=fmax(hi,v);
            }
            lower[3*node+c]=lo; upper[3*node+c]=hi; candidate[3*node+c]=1;
        }
    }
};
struct SimpleBlendSamples {
    SimpleMesh mesh; SimpleState state;
    const double *lower,*upper; double* candidate;
    const unsigned char* owned=nullptr;
    MARS_BLEND_HD void sample(int node,const double* point) const {
        if (owned && !owned[node]) return;
        const double delta[3]={point[0]-mesh.x[node],point[1]-mesh.y[node],point[2]-mesh.z[node]};
        for (int c=0;c<3;++c) {
            double r=0;
            for (int j=0;j<3;++j) r+=state.velocity_gradient[9*node+3*c+j]*delta[j];
            if (!geometry_finite(r)) { simple_error(state.error); continue; }
            blend_min(candidate+3*node+c,velocity_blend_bound(state.velocity[3*node+c],lower[3*node+c],upper[3*node+c],r));
        }
    }
};
struct SimpleBlendInterior : SimpleBlendSamples {
    MARS_BLEND_HD void operator()(int e) const {
        int nodes[4]; double xyz[12]; mesh.cell(e,nodes,xyz);
        for (int s=0;s<6;++s) {
            double shape[4],point[3]{}; tet_sample_shape(s,false,shape);
            for (int k=0;k<4;++k) for (int j=0;j<3;++j) point[j]+=shape[k]*xyz[3*k+j];
            for (int end=0;end<2;++end) sample(nodes[tet_edge_node(s,end)],point);
        }
    }
};
struct SimpleBlendBoundary : SimpleBlendSamples {
    MARS_BLEND_HD void operator()(int i) const {
        const auto f=mesh.faces[i];
        int nodes[4]; double xyz[12]; mesh.cell(f.element,nodes,xyz);
        // Every boundary sample constrains the limiter, including closed outlets.
        for (int s=0;s<3;++s) {
            double shape[3],point[3]{}; tri_sample_shape(s,false,shape);
            for (int k=0;k<3;++k) for (int j=0;j<3;++j) point[j]+=shape[k]*xyz[3*tet_face_node(f.ordinal,k)+j];
            sample(nodes[tet_face_node(f.ordinal,s)],point);
        }
    }
};
struct SimpleBlendFinish {
    const double* candidate; double* blend;
    MARS_BLEND_HD void operator()(int node) const {
        for (int c=0;c<3;++c) blend[3*node+c]=relax_velocity_blend(candidate[3*node+c],blend[3*node+c]);
    }
};
} // namespace mars::segregated

#undef MARS_BLEND_HD
