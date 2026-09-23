#pragma once
#include "mars_segregated_assembly.hpp"
#include "mars_segregated_update.hpp"

#if defined(__CUDACC__)
#define MARS_SIMPLE_HD __host__ __device__
#else
#define MARS_SIMPLE_HD
#endif
namespace mars::segregated {

// Restricted steady, laminar SIMPLE profile: fixed frame, upwind, one open outlet.
struct SimpleControls {
    double density=1, viscosity=.1, pseudo_dt=.01, inlet_speed=.1;
    double alpha_u=.3, alpha_p=.3, alpha_mass=.75, beta=.05, pressure_reference=0;
};
struct SimpleFace { int element, ordinal, kind; }; // inlet=0, outlet=1, no-slip wall=2

MARS_SIMPLE_HD inline TetInteriorInput simple_interior(int stage, const int* nodes,
    const double* coordinates, const double* velocity, const double* pressure,
    const double* flux, const SimpleControls& c)
{
    TetInteriorInput x{}; x.stage=stage;
    for (int i=0;i<12;++i) x.coordinates[i]=coordinates[i];
    for (int s=0;s<6;++s) x.stored_flux[s]=flux[s];
    for (int k=0;k<4;++k) {
        x.density[k]=c.density; x.viscosity[k]=c.viscosity; x.pressure[k]=pressure[nodes[k]];
        for (int j=0;j<3;++j) x.velocity[3*k+j]=velocity[3*nodes[k]+j];
    }
    return x;
}

MARS_SIMPLE_HD inline BoundaryAssemblyInput simple_boundary(bool momentum, SimpleFace face,
    const int* nodes, const TetGeometry<double>& g, const double* velocity, const double* pressure,
    const double* trace, const double* flux, const SimpleControls& c, bool wall_initialized, const int* reversal=nullptr)
{
    BoundaryAssemblyInput input{}; input.element=face.element; input.face=face.ordinal;
    auto& x=input.values; x.stage=face.kind+(momentum?3:0);
    double area[3]; tet_boundary_area(g,face.ordinal,area);
    double magnitude=0; for (int j=0;j<3;++j) magnitude+=area[j]*area[j]; magnitude=sqrt(magnitude);
    const int opposite=tet_opposite_node(face.ordinal);
    for (int k=0;k<4;++k) {
        input.nodes[k]=nodes[k]; x.pressure[k]=pressure[nodes[k]];
        for (int j=0;j<3;++j) x.velocity[3*k+j]=velocity[3*nodes[k]+j];
    }
    for (int f=0;f<3;++f) {
        const int local=tet_face_node(face.ordinal,f);
        x.face_nodes[f]=x.nearest[f]=local; x.opposing[f]=opposite;
        x.reversal[f]=face.kind==1 && reversal?reversal[f]:0;
        x.density[f]=c.density; x.viscosity[f]=c.viscosity; x.stored_flux[f]=flux[f];
        for (int j=0;j<3;++j) {
            x.boundary_velocity[3*f+j]=face.kind==0?-c.inlet_speed*area[j]/magnitude:0;
            if (x.stage==3) x.velocity[3*local+j]=x.boundary_velocity[3*f+j];
        }
        if (x.stage==1) x.pressure[local]=trace[f];
    }
    if (x.stage==5) {
        // The full normal edge length is 3V/|A_face| = V/|A_sample|.
        // OpenAccel initializes this laminar coefficient after its first outer solve.
        for (int f=0;f<3;++f) {
            const int node=nodes[tet_face_node(face.ordinal,f)];
            input.nodes[f]=node; x.face_nodes[f]=x.nearest[f]=f;
            for (int j=0;j<3;++j) x.velocity[3*f+j]=velocity[3*node+j];
            x.wall_coefficient[f]=wall_initialized?4*c.viscosity*magnitude*magnitude/g.volume:0;
        }
    }
    return input;
}

// Expand sorted node blocks to sorted scalar rows without transposing components.
template<int Components>
MARS_SIMPLE_HD inline void simple_scalar_row(BlockCsrView<Components> a, int row,
    int* offsets, int* columns, double* values)
{
    const int node=row/Components, component=row%Components;
    const int begin=a.offsets[node], end=a.offsets[node+1];
    const int start=Components*Components*begin+component*Components*(end-begin);
    offsets[row]=start;
    for (int k=begin;k<end;++k) for (int j=0;j<Components;++j) {
        const int out=start+(k-begin)*Components+j;
        columns[out]=Components*a.columns[k]+j;
        values[out]=a.values[Components*Components*k+Components*component+j];
    }
    if (row+1==Components*a.nodes) offsets[row+1]=Components*Components*end;
}

MARS_SIMPLE_HD inline void simple_error(int* error) {
#if defined(__CUDA_ARCH__)
    atomicExch(error,1);
#else
    *error=1;
#endif
}
struct SimpleMesh {
    int node_count, element_count, face_count;
    const int* nodes[4];
    const double *x,*y,*z;
    const SimpleFace* faces;
    TetGeometry<double>* geometry;
    MARS_SIMPLE_HD void cell(int e,int* n,double* xyz) const {
        for (int k=0;k<4;++k) { n[k]=nodes[k][e]; xyz[3*k]=x[n[k]]; xyz[3*k+1]=y[n[k]]; xyz[3*k+2]=z[n[k]]; }
    }
};
struct SimpleState {
    double *velocity,*pressure,*velocity_gradient,*pressure_gradient,*influence;
    double *volume,*mass_divergence,*interior_flux,*boundary_flux,*trace,*boundary_factor;
    int* error;
    int* reversal=nullptr;
};
struct SimpleGeometry {
    SimpleMesh mesh; SimpleState state;
    MARS_SIMPLE_HD void operator()(int e) const {
        int nodes[4]; double xyz[12]; mesh.cell(e,nodes,xyz);
        if (!tet_geometry(xyz,mesh.geometry[e])) { simple_error(state.error); return; }
        for (int k=0;k<4;++k) assembly_add(state.volume+nodes[k],mesh.geometry[e].volume/4);
    }
};
struct SimpleBoundaryFactor {
    SimpleMesh mesh; double* factor;
    MARS_SIMPLE_HD void operator()(int i) const {
        const auto f=mesh.faces[i];
        for (int j=0;j<3;++j)
            assembly_add(factor+mesh.nodes[tet_face_node(f.ordinal,j)][f.element],1.);
    }
};
template<int Components> struct SimpleGradientInterior {
    SimpleMesh mesh; const double* field; double* sum;
    MARS_SIMPLE_HD void operator()(int e) const {
        double values[4*Components],local[12*Components];
        for (int k=0;k<4;++k) for (int c=0;c<Components;++c) values[k*Components+c]=field[mesh.nodes[k][e]*Components+c];
        tet_gradient_numerator<Components>(mesh.geometry[e],values,Components==1,true,local);
        for (int k=0;k<4;++k) for (int c=0;c<3*Components;++c) assembly_add(sum+mesh.nodes[k][e]*3*Components+c,local[k*3*Components+c]);
    }
};
template<int Components> struct SimpleGradientBoundary {
    SimpleMesh mesh; const double* field; double* sum;
    MARS_SIMPLE_HD void operator()(int i) const {
        auto f=mesh.faces[i]; double values[3*Components],local[9*Components],area[3];
        tet_boundary_area(mesh.geometry[f.element],f.ordinal,area);
        for (int k=0;k<3;++k) for (int c=0;c<Components;++c) values[k*Components+c]=field[mesh.nodes[tet_face_node(f.ordinal,k)][f.element]*Components+c];
        tri_gradient_numerator<Components>(area,values,Components==1,true,local);
        for (int k=0;k<3;++k) for (int c=0;c<3*Components;++c) assembly_add(sum+mesh.nodes[tet_face_node(f.ordinal,k)][f.element]*3*Components+c,local[k*3*Components+c]);
    }
};
template<int Components> struct SimpleGradientFinish {
    const double *volume,*sum; double* gradient; int* error;
    MARS_SIMPLE_HD void operator()(int i) const {
        if (!(volume[i/(3*Components)]>0)) { simple_error(error); return; }
        gradient[i]=sum[i]/volume[i/(3*Components)];
    }
};
template<int Components> struct SimpleInterior {
    SimpleMesh mesh; SimpleState state; SimpleControls controls; BlockCsrView<Components> matrix;
    bool update_flux=false;
    MARS_SIMPLE_HD void operator()(int e) const {
        int nodes[4]; double xyz[12]; mesh.cell(e,nodes,xyz);
        auto x=simple_interior(Components==3?1:0,nodes,xyz,state.velocity,state.pressure,state.interior_flux+6*e,controls);
        native_interior(x,mesh.geometry[e],nodes,state.velocity_gradient,state.pressure_gradient,state.influence);
        TetInteriorOutput y; tet_interior(x,y);
        if (!update_flux) {
            if (!scatter_block(matrix,nodes,4,y.lhs,y.rhs)) simple_error(state.error);
        } else for (int s=0;s<6;++s) {
            const double q=controls.alpha_mass*y.flux[s]+(1-controls.alpha_mass)*state.interior_flux[6*e+s];
            state.interior_flux[6*e+s]=q;
            assembly_add(state.mass_divergence+nodes[tet_edge_node(s,0)],q);
            assembly_add(state.mass_divergence+nodes[tet_edge_node(s,1)],-q);
        }
    }
};
template<int Components> struct SimpleBoundary {
    SimpleMesh mesh; SimpleState state; SimpleControls controls; BlockCsrView<Components> matrix;
    bool wall_initialized=false, update_flux=false;
    MARS_SIMPLE_HD void operator()(int i) const {
        const auto face=mesh.faces[i]; int nodes[4]; double xyz[12]; mesh.cell(face.element,nodes,xyz);
        const auto& g=mesh.geometry[face.element];
        auto input=simple_boundary(Components==3,face,nodes,g,state.velocity,state.pressure,state.trace+3*i,state.boundary_flux+3*i,controls,wall_initialized,state.reversal?state.reversal+3*i:nullptr);
        auto& x=input.values;
        if (!native_boundary(x,input,g,nodes,state.pressure_gradient,state.influence)) { simple_error(state.error); return; }
        BoundaryOutput y; boundary_block(x,y);
        if (!update_flux) {
            if (!scatter_block(matrix,input.nodes,x.stage==5?3:4,y.lhs,y.rhs)) simple_error(state.error);
        } else {
            double flux[3];
            for (int j=0;j<3;++j) flux[j]=x.reversal[j]?0:controls.alpha_mass*y.flux[j]+(1-controls.alpha_mass)*state.boundary_flux[3*i+j];
            if (face.kind==1) {
                int next[3]; double velocity[9],pressure[3],filtered[3];
                for (int f=0;f<3;++f) {
                    const int n=nodes[tet_face_node(face.ordinal,f)]; pressure[f]=state.pressure[n];
                    for (int j=0;j<3;++j) velocity[3*f+j]=state.velocity[3*n+j];
                }
                outlet_reversal_update(flux,x.reversal,velocity,pressure,state.trace+3*i,x.area,false,filtered,next);
                for (int j=0;j<3;++j) {
                    if (state.reversal) state.reversal[3*i+j]=next[j];
                    else if (next[j]) simple_error(state.error);
                    flux[j]=filtered[j];
                }
            }
            for (int j=0;j<3;++j) {
                state.boundary_flux[3*i+j]=flux[j];
                assembly_add(state.mass_divergence+nodes[tet_face_node(face.ordinal,j)],flux[j]);
            }
        }
    }
};
struct SimpleMomentumNode {
    SimpleState state; SimpleControls controls; BlockCsrView<3> matrix;
    MARS_SIMPLE_HD void operator()(int n) const {
        SteadyMomentumNode x{}; x.density=controls.density; x.pseudo_dt=controls.pseudo_dt;
        x.volume=state.volume[n]; x.mass_divergence=state.mass_divergence[n];
        for (int j=0;j<3;++j) { x.velocity[j]=state.velocity[3*n+j]; x.pressure_gradient[j]=state.pressure_gradient[3*n+j]; }
        double a[9],b[3],dt[3]; steady_momentum_node(x,a,b);
        if (!scatter_block(matrix,&n,1,a,b) || !finish_momentum_row(matrix,n,x.volume,controls.alpha_u,state.boundary_factor[n]>0?.75:1.,false,state.influence+3*n,dt)) simple_error(state.error);
    }
};
struct SimpleTraceMoment {
    SimpleMesh mesh; SimpleState state; double* moment;
    MARS_SIMPLE_HD void operator()(int i) const {
        const auto f=mesh.faces[i]; if (f.kind!=1 || (state.reversal && state.reversal[3*i])) return;
        double area[3]; tet_boundary_area(mesh.geometry[f.element],f.ordinal,area);
        double magnitude=0; for (int j=0;j<3;++j) magnitude+=area[j]*area[j]; magnitude=sqrt(magnitude);
        for (int j=0;j<3;++j) {
            assembly_add(moment,state.pressure[mesh.nodes[tet_face_node(f.ordinal,j)][f.element]]*magnitude);
            assembly_add(moment+1,magnitude);
        }
    }
};
struct SimpleTrace {
    SimpleMesh mesh; SimpleState state; SimpleControls controls; const double* moment;
    MARS_SIMPLE_HD void operator()(int i) const {
        auto f=mesh.faces[i]; if (f.kind!=1 || (state.reversal && state.reversal[3*i])) return;
        if (!(moment[1]>0)) { simple_error(state.error); return; }
        for (int j=0;j<3;++j) state.trace[3*i+j]=outlet_trace_update(state.pressure[mesh.nodes[tet_face_node(f.ordinal,j)][f.element]],controls.pressure_reference,moment[0]/moment[1],controls.beta,state.trace[3*i+j],false);
    }
};
struct SimpleAddIncrement {
    double* field; const double* increment; double alpha;
    MARS_SIMPLE_HD void operator()(int i) const { field[i]+=alpha*increment[i]; }
};
struct SimpleCorrectVelocity {
    SimpleState state; const double* increment_gradient;
    MARS_SIMPLE_HD void operator()(int n) const { velocity_update(state.velocity+3*n,state.influence+3*n,increment_gradient+3*n,state.velocity+3*n); }
};
template<int Components> struct SimpleScalarRows {
    BlockCsrView<Components> matrix; int *offsets,*columns; double* values;
    MARS_SIMPLE_HD void operator()(int row) const { simple_scalar_row(matrix,row,offsets,columns,values); }
};
struct SimpleResidual {
    const int *offsets,*columns; const double *values,*rhs,*solution; double* squares;
    MARS_SIMPLE_HD void operator()(int row) const {
        double r=-rhs[row]; for (int k=offsets[row];k<offsets[row+1];++k) r+=values[k]*solution[columns[k]];
        squares[2*row]=r*r; squares[2*row+1]=rhs[row]*rhs[row];
    }
};
} // namespace mars::segregated
#undef MARS_SIMPLE_HD
