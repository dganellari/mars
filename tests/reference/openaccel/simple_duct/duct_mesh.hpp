#pragma once
// C++ mirror of duct_mesh.py (same lattice, numbering, orientation and side sets; the
// canonical() text must match `duct_mesh.py --dump` byte for byte): x in [0,L] (inlet x=0,
// outlet x=L), y in [-W/2,W/2], z in [-H/2,H/2], no-slip walls on the four sides. nx*ny*nz
// hexes, six positively oriented Kuhn tets per hex around the hex diagonal, so the mesh is
// invariant under x -> x+hx and a discrete fully developed state exists.
//
// slab_view(): a test-only partition in the shape ElementDomain presents (keys = lattice ids,
// sorted local nodes, lower-halo/own/upper-halo element bands, lowest-rank node ownership,
// complete element stars, NodeHaloTopology-style lists), for the production simple_partition.
#include "mars_segregated_simple_partition.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <cstdio>
#include <string>
#include <vector>

namespace duct {

struct Lattice {
    int nx=0, ny=0, nz=0;
    double length=7, width=2, height=1;

    // cells across the height; hx = stretch*hz. Sizes must come out as integers, ny and nz even
    // (the centerline is a node line).
    static Lattice make(int cells,double length,double width,double height,double stretch) {
        Lattice l; l.length=length; l.width=width; l.height=height; l.nz=cells;
        const double ny=cells*width/height, nx=cells*length/(stretch*height);
        l.ny=int(std::lround(ny)); l.nx=int(std::lround(nx));
        if (cells<2 || cells%2 || std::abs(ny-l.ny)>1e-9 || l.ny%2 || std::abs(nx-l.nx)>1e-9 || l.nx<1)
            throw std::runtime_error("duct lattice: need even cells, W/H*cells even, L*cells/(stretch*H) integer");
        return l;
    }
    double hx() const { return length/nx; }
    double hy() const { return width/ny; }
    double hz() const { return height/nz; }
    long long nodes() const { return (long long)(nx+1)*(ny+1)*(nz+1); }
    long long elements() const { return 6LL*nx*ny*nz; }
    long long id(int i,int j,int k) const { return i+(long long)(nx+1)*(j+(long long)(ny+1)*k); }
    std::array<int,3> ijk(long long g) const {
        const int i=int(g%(nx+1)); g/=nx+1; return {i,int(g%(ny+1)),int(g/(ny+1))};
    }
    double x(int i) const { return i*hx(); }
    double y(int j) const { return -width/2+j*hy(); }
    double z(int k) const { return -height/2+k*hz(); }
    std::array<double,3> coordinates(long long g) const { const auto c=ijk(g); return {x(c[0]),y(c[1]),z(c[2])}; }

    // Element el = t + 6*(hex i + nx*(j + ny*k)), t the Kuhn path permutation.
    std::array<long long,4> element(long long el) const {
        static constexpr int axes[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};
        const int t=int(el%6); long long h=el/6;
        const int i=int(h%nx); h/=nx; const int j=int(h%ny), k=int(h/ny);
        int c[3]={i,j,k}; std::array<long long,4> n;
        n[0]=id(i,j,k);
        for (int s=0;s<2;++s) { ++c[axes[t][s]]; n[s+1]=id(c[0],c[1],c[2]); }
        n[3]=id(i+1,j+1,k+1);
        // Orientation from the lattice steps (the same for every hex).
        double m[3][3];
        for (int r=0;r<3;++r) { const auto a=coordinates(n[r+1]), b=coordinates(n[0]); for (int q=0;q<3;++q) m[r][q]=a[q]-b[q]; }
        const double det=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])
                        +m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
        if (det<0) std::swap(n[2],n[3]);
        return n;
    }
    // 0 inlet, 1 outlet, 2 wall, -1 not on the duct boundary.
    int kind(const long long* face) const {
        std::array<int,3> c[3]; for (int j=0;j<3;++j) c[j]=ijk(face[j]);
        auto all=[&](int axis,int value) { return c[0][axis]==value && c[1][axis]==value && c[2][axis]==value; };
        if (all(0,0)) return 0;
        if (all(0,nx)) return 1;
        return all(1,0) || all(1,ny) || all(2,0) || all(2,nz) ? 2 : -1;
    }
    // Boundary triangles per kind: 2 per hex face.
    std::array<long long,3> boundary_faces() const {
        return {2LL*ny*nz,2LL*ny*nz,4LL*nx*nz+4LL*nx*ny};
    }
    bool on_wall(long long g) const { const auto c=ijk(g); return c[1]==0 || c[1]==ny || c[2]==0 || c[2]==nz; }
};

// The whole mesh as a SimpleInput in lattice numbering (one-rank production SimpleRunner).
inline mars::segregated::SimpleInput global_input(const Lattice& l) {
    mars::segregated::SimpleInput f;
    for (long long g=0;g<l.nodes();++g) { const auto c=l.coordinates(g); f.x.push_back(c[0]); f.y.push_back(c[1]); f.z.push_back(c[2]); }
    for (long long el=0;el<l.elements();++el) { const auto n=l.element(el); for (int k=0;k<4;++k) f.nodes[k].push_back(int(n[k])); }
    for (long long el=0;el<l.elements();++el) {
        const auto n=l.element(el);
        for (int o=0;o<4;++o) {
            long long face[3]; for (int j=0;j<3;++j) face[j]=n[mars::segregated::tet_face_node(o,j)];
            const int kind=l.kind(face);
            if (kind>=0) f.faces.push_back({int(el),o,kind});
        }
    }
    return f;
}

// Canonical text shared with duct_mesh.py --dump.
inline std::string canonical(const Lattice& l) {
    std::string out; char line[160];
    auto add=[&](int n) { out.append(line,std::size_t(n)); };
    add(std::snprintf(line,sizeof line,"duct %d %d %d %.17g %.17g %.17g\n",l.nx,l.ny,l.nz,l.length,l.width,l.height));
    for (long long g=0;g<l.nodes();++g) { const auto c=l.coordinates(g); add(std::snprintf(line,sizeof line,"n %.17g %.17g %.17g\n",c[0],c[1],c[2])); }
    for (long long el=0;el<l.elements();++el) { const auto n=l.element(el); add(std::snprintf(line,sizeof line,"e %lld %lld %lld %lld\n",n[0],n[1],n[2],n[3])); }
    const auto input=global_input(l);
    for (int kind=0;kind<3;++kind) {
        const char* name[]={"inlet","outlet","walls"};
        std::vector<std::pair<int,int>> faces;
        for (const auto& f:input.faces) if (f.kind==kind) faces.push_back({f.element,f.ordinal});
        std::sort(faces.begin(),faces.end());
        for (const auto& f:faces) add(std::snprintf(line,sizeof line,"s %s %d %d\n",name[kind],f.first,f.second));
    }
    return out;
}

// Element owner of the test partition: x slabs at uneven cuts; four ranks are a 2x2 split
// (x and y), so partition boundaries meet the walls, the inlet or outlet planes and each other.
inline int slab_owner(const Lattice& l,long long el,int ranks) {
    const auto n=l.element(el); double cx=0, cy=0;
    for (long long g:n) { const auto c=l.coordinates(g); cx+=c[0]/4; cy+=c[1]/4; }
    const double fx=cx/l.length, fy=(cy+l.width/2)/l.width;
    if (ranks==4) return 2*(fx>.41)+(fy>.37);
    int r=0;
    for (int q=1;q<ranks;++q) r+=fx>(q+.15*((q%2)?1:-1))/ranks;
    return r;
}

// Rank r's state as ElementDomain presents it. Test-only: every rank scans the whole lattice.
inline mars::segregated::runtime::DomainView<long long> slab_view(const Lattice& l,int rank,int ranks) {
    const long long N=l.nodes(), E=l.elements();
    std::vector<int> element_owner(std::size_t(E),0), node_owner(std::size_t(N),ranks);
    for (long long el=0;el<E;++el) {
        element_owner[el]=slab_owner(l,el,ranks);
        for (long long g:l.element(el)) node_owner[g]=std::min(node_owner[g],element_owner[el]);
    }
    // Present on q: owned by q, or touching a node q owns (complete stars). Held: their nodes.
    std::vector<std::vector<char>> held(std::size_t(ranks),std::vector<char>(std::size_t(N),0));
    std::vector<long long> mine;
    for (long long el=0;el<E;++el) {
        const auto n=l.element(el);
        for (int q=0;q<ranks;++q) {
            bool present=element_owner[el]==q;
            for (long long g:n) present=present || node_owner[g]==q;
            if (!present) continue;
            for (long long g:n) held[q][g]=1;
            if (q==rank) mine.push_back(el);
        }
    }
    mars::segregated::runtime::DomainView<long long> v;
    std::vector<int> local(std::size_t(N),-1);
    for (long long g=0;g<N;++g) if (held[rank][g]) {
        local[g]=int(v.key.size()); v.key.push_back(g);
        const auto c=l.coordinates(g); v.x.push_back(c[0]); v.y.push_back(c[1]); v.z.push_back(c[2]);
        v.owned.push_back(node_owner[g]==rank);
    }
    auto band=[&](long long el) { return element_owner[el]<rank?0:element_owner[el]==rank?1:2; };
    std::stable_sort(mine.begin(),mine.end(),[&](long long a,long long b) { return band(a)<band(b); });
    for (long long el:mine) { const auto n=l.element(el); for (int k=0;k<4;++k) v.nodes[k].push_back(local[n[k]]); }
    v.element_begin=int(std::count_if(mine.begin(),mine.end(),[&](long long el) { return band(el)==0; }));
    v.element_end=v.element_begin+int(std::count_if(mine.begin(),mine.end(),[&](long long el) { return band(el)==1; }));
    // Both sides list a peer's nodes in lattice order, so the slots line up.
    for (int q=0;q<ranks;++q) {
        if (q==rank) continue;
        std::vector<int> send, recv;
        for (std::size_t u=0;u<v.key.size();++u) {
            const long long g=v.key[u];
            if (node_owner[g]==q) recv.push_back(int(u));
            if (node_owner[g]==rank && held[q][g]) send.push_back(int(u));
        }
        if (send.empty() && recv.empty()) continue;
        v.peers.push_back(q);
        v.send_nodes.insert(v.send_nodes.end(),send.begin(),send.end()); v.send_offsets.push_back(int(v.send_nodes.size()));
        v.recv_nodes.insert(v.recv_nodes.end(),recv.begin(),recv.end()); v.recv_offsets.push_back(int(v.recv_nodes.size()));
    }
    return v;
}

} // namespace duct
