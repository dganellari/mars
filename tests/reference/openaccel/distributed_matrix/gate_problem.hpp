#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <vector>

namespace dmatrix_gate {
inline std::uint64_t mix(std::uint64_t x) {
    x+=0x9e3779b97f4a7c15ull; x=(x^(x>>30))*0xbf58476d1ce4e5b9ull; x=(x^(x>>27))*0x94d049bb133111ebull;
    return x^(x>>31);
}
inline std::uint64_t key(std::uint64_t a,std::uint64_t b,std::uint64_t c,std::uint64_t d) { return mix(mix(mix(mix(a)^b)^c)^d); }
inline double symmetric(std::uint64_t k) { return double(mix(k)>>11)*0x1.0p-52-1.0; } // [-1,1)
template<class T> void shuffle(std::vector<T>& v,std::uint64_t seed) {
    for (std::size_t i=v.size();i>1;--i) std::swap(v[i-1],v[mix(seed^(i*0x100000001b3ull))%i]);
}

struct Options {
    int nx=6, ny=5, nz=4;          // 120 nodes, 27-point couplings
    bool empty_rank=false;         // rank 1 owns no rows (needs >= 2 ranks)
    bool empty_rank_holds_nothing=false;
    bool far_ghosts=true;          // partial second ghost layer, never used by owned rows
    bool far_ghosts_without_id=false;
    bool omit_coupling=false;      // last rank drops its largest ghost coupling from one owned row
    bool diffusion=false;         // connected seven-point operator for multilevel AMG coverage
};

// Public systems with a source-keyed solution and an uneven, shuffled partition.
struct Problem {
    int C, ranks, nodes;
    Options options;
    std::vector<std::vector<int>> neighbors, owned_order;
    std::vector<int> owner, weight, solver_to_global;
    std::vector<long long> first, solver_node, source_id;
    Problem(int components,int rank_count,Options o):C(components),ranks(rank_count),nodes(o.nx*o.ny*o.nz),options(o) {
        neighbors.resize(nodes);
        for (int z=0;z<o.nz;++z) for (int y=0;y<o.ny;++y) for (int x=0;x<o.nx;++x)
            for (int dz=-1;dz<=1;++dz) for (int dy=-1;dy<=1;++dy) for (int dx=-1;dx<=1;++dx) {
                if (o.diffusion && std::abs(dx)+std::abs(dy)+std::abs(dz)>1) continue;
                const int X=x+dx, Y=y+dy, Z=z+dz;
                if (X>=0 && Y>=0 && Z>=0 && X<o.nx && Y<o.ny && Z<o.nz) neighbors[x+o.nx*(y+o.ny*z)].push_back(X+o.nx*(Y+o.ny*Z));
            }
        for (auto& n:neighbors) std::sort(n.begin(),n.end());
        std::vector<int> permutation(nodes); for (int g=0;g<nodes;++g) permutation[g]=g;
        shuffle(permutation,7); source_id.resize(nodes);
        for (int g=0;g<nodes;++g) source_id[g]=1000003+17LL*permutation[g];
        // Uneven weights, then irregular seams: every fifth node moves to the next owning rank.
        weight.resize(ranks); int total=0;
        for (int r=0;r<ranks;++r) total+=weight[r]=(o.empty_rank && r==1)?0:r+1;
        owner.resize(nodes);
        for (int g=0,r=0,acc=0;g<nodes;++g) {
            while (r<ranks-1 && (long long)g*total>=(long long)(acc+weight[r])*nodes) acc+=weight[r++];
            while (weight[r]==0) r=(r+1)%ranks;
            owner[g]=r;
            if (ranks>1 && mix(31ull*g+7)%5==0) { int q=(r+1)%ranks; while (weight[q]==0) q=(q+1)%ranks; owner[g]=q; }
        }
        owned_order.resize(ranks); first.assign(ranks+1,0); solver_node.resize(nodes); solver_to_global.resize(nodes);
        for (int g=0;g<nodes;++g) owned_order[owner[g]].push_back(g);
        for (int r=0;r<ranks;++r) {
            shuffle(owned_order[r],1000+r); first[r+1]=first[r]+(long long)owned_order[r].size();
            for (std::size_t k=0;k<owned_order[r].size();++k) {
                solver_node[owned_order[r][k]]=first[r]+(long long)k; solver_to_global[first[r]+k]=owned_order[r][k];
            }
        }
    }
    // The default random entries expose source-id mapping mistakes.
    double value(int gi,int gj,int c,int j,unsigned seed) const {
        if (options.diffusion) {
            if (c!=j) return 0;
            const double wx=1+.125*seed, wy=1+.25*seed, wz=1+.375*seed;
            // Positive reaction and zero exterior Dirichlet values remove the nullspace.
            if (gi==gj) return 2*(wx+wy+wz)+.125;
            return gi%options.nx!=gj%options.nx?-wx:
                (gi/options.nx)%options.ny!=(gj/options.nx)%options.ny?-wy:-wz;
        }
        const double v=symmetric(key(std::uint64_t(source_id[gi]),std::uint64_t(source_id[gj]),std::uint64_t(8*c+j),seed));
        return gi==gj && c==j?40.0*C+1.0+0.5*v:0.5*v;
    }
    double rhs(int g,int c,unsigned seed) const { return symmetric(key(std::uint64_t(source_id[g]),77,std::uint64_t(c),seed)); }
    double solution(int g,int c,unsigned seed) const { return symmetric(key(std::uint64_t(source_id[g]),99,std::uint64_t(c),seed)); }
    double product(int g,int c,unsigned value_seed,unsigned x_seed) const {
        long double s=0;
        for (int u:neighbors[g]) for (int j=0;j<C;++j) s+=(long double)value(g,u,c,j,value_seed)*solution(u,j,x_seed);
        return double(s);
    }
    bool coupled(int g,int u) const { return std::binary_search(neighbors[g].begin(),neighbors[g].end(),u); }
    // Nodes present on rank r: owned, every coupled neighbour (complete owned rows), and far ghosts.
    std::vector<int> holds(int r,std::set<int>* far=nullptr) const {
        std::set<int> present(owned_order[r].begin(),owned_order[r].end()), ghosts, extra;
        for (int g:owned_order[r]) for (int u:neighbors[g]) if (!present.count(u)) ghosts.insert(u);
        if (options.empty_rank && weight[r]==0 && !options.empty_rank_holds_nothing)
            for (int g=0;g<nodes;++g) if (mix(std::uint64_t(g))%11==0) ghosts.insert(g);
        if (options.far_ghosts)
            for (int g:ghosts) for (int u:neighbors[g])
                if (!present.count(u) && !ghosts.count(u) && mix(std::uint64_t(u)+5)%3==0) extra.insert(u);
        present.insert(ghosts.begin(),ghosts.end()); present.insert(extra.begin(),extra.end());
        if (far) *far=extra;
        return {present.begin(),present.end()};
    }
};

}
