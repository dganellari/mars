"""Exercise the actual host inspection method with lightweight coordinate storage."""
import argparse
from pathlib import Path
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--compiler', required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--sysroot', default='')
    args = parser.parse_args()
    source = (args.source / 'backend/distributed/unstructured/domain.hpp').read_text()
    start = source.index('    auto getElementNodeCoordinates(size_t elemIdx)')
    end = source.index('\n    // Helper to get local node ID', start)
    method = source[start:end]
    fixture = r'''
#include <array>
#include <tuple>
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <utility>
using std::size_t;
template<size_t... I> auto connectivity(std::index_sequence<I...>) {
    return std::make_tuple(std::array<uint64_t,2>{sizeof...(I)-1-I,I}...);
}
template<int N> struct Domain {
    using KeyType=uint64_t; using RealType=double;
    static constexpr int NodesPerElement=N;
    struct Cache { std::array<double,N> d_node_x_,d_node_y_,d_node_z_; } cache, original;
    Cache *coordCache_=nullptr, *originalCoords_=nullptr;
    decltype(connectivity(std::make_index_sequence<N>{})) d_conn_keys_=connectivity(std::make_index_sequence<N>{});
    Domain(int mode) {
        for(int i=0;i<N;++i) {
            cache.d_node_x_[i]=10+i; cache.d_node_y_[i]=20+i; cache.d_node_z_[i]=30+i;
            original.d_node_x_[i]=100+i; original.d_node_y_[i]=200+i; original.d_node_z_[i]=300+i;
        }
        if(mode==1 || mode==3) coordCache_=&cache;
        if(mode==2 || mode==3) originalCoords_=&original;
    }
    const auto& getElementToNodeConnectivity() { return d_conn_keys_; }
    const auto& getNodeX() { return (originalCoords_ ? originalCoords_ : coordCache_)->d_node_x_; }
    const auto& getNodeY() { return (originalCoords_ ? originalCoords_ : coordCache_)->d_node_y_; }
    const auto& getNodeZ() { return (originalCoords_ ? originalCoords_ : coordCache_)->d_node_z_; }
    auto sfcToPhysicalCoordinate(uint64_t k) { return std::make_tuple(double(k),double(k+1),double(k+2)); }
METHOD
};
template<int N> void check() {
    for(int mode=0;mode<4;++mode) for(int element=0;element<2;++element) {
        Domain<N> domain(mode);
        auto [x,y,z]=domain.getElementNodeCoordinates(element);
        for(int i=0;i<N;++i) {
            int k=element ? i : N-1-i;
            assert(x[i]==k+(mode==0 ? 0 : mode==1 ? 10 : 100));
            assert(y[i]==k+(mode==0 ? 1 : mode==1 ? 20 : 200));
            assert(z[i]==k+(mode==0 ? 2 : mode==1 ? 30 : 300));
        }
    }
}
int main() { check<4>(); check<8>(); }
'''.replace('METHOD', method)
    with tempfile.TemporaryDirectory(dir=args.work) as tmp:
        work = Path(tmp)
        (work/'test.cpp').write_text(fixture)
        subprocess.run([args.compiler, '-std=c++20'] + (['-isysroot', args.sysroot] if args.sysroot else []) + [str(work/'test.cpp'), '-o', str(work/'test')], check=True)
        subprocess.run([str(work/'test')], check=True)
    print('PASS: Tet4/Hex8, every corner, two elements, decoded/cached/original/both modes')


if __name__ == '__main__':
    main()
