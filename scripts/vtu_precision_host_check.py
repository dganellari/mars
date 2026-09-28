#!/usr/bin/env python3
"""Exercise the production VTU formatter with invented data and host I/O stubs."""
from pathlib import Path
import subprocess
import tempfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
STUBS = r'''
#include <array>
#include <cstring>
#include <tuple>
#include <vector>
namespace cstone { namespace execution { struct Gpu {}; } template<class T> using DeviceVector=std::vector<T>; }
namespace thrust { template<class T> T* raw_pointer_cast(T* p){return p;} }
constexpr int MPI_COMM_WORLD=0, cudaMemcpyDeviceToHost=0;
inline int MPI_Comm_rank(int,int* r){*r=0;return 0;}
inline int MPI_Comm_size(int,int* n){*n=1;return 0;}
inline int cudaMemcpy(void* dst,const void* src,size_t n,int){std::memcpy(dst,src,n);return 0;}
struct HexTag {};
template<class T> struct ElemTraits {static constexpr int NodesPerElem=8;};
namespace mars {
template<class Tag,class R,class K,class Accel> struct ElementDomain {
    std::vector<R> x=std::vector<R>(8,1.000000000000003),y=x,z=x;
    std::array<std::vector<K>,8> conn;
    ElementDomain(){for(int i=0;i<8;++i)conn[i]={K(i)};}
    size_t startIndex()const{return 0;} size_t endIndex()const{return 1;}
    size_t getNodeCount()const{return 8;}
    const auto& getDomain()const{return *this;}
    const auto& getNodeX()const{return x;}
    const auto& getNodeY()const{return y;}
    const auto& getNodeZ()const{return z;}
    const auto& getElementToNodeConnectivity()const{return conn;}
};
}
'''
MAIN = r'''
int main(int argc,char** argv){
    if(argc!=2)return 1;
    mars::ElementDomain<HexTag,double,unsigned,cstone::execution::Gpu> domain;
    std::vector<double> u(8,0.12345678901234567),p(8,123456789.01234567),cell(1,1.2345678901234567);
    using Writer=mars::fem::VTUParallelWriter<unsigned,double>;
    using FD=Writer::FieldDesc;
    const std::vector<FD> fields{{"u",FD::Kind::PointScalar,&u,nullptr,nullptr},
        {"p",FD::Kind::PointScalar,&p,nullptr,nullptr},
        {"velocity",FD::Kind::PointVector3,&u,&u,&u},
        {"cell",FD::Kind::CellScalar,&cell,nullptr,nullptr}};
    Writer exact(std::string(argv[1])+"/exact",true);
    exact.writeMultiFieldFrame(1,0.12345678901234567,domain,fields);
    exact.writeMultiFieldFrame(2,0.24691357802469134,domain,fields);
    Writer legacy(std::string(argv[1])+"/legacy");
    legacy.writeMultiFieldFrame(1,1,domain,fields);
}
'''


def main():
    header = (ROOT / 'backend/distributed/unstructured/utils/mars_vtu_parallel_writer.hpp').read_text()
    body = '\n'.join(line for line in header.splitlines()
                     if not (line.startswith('#pragma once') or
                             line.startswith('#include "') or
                             line.startswith('#include <thrust/') or
                             line in ('#include <mpi.h>', '#include <cuda_runtime.h>')))
    scratch = ROOT / '.local-worktrees' / 'poiseuille-host'
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='vtu-', dir=str(scratch)) as temporary:
        work = Path(temporary)
        source, executable = work / 'writer.cpp', work / 'writer'
        source.write_text(STUBS + body + MAIN)
        subprocess.run(['c++', '-std=c++17', '-Wall', '-Wextra', '-Werror',
                        str(source), '-o', str(executable)], check=True)
        subprocess.run([str(executable), str(work)], check=True)
        collection = ET.parse(str(work / 'exact.pvd')).findall('.//DataSet')
        assert [float(item.get('timestep')) for item in collection] == [
            0.12345678901234567, 0.24691357802469134]
        for item in collection:
            master = ET.parse(str(work / item.get('file')))
            for array in master.findall('.//PDataArray'):
                assert array.get('type') == 'Float64'
            piece = ET.parse(str(work / master.find('.//Piece').get('Source')))
            expected = {'u': 0.12345678901234567, 'p': 123456789.01234567,
                        'velocity': 0.12345678901234567, 'cell': 1.2345678901234567}
            for group in ('Points', 'PointData', 'CellData'):
                for array in piece.findall('.//' + group + '/DataArray'):
                    assert array.get('type') == 'Float64'
                    value = expected.get(array.get('Name'), 1.000000000000003)
                    assert all(float(x) == value for x in array.text.split())
        master_path = ET.parse(str(work / 'legacy.pvd')).find('.//DataSet').get('file')
        assert ET.parse(str(work / master_path)).find('.//PPoints/PDataArray').get('type') == 'Float32'
    print('PASS: Float64 points, scalars, vectors, cell data and time round-trip; CUDA not tested')


if __name__ == '__main__':
    main()
