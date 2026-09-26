#pragma once
#include "mars_segregated_simple_input.hpp"
#include "mars_segregated_native_mapping.hpp"
#include "mars_segregated_simple_runtime.hpp"

#ifdef MARS_REPLAY_CUDA
#include <thrust/adjacent_find.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
namespace mars::segregated::runtime {
inline SimpleInput load_simple_exodus(const std::string& path) {
#ifndef MARS_HAVE_NETCDF
    (void)path;
    throw std::runtime_error("native SIMPLE Exodus input requires netCDF support");
#else
    // Bound the existing reader to this public, single-block Tet4 profile before allocation.
    auto check_netcdf=[](int status) {
        if (status!=NC_NOERR) throw std::runtime_error(std::string("NetCDF error: ")+nc_strerror(status));
    };
    int file=-1;
    check_netcdf(nc_open(path.c_str(),NC_NOWRITE,&file));
    try {
        auto dimension=[&](const char* name,size_t expected) {
            int id; size_t size;
            check_netcdf(nc_inq_dimid(file,name,&id)); check_netcdf(nc_inq_dimlen(file,id,&size));
            ensure(size==expected,"unsupported native SIMPLE public mesh dimensions");
        };
        dimension("num_nodes",425); dimension("num_elem",1536); dimension("num_el_blk",1);
        dimension("num_el_in_blk1",1536); dimension("num_nod_per_el1",4);
        int var,rank,dimensions[NC_MAX_VAR_DIMS];
        check_netcdf(nc_inq_varid(file,"connect1",&var)); check_netcdf(nc_inq_varndims(file,var,&rank));
        ensure(rank==2,"invalid Tet4 connectivity dimensions");
        check_netcdf(nc_inq_vardimid(file,var,dimensions));
        size_t rows,columns;
        check_netcdf(nc_inq_dimlen(file,dimensions[0],&rows)); check_netcdf(nc_inq_dimlen(file,dimensions[1],&columns));
        ensure(rows==1536 && columns==4,"invalid Tet4 connectivity shape");
        std::vector<int> connectivity(rows*columns);
        check_netcdf(nc_get_var_int(file,var,connectivity.data()));
        for (int node:connectivity) ensure(node>=1 && node<=425,"invalid Exodus node index");
        for (const char* name:{"coordx","coordy","coordz"}) {
            check_netcdf(nc_inq_varid(file,name,&var)); check_netcdf(nc_inq_varndims(file,var,&rank));
            ensure(rank==1,"native input requires split Exodus coordinates");
            check_netcdf(nc_inq_vardimid(file,var,dimensions)); check_netcdf(nc_inq_dimlen(file,dimensions[0],&rows));
            ensure(rows==425,"invalid Exodus coordinate shape");
        }
    } catch (...) { nc_close(file); throw; }
    check_netcdf(nc_close(file));
    auto [n,e,x,y,z,conn,ids,unused]=mars::readExodusMeshWithElementPartitioning<4,double,uint64_t>(path,0,1);
    ensure(n==425 && e==1536,"native SIMPLE currently accepts only the public 425-node channel");
    SimpleInput input; input.x=std::move(x); input.y=std::move(y); input.z=std::move(z);
    const std::vector<uint64_t>* columns[]={&std::get<0>(conn),&std::get<1>(conn),&std::get<2>(conn),&std::get<3>(conn)};
    for (int j=0;j<4;++j) for (auto node:*columns[j]) {
        ensure(node<n,"invalid Exodus connectivity"); input.nodes[j].push_back(int(node));
    }
    std::map<std::array<double,3>,int> coordinates;
    for (int i=0;i<int(n);++i) {
        ensure(ids[i]==uint64_t(i),"public input must use every Exodus node");
        ensure(std::isfinite(input.x[i]) && std::isfinite(input.y[i]) && std::isfinite(input.z[i]),"nonfinite Exodus coordinates");
        ensure(coordinates.emplace(std::array<double,3>{input.x[i],input.y[i],input.z[i]},i).second,"coincident source nodes");
    }
    using FaceKey=std::array<int,3>;
    std::map<FaceKey,std::vector<SimpleFace>> incidence;
    for (int i=0;i<int(e);++i) {
        int cell[4]; for (int j=0;j<4;++j) cell[j]=input.nodes[j][i];
        auto key=simple_cell_key(cell);
        for (int j=1;j<4;++j) ensure(key.nodes[j]!=key.nodes[j-1],"degenerate source cell");
        for (int f=0;f<4;++f) {
            FaceKey face; for (int j=0;j<3;++j) face[j]=cell[tet_face_node(f,j)];
            std::sort(face.begin(),face.end()); incidence[face].push_back({i,f,-1});
        }
    }
    int exterior=0;
    for (const auto& item:incidence) { ensure(item.second.size()<=2,"nonmanifold source topology"); exterior+=item.second.size()==1; }
    ensure(exterior==576,"wrong public exterior coverage");
    const auto side_sets=mars::readExodusSideSetsTet4(path,0);
    ensure(side_sets.triangleCoordsByName.size()==3,"expected inlet, outlet and walls side sets");
    const char* names[]={"inlet","outlet","walls"}; const int counts[]={32,32,512};
    std::set<FaceKey> tagged;
    for (int kind=0;kind<3;++kind) {
        const auto found=side_sets.triangleCoordsByName.find(names[kind]);
        ensure(found!=side_sets.triangleCoordsByName.end(),"missing public side set");
        const auto& triangles=found->second;
        ensure(triangles.size()==size_t(3*counts[kind]),"wrong public side-set coverage");
        for (size_t i=0;i<triangles.size();i+=3) {
            FaceKey face;
            for (int j=0;j<3;++j) {
                auto node=coordinates.find(triangles[i+j]); ensure(node!=coordinates.end(),"unresolved side-set coordinate");
                face[j]=node->second;
            }
            std::sort(face.begin(),face.end()); auto cell=incidence.find(face);
            ensure(cell!=incidence.end() && cell->second.size()==1,"side set is not an exterior face");
            ensure(tagged.insert(face).second,"duplicate side-set face");
            auto boundary=cell->second[0]; boundary.kind=kind; input.faces.push_back(boundary);
        }
    }
    return input;
#endif
}

using SimpleDomain=mars::ElementDomain<mars::TetTag,double,uint64_t,cstone::GpuTag>;
struct NativeNodeMap {
    const double *sx,*sy,*sz; const uint64_t* keys; int count;
    cstone::Box<double> box;
    double *x,*y,*z; int *source_to_local,*source_node,*hits,*error;
    __device__ void operator()(int i) const {
        const auto key=cstone::sfc3D<cstone::HilbertKey<uint64_t>>(sx[i],sy[i],sz[i],box).value();
        int local=simple_find_key(keys,count,key);
        if (local<0) { atomicExch(error,1); return; }
        if (atomicAdd(hits+local,1)!=0) { atomicExch(error,1); return; }
        source_to_local[i]=local; source_node[local]=i;
        x[local]=sx[i]; y[local]=sy[i]; z[local]=sz[i];
    }
};
struct NativeCellKeys {
    const int *nodes[4]; const int* map; SimpleCellKey* keys;
    __device__ void operator()(int e) const {
        int cell[4]; for (int j=0;j<4;++j) cell[j]=map[nodes[j][e]];
        keys[e]=simple_cell_key(cell);
    }
};
struct NativeCellMap {
    const uint64_t *nodes[4]; int *out[4];
    const SimpleCellKey* keys; const int* source_elements;
    int node_count,element_count; int *source_to_local,*hits,*error;
    __device__ void operator()(int e) const {
        int cell[4];
        for (int j=0;j<4;++j) {
            if (nodes[j][e]>=uint64_t(node_count)) { atomicExch(error,1); return; }
            cell[j]=int(nodes[j][e]); out[j][e]=cell[j];
        }
        int found=simple_find_key(keys,element_count,simple_cell_key(cell));
        if (found<0) { atomicExch(error,1); return; }
        int source=source_elements[found];
        if (atomicAdd(hits+source,1)!=0) { atomicExch(error,1); return; }
        source_to_local[source]=e;
    }
};
struct NativeFaceMap {
    const SimpleFace* source; SimpleFace* faces;
    const int *source_nodes[4],*nodes[4],*node_map,*element_map;
    int* error;
    __device__ void operator()(int i) const {
        const auto face=source[i]; const int e=element_map[face.element];
        int triangle[3],cell[4];
        for (int j=0;j<3;++j) triangle[j]=node_map[source_nodes[tet_face_node(face.ordinal,j)][face.element]];
        for (int j=0;j<4;++j) cell[j]=nodes[j][e];
        const int ordinal=simple_native_face(triangle,cell);
        if (ordinal<0) { atomicExch(error,1); return; }
        faces[i]={e,ordinal,face.kind};
    }
};
struct NativeSimpleInput {
    thrust::device_vector<double> x,y,z;
    std::array<thrust::device_vector<int>,4> nodes;
    thrust::device_vector<SimpleFace> faces;
    thrust::device_vector<int> source_node;
    explicit NativeSimpleInput(const std::string& path) {
        const auto source=load_simple_exodus(path);
        SimpleDomain::HostCoordsTuple coords{source.x,source.y,source.z};
        SimpleDomain::HostConnectivityTuple connectivity;
        std::get<0>(connectivity).assign(source.nodes[0].begin(),source.nodes[0].end());
        std::get<1>(connectivity).assign(source.nodes[1].begin(),source.nodes[1].end());
        std::get<2>(connectivity).assign(source.nodes[2].begin(),source.nodes[2].end());
        std::get<3>(connectivity).assign(source.nodes[3].begin(),source.nodes[3].end());
        SimpleDomain domain(coords,connectivity,0,1);
        const auto& conn=domain.getElementToNodeConnectivity();
        const auto& keys=domain.getLocalToGlobalSfcMap();
        const int n=int(source.x.size()),e=int(source.nodes[0].size());
        ensure(domain.getNodeCount()==size_t(n) && domain.getElementCount()==size_t(e)
               && domain.localElementCount()==size_t(e) && domain.startIndex()==0 && keys.size()==size_t(n),
               "ElementDomain changed public topology counts or ownership");
        x.resize(n); y.resize(n); z.resize(n); source_node.resize(n); faces.resize(source.faces.size());
        for (auto& column:nodes) column.resize(e);
        Array<double> sx(source.x),sy(source.y),sz(source.z);
        Array<int> node_map(n),node_hits(n),element_map(e),element_hits(e),error(1);
        // The original doubles, not decoded SFC cell centres, define the PDE geometry.
        launch(n,NativeNodeMap{sx.data(),sy.data(),sz.data(),keys.data(),n,domain.getBoundingBox(),
            device_data(x),device_data(y),device_data(z),node_map.data(),device_data(source_node),node_hits.data(),error.data()});
        ensure(error.host()[0]==0,"native node mapping is not a bijection");
        Array<int> s0(source.nodes[0]),s1(source.nodes[1]),s2(source.nodes[2]),s3(source.nodes[3]);
        thrust::device_vector<SimpleCellKey> cell_keys(e); thrust::device_vector<int> cell_ids(e);
        thrust::sequence(cell_ids.begin(),cell_ids.end());
        launch(e,NativeCellKeys{{s0.data(),s1.data(),s2.data(),s3.data()},node_map.data(),device_data(cell_keys)});
        thrust::sort_by_key(cell_keys.begin(),cell_keys.end(),cell_ids.begin());
        ensure(thrust::adjacent_find(cell_keys.begin(),cell_keys.end())==cell_keys.end(),"duplicate source cells");
        launch(e,NativeCellMap{{std::get<0>(conn).data(),std::get<1>(conn).data(),std::get<2>(conn).data(),std::get<3>(conn).data()},
            {device_data(nodes[0]),device_data(nodes[1]),device_data(nodes[2]),device_data(nodes[3])},
            device_data(cell_keys),device_data(cell_ids),n,e,element_map.data(),element_hits.data(),error.data()});
        ensure(error.host()[0]==0,"native cell mapping is not a bijection");
        Array<SimpleFace> source_faces(source.faces);
        launch(int(faces.size()),NativeFaceMap{source_faces.data(),device_data(faces),
            {s0.data(),s1.data(),s2.data(),s3.data()},
            {device_data(nodes[0]),device_data(nodes[1]),device_data(nodes[2]),device_data(nodes[3])},
            node_map.data(),element_map.data(),error.data()});
        ensure(error.host()[0]==0,"native boundary face mapping failed");
        assembly_cuda_check(cudaDeviceSynchronize());
        std::cout<<"PASS: ElementDomain node/cell bijections and all 576 tagged exterior faces; original coordinates retained\n";
    }
};
}
#endif
