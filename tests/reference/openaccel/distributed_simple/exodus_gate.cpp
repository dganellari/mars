#include "mars_segregated_simple_native_mesh.hpp"
#include <filesystem>
#include <iostream>
using namespace mars::segregated::runtime;

void nc_check(int code) { if (code!=NC_NOERR) throw std::runtime_error(nc_strerror(code)); }
void write_fixture(const std::string& path,int variant) {
    int id; nc_check(nc_create(path.c_str(),NC_CLOBBER,&id));
    auto dimension=[&](const char* name,size_t count) { int d; nc_check(nc_def_dim(id,name,count,&d)); return d; };
    auto variable=[&](const char* name,nc_type type,std::initializer_list<int> dimensions) {
        int v; nc_check(nc_def_var(id,name,type,int(dimensions.size()),dimensions.begin(),&v)); return v;
    };
    const int xyz=dimension("num_dim",3),nodes=dimension("num_nodes",4),elements=dimension("num_elem",1);
    dimension("num_el_blk",1); dimension("num_el_in_blk1",1);
    const int corners=dimension("num_nod_per_el1",4),sets=dimension("num_side_sets",3),width=dimension("len_name",8);
    int coordinates[3]={-1,-1,-1};
    if (variant==1) coordinates[0]=variable("coord",NC_DOUBLE,{xyz,nodes});
    else {
        coordinates[0]=variable("coordx",NC_DOUBLE,{nodes}); coordinates[1]=variable("coordy",NC_DOUBLE,{nodes});
        if (variant!=4) coordinates[2]=variable("coordz",NC_DOUBLE,{nodes});
    }
    const int connectivity=variable("connect1",NC_INT,variant==2?std::initializer_list<int>{corners}:std::initializer_list<int>{elements,corners});
    nc_check(nc_put_att_text(id,connectivity,"elem_type",variant==3?0:4,"TET4"));
    const int names=variable("ss_names",NC_CHAR,{sets,width});
    int side_elements[3],sides[3];
    for (int i=0;i<3;++i) {
        const auto suffix=std::to_string(i+1);
        const int d=dimension(("num_side_ss"+suffix).c_str(),i==2?2:1);
        side_elements[i]=variable(("elem_ss"+suffix).c_str(),NC_INT,{d});
        sides[i]=variable(("side_ss"+suffix).c_str(),NC_INT,{d});
    }
    nc_check(nc_enddef(id));
    double coords[12]={0,1,0,0, 0,0,1,0, 0,0,0,1};
    if (variant==7) coords[0]=std::numeric_limits<double>::quiet_NaN();
    if (variant==1) nc_check(nc_put_var_double(id,coordinates[0],coords));
    else for (int j=0;j<3;++j) if (coordinates[j]>=0) nc_check(nc_put_var_double(id,coordinates[j],coords+4*j));
    int conn[4]={1,2,3,4}; if (variant==8) conn[3]=5; if (variant==9) conn[3]=3;
    nc_check(nc_put_var_int(id,connectivity,conn));
    const char labels[3][8]={"inlet","outlet","walls"}; nc_check(nc_put_var_text(id,names,&labels[0][0]));
    for (int j=0;j<3;++j) {
        int e[2]={1,1},f[2]={j+1,4};
        if (variant==5 && j==2) f[1]=1;
        if (variant==6 && j==2) f[1]=5;
        nc_check(nc_put_var_int(id,side_elements[j],e)); nc_check(nc_put_var_int(id,sides[j],f));
    }
    nc_check(nc_close(id));
}
template<class T> std::vector<T> file_gate_host(const Buffer<T>& values) {
#ifdef MARS_REPLAY_CUDA
    std::vector<T> out(values.size()); thrust::copy(values.begin(),values.end(),out.begin()); return out;
#else
    return values;
#endif
}
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv); int rank,ranks; MPI_Comm_rank(MPI_COMM_WORLD,&rank); MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try {
#ifdef MARS_REPLAY_CUDA
        int devices,local; MPI_Comm shared; MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&shared);
        MPI_Comm_rank(shared,&local); MPI_Comm_free(&shared);
        assembly_cuda_check(cudaGetDeviceCount(&devices)); ensure(devices>0,"no CUDA device"); assembly_cuda_check(cudaSetDevice(local%devices));
#endif
        ensure(argc>=2 && argc<=4,"supply a fixture directory and optional public channel mesh");
        if (!rank) std::filesystem::create_directories(argv[1]);
        MPI_Barrier(MPI_COMM_WORLD);
        const char* expected[]={"","","native Exodus read failed","native Exodus read failed","native Exodus read failed",
            "invalid or repeated boundary face","invalid or repeated boundary face","invalid source coordinates or connectivity",
            "invalid source coordinates or connectivity","degenerate source connectivity"};
        for (int variant=0;variant<10;++variant) {
            const auto path=std::string(argv[1])+"/fixture-"+std::to_string(variant)+".exo";
            if (!rank) write_fixture(path,variant);
            MPI_Barrier(MPI_COMM_WORLD);
            bool rejected=false;
            try {
                const auto mesh=read_simple_mesh(MPI_COMM_WORLD,path);
                bool correct=file_gate_host(mesh.x)==std::vector<double>{0,1,0,0}
                    && file_gate_host(mesh.y)==std::vector<double>{0,0,1,0}
                    && file_gate_host(mesh.z)==std::vector<double>{0,0,0,1};
                for (int j=0;j<4;++j) correct=correct && file_gate_host(mesh.nodes[j])==std::vector<int>{j};
                const auto faces=file_gate_host(mesh.faces); correct=correct && faces.size()==4;
                for (size_t j=0;j<faces.size();++j) correct=correct && faces[j].element==0 && faces[j].ordinal==int(j) && faces[j].kind==std::min(int(j),2);
                simple_collective(MPI_COMM_WORLD,correct,"Exodus arrays differ from file data");
            } catch (const std::exception& e) {
                rejected=variant>=2 && std::string(e.what()).find(expected[variant])!=std::string::npos;
                if (!rejected) throw;
            }
            simple_collective(MPI_COMM_WORLD,rejected==(variant>=2),"malformed Exodus fixture accepted");
        }
        if (argc>=3) {
            const auto channel=read_simple_mesh(MPI_COMM_WORLD,argv[2]);
            simple_collective(MPI_COMM_WORLD,channel.x.size()==425 && channel.nodes[0].size()==1536 && channel.faces.size()==576,
                              "public channel dimensions changed");
        }
        if (argc==4) {
            const mars::segregated::SimpleBoundaryNames selection{"feed","exit",{"casing","cover"}};
            const auto mesh=read_simple_mesh(MPI_COMM_WORLD,argv[3],selection);
            simple_collective(MPI_COMM_WORLD,mesh.x.size()==81 && mesh.nodes[0].size()==192 && mesh.faces.size()==144,
                              "oblique public mesh dimensions changed");
            for (auto bad:{mars::segregated::SimpleBoundaryNames{},
                          mars::segregated::SimpleBoundaryNames{"feed","exit",{"casing","missing"}},
                          mars::segregated::SimpleBoundaryNames{"feed","exit",{"casing","casing"}}}) {
                bool rejected=false;
                try { read_simple_mesh(MPI_COMM_WORLD,argv[3],bad); } catch (const std::exception&) { rejected=true; }
                simple_collective(MPI_COMM_WORLD,rejected,"invalid side-set mapping accepted");
            }
        }
        if (!rank) std::cout<<"PASS: native Exodus input, packed/split coordinates and 8 malformed cases ranks="<<ranks<<'\n';
    } catch (const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<'\n'; MPI_Abort(MPI_COMM_WORLD,1); }
    MPI_Finalize();
}
