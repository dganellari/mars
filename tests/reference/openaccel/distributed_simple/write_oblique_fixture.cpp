// Invented public geometry; no input mesh or case data is read.
#include "configured_case.hpp"
#include <netcdf.h>
#include <filesystem>
#include <iostream>

void nc_check(int status) { if (status!=NC_NOERR) throw std::runtime_error(nc_strerror(status)); }
int main(int argc,char** argv) {
    try {
        if (argc!=2) throw std::runtime_error("supply a new output .exo path");
        if (std::filesystem::exists(argv[1])) throw std::runtime_error("output already exists");
        auto mesh=dsimple_gate::channel(8,2,2);
        std::array<std::vector<int>,4> elements,sides;
        for (auto f:mesh.faces) {
            int patch=f.kind;
            if (patch==2) {
                bool zface=true;
                const double first=mesh.z[mesh.nodes[mars::segregated::tet_face_node(f.ordinal,0)][f.element]];
                for (int j=1;j<3;++j) zface=zface && mesh.z[mesh.nodes[mars::segregated::tet_face_node(f.ordinal,j)][f.element]]==first;
                if (zface) patch=3;
            }
            elements[patch].push_back(f.element+1); sides[patch].push_back(f.ordinal+1);
        }
        dsimple_gate::rotate_channel(mesh);
        int file; nc_check(nc_create(argv[1],NC_NOCLOBBER,&file));
        const float version=7.22f; const int word=8,large=1;
        const std::string title="Public synthetic oblique SIMPLE channel";
        nc_check(nc_put_att_text(file,NC_GLOBAL,"title",title.size(),title.c_str()));
        nc_check(nc_put_att_float(file,NC_GLOBAL,"api_version",NC_FLOAT,1,&version));
        nc_check(nc_put_att_float(file,NC_GLOBAL,"version",NC_FLOAT,1,&version));
        nc_check(nc_put_att_int(file,NC_GLOBAL,"floating_point_word_size",NC_INT,1,&word));
        nc_check(nc_put_att_int(file,NC_GLOBAL,"file_size",NC_INT,1,&large));
        auto dim=[&](const char* name,size_t n) { int d; nc_check(nc_def_dim(file,name,n,&d)); return d; };
        auto var=[&](const std::string& name,nc_type type,std::initializer_list<int> dims) {
            int v; nc_check(nc_def_var(file,name.c_str(),type,int(dims.size()),dims.begin(),&v)); return v;
        };
        dim("num_dim",3); const int n=dim("num_nodes",mesh.x.size());
        dim("num_elem",mesh.nodes[0].size()); const int block=dim("num_el_blk",1);
        const int e=dim("num_el_in_blk1",mesh.nodes[0].size()),four=dim("num_nod_per_el1",4);
        const int sets=dim("num_side_sets",4),width=dim("len_name",33);
        const int coordinates[3]={var("coordx",NC_DOUBLE,{n}),var("coordy",NC_DOUBLE,{n}),var("coordz",NC_DOUBLE,{n})};
        const int conn=var("connect1",NC_INT,{e,four}); nc_check(nc_put_att_text(file,conn,"elem_type",6,"TETRA4"));
        const int names=var("ss_names",NC_CHAR,{sets,width});
        const int bprop=var("eb_prop1",NC_INT,{block}),sprop=var("ss_prop1",NC_INT,{sets});
        nc_check(nc_put_att_text(file,bprop,"name",2,"ID")); nc_check(nc_put_att_text(file,sprop,"name",2,"ID"));
        const int bs=var("eb_status",NC_INT,{block}),ss=var("ss_status",NC_INT,{sets});
        int se[4],sf[4];
        for (int i=0;i<4;++i) {
            const auto suffix=std::to_string(i+1);
            const int count=dim(("num_side_ss"+suffix).c_str(),elements[i].size());
            se[i]=var("elem_ss"+suffix,NC_INT,{count}); sf[i]=var("side_ss"+suffix,NC_INT,{count});
        }
        nc_check(nc_enddef(file));
        nc_check(nc_put_var_double(file,coordinates[0],mesh.x.data()));
        nc_check(nc_put_var_double(file,coordinates[1],mesh.y.data()));
        nc_check(nc_put_var_double(file,coordinates[2],mesh.z.data()));
        std::vector<int> connectivity(4*mesh.nodes[0].size());
        for (size_t i=0;i<mesh.nodes[0].size();++i) for (int j=0;j<4;++j) connectivity[4*i+j]=mesh.nodes[j][i]+1;
        nc_check(nc_put_var_int(file,conn,connectivity.data()));
        const char labels[4][33]={"feed","exit","casing","cover"}; nc_check(nc_put_var_text(file,names,&labels[0][0]));
        const int ids[4]={1,2,3,4},ones[4]={1,1,1,1};
        nc_check(nc_put_var_int(file,bprop,ids)); nc_check(nc_put_var_int(file,sprop,ids));
        nc_check(nc_put_var_int(file,bs,ones)); nc_check(nc_put_var_int(file,ss,ones));
        for (int i=0;i<4;++i) { nc_check(nc_put_var_int(file,se[i],elements[i].data())); nc_check(nc_put_var_int(file,sf[i],sides[i].data())); }
        nc_check(nc_close(file));
        std::cout<<"Wrote public oblique channel: "<<argv[1]<<'\n';
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
