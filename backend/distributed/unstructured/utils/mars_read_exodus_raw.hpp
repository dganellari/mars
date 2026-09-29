#pragma once
// Read file arrays without constructing CPU topology or local/global node maps.
#include <array>
#include <climits>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef MARS_HAVE_NETCDF
#include <netcdf.h>
#endif
namespace mars {
struct ExodusRawSideSet {
    std::string name;
    std::vector<long long> elements,sides; // Exodus one-based file indices
    long long id=0; // Zero means the optional ID array was absent; never infer an ID from sequence.
};
struct ExodusRawTet4 {
    std::array<std::vector<double>,3> coordinates;
    std::vector<long long> connectivity; // element-major, one-based file indices
    std::vector<ExodusRawSideSet> side_sets;
};
inline ExodusRawTet4 readExodusTet4Raw(const std::string& path) {
#ifndef MARS_HAVE_NETCDF
    (void)path;
    throw std::runtime_error("native SIMPLE Exodus input requires netCDF support");
#else
    auto check=[](int result) {
        if (result!=NC_NOERR) throw std::runtime_error(std::string("NetCDF: ")+nc_strerror(result));
    };
    struct File {
        int id=-1;
        ~File() { if (id>=0) nc_close(id); }
    } file;
    check(nc_open(path.c_str(),NC_NOWRITE,&file.id));
    auto dim=[&](const std::string& name) {
        int id; size_t size;
        check(nc_inq_dimid(file.id,name.c_str(),&id)); check(nc_inq_dimlen(file.id,id,&size));
        return size;
    };
    auto variable=[&](const std::string& name,const std::vector<size_t>& shape) {
        int id,rank; int dims[NC_MAX_VAR_DIMS];
        check(nc_inq_varid(file.id,name.c_str(),&id)); check(nc_inq_varndims(file.id,id,&rank));
        if (size_t(rank)!=shape.size()) throw std::runtime_error("invalid Exodus variable rank: "+name);
        check(nc_inq_vardimid(file.id,id,dims));
        for (int j=0;j<rank;++j) { size_t size; check(nc_inq_dimlen(file.id,dims[j],&size));
            if (size!=shape[size_t(j)]) throw std::runtime_error("invalid Exodus variable shape: "+name);
        }
        return id;
    };
    const size_t n=dim("num_nodes"),e=dim("num_elem");
    if (dim("num_dim")!=3 || dim("num_el_blk")!=1 || dim("num_nod_per_el1")!=4 || dim("num_el_in_blk1")!=e)
        throw std::runtime_error("native SIMPLE currently requires one three-dimensional Tet4 block");
    if (!n || !e || n>size_t(INT_MAX/9) || e>size_t(INT_MAX/16))
        throw std::runtime_error("Exodus mesh exceeds SIMPLE local index capacity");
    const int conn=variable("connect1",{e,4});
    size_t type_length=0; check(nc_inq_attlen(file.id,conn,"elem_type",&type_length));
    if (!type_length || type_length>64) throw std::runtime_error("invalid Exodus element type");
    std::string type(type_length,'\0'); check(nc_get_att_text(file.id,conn,"elem_type",&type[0]));
    while (!type.empty() && (type.back()=='\0' || type.back()==' ')) type.pop_back();
    if (type!="TETRA" && type!="TETRA4" && type!="TET4") throw std::runtime_error("Exodus element block is not Tet4");
    ExodusRawTet4 out; out.connectivity.resize(4*e);
    check(nc_get_var_longlong(file.id,conn,out.connectivity.data()));
    int split;
    if (nc_inq_varid(file.id,"coordx",&split)==NC_NOERR) {
        const char* names[]={"coordx","coordy","coordz"};
        for (int j=0;j<3;++j) { out.coordinates[j].resize(n);
            check(nc_get_var_double(file.id,variable(names[j],{n}),out.coordinates[j].data()));
        }
    } else {
        const int coord=variable("coord",{3,n}); const size_t count[]={1,n};
        for (int j=0;j<3;++j) { out.coordinates[j].resize(n); const size_t start[]={size_t(j),0};
            check(nc_get_vara_double(file.id,coord,start,count,out.coordinates[j].data()));
        }
    }
    const size_t sets=dim("num_side_sets");
    if (sets>size_t(INT_MAX)) throw std::runtime_error("too many Exodus side sets");
    int name_id,rank,dims[NC_MAX_VAR_DIMS]; size_t width=0;
    const int name_status=nc_inq_varid(file.id,"ss_names",&name_id);
    std::vector<char> names;
    if (name_status!=NC_ENOTVAR) {
        check(name_status); check(nc_inq_varndims(file.id,name_id,&rank));
        if (rank!=2) throw std::runtime_error("invalid Exodus side-set name array");
        check(nc_inq_vardimid(file.id,name_id,dims)); check(nc_inq_dimlen(file.id,dims[1],&width));
        if (!width || width>4096 || sets>size_t(INT_MAX)/width) throw std::runtime_error("invalid Exodus side-set names");
        names.resize(sets*width); check(nc_get_var_text(file.id,variable("ss_names",{sets,width}),names.data()));
    }
    std::vector<long long> ids(sets,0);
    int id_variable; const int id_status=nc_inq_varid(file.id,"ss_prop1",&id_variable);
    if (id_status!=NC_ENOTVAR) {
        check(id_status); check(nc_get_var_longlong(file.id,variable("ss_prop1",{sets}),ids.data()));
        for (auto id:ids) if (id<=0) throw std::runtime_error("invalid Exodus side-set ID");
    }
    out.side_sets.resize(sets);
    for (size_t i=0;i<sets;++i) {
        auto& set=out.side_sets[i]; size_t len=0;
        while (len<width && names[i*width+len]) ++len;
        if (width) set.name.assign(names.data()+i*width,len);
        while (!set.name.empty() && set.name.back()==' ') set.name.pop_back();
        set.id=ids[i];
        const std::string index=std::to_string(i+1);
        const size_t count=dim("num_side_ss"+index);
        if (count>size_t(INT_MAX)) throw std::runtime_error("side set exceeds SIMPLE index capacity");
        set.elements.resize(count); set.sides.resize(count);
        check(nc_get_var_longlong(file.id,variable("elem_ss"+index,{count}),set.elements.data()));
        check(nc_get_var_longlong(file.id,variable("side_ss"+index,{count}),set.sides.data()));
    }
    return out;
#endif
}
} // namespace mars
