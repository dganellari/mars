#pragma once
#include <HYPRE_krylov.h>
#include <_hypre_parcsr_ls.h>
#if HYPRE_RELEASE_NUMBER >= 30000
#include <_hypre_krylov.h>
#else
#include <krylov.h>
#endif
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <istream>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>

namespace mars::fem::pressure_settings {
using Values=std::map<std::string,double>;
inline void checked(HYPRE_Int error) { if (error) throw std::runtime_error("private Hypre settings query failed"); }
// Matching installed internal headers are required, as in the defaults probe.
inline Values snapshot(HYPRE_Solver solver,HYPRE_Solver preconditioner,bool flex) {
    if (!solver || !preconditioner) throw std::runtime_error("pressure replay requires BoomerAMG");
    HYPRE_Int major,minor,patch,release;
    checked(HYPRE_VersionNumber(&major,&minor,&patch,&release));
    if (release!=HYPRE_RELEASE_NUMBER) throw std::runtime_error("Hypre header version mismatch");
    Values v{{"method",flex?1.:0.},{"hypre_release",double(release)}};
    auto integer=[&](const char* name,auto get) { HYPRE_Int value; checked(get(solver,&value)); v[name]=value; };
    auto real=[&](const char* name,auto get) { HYPRE_Real value; checked(get(solver,&value)); v[name]=value; };
    integer("kdim",flex?HYPRE_FlexGMRESGetKDim:HYPRE_GMRESGetKDim);
    integer("miniter",flex?HYPRE_FlexGMRESGetMinIter:HYPRE_GMRESGetMinIter);
    integer("maxiter",flex?HYPRE_FlexGMRESGetMaxIter:HYPRE_GMRESGetMaxIter);
    real("rtol",flex?HYPRE_FlexGMRESGetTol:HYPRE_GMRESGetTol);
    HYPRE_Real absolute;
    if (flex) checked(hypre_FlexGMRESGetAbsoluteTol(solver,&absolute));
    else checked(HYPRE_GMRESGetAbsoluteTol(solver,&absolute));
    v["atol"]=absolute;
    auto* a=reinterpret_cast<hypre_ParAMGData*>(preconditioner);
#define MARS_PRESSURE_SETTING(key,field) v[key]=hypre_ParAMGData##field(a)
    MARS_PRESSURE_SETTING("coarsentype",CoarsenType);
    MARS_PRESSURE_SETTING("interptype",InterpType);
    MARS_PRESSURE_SETTING("relaxorder",RelaxOrder);
    MARS_PRESSURE_SETTING("pmax",PMaxElmts);
    MARS_PRESSURE_SETTING("maxlevels",MaxLevels);
    MARS_PRESSURE_SETTING("mincoarsesize",MinCoarseSize);
    MARS_PRESSURE_SETTING("maxcoarsesize",MaxCoarseSize);
    MARS_PRESSURE_SETTING("numfunctions",NumFunctions);
    MARS_PRESSURE_SETTING("coarsencutfactor",CoarsenCutFactor);
    MARS_PRESSURE_SETTING("cycletype",CycleType);
    MARS_PRESSURE_SETTING("fcycle",FCycle);
    MARS_PRESSURE_SETTING("strongthreshold",StrongThreshold);
    MARS_PRESSURE_SETTING("truncfactor",TruncFactor);
    MARS_PRESSURE_SETTING("jacobitruncthreshold",JacobiTruncThreshold);
    MARS_PRESSURE_SETTING("maxrowsum",MaxRowSum);
    MARS_PRESSURE_SETTING("aggnumlevels",AggNumLevels);
    MARS_PRESSURE_SETTING("agginterptype",AggInterpType);
    MARS_PRESSURE_SETTING("aggtruncfactor",AggTruncFactor);
    MARS_PRESSURE_SETTING("numpaths",NumPaths);
    MARS_PRESSURE_SETTING("keeptranspose",KeepTranspose);
    MARS_PRESSURE_SETTING("nodal",Nodal);
    MARS_PRESSURE_SETTING("nodaldiag",NodalDiag);
    MARS_PRESSURE_SETTING("relaxtype",UserRelaxType);
    MARS_PRESSURE_SETTING("coarserelax",UserCoarseRelaxType);
    MARS_PRESSURE_SETTING("effective_levels",NumLevels);
#undef MARS_PRESSURE_SETTING
    v["numsweeps"]=hypre_ParAMGDataNumGridSweeps(a)[0];
    for (int cycle=1;cycle<=3;++cycle) {
        HYPRE_Int relax,sweeps;
        checked(hypre_BoomerAMGGetCycleRelaxType(preconditioner,&relax,cycle));
        checked(hypre_BoomerAMGGetCycleNumSweeps(preconditioner,&sweeps,cycle));
        if (relax!=hypre_ParAMGDataGridRelaxType(a)[cycle] || sweeps!=hypre_ParAMGDataNumGridSweeps(a)[cycle])
            throw std::runtime_error("Hypre settings layout mismatch");
        const auto suffix=std::to_string(cycle);
        v["effective_relax_"+suffix]=relax; v["sweeps_"+suffix]=sweeps;
    }
    return v;
}
inline void write(std::ostream& out,const Values& values) {
    out<<std::setprecision(17);
    for (const auto& value:values) out<<value.first<<' '<<value.second<<'\n';
}
inline Values read(std::istream& in) {
    Values values; std::string key; double value;
    while (in>>key) {
        if (!(in>>value) || !std::isfinite(value) || !values.emplace(key,value).second)
            throw std::runtime_error("invalid private solver settings");
    }
    if (!in.eof()) throw std::runtime_error("invalid private solver settings");
    return values;
}
inline void apply(const Values& v,HYPRE_Solver solver,HYPRE_Solver amg,bool flex) {
    auto i=[](double value) { if (!std::isfinite(value) || value!=std::floor(value) || value<INT32_MIN || value>INT32_MAX)
        throw std::runtime_error("invalid integer setting"); return HYPRE_Int(value); };
    for (const auto& [key,value]:v) {
        if (key=="method" || key=="hypre_release" || key.rfind("effective_",0)==0) continue;
        if (key=="kdim") checked((flex?HYPRE_FlexGMRESSetKDim:HYPRE_GMRESSetKDim)(solver,i(value)));
        else if (key=="miniter") checked((flex?HYPRE_FlexGMRESSetMinIter:HYPRE_GMRESSetMinIter)(solver,i(value)));
        else if (key=="maxiter") checked((flex?HYPRE_FlexGMRESSetMaxIter:HYPRE_GMRESSetMaxIter)(solver,i(value)));
        else if (key=="rtol") checked((flex?HYPRE_FlexGMRESSetTol:HYPRE_GMRESSetTol)(solver,value));
        else if (key=="atol") checked((flex?HYPRE_FlexGMRESSetAbsoluteTol:HYPRE_GMRESSetAbsoluteTol)(solver,value));
#define MARS_PRESSURE_INT(name,setter) else if (key==name) checked(HYPRE_BoomerAMGSet##setter(amg,i(value)))
#define MARS_PRESSURE_REAL(name,setter) else if (key==name) checked(HYPRE_BoomerAMGSet##setter(amg,value))
        MARS_PRESSURE_INT("coarsentype",CoarsenType); MARS_PRESSURE_INT("interptype",InterpType);
        MARS_PRESSURE_INT("relaxorder",RelaxOrder); MARS_PRESSURE_INT("pmax",PMaxElmts);
        MARS_PRESSURE_INT("maxlevels",MaxLevels); MARS_PRESSURE_INT("mincoarsesize",MinCoarseSize);
        MARS_PRESSURE_INT("maxcoarsesize",MaxCoarseSize); MARS_PRESSURE_INT("numfunctions",NumFunctions);
        MARS_PRESSURE_INT("coarsencutfactor",CoarsenCutFactor); MARS_PRESSURE_INT("cycletype",CycleType);
        MARS_PRESSURE_INT("fcycle",FCycle); MARS_PRESSURE_INT("aggnumlevels",AggNumLevels);
        MARS_PRESSURE_INT("agginterptype",AggInterpType); MARS_PRESSURE_INT("numpaths",NumPaths);
        MARS_PRESSURE_INT("keeptranspose",KeepTranspose); MARS_PRESSURE_INT("nodal",Nodal);
        MARS_PRESSURE_INT("nodaldiag",NodalDiag);
        MARS_PRESSURE_REAL("strongthreshold",StrongThreshold); MARS_PRESSURE_REAL("truncfactor",TruncFactor);
        MARS_PRESSURE_REAL("jacobitruncthreshold",JacobiTruncThreshold); MARS_PRESSURE_REAL("maxrowsum",MaxRowSum);
        MARS_PRESSURE_REAL("aggtruncfactor",AggTruncFactor);
#undef MARS_PRESSURE_INT
#undef MARS_PRESSURE_REAL
        // Apply these in a fixed order below: generic setters overwrite cycle settings.
        else if (key!="relaxtype" && key!="coarserelax" && key!="numsweeps" && key!="sweeps_1" && key!="sweeps_2" && key!="sweeps_3")
            throw std::runtime_error("unsupported private solver setting");
    }
    if (v.count("relaxtype") && v.at("relaxtype")>=0) checked(HYPRE_BoomerAMGSetRelaxType(amg,i(v.at("relaxtype"))));
    if (v.count("coarserelax") && v.at("coarserelax")>=0) checked(HYPRE_BoomerAMGSetCycleRelaxType(amg,i(v.at("coarserelax")),3));
    if (v.count("numsweeps")) checked(HYPRE_BoomerAMGSetNumSweeps(amg,i(v.at("numsweeps"))));
    for (int cycle=1;cycle<=3;++cycle) {
        const auto key="sweeps_"+std::to_string(cycle);
        if (v.count(key)) checked(HYPRE_BoomerAMGSetCycleNumSweeps(amg,i(v.at(key)),cycle));
    }
    checked(HYPRE_BoomerAMGSetMaxIter(amg,1)); checked(HYPRE_BoomerAMGSetTol(amg,0));
}
}
