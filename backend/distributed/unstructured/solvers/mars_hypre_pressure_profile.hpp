#pragma once
#include "mars_hypre_pressure_settings.hpp"
#include <set>

namespace mars::fem::pressure_settings {
// Reference AMG permits Hypre's host staging for HMIS and the direct coarse solve.
inline void validate_gpu_profile(const Values& v,bool reference_amg=false) {
#if HYPRE_RELEASE_NUMBER < 23300
    if (reference_amg) throw std::runtime_error("reference AMG requires Hypre 2.33 or newer");
#endif
    const std::set<std::string> integers={"method","kdim","miniter","maxiter","coarsentype","interptype",
        "relaxorder","pmax","maxlevels","mincoarsesize","maxcoarsesize","numfunctions","coarsencutfactor",
        "cycletype","fcycle","aggnumlevels","agginterptype","numpaths","keeptranspose","nodal","nodaldiag",
        "relaxtype","coarserelax","numsweeps","relax_down","relax_up","sweeps_1","sweeps_2","sweeps_3"};
    auto keys=integers;
    for (const auto* key:{"rtol","atol","strongthreshold","truncfactor","jacobitruncthreshold","maxrowsum","aggtruncfactor"}) keys.insert(key);
    bool valid=v.size()==keys.size();
    for (const auto& [key,value]:v)
        valid=valid && keys.count(key) && std::isfinite(value)
            && (!integers.count(key) || (value==std::floor(value) && value>=0 && value<=INT32_MAX));
    if (!valid) throw std::runtime_error("invalid pressure profile fields");
    auto get=[&](const char* key) { return v.at(key); };
    auto smoother=[](double value) { return value==6 || value==13 || value==14 || value==18; };
    valid=(get("method")==0 || get("method")==1) && get("kdim")>0 && get("maxiter")>0
        && get("miniter")<=get("maxiter") && get("rtol")>0 && get("rtol")<1 && get("atol")>=0
        && get("coarsentype")==(reference_amg?10:8) && get("interptype")==6 && get("relaxorder")==0
        && get("aggnumlevels")==0 && get("nodal")==0 && get("nodaldiag")==0 && get("numfunctions")==1
        && get("cycletype")==1 && get("fcycle")==0 && get("keeptranspose")<=1
        && get("maxlevels")>0 && get("maxcoarsesize")>0 && get("mincoarsesize")<=get("maxcoarsesize")
        && get("numpaths")>0 && get("coarserelax")==(reference_amg?9:18)
        && smoother(get("relaxtype")) && smoother(get("relax_down")) && smoother(get("relax_up"))
        && get("numsweeps")>0 && get("sweeps_1")>0 && get("sweeps_2")>0 && get("sweeps_3")>0;
    for (const auto* key:{"strongthreshold","truncfactor","jacobitruncthreshold","maxrowsum","aggtruncfactor"})
        valid=valid && get(key)>=0 && get(key)<=1;
    if (!valid) throw std::runtime_error("unsupported pressure profile for the selected AMG mode");
}

inline void apply_gpu_amg_profile(const Values& profile,HYPRE_Solver amg,bool reference_amg=false) {
    validate_gpu_profile(profile,reference_amg);
    auto v=profile;
    for (const auto* key:{"method","kdim","miniter","maxiter","rtol","atol","relax_down","relax_up"}) v.erase(key);
    apply(v,nullptr,amg,false);
    checked(HYPRE_BoomerAMGSetCycleRelaxType(amg,int(profile.at("relax_down")),1));
    checked(HYPRE_BoomerAMGSetCycleRelaxType(amg,int(profile.at("relax_up")),2));
}

template<class Solver> void configure_gpu_profile(Solver& solver,const Values& profile,bool reference_amg=false) {
    validate_gpu_profile(profile,reference_amg);
    solver.set_krylov_controls(profile.at("method")!=0,int(profile.at("kdim")),
                               int(profile.at("miniter")),int(profile.at("maxiter")));
    solver.set_stopping_tolerances(profile.at("rtol"),profile.at("atol"));
    solver.enable_true_residual_check(profile.at("atol"),profile.at("rtol"),true);
    solver.set_amg_controls([profile,reference_amg](HYPRE_Solver amg) { apply_gpu_amg_profile(profile,amg,reference_amg); });
}
}
