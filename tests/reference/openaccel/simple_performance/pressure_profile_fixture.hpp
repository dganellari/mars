#pragma once
#include "../../../../backend/distributed/unstructured/solvers/mars_hypre_pressure_profile.hpp"

inline mars::fem::pressure_settings::Values pressure_profile_fixture(bool flex=false,bool one_level=false) {
    return {{"method",flex?1.:0.},{"kdim",5},{"miniter",0},{"maxiter",200},{"rtol",1e-10},{"atol",1e-13},
        {"coarsentype",8},{"interptype",6},{"relaxorder",0},{"pmax",4},{"maxlevels",one_level?1.:25.},
        {"mincoarsesize",0},{"maxcoarsesize",9},{"numfunctions",1},{"coarsencutfactor",0},{"cycletype",1},
        {"fcycle",0},{"strongthreshold",.25},{"truncfactor",0},{"jacobitruncthreshold",.01},{"maxrowsum",.9},
        {"aggnumlevels",0},{"agginterptype",4},{"aggtruncfactor",0},{"numpaths",1},{"keeptranspose",0},
        {"nodal",0},{"nodaldiag",0},{"relaxtype",6},{"coarserelax",18},{"numsweeps",1},
        {"relax_down",13},{"relax_up",14},{"sweeps_1",1},{"sweeps_2",1},{"sweeps_3",1}};
}

inline bool pressure_profile_controls_match(const mars::fem::pressure_settings::Values& expected,
                                           const mars::fem::pressure_settings::Values& actual,
                                           std::ostream& detail) {
    bool matches=true;
    auto compare=[&](const std::string& name,double value) {
        const auto found=actual.find(name);
        if (found!=actual.end() && found->second==value) return;
        matches=false;
        detail<<name<<" expected="<<value<<" actual=";
        if (found==actual.end()) detail<<"missing";
        else detail<<found->second;
        detail<<"; ";
    };
    for (const auto& [key,value]:expected) {
        const auto name=key=="relax_down"?"effective_relax_1":key=="relax_up"?"effective_relax_2":key;
        compare(name,value);
    }
    compare("effective_relax_3",expected.at("coarserelax"));
    return matches;
}

inline bool pressure_profile_levels_match(const mars::fem::pressure_settings::Values& actual,bool one_level) {
    const auto found=actual.find("effective_levels");
    if (found==actual.end() || !std::isfinite(found->second) || found->second!=std::floor(found->second)) return false;
    return one_level?found->second==1:found->second>1;
}
