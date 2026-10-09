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

inline bool pressure_profile_matches(const mars::fem::pressure_settings::Values& expected,
                                    const mars::fem::pressure_settings::Values& actual,bool one_level) {
    for (const auto& [key,value]:expected) {
        const auto name=key=="relax_down"?"effective_relax_1":key=="relax_up"?"effective_relax_2":key;
        if (actual.at(name)!=value) return false;
    }
    return actual.at("effective_relax_3")==expected.at("coarserelax")
        && (one_level?actual.at("effective_levels")==1:actual.at("effective_levels")>1);
}
