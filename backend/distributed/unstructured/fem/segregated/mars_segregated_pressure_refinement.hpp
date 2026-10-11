#pragma once
#include <cmath>

namespace mars::segregated {
struct PressureRefinementResult {
    bool accepted=false;
    int rounds=0,iterations=0;
};
struct PressureCorrectionResult { bool accepted=false; int iterations=0; };

// Operations publish every candidate before evaluating owned rows and return
// collective decisions. The original equation and tolerance never change.
template<class Operations>
PressureRefinementResult refine_pressure(Operations& op,int used,int budget) {
    PressureRefinementResult result{false,0,used};
    auto current=op.defect();
    for (int round=0;round<3 && result.iterations>=0 && result.iterations<budget;++round) {
        if (!current.finite || !(current.residual2>0) || !std::isfinite(current.residual2)) break;
        const int remaining=budget-result.iterations;
        const auto correction=op.correct(remaining);
        ++result.rounds;
        if (correction.iterations<0 || correction.iterations>remaining) break;
        result.iterations+=correction.iterations;
        if (!correction.accepted || correction.iterations==0) break;
        const auto next=op.candidate();
        if (!next.finite || !std::isfinite(next.residual2) || !(next.residual2<current.residual2)) break;
        const bool accepted=next.passed && op.verify();
        op.keep();
        current=next;
        if (accepted) { result.accepted=true; return result; }
    }
    op.restore();
    return result;
}
} // namespace mars::segregated
