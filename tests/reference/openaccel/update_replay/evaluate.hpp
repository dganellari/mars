#pragma once
#include "mars_segregated_update.hpp"

#if defined(__CUDACC__)
#define UPDATE_HD __host__ __device__
#else
#define UPDATE_HD
#endif
namespace update_replay {
inline constexpr int input_width[] = {3,1,10,24,9,25,31,10,6,2};
inline constexpr int output_width[] = {1,1,3,1,1,1,6,2,1,1};
struct Input { int stage = 0; double values[31]{}; };
struct Output { double values[6]{}; };
inline bool valid(const Input& in)
{
    const double* x=in.values;
    auto flag=[](double v){return v==0 || v==1;};
    auto alpha=[](double v){return v>0 && v<=1;};
    switch(in.stage) {
    case 0: return alpha(x[2]);
    case 1: return true;
    case 2: return flag(x[9]);
    case 3: return x[0]>0 && alpha(x[23]);
    case 4: return x[0]>0 && alpha(x[8]);
    case 5: return x[0]>0 && alpha(x[23]) && flag(x[24]);
    case 6: {
        double squared=0;
        for(int j=0;j<3;++j) {const double sum=x[21+j]+x[24+j]+x[27+j];squared+=sum*sum;}
        return flag(x[3]) && x[3]==x[4] && x[3]==x[5] && flag(x[30]) && squared>0;
    }
    case 7: return flag(x[9]) && (x[9]==1 || (x[6]*x[6]+x[7]*x[7]+x[8]*x[8]>0));
    case 8: return x[3]>=0 && x[3]<=1 && flag(x[5]);
    case 9: return x[1]>0;
    default: return false;
    }
}
UPDATE_HD inline void evaluate(const Input& in, Output& out)
{
    using namespace mars::segregated;
    const double* x = in.values;
    out = Output{};
    switch (in.stage) {
    case 0: out.values[0] = pressure_update(x[0],x[1],x[2]); break;
    case 1: out.values[0] = x[0]; break;
    case 2: velocity_update(x,x+3,x+6,out.values); break;
    case 3: case 5: {
        FluxUpdate f{};
        f.density=x[0]; f.old=x[22]; f.alpha=x[23];
        for(int j=0;j<3;++j) {
            f.velocity[j]=x[1+j]; f.influence[j]=x[4+j]; f.compact_gradient[j]=x[7+j];
            f.reconstructed_gradient[j]=x[10+j]; f.original_force[j]=x[13+j];
            f.reconstructed_force[j]=x[16+j]; f.area[j]=x[19+j];
        }
        out.values[0]=mass_flux_update(f,in.stage==5,in.stage==5 && x[24]!=0); break;
    }
    case 4: out.values[0]=inlet_flux_update(x[0],x+1,x+4,x[7],x[8]); break;
    case 6: {
        int flags[3], updated[3];
        for(int s=0;s<3;++s) flags[s]=int(x[3+s]);
        outlet_reversal_update(x,flags,x+6,x+15,x+18,x+21,x[30]!=0,out.values,updated);
        for(int s=0;s<3;++s) out.values[3+s]=updated[s];
        break;
    }
    case 7: outlet_trace_moment(x,x+3,x+6,x[9]!=0,out.values); break;
    case 8: out.values[0]=outlet_trace_update(x[0],x[1],x[2],x[3],x[4],x[5]!=0); break;
    case 9: out.values[0]=x[0]/x[1]; break;
    }
}
} // namespace update_replay
#undef UPDATE_HD
