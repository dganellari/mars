#pragma once
#include "mars_nonconformal_overlap.hpp"

namespace mars::fem
{
// Physical pressure [Pa], velocity [m/s], rho [kg/m^3], dynamic mu [Pa s].
// SIMPLE influence is volume / momentum diagonal [m^3 s/kg].
template<class T>
struct InterfaceFlowState
{
    Vec3<T> velocity{}, pressureGradient{}, influence{};
    T pressure{};
};
template<class T>
struct InterfaceFlowFlux
{
    bool valid = false;
    T mass{}, pressureCoefficient[8]{}, velocityCoefficient[8][3]{};
    T momentum[3]{}, momentumCoefficient[8][3][3]{};
    T pressureForce[3]{}, pressureForceCoefficient[8][3]{};
};
namespace nonconformal_flow_detail
{
template<class T>
MARS_FS_FN T component(Vec3<T> v, int i)
{
    return i == 0 ? v.x : (i == 1 ? v.y : v.z);
}
template<class T>
MARS_FS_FN bool finite(Vec3<T> v)
{
    return nonconformal_detail::finite(v.x) && nonconformal_detail::finite(v.y) && nonconformal_detail::finite(v.z);
}
} // namespace nonconformal_flow_detail

// A single outgoing mass flux enters continuity and the lagged momentum flux.
// Pressure coefficients differentiate this flux with reconstructed gradients
// and influence held fixed, exactly as in SIMPLE's pressure-correction step.
// The caller adds every output to row[0] and subtracts it from row[1].
template<class T>
MARS_FS_FN InterfaceFlowFlux<T> nonconformalFlowFlux(
    const InterfacePatch<T>& p, const InterfaceFlowState<T>* state, T density, T muLeft, T muRight, T advectingMass)
{
    using namespace nonconformal_detail;
    using nonconformal_flow_detail::component;
    InterfaceFlowFlux<T> out;
    if (!state || !finite(density) || density <= 0 || !finite(muLeft) || muLeft <= 0 || !finite(muRight) ||
        muRight <= 0 || !finite(advectingMass) || !finite(p.area) || p.area <= 0 ||
        !nonconformal_flow_detail::finite(p.normal) || std::abs(dot(p.normal, p.normal) - T(1)) > T(1e-10))
        return out;
    for (int side = 0; side < 2; ++side)
        if (!finite(p.inverseHeight[side]) || p.inverseHeight[side] <= 0) return out;
    T d[3]{}, sumShape[2]{}, reconstructed[3]{};
    for (int k = 0; k < 8; ++k)
    {
        if (!finite(p.shape[k]) || !nonconformal_flow_detail::finite(p.gradient[k]) ||
            !nonconformal_flow_detail::finite(state[k].velocity) || !finite(state[k].pressure) ||
            !nonconformal_flow_detail::finite(state[k].pressureGradient) ||
            !nonconformal_flow_detail::finite(state[k].influence))
            return out;
        sumShape[k / 4] += p.shape[k];
        for (int j = 0; j < 3; ++j)
        {
            T influence = component(state[k].influence, j);
            if (influence < 0) return out;
            d[j] += T(.5) * p.shape[k] * influence;
            reconstructed[j] += T(.5) * p.shape[k] * component(state[k].pressureGradient, j);
        }
    }
    if (std::abs(sumShape[0] - T(1)) > T(1e-10) || std::abs(sumShape[1] - T(1)) > T(1e-10)) return out;
    T dn = 0;
    for (int j = 0; j < 3; ++j)
    {
        T normal = component(p.normal, j);
        if (d[j] < 0) return out;
        dn += normal * normal * d[j];
        out.mass += density * p.area * d[j] * reconstructed[j] * normal;
    }
    const T inverseHeight  = T(.5) * (p.inverseHeight[0] + p.inverseHeight[1]);
    const T viscousPenalty = T(.5) * (muLeft * p.inverseHeight[0] + muRight * p.inverseHeight[1]);
    for (int k = 0; k < 8; ++k)
    {
        const T sign = k < 4 ? T(1) : T(-1);
        const T mu   = k < 4 ? muLeft : muRight;
        T pressure   = dn * inverseHeight * sign * p.shape[k];
        for (int i = 0; i < 3; ++i)
        {
            const T ni = component(p.normal, i);
            pressure -= T(.5) * d[i] * component(p.gradient[k], i) * ni;
            out.velocityCoefficient[k][i] = density * p.area * T(.5) * p.shape[k] * ni;
            out.mass += out.velocityCoefficient[k][i] * component(state[k].velocity, i);
            out.pressureForceCoefficient[k][i] = p.area * T(.5) * p.shape[k] * ni;
            out.pressureForce[i] += out.pressureForceCoefficient[k][i] * state[k].pressure;
            for (int j = 0; j < 3; ++j)
            {
                const T nj = component(p.normal, j);
                T value    = -T(.5) * mu * component(p.gradient[k], i) * nj;
                if (i == j) value -= T(.5) * mu * dot(p.gradient[k], p.normal);
                value += viscousPenalty * sign * p.shape[k] * (T(i == j) + ni * nj);
                value *= p.area;
                if (i == j && ((advectingMass >= 0) == (k < 4))) value += advectingMass * p.shape[k];
                out.momentumCoefficient[k][i][j] = value;
                out.momentum[i] += value * component(state[k].velocity, j);
            }
        }
        out.pressureCoefficient[k] = density * p.area * pressure;
        out.mass += out.pressureCoefficient[k] * state[k].pressure;
    }
    out.valid = finite(out.mass);
    for (int i = 0; i < 3; ++i)
    {
        out.valid = out.valid && finite(out.momentum[i]) && finite(out.pressureForce[i]);
        for (int k = 0; k < 8; ++k)
        {
            out.valid = out.valid && finite(out.pressureCoefficient[k]) && finite(out.velocityCoefficient[k][i]) &&
                        finite(out.pressureForceCoefficient[k][i]);
            for (int j = 0; j < 3; ++j)
                out.valid = out.valid && finite(out.momentumCoefficient[k][i][j]);
        }
    }
    return out;
}
} // namespace mars::fem
