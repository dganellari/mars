#include "fixture.hpp"
#include "backend/distributed/unstructured/fem/mars_nonconformal_flow.hpp"
#include <limits>

int main()
{
    using nonconformal_flow_detail::component;
    auto faces = fixture();
    for (int b = 1; b < int(faces.size()); ++b)
    {
        InterfaceStencil<double> rows[nonconformal_pair_capacity];
        InterfacePatch<double> patch[nonconformal_pair_capacity / 2];
        int count = 0;
        require(nonconformalConservativePair(faces[0], faces[b], {1e-12, 1e-12}, rows, nonconformal_pair_capacity,
                                             count, patch) == InterfaceStatus::valid,
                "flow overlap geometry");
        for (int q = 0; q < count / 2; ++q)
        {
            auto& p = patch[q];
            InterfaceFlowState<double> s[8];
            double mu = .3, rho = 1.7;
            for (int k = 0; k < 8; ++k)
            {
                auto x                = (k < 4 ? faces[0] : faces[b]).x[k % 4];
                s[k].velocity         = {1. + 2 * x.x + 3 * x.y, -.5 + 4 * x.z, .7 - 2 * x.y};
                s[k].pressure         = 3. + x.x - 2 * x.y + 5 * x.z;
                s[k].pressureGradient = {1., -2., 5.};
                s[k].influence        = {.11, .17, .23};
            }
            auto f = nonconformalFlowFlux(p, s, rho, mu, mu, 0.);
            require(f.valid, "finite affine flow flux");
            double meanVelocity[3]{}, meanPressure = 0;
            for (int k = 0; k < 8; ++k)
            {
                for (int j = 0; j < 3; ++j)
                    meanVelocity[j] += .5 * p.shape[k] * component(s[k].velocity, j);
                meanPressure += .5 * p.shape[k] * s[k].pressure;
            }
            double affineMass = 0;
            const double gradient[3][3]{{2, 3, 0}, {0, 0, 4}, {0, -2, 0}};
            for (int i = 0; i < 3; ++i)
            {
                affineMass += rho * p.area * meanVelocity[i] * component(p.normal, i);
                double stress = 0;
                for (int j = 0; j < 3; ++j)
                    stress -= mu * p.area * (gradient[i][j] + gradient[j][i]) * component(p.normal, j);
                near(f.momentum[i], stress, "full affine symmetric viscous stress");
                near(f.pressureForce[i], p.area * meanPressure * component(p.normal, i), "physical pressure traction");
            }
            near(f.mass, affineMass, "Rhie Chow correction vanishes for common affine pressure");
            double constantPressure = 0;
            for (int k = 0; k < 8; ++k)
                constantPressure += f.pressureCoefficient[k];
            near(constantPressure, 0., "pressure gauge nullspace");
            for (int k = 0; k < 8; ++k)
            {
                auto saved = s[k];
                s[k].pressure += .125;
                auto changed = nonconformalFlowFlux(p, s, rho, mu, mu, 0.);
                near((changed.mass - f.mass) / .125, f.pressureCoefficient[k], "pressure flux derivative");
                s[k] = saved;
                for (int j = 0; j < 3; ++j)
                {
                    if (j == 0) s[k].velocity.x += .125;
                    if (j == 1) s[k].velocity.y += .125;
                    if (j == 2) s[k].velocity.z += .125;
                    changed = nonconformalFlowFlux(p, s, rho, mu, mu, 0.);
                    near((changed.mass - f.mass) / .125, f.velocityCoefficient[k][j], "velocity mass derivative");
                    for (int i = 0; i < 3; ++i)
                        near((changed.momentum[i] - f.momentum[i]) / .125, f.momentumCoefficient[k][i][j],
                             "coupled viscous block derivative");
                    s[k] = saved;
                }
            }
            // A pressure discontinuity must drive a positive flux from high to low.
            for (int k = 0; k < 8; ++k)
            {
                s[k].velocity         = {};
                s[k].pressure         = k < 4 ? 2. : 1.;
                s[k].pressureGradient = {};
            }
            auto jump = nonconformalFlowFlux(p, s, rho, mu, mu, 0.);
            require(jump.valid && jump.mass > 0, "pressure jump drives flux to low pressure");
            for (double advecting : {-2., 3.})
            {
                for (int k = 0; k < 8; ++k)
                    s[k].velocity = {1.2, -.4, .7};
                auto adv = nonconformalFlowFlux(p, s, rho, mu, mu, advecting);
                require(adv.valid, "both advection directions valid");
                near(adv.momentum[0], 1.2 * advecting, "shared mass flux transports momentum x");
                near(adv.momentum[1], -.4 * advecting, "shared mass flux transports momentum y");
                near(adv.momentum[2], .7 * advecting, "shared mass flux transports momentum z");
            }
            for (int k = 0; k < 8; ++k)
            {
                auto x        = (k < 4 ? faces[0] : faces[b]).x[k % 4];
                s[k].velocity = {x.z, 0., -x.x};
            }
            auto rotation = nonconformalFlowFlux(p, s, rho, mu, mu, 0.);
            require(rotation.valid, "rigid rotation valid");
            for (int i = 0; i < 3; ++i)
                near(rotation.momentum[i], 0., "rigid rotation has no viscous stress");
            s[0].influence.x = -1.;
            require(!nonconformalFlowFlux(p, s, rho, mu, mu, 0.).valid, "negative influence refused");
            s[0].influence.x = .1;
            s[0].pressure    = std::numeric_limits<double>::quiet_NaN();
            require(!nonconformalFlowFlux(p, s, rho, mu, mu, 0.).valid, "nonfinite state refused");
        }
    }
    std::printf("PASS: stationary non-conformal flow flux checks\n");
}
