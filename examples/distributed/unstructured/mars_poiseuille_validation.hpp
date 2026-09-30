#pragma once

// Validation of the Poiseuille example, kept apart from the tutorial code in
// mars_poiseuille_flow.cu.
//
// After the last step, against the developed plane-Poiseuille solution:
//   profile    RMS of u - 1.5 U (1 - eta^2) on the plane 90% down the channel,
//              and the profile of one node plane as CSV for plotting
//   flux       Q(x) / Q(inlet) at 25, 50 and 75% of the length, from the computed u
//   pressure   -dp/dx between 60% and 90% against 12 rho nu U / H^2
//
// --check turns this into the release gate. It also checks the projection of
// the first three and the last step (NavierStokes::projectionReport), requires the
// velocity to be steady over the last 20 steps, and exits 1 if anything fails.
// --comparison-output writes u, v, w, p at full precision for comparisons with
// other codes.

#include "backend/distributed/unstructured/fem/mars_navier_stokes.hpp"
#include "backend/distributed/unstructured/utils/mars_vtu_parallel_writer.hpp"

#include <thrust/copy.h>
#include <thrust/for_each.h>
#include <thrust/pair.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace poiseuille
{

using KeyType  = uint64_t;
using RealType = double;
using Solver   = mars::fem::NavierStokes<KeyType, RealType>;
using Domain   = Solver::Domain;
using Vector   = Solver::Vector;

struct Options
{
    bool check           = false;
    double rmsTol        = 6e-3; // the reference tolerance
    double fluxTol       = 0.10;
    double steadyTol     = 1e-6; // velocity change over the last 20 steps / U
    double continuityTol = 1e-6; // continuity RMS * H / U and boundary flux balance
    double profileX      = -1;   // default: 90% down the channel
    double profileXTol   = -1;   // default: 2% of the length
    std::string comparisonPrefix;
    int comparisonEvery = 50;
};

// Takes a validation option. Returns false if arg is not one; sets bad for a malformed value.
inline bool parseOption(const std::string& arg, Options& o, bool& bad)
{
    auto take = [&arg, &bad](const char* key, auto& out) {
        const std::string prefix = std::string("--") + key + "=";
        if (arg.rfind(prefix, 0) != 0) return false;
        std::istringstream in(arg.substr(prefix.size()));
        bad = bad || !(in >> out) || !in.eof();
        return true;
    };
    if (arg == "--check") return o.check = true;
    bool known = take("rms-tol", o.rmsTol) || take("flux-tol", o.fluxTol) || take("steady-tol", o.steadyTol) ||
                 take("continuity-tol", o.continuityTol) || take("profile-x", o.profileX) ||
                 take("profile-xtol", o.profileXTol) || take("comparison-output", o.comparisonPrefix) ||
                 take("comparison-every", o.comparisonEvery);
    bad = bad || !(o.rmsTol > 0 && o.fluxTol > 0 && o.steadyTol > 0 && o.continuityTol > 0 && o.comparisonEvery > 0);
    return known;
}

inline void printOptions()
{
    std::cout << "Validation:\n"
                 "  --check               release gate: exit 1 unless every check below passes\n"
                 "  --rms-tol=X           profile RMS tolerance (default 6e-3)\n"
                 "  --flux-tol=X          flux ratio tolerance around 1 (default 0.10)\n"
                 "  --steady-tol=X        velocity change over the last 20 steps / U (default 1e-6)\n"
                 "  --continuity-tol=X    continuity RMS * H / U and flux balance (default 1e-6)\n"
                 "  --profile-x=X         probe plane (default 90% down the channel)\n"
                 "  --profile-xtol=X      probe half width (default 2% of the length)\n"
                 "  --comparison-output=PREFIX  u, v, w, p at full precision; --comparison-every=N (default 50)\n";
}

template<class Function>
struct OwnedTerm
{
    const uint8_t* owned;
    Function f;
    __device__ double operator()(size_t i) const { return owned[i] == 1 ? f(i) : 0.0; }
};

// Global sum over the owned nodes of f(i).
template<class Function>
double ownedSum(const Domain& domain, Function f)
{
    double local = thrust::transform_reduce(
        thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(domain.getNodeCount()),
        OwnedTerm<Function>{domain.getNodeOwnershipMap().data(), f}, 0.0, thrust::plus<double>());
    double global = 0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
}

// x of the node plane nearest target, the same on every rank; a tie goes to the smaller x.
inline double nearestNodePlane(const Domain& domain, double target)
{
    using Candidate      = thrust::pair<double, double>;
    const double inf     = std::numeric_limits<double>::infinity();
    const RealType* x    = domain.getNodeX().data();
    const Candidate best = thrust::transform_reduce(
        thrust::device, thrust::counting_iterator<size_t>(0), thrust::counting_iterator<size_t>(domain.getNodeCount()),
        [x, target] __device__(size_t i) -> Candidate { return Candidate(fabs(double(x[i]) - target), double(x[i])); },
        Candidate(inf, inf),
        [] __device__(const Candidate& a, const Candidate& b) -> Candidate { return b < a ? b : a; });
    double distance = inf, plane = inf;
    MPI_Allreduce(&best.first, &distance, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    const double candidate = best.first == distance ? best.second : inf;
    MPI_Allreduce(&candidate, &plane, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return plane;
}

// Mass-weighted change of (u, v) against a snapshot, relative to U.
inline double velocityChange(const Solver& s, const Domain& domain, const Vector& u0, const Vector& v0, double U)
{
    double change = ownedSum(domain, [u = s.u.data(), v = s.v.data(), u0 = u0.data(), v0 = v0.data(),
                                      m = s.mass().data()] __device__(size_t i) -> double {
        double du = u[i] - u0[i], dv = v[i] - v0[i];
        return m[i] * (du * du + dv * dv);
    });
    double volume = ownedSum(domain, [m = s.mass().data()] __device__(size_t i) -> double { return m[i]; });
    return std::sqrt(change / volume) / U;
}

// (y, u) of the owned nodes on the plane x = xPlane, gathered on rank 0 and sorted by y.
inline std::vector<std::pair<double, double>> gatherPlane(const Solver& s, const Domain& domain, double xPlane,
                                                          double tolerance)
{
    const size_t n     = domain.getNodeCount();
    const uint8_t* own = domain.getNodeOwnershipMap().data();
    const RealType* x  = domain.getNodeX().data();
    cstone::DeviceVector<int> ids(n);
    int count = int(thrust::copy_if(thrust::device, thrust::counting_iterator<int>(0),
                                    thrust::counting_iterator<int>(int(n)), ids.data(),
                                    [own, x, xPlane, tolerance] __device__(int i) -> bool {
                                        return own[i] == 1 && fabs(double(x[i]) - xPlane) < tolerance;
                                    }) -
                    ids.data());
    cstone::DeviceVector<double> pairs(2 * size_t(count));
    thrust::for_each(thrust::device, thrust::counting_iterator<int>(0), thrust::counting_iterator<int>(count),
                     [ids = ids.data(), y = domain.getNodeY().data(), u = s.u.data(), out = pairs.data()] __device__(
                         int k) {
                         out[2 * k]     = y[ids[k]];
                         out[2 * k + 1] = u[ids[k]];
                     });
    std::vector<double> local(pairs.size());
    cudaMemcpy(local.data(), pairs.data(), local.size() * sizeof(double), cudaMemcpyDeviceToHost);

    int rank = 0, numRanks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
    int sendCount = int(local.size());
    std::vector<int> counts(numRanks, 0), displs(numRanks, 0);
    MPI_Gather(&sendCount, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    for (int r = 1; r < numRanks; ++r) displs[r] = displs[r - 1] + counts[r - 1];
    std::vector<double> all(rank == 0 ? size_t(displs.back() + counts.back()) : 0);
    MPI_Gatherv(local.data(), sendCount, MPI_DOUBLE, all.data(), counts.data(), displs.data(), MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    std::vector<std::pair<double, double>> points;
    for (size_t k = 0; k + 1 < all.size(); k += 2) points.emplace_back(all[k], all[k + 1]);
    std::sort(points.begin(), points.end());
    return points;
}

// -dp/dx from a quadratic fit u = a y^2 + b y + c over the core (u > Umax / 2): G = -2 a rho nu.
inline double pressureGradientFromProfile(const std::vector<std::pair<double, double>>& points, double uMax,
                                          double rho, double nu)
{
    double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, b0 = 0, b1 = 0, b2 = 0;
    for (auto [y, u] : points)
    {
        if (u <= 0.5 * uMax) continue;
        double y2 = y * y;
        s0 += 1, s1 += y, s2 += y2, s3 += y2 * y, s4 += y2 * y2;
        b0 += u, b1 += y * u, b2 += y2 * u;
    }
    if (s0 < 5) return std::numeric_limits<double>::quiet_NaN();
    auto det3 = [](double a11, double a12, double a13, double a21, double a22, double a23, double a31, double a32,
                   double a33) {
        return a11 * (a22 * a33 - a23 * a32) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31);
    };
    double D  = det3(s4, s3, s2, s3, s2, s1, s2, s1, s0);
    double Da = det3(b2, s3, s2, b1, s2, s1, b0, s1, s0);
    return D != 0 ? -2.0 * Da / D * rho * nu : std::numeric_limits<double>::quiet_NaN();
}

// Mean pressure over the core (u > Umax / 2) of the node plane nearest x.
inline double corePressure(const Solver& s, const Domain& domain, double x, double tolerance, double uMax)
{
    const double plane = nearestNodePlane(domain, x);
    const RealType* xs = domain.getNodeX().data();
    const RealType* u  = s.u.data();
    double sum = ownedSum(domain, [xs, u, p = s.p.data(), plane, tolerance, uMax] __device__(size_t i) -> double {
        return fabs(xs[i] - plane) < tolerance && u[i] > 0.5 * uMax ? p[i] : 0.0;
    });
    double count = ownedSum(domain, [xs, u, plane, tolerance, uMax] __device__(size_t i) -> double {
        return fabs(xs[i] - plane) < tolerance && u[i] > 0.5 * uMax ? 1.0 : 0.0;
    });
    return count > 0 ? sum / count : std::numeric_limits<double>::quiet_NaN();
}

// Flux of u in +x through the node plane nearest x.
inline double planeFlux(const Solver& s, const Domain& domain, double x, double tolerance, double xMax)
{
    const double plane = nearestNodePlane(domain, x);
    const int side     = plane < xMax - tolerance ? +1 : -1; // elements on the +x side, except at the outlet
    Vector area;
    mars::fem::openingFaceAreas<KeyType>(domain, mars::fem::Opening<RealType>{0, plane, side}, RealType(tolerance),
                                         area);
    // The outward normal points to -x on side +1.
    return -side * ownedSum(domain, [a = area.data(), u = s.u.data()] __device__(size_t i) -> double {
               return u[i] * a[i];
           });
}

// The per-step part: comparison frames, the steady-state snapshot and the
// projection checks. finish() runs the final checks and returns the exit code.
class Monitor
{
public:
    Monitor(const Options& o, const Solver& s, int numSteps, double inflow)
        : opt_(o)
        , numSteps_(numSteps)
        , inflow_(inflow)
        , steadyWindow_(std::min(20, numSteps - 1))
    {
        if (!o.comparisonPrefix.empty())
            comparison_ = std::make_unique<mars::fem::VTUParallelWriter<KeyType, RealType>>(o.comparisonPrefix, true);
        zeroW_.resize(s.u.size());
        cudaMemset(zeroW_.data(), 0, zeroW_.size() * sizeof(RealType));
    }

    void afterStep(Solver& s, const Domain& domain, int step, double t)
    {
        if (comparison_ && (step == 0 || step % opt_.comparisonEvery == 0 || step == numSteps_))
        {
            using FD = mars::fem::VTUParallelWriter<KeyType, RealType>::FieldDesc;
            std::vector<FD> fields{{"u", FD::Kind::PointScalar, &s.u, nullptr, nullptr},
                                   {"v", FD::Kind::PointScalar, &s.v, nullptr, nullptr},
                                   {"w", FD::Kind::PointScalar, &zeroW_, nullptr, nullptr},
                                   {"p", FD::Kind::PointScalar, &s.p, nullptr, nullptr}};
            comparison_->writeMultiFieldFrame(step, t, domain, fields);
        }
        if (!opt_.check) return;
        if (step == numSteps_ - steadyWindow_)
        {
            u0_ = s.u;
            v0_ = s.v;
        }
        if (step >= 1 && (step <= 3 || step == numSteps_))
        {
            last_ = s.projectionReport();
            // When both divergences approach zero the relative identity loses its meaning; then
            // the absolute one must be at roundoff of the velocity gradient scale U / H.
            const auto& box      = s.box();
            double identityFloor = 1e-10 * inflow_ / (box.hi[1] - box.lo[1]);
            bool ok = last_.finite && !(last_.identity > 1e-7 && last_.identityRms > identityFloor) &&
                      last_.balanceIdentity <= 1e-10 && last_.unreachedMax <= 1e-8;
            projectionOk_ = projectionOk_ && ok;
            if (domain.rank() == 0)
                std::cout << std::scientific << std::setprecision(6) << "[channel-projection] step=" << step
                          << " identity=" << last_.identity << " identity_rms=" << last_.identityRms
                          << " continuity_rms=" << last_.continuityRms << " continuity_max=" << last_.continuityMax
                          << " unreached_max=" << last_.unreachedMax << " Qin=" << last_.inflow
                          << " Qout=" << last_.outflow
                          << " balance=" << last_.balance << " balance_identity=" << last_.balanceIdentity
                          << " gate=" << (ok ? "PASS" : "FAIL") << "\n"
                          << std::defaultfloat;
        }
    }

    int finish(const Solver& s, const Domain& domain, const std::string& csvPrefix, double rho, double nu)
    {
        const bool root  = domain.rank() == 0;
        const auto& box  = s.box();
        const double L   = box.hi[0] - box.lo[0];
        const double H   = box.hi[1] - box.lo[1];
        const double uMax = 1.5 * inflow_;
        const double G    = 12 * rho * nu * inflow_ / (H * H);
        const double snap = 1e-9 * std::max(1.0, L);

        // Profile RMS over a slab around the probe plane.
        const double xProbe = opt_.profileX >= 0 ? opt_.profileX : box.lo[0] + 0.9 * L;
        const double xTol   = opt_.profileXTol > 0 ? opt_.profileXTol : 0.02 * L;
        const double yMid = 0.5 * (box.lo[1] + box.hi[1]), halfH = 0.5 * H;
        const RealType* x = domain.getNodeX().data();
        double sumSq      = ownedSum(domain, [x, y = domain.getNodeY().data(), u = s.u.data(), xProbe, xTol, uMax, yMid,
                                         halfH] __device__(size_t i) -> double {
            if (!(fabs(x[i] - xProbe) < xTol)) return 0.0;
            double eta = (y[i] - yMid) / halfH, err = u[i] - uMax * (1 - eta * eta);
            return err * err;
        });
        double count = ownedSum(domain, [x, xProbe, xTol] __device__(size_t i) -> double {
            return fabs(x[i] - xProbe) < xTol ? 1.0 : 0.0;
        });
        double rms   = count > 0 ? std::sqrt(sumSq / count) : -1.0;

        // One node plane for plotting (scripts/plot_poiseuille_profile.py).
        const double xPlane = nearestNodePlane(domain, xProbe);
        auto points         = gatherPlane(s, domain, xPlane, snap);
        if (root)
        {
            std::string name = (csvPrefix.empty() ? std::string("poiseuille") : csvPrefix) + "_profile.csv";
            std::ofstream csv(name);
            csv << "y,u_solved,u_analytic\n";
            for (auto [y, u] : points)
            {
                double eta = (y - yMid) / halfH;
                csv << y << "," << u << "," << uMax * (1 - eta * eta) << "\n";
            }
            std::cout << "Profile CSV: " << name << " (plane x=" << xPlane << ", " << points.size() << " nodes)\n";
        }

        double qIn = planeFlux(s, domain, box.lo[0], snap, box.hi[0]);
        const double fractions[3] = {0.25, 0.50, 0.75};
        double ratio[3];
        for (int k = 0; k < 3; ++k)
            ratio[k] = qIn != 0 ? planeFlux(s, domain, box.lo[0] + fractions[k] * L, snap, box.hi[0]) / qIn : 0.0;

        double pA = corePressure(s, domain, box.lo[0] + 0.6 * L, snap, uMax);
        double pB = corePressure(s, domain, box.lo[0] + 0.9 * L, snap, uMax);
        double gradP = (pA - pB) / (0.3 * L);

        double steady = std::numeric_limits<double>::infinity();
        if (opt_.check && steadyWindow_ > 0) steady = velocityChange(s, domain, u0_, v0_, inflow_);

        if (root)
        {
            std::cout << std::scientific << std::setprecision(6) << "\nPoiseuille validation\n"
                      << "  profile RMS at x=" << std::fixed << std::setprecision(4) << xProbe << " +/- " << xTol
                      << ": " << std::scientific << std::setprecision(6) << rms << " (U_max=" << uMax
                      << ", " << long(count) << " nodes)\n"
                      << "  flux Q(x)/Q(inlet) at 25/50/75%: " << std::fixed << std::setprecision(4) << ratio[0]
                      << " / " << ratio[1] << " / " << ratio[2] << "\n"
                      << std::scientific << std::setprecision(4) << "  -dp/dx from p: " << gradP
                      << "   from the u profile: " << pressureGradientFromProfile(points, uMax, rho, nu)
                      << "   exact: " << G << "\n"
                      << std::defaultfloat;
        }
        if (!opt_.check) return 0;

        bool rmsOk        = rms >= 0 && rms < opt_.rmsTol;
        bool fluxOk       = std::all_of(ratio, ratio + 3, [&](double r) { return std::abs(r - 1) <= opt_.fluxTol; });
        bool steadyOk     = std::isfinite(steady) && steady <= opt_.steadyTol;
        double continuity = last_.continuityRms * H / inflow_;
        bool continuityOk = std::isfinite(continuity) && continuity <= opt_.continuityTol &&
                            std::isfinite(last_.balance) && last_.balance <= opt_.continuityTol;
        bool pass         = rmsOk && fluxOk && steadyOk && continuityOk && projectionOk_;
        if (root)
            std::cout << std::scientific << std::setprecision(3) << "VALIDATION " << (pass ? "PASS" : "FAIL")
                      << ": RMS=" << rms << (rmsOk ? " < " : " >= ") << opt_.rmsTol << ", flux "
                      << (fluxOk ? "PASS" : "FAIL") << ", steady=" << steady << (steadyOk ? " PASS" : " FAIL")
                      << ", continuity*H/U=" << continuity << " balance=" << last_.balance
                      << (continuityOk ? " PASS" : " FAIL") << ", projection " << (projectionOk_ ? "PASS" : "FAIL")
                      << "\n"
                      << std::defaultfloat;
        return pass ? 0 : 1;
    }

private:
    Options opt_;
    int numSteps_;
    double inflow_;
    int steadyWindow_;
    std::unique_ptr<mars::fem::VTUParallelWriter<KeyType, RealType>> comparison_;
    Vector zeroW_, u0_, v0_;
    Solver::ProjectionReport last_{};
    bool projectionOk_ = true;
};

} // namespace poiseuille
