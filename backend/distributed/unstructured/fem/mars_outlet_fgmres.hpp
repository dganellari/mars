#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// GPT/Codex, 2026-09-10. The controller is shared by device operations and host gates.
// Ops owns persistent vectors; only the small Hessenberg problem lives on the host.
struct OutletKrylovResult
{
    bool converged;
    int iterations;
    double relative_residual;
};

template<class Ops>
OutletKrylovResult outlet_fgmres(Ops& ops, int restart, int max_iterations,
                                double tolerance, double epsilon)
{
    if (restart < 1 || max_iterations < 1 || !(tolerance > 0) || !std::isfinite(tolerance))
        return {false, 0, INFINITY};
    auto x = ops.solution();
    auto r = ops.residual();
    auto w = ops.work();
    ops.zero(x);
    ops.copy(ops.rhs(), r);
    const double b_norm = ops.norm(r);
    if (!std::isfinite(b_norm)) return {false, 0, INFINITY};
    if (b_norm == 0) return {true, 0, 0};
    std::vector<double> h(size_t(restart + 1)*restart), cs(restart), sn(restart);
    std::vector<double> g(restart + 1), y(restart), column(restart + 1);
    auto entry = [&](int i, int j) -> double& { return h[size_t(j)*(restart + 1) + i]; };
    int iterations = 0;
    double beta = b_norm;
    while (iterations < max_iterations)
    {
        std::fill(h.begin(), h.end(), 0.0);
        std::fill(g.begin(), g.end(), 0.0);
        g[0] = beta;
        ops.copy(r, ops.basis(0));
        ops.scale(1/beta, ops.basis(0));
        const int count = std::min(restart, max_iterations - iterations);
        int used = 0;
        bool breakdown = false;
        for (int j = 0; j < count; ++j)
        {
            if (!ops.precondition(ops.basis(j), ops.direction(j)))
                return {false, iterations, beta/b_norm};
            ops.apply(ops.direction(j), w);
            const double original_norm = ops.norm(w);
            if (!std::isfinite(original_norm)) return {false, iterations, INFINITY};
            // Two global classical Gram-Schmidt passes; no per-column host reduction.
            ops.orthogonalize(w, j + 1, column.data());
            for (int i = 0; i <= j; ++i) entry(i, j) = column[i];
            const double next_norm = ops.norm(w);
            if (!std::isfinite(next_norm)) return {false, iterations, INFINITY};
            entry(j + 1, j) = next_norm;
            breakdown = next_norm <= 32*epsilon*original_norm;
            if (!breakdown)
            {
                ops.copy(w, ops.basis(j + 1));
                ops.scale(1/next_norm, ops.basis(j + 1));
            }
            for (int i = 0; i < j; ++i)
            {
                const double a = entry(i, j), b = entry(i + 1, j);
                entry(i, j) = cs[i]*a + sn[i]*b;
                entry(i + 1, j) = -sn[i]*a + cs[i]*b;
            }
            const double diagonal = std::hypot(entry(j, j), entry(j + 1, j));
            ++iterations;
            if (!(diagonal > 0) || !std::isfinite(diagonal))
                return {false, iterations, beta/b_norm};
            cs[j] = entry(j, j)/diagonal;
            sn[j] = entry(j + 1, j)/diagonal;
            entry(j, j) = diagonal;
            entry(j + 1, j) = 0;
            g[j + 1] = -sn[j]*g[j];
            g[j] *= cs[j];
            used = j + 1;
            if (breakdown || std::abs(g[j + 1])/b_norm <= tolerance) break;
        }
        for (int i = used - 1; i >= 0; --i)
        {
            double value = g[i];
            for (int j = i + 1; j < used; ++j) value -= entry(i, j)*y[j];
            y[i] = value/entry(i, i);
            if (!std::isfinite(y[i])) return {false, iterations, INFINITY};
        }
        // Z contains the actual preconditioner output from each iteration.
        for (int i = 0; i < used; ++i) ops.axpy(y[i], ops.direction(i), x);
        ops.apply(x, w);
        ops.copy(ops.rhs(), r);
        ops.axpy(-1, w, r);
        beta = ops.norm(r);
        if (!std::isfinite(beta)) return {false, iterations, INFINITY};
        ops.report(iterations, beta/b_norm);
        // A small Hessenberg residual or happy breakdown alone never certifies the solve.
        if (beta/b_norm <= tolerance) return {true, iterations, beta/b_norm};
        if (breakdown) return {false, iterations, beta/b_norm};
    }
    return {false, iterations, beta/b_norm};
}
