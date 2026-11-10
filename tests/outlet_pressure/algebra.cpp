#include "pressure_kernels.hpp"
#include <array>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using Vector = std::array<double, 3>;
using Matrix = std::array<double, 9>;
Vector solve(Matrix a, Vector b)
{
    for (int k = 0; k < 3; ++k)
    {
        if (std::abs(a[3 * k + k]) < 1e-12) throw std::runtime_error("singular toy matrix");
        for (int i = k + 1; i < 3; ++i)
        {
            double factor = a[3 * i + k] / a[3 * k + k];
            for (int j = k; j < 3; ++j)
                a[3 * i + j] -= factor * a[3 * k + j];
            b[i] -= factor * b[k];
        }
    }
    Vector x{};
    for (int i = 2; i >= 0; --i)
    {
        x[i] = b[i];
        for (int j = i + 1; j < 3; ++j)
            x[i] -= a[3 * i + j] * x[j];
        x[i] /= a[3 * i + i];
    }
    return x;
}
int main()
{
    const int rows[4] = {0, 3, 6, 9}, cols[9] = {0, 1, 2, 0, 1, 2, 0, 1, 2}, diag[3] = {0, 4, 8}, map[3] = {0, 1, 2};
    const uint8_t mask[3] = {0, 0, 1};
    const Matrix original{1, -1, 0, -1, 2, -1, 0, -1, 1};
    std::cout << "[";
    for (int pass = 0; pass < 2; ++pass)
    {
        const double old = pass ? 2 : 0, target = pass ? 0 : 2;
        Vector p{old, old, old}, targets{0, 0, target}, rhs{};
        Matrix row = original;
        for (int i = 0; i < 3; ++i)
        {
            threadIdx.x = i;
            enforceBcMatrixKernel(mask, rows, cols, diag, row.data(), 3);
            enforcePressureBcRhsLiftKernel(mask, targets.data(), p.data(), map, rhs.data(), 3);
        }
        Matrix symmetric = row;
        for (int i = 0; i < 3; ++i)
        {
            threadIdx.x = i;
            enforceBcColMatrixKernelByNode(mask, map, map, 3, rows, cols, symmetric.data(), 3);
        }
        Vector lifted = rhs;
        for (int i = 0; i < 2; ++i)
            lifted[i] -= original[3 * i + 2] * (target - old);
        const auto x = solve(row, rhs), wrong = solve(symmetric, rhs), right = solve(symmetric, lifted);
        for (int i = 0; i < 3; ++i)
            if (std::abs(x[i] + old - target) > 1e-12 || std::abs(right[i] + old - target) > 1e-12) return 1;
        if (std::abs(wrong[0] + old - target) < 1) return 1;
        double skew = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                skew = std::max(skew, std::abs(row[3 * i + j] - row[3 * j + i]));
        if (pass) std::cout << ",";
        std::cout << "{\"target\":" << target << ",\"old_pressure\":" << old << ",\"row_only_skew\":" << skew
                  << ",\"correct_interior_pressure\":" << x[0] + old
                  << ",\"column_zero_without_lift_interior\":" << wrong[0] + old
                  << ",\"column_zero_with_lift_interior\":" << right[0] + old << "}";
    }
    // Distinct outlet pressures give a linear interior value, not either boundary value.
    const uint8_t ends[3] = {1, 0, 1};
    Vector old_pressure{0, 0, 0}, targets{-2, 99, 4}, rhs{};
    Matrix row = original;
    for (int i = 0; i < 3; ++i)
    {
        threadIdx.x = i;
        enforceBcMatrixKernel(ends, rows, cols, diag, row.data(), 3);
        enforcePressureBcRhsLiftKernel(ends, targets.data(), old_pressure.data(), map, rhs.data(), 3);
    }
    const auto unequal = solve(row, rhs);
    if (unequal != Vector{-2, 1, 4} || rhs[1] != 0) return 1;
    std::cout << "]\n";
}
