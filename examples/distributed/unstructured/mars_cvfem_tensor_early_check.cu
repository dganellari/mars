#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <vector>
#ifndef MARS_TENSOR_EARLY_HOST_CHECK
#include "backend/distributed/unstructured/fem/mars_cvfem_hex_kernel_tensor.hpp"
#include "backend/distributed/unstructured/fem/mars_cvfem_hex_kernel_tensor_early.hpp"
#endif

namespace {
#ifdef MARS_TENSOR_EARLY_HOST_CHECK
using Stream = int;
#define CHECK_GPU(call) ((void)0)
#define LAUNCH_KERNEL(name, count, ...) \
    for (size_t work = 0; work < ((count + block - 1) / block) * block; ++work) { \
        blockIdx.x = unsigned(work / block); blockDim.x = block; threadIdx.x = unsigned(work % block); \
        mars::fem::name<uint64_t, double>(__VA_ARGS__); \
    }
#else
using Stream = cudaStream_t;
void check_gpu(cudaError_t status) {
    if (status != cudaSuccess) throw std::runtime_error(cudaGetErrorString(status));
}
#define CHECK_GPU(call) check_gpu(call)
#define LAUNCH_KERNEL(name, count, ...) \
    mars::fem::name<uint64_t, double><<<(count + block - 1) / block, block, 0, stream>>>(__VA_ARGS__)
#endif

template<class T> struct Buffer {
    T* ptr = nullptr;
    size_t size = 0;
#ifdef MARS_TENSOR_EARLY_HOST_CHECK
    std::vector<T> storage;
#endif
    Buffer() = default;
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;
    void set(const std::vector<T>& input) {
        size = input.size();
#ifdef MARS_TENSOR_EARLY_HOST_CHECK
        storage = input; ptr = storage.data();
#else
        CHECK_GPU(cudaMalloc(reinterpret_cast<void**>(&ptr), size * sizeof(T)));
        CHECK_GPU(cudaMemcpy(ptr, input.data(), size * sizeof(T), cudaMemcpyHostToDevice));
#endif
    }
    std::vector<T> get() const {
#ifdef MARS_TENSOR_EARLY_HOST_CHECK
        return storage;
#else
        std::vector<T> result(size);
        CHECK_GPU(cudaMemcpy(result.data(), ptr, size * sizeof(T), cudaMemcpyDeviceToHost));
        return result;
#endif
    }
    ~Buffer() {
#ifndef MARS_TENSOR_EARLY_HOST_CHECK
        if (ptr) cudaFree(ptr);
#endif
    }
};

struct Fixture {
    int elements, nodes;
    std::array<std::vector<uint64_t>, 8> conn;
    std::array<std::vector<double>, 9> field;
    std::array<std::vector<double>, 4> scs;
    std::vector<int> dofs, rowptr, columns, diagonal;
    std::vector<uint8_t> ownership;
    std::vector<double> values, rhs;
    int owned_rows;

    Fixture(int count, int variant) : elements(count), nodes(4 * (count + 1)), owned_rows(nodes) {
        for (auto& v : field) v.resize(nodes);
        for (auto& v : scs) v.resize(12 * count);
        dofs.resize(nodes); ownership.assign(nodes, 1);
        for (int n = 0; n < nodes; ++n) {
            const double x = n / 4, y = (n % 4) % 2, z = (n % 4) / 2;
            const double warp = variant == 0 ? 0 : 0.013;
            field[0][n] = x + .13 * y + warp * y * z;
            field[1][n] = y + .07 * z + warp * x * z;
            field[2][n] = z + .03 * x + warp * x * y;
            field[3][n] = .8 + .003 * n;
            field[4][n] = .2 + .017 * n + .03 * std::sin(double(n));
            field[5][n] = .1 + .04 * (n % 7);
            field[6][n] = .2 + .001 * n;
            field[7][n] = -.3 + .002 * n;
            field[8][n] = .05 - .003 * n;
            dofs[n] = nodes - 1 - n;
        }
        for (int e = 0; e < count; ++e) {
            const uint64_t corners[8] = {unsigned(4*e), unsigned(4*e+4), unsigned(4*e+5), unsigned(4*e+1),
                                          unsigned(4*e+2), unsigned(4*e+6), unsigned(4*e+7), unsigned(4*e+3)};
            for (int i = 0; i < 8; ++i) conn[i].push_back(corners[i]);
            for (int ip = 0; ip < 12; ++ip) {
                scs[0][12*e+ip] = .1 * ((ip % 3) - 1);
                // Arbitrary oriented areas isolate assembly algebra from geometry preprocessing.
                scs[1][12*e+ip] = .07 + .013 * ip;
                scs[2][12*e+ip] = -.11 + .008 * ip;
                scs[3][12*e+ip] = .21 - .007 * ip;
            }
        }
        if (variant == 2) {
            dofs[nodes-1] = dofs[0]; // Different corners alias the same diagonal.
            dofs[2] = -1;
            ownership[3] = 0;
            owned_rows = nodes - 3; // Owned nodes can still map to ghost rows.
        }
        if (variant == 3) owned_rows = -1;
        for (int row = 0; row < nodes; ++row) {
            rowptr.push_back(int(columns.size()));
            for (int col = 0; col < nodes; ++col) {
                if (variant == 3 && row != col && std::abs(row-col) > 5) continue;
                if (row == col) diagonal.push_back(int(columns.size()));
                columns.push_back(col);
            }
        }
        rowptr.push_back(int(columns.size()));
        values.assign(columns.size(), .03125);
        rhs.assign(nodes, .125);
    }
};

struct DeviceFixture {
    std::array<Buffer<uint64_t>, 8> conn;
    std::array<Buffer<double>, 9> field;
    std::array<Buffer<double>, 4> scs;
    Buffer<int> dofs, rowptr, columns, diagonal;
    Buffer<uint8_t> ownership;
    Buffer<double> values, rhs;
    Buffer<mars::fem::CSRMatrix<double>> matrix;
    explicit DeviceFixture(const Fixture& f) {
        for (int i = 0; i < 8; ++i) conn[i].set(f.conn[i]);
        for (int i = 0; i < 9; ++i) field[i].set(f.field[i]);
        for (int i = 0; i < 4; ++i) scs[i].set(f.scs[i]);
        dofs.set(f.dofs); rowptr.set(f.rowptr); columns.set(f.columns); diagonal.set(f.diagonal);
        ownership.set(f.ownership); values.set(f.values); rhs.set(f.rhs);
        matrix.set({{rowptr.ptr, columns.ptr, values.ptr, diagonal.ptr, f.nodes, int(f.values.size()), f.owned_rows}});
    }
};

void launch(DeviceFixture& d, size_t elements, unsigned block, Stream stream, bool rows) {
#define INPUTS d.conn[0].ptr, d.conn[1].ptr, d.conn[2].ptr, d.conn[3].ptr, \
    d.conn[4].ptr, d.conn[5].ptr, d.conn[6].ptr, d.conn[7].ptr, elements, \
    d.field[0].ptr, d.field[1].ptr, d.field[2].ptr, d.field[3].ptr, d.field[4].ptr, \
    d.field[5].ptr, d.field[6].ptr, d.field[7].ptr, d.field[8].ptr, \
    d.scs[0].ptr, d.scs[1].ptr, d.scs[2].ptr, d.scs[3].ptr, d.dofs.ptr, d.ownership.ptr, d.matrix.ptr, d.rhs.ptr
    // One block also exercises the kernels' empty-input return without dereferencing data.
    if (rows) {
        LAUNCH_KERNEL(cvfem_hex_assembly_kernel_tensor_early, std::max(size_t(1), elements), INPUTS);
    }
    else      { LAUNCH_KERNEL(cvfem_hex_assembly_kernel_tensor, std::max(size_t(1), elements), INPUTS); }
    CHECK_GPU(cudaGetLastError());
#undef INPUTS
}

void compare(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) throw std::runtime_error("size mismatch");
    for (size_t i = 0; i < a.size(); ++i)
        if (!std::isfinite(a[i]) || !std::isfinite(b[i]) ||
            std::abs(a[i]-b[i]) > 2e-12 * std::max({1., std::abs(a[i]), std::abs(b[i])})) {
            std::fprintf(stderr, "entry %zu: %.17g versus %.17g\n", i, a[i], b[i]);
            throw std::runtime_error("entry mismatch");
        }
}

void conservation(const Fixture& f, const std::vector<double>& values, const std::vector<double>& rhs) {
    double balance = 0, scale = 1;
    std::vector<double> columns(f.nodes, 0);
    for (size_t i = 0; i < values.size(); ++i) {
        columns[f.columns[i]] += values[i] - f.values[i];
        scale += std::abs(values[i] - f.values[i]);
    }
    for (int n = 0; n < f.nodes; ++n) { balance += rhs[n] - f.rhs[n]; scale += std::abs(rhs[n]-f.rhs[n]); }
    if (std::abs(balance) > 2e-12*scale) throw std::runtime_error("internal flux imbalance");
    for (double sum : columns) if (std::abs(sum) > 2e-12*scale) throw std::runtime_error("matrix column imbalance");
}
} // namespace

int main() {
    try {
        Stream stream{};
#if defined(__HIPCC__)
        CHECK_GPU(hipStreamCreateWithFlags(&stream, hipStreamNonBlocking));
#else
        CHECK_GPU(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
#endif
        int cases = 0;
        for (int elements : {1, 3, 17}) for (int variant = 0; variant < 4; ++variant)
        for (unsigned block : {128u, 256u, 512u}) {
            Fixture f(elements, variant);
            DeviceFixture baseline(f), rows(f);
            for (int repeat = 0; repeat < 2; ++repeat) {
                launch(baseline, elements, block, stream, false);
                launch(rows, elements, block, stream, true);
            }
            CHECK_GPU(cudaStreamSynchronize(stream));
            auto values = rows.values.get(), rhs = rows.rhs.get();
            compare(baseline.values.get(), values); compare(baseline.rhs.get(), rhs);
            if (variant < 2) conservation(f, values, rhs);
            if (variant == 2) {
                for (int row = f.owned_rows; row < f.nodes; ++row) {
                    if (rhs[row] != f.rhs[row]) throw std::runtime_error("ghost RHS modified");
                    for (int i = f.rowptr[row]; i < f.rowptr[row+1]; ++i)
                        if (values[i] != f.values[i]) throw std::runtime_error("ghost row modified");
                }
            }
            ++cases;
        }
        DeviceFixture empty(Fixture(1, 0));
        auto values = empty.values.get(), rhs = empty.rhs.get();
        launch(empty, 0, 128, stream, false); launch(empty, 0, 128, stream, true);
        CHECK_GPU(cudaStreamSynchronize(stream));
        compare(values, empty.values.get()); compare(rhs, empty.rhs.get());
        CHECK_GPU(cudaStreamDestroy(stream));
        std::printf("PASS tensor early: %d entrywise cases, conservation, ghost rows and empty input\n", cases);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "FAIL tensor assembly: %s\n", e.what()); return 1;
    }
}
