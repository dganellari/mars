// Opt-in adapter around actual OpenAccel local assembly.
#pragma once
#include "export_writer.hpp"
#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <type_traits>

namespace mars_reference {
template<class Bulk> class ExportScope {
    Bulk& bulk_;
    std::unique_ptr<ExportWriter> writer_;

    template<class F> void checked(F&& operation)
    {
        try { operation(); }
        catch (const std::exception& error) {
            std::cerr << "Reference export failed: " << error.what() << '\n';
            MPI_Abort(bulk_.parallel(), 1);
            std::abort();
        }
    }

public:
    ExportScope(Bulk& bulk, const char* stage, int spatial_dimension) : bulk_(bulk)
    {
        const char* directory = std::getenv("MARS_OPENACCEL_EXPORT_DIR");
        if (!directory || !*directory) return;
        checked([&] {
            require(spatial_dimension == 3, "reference export requires a 3D build");
            const char* fixture = std::getenv("MARS_OPENACCEL_PUBLIC_FIXTURE");
            require(fixture && *fixture, "set MARS_OPENACCEL_PUBLIC_FIXTURE for a public fixture only");
            static std::map<std::string, std::uint64_t> calls;
            const auto call = ++calls[stage];
            int ranks = 0;
            require(MPI_Comm_size(bulk_.parallel(), &ranks) == MPI_SUCCESS, "MPI_Comm_size failed");
            const int rank = bulk_.parallel_rank();
            const auto path = std::filesystem::path(directory);
            // Every rank emits a header/footer, including an empty element partition.
            require(std::filesystem::is_directory(path), "create a fresh export directory first");
            writer_ = std::make_unique<ExportWriter>(
                path/(std::string(stage)+".call"+std::to_string(call)+".rank"+std::to_string(rank)+".jsonl"),
                fixture, stage, call, rank, ranks, "openaccel", true);
        });
    }

    ~ExportScope() { if (writer_) checked([&] { writer_->finish(); }); }

    bool active() const { return bool(writer_); }

    template<class Entity, class Nodes, class Index>
    void inputs(Entity parent, const Nodes& nodes, const Index* adjacent, int samples,
                std::initializer_list<InputField> fields)
    {
        if (!writer_ || bulk_.parallel_owner_rank(parent) != bulk_.parallel_rank()) return;
        checked([&] {
            require(nodes.size() == 4 && samples == 6, "frozen inputs require Tet4");
            std::vector<std::uint64_t> ids, edges;
            for (auto node : nodes) ids.push_back(bulk_.identifier(node));
            for (int i = 0; i < 2*samples; ++i) {
                require(adjacent[i] >= 0 && adjacent[i] < 4, "invalid sample endpoint");
                edges.push_back(ids[adjacent[i]]);
            }
            writer_->inputs(bulk_.identifier(parent), ids, edges, fields);
        });
    }

    template<class Entity, class Nodes, class Values>
    void block(Entity parent, const Nodes& nodes, int components, const Values& lhs, const Values& rhs)
    {
        static_assert(std::is_same_v<typename Values::value_type, double>, "reference export requires FP64");
        if (!writer_ || bulk_.parallel_owner_rank(parent) != bulk_.parallel_rank()) return;
        checked([&] {
            std::vector<std::uint64_t> ids;
            for (auto node : nodes) ids.push_back(bulk_.identifier(node));
            writer_->block(bulk_.identifier(parent), ids, components, lhs, rhs);
        });
    }

    template<class Entity, class Scalar>
    void sample(Entity parent, Entity left, Entity right, Scalar flux, const Scalar* area)
    {
        static_assert(std::is_same_v<Scalar, double>, "reference export requires FP64");
        if (!writer_ || bulk_.parallel_owner_rank(parent) != bulk_.parallel_rank()) return;
        checked([&] {
            writer_->sample(bulk_.identifier(parent), bulk_.identifier(left), bulk_.identifier(right),
                            flux, std::vector<double>(area, area+3));
        });
    }
};
} // namespace mars_reference
