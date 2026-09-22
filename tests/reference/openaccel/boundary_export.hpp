#pragma once
#include "export_writer.hpp"
#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <map>
#include <utility>

namespace mars_reference {
template<class Bulk> class BoundaryExport {
    Bulk& bulk_;
    std::ofstream out_;
    std::size_t count_ = 0;
    template<class F> void checked(F&& f) {
        try { f(); }
        catch (const std::exception& e) {
            std::cerr << "Boundary export failed: " << e.what() << '\n';
            MPI_Abort(bulk_.parallel(), 1); std::abort();
        }
    }
    template<class Values> void array(const Values& values) {
        out_ << '['; bool first = true;
        for (auto value : values) {
            require(std::isfinite(double(value)), "nonfinite boundary export");
            if (!first) out_ << ',';
            first = false; out_ << value;
        }
        out_ << ']';
    }
public:
    BoundaryExport(Bulk& bulk, const char* stage, int dimension) : bulk_(bulk) {
        const char* root = std::getenv("MARS_OPENACCEL_EXPORT_DIR");
        if (!root || !*root) return;
        checked([&] {
            const char* fixture = std::getenv("MARS_OPENACCEL_PUBLIC_FIXTURE");
            require(dimension == 3 && fixture && std::string(fixture) == "public_channel", "boundary export requires public 3D channel");
            int ranks = 0;
            require(MPI_Comm_size(bulk.parallel(), &ranks) == MPI_SUCCESS, "MPI_Comm_size failed");
            static std::map<std::string, int> calls;
            const int call = ++calls[stage], rank = bulk.parallel_rank();
            require(std::filesystem::is_directory(root), "create fresh export root first");
            const auto directory = std::filesystem::path(root)/"boundary";
            std::filesystem::create_directory(directory);
            const auto path = directory/(std::string(stage)+".call"+std::to_string(call)+".rank"+std::to_string(rank)+".jsonl");
            require(!std::filesystem::exists(path), "boundary export already exists");
            out_.exceptions(std::ios::failbit | std::ios::badbit);
            out_.open(path); out_.imbue(std::locale::classic()); out_ << std::setprecision(17);
            out_ << "{\"kind\":\"header\",\"schema\":1,\"producer\":\"openaccel\",\"fixture\":\"public_channel\","
                 << "\"reference_revision\":\"" << reference_revision << "\",\"solver_revision\":\"" << solver_revision
                 << "\",\"stage\":\"" << stage << "\",\"call\":" << call << ",\"rank\":" << rank << ",\"ranks\":" << ranks << "}\n";
        });
    }
    bool active() const { return out_.is_open(); }
    template<class F> void capture(F&& f) { if (active()) checked(std::forward<F>(f)); }
    template<class Entity, class Nodes> void block(Entity side, const Nodes& nodes, int components,
            std::initializer_list<InputField> fields, const std::vector<double>& lhs, const std::vector<double>& rhs) {
        if (!active() || bulk_.parallel_owner_rank(side) != bulk_.parallel_rank()) return;
        checked([&] {
            require((nodes.size() == 3 || nodes.size() == 4) && (components == 1 || components == 3), "invalid boundary dimensions");
            const auto width = nodes.size()*components;
            require(lhs.size() == width*width && rhs.size() == width, "invalid boundary block size");
            std::vector<std::uint64_t> ids;
            for (auto node : nodes) ids.push_back(bulk_.identifier(node));
            out_ << "{\"kind\":\"face\",\"id\":" << bulk_.identifier(side) << ",\"nodes\":";
            array(ids); out_ << ",\"components\":" << components << ",\"inputs\":{";
            bool first = true; std::set<std::string> names;
            for (const auto& field : fields) {
                const std::string name(field.name);
                require(!name.empty() && name.find_first_not_of("abcdefghijklmnopqrstuvwxyz_") == std::string::npos
                        && names.insert(name).second && field.data && field.size, "invalid boundary field");
                if (!first) out_ << ',';
                first = false; out_ << '"' << name << "\":";
                array(std::vector<double>(field.data, field.data+field.size));
            }
            out_ << "},\"outputs\":{\"lhs\":"; array(lhs); out_ << ",\"rhs\":"; array(rhs);
            out_ << "}}\n"; ++count_;
        });
    }
    ~BoundaryExport() {
        if (active()) checked([&] { out_ << "{\"kind\":\"end\",\"records\":" << count_ << "}\n"; out_.close(); });
    }
};
} // namespace mars_reference
