#pragma once
#include "export_writer.hpp"
#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <map>

namespace mars_reference {
class NodeExport {
    MPI_Comm comm_;
    std::ofstream out_;
    std::size_t count_ = 0;

    template<class F> void checked(F&& action) {
        try { action(); }
        catch (const std::exception& e) {
            std::cerr << "Node export failed: " << e.what() << '\n';
            MPI_Abort(comm_, 1);
            std::abort();
        }
    }
    void fields(std::initializer_list<InputField> values) {
        out_ << '{';
        bool first = true;
        std::set<std::string> names;
        for (const auto& f : values) {
            const std::string name(f.name);
            require(!name.empty() && name.find_first_not_of("abcdefghijklmnopqrstuvwxyz_") == std::string::npos
                    && names.insert(name).second && f.data && f.size, "invalid node field");
            if (!first) out_ << ',';
            first = false;
            out_ << '"' << name << "\":[";
            for (std::size_t i = 0; i < f.size; ++i) {
                require(std::isfinite(f.data[i]), "nonfinite node value");
                if (i) out_ << ',';
                out_ << f.data[i];
            }
            out_ << ']';
        }
        out_ << '}';
    }
public:
    NodeExport(MPI_Comm comm, const char* stage, bool enabled = true) : comm_(comm) {
        const char* root = std::getenv("MARS_OPENACCEL_EXPORT_DIR");
        if (!enabled || !root || !*root) return;
        checked([&] {
            const char* fixture = std::getenv("MARS_OPENACCEL_PUBLIC_FIXTURE");
            require(fixture && std::string(fixture) == "public_channel", "node export requires the public channel");
            int rank = 0, ranks = 0;
            require(MPI_Comm_rank(comm, &rank) == MPI_SUCCESS && MPI_Comm_size(comm, &ranks) == MPI_SUCCESS,
                    "node export communicator failure");
            static std::map<std::string, int> calls;
            const int call = ++calls[stage];
            const auto directory = std::filesystem::path(root)/"nodes";
            require(std::filesystem::is_directory(root), "create fresh export root first");
            std::filesystem::create_directory(directory);
            const auto path = directory/(std::string(stage)+".call"+std::to_string(call)+".rank"+std::to_string(rank)+".jsonl");
            require(!std::filesystem::exists(path), "node export exists");
            out_.exceptions(std::ios::failbit | std::ios::badbit);
            out_.open(path);
            out_.imbue(std::locale::classic());
            out_ << std::setprecision(17) << "{\"schema\":1,\"kind\":\"header\",\"producer\":\"openaccel\","
                 << "\"fixture\":\"public_channel\",\"reference_revision\":\"" << reference_revision
                 << "\",\"solver_revision\":\"" << solver_revision << "\",\"stage\":\"" << stage
                 << "\",\"call\":" << call << ",\"rank\":" << rank << ",\"ranks\":" << ranks << "}\n";
        });
    }
    bool active() const { return out_.is_open(); }
    void node(std::uint64_t id, std::initializer_list<InputField> inputs,
              std::initializer_list<InputField> outputs) {
        if (!active()) return;
        checked([&] {
            require(id > 0, "invalid node id");
            out_ << "{\"kind\":\"node\",\"id\":" << id << ",\"inputs\":";
            fields(inputs); out_ << ",\"outputs\":"; fields(outputs); out_ << "}\n";
            ++count_;
        });
    }
    ~NodeExport() {
        if (active()) checked([&] { out_ << "{\"kind\":\"end\",\"records\":" << count_ << "}\n"; out_.close(); });
    }
};
} // namespace mars_reference
