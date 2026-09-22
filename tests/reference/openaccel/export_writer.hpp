// Public reference instrumentation; no solver equations.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <initializer_list>
#include <limits>
#include <locale>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace mars_reference {
inline constexpr const char* reference_revision = "0d69041ba1afda63e9e4328d9e0d9834bba37756";
inline constexpr const char* solver_revision = "e351ba5eeaf3537dcc53d0aba09a8347f0a44cd0";

inline void require(bool ok, const char* message)
{
    if (!ok) throw std::runtime_error(message);
}

struct InputField {
    const char* name;
    const double* data;
    std::size_t size;
    InputField(const char* name_, const std::vector<double>& values)
        : name(name_), data(values.data()), size(values.size()) {}
    InputField(const char* name_, const double* data_, std::size_t size_)
        : name(name_), data(data_), size(size_) {}
};

class ExportWriter {
    std::ofstream out_;
    std::size_t records_ = 0;
    bool finished_ = false;
    bool frozen_inputs_ = false;

    template<class Values> void array(const Values& values)
    {
        out_ << '[';
        bool first = true;
        for (auto value : values) {
            require(std::isfinite(double(value)), "nonfinite reference export");
            if (!first) out_ << ',';
            out_ << value;
            first = false;
        }
        out_ << ']';
    }

public:
    ExportWriter(const std::filesystem::path& path, const std::string& fixture,
                 const std::string& stage, std::uint64_t call, int rank, int ranks,
                 const std::string& producer = "openaccel", bool frozen_inputs = false)
        : frozen_inputs_(frozen_inputs)
    {
        require(fixture == "unit_tet" || fixture == "skew_tet"
                || fixture == "two_tets" || fixture == "public_channel",
                "export requires a named public contract fixture");
        require(stage == "momentum.interior" || stage == "pressure.interior",
                "unsupported reference export stage");
        require(producer == "openaccel" || producer == "mars" || producer == "harness-test",
                "invalid export producer");
        require(call > 0 && ranks > 0 && rank >= 0 && rank < ranks, "invalid export rank/call");
        require(!std::filesystem::exists(path), "refusing to overwrite a reference export");
        out_.exceptions(std::ios::badbit | std::ios::failbit);
        out_.open(path);
        out_.imbue(std::locale::classic());
        out_ << std::setprecision(std::numeric_limits<double>::max_digits10);
        out_ << "{\"kind\":\"header\",\"schema\":" << (frozen_inputs_ ? 2 : 1)
             << ",\"producer\":\"" << producer
             << "\",\"reference_revision\":\"" << reference_revision
             << "\",\"solver_revision\":\"" << solver_revision
             << "\",\"fixture\":\"" << fixture << "\",\"stage\":\"" << stage
             << "\",\"call\":" << call << ",\"rank\":" << rank << ",\"ranks\":" << ranks
             << ",\"precision\":\"float64\",\"coverage\":\"local-interior-only\"}\n";
    }

    void inputs(std::uint64_t parent, const std::vector<std::uint64_t>& nodes,
                const std::vector<std::uint64_t>& edges, std::initializer_list<InputField> fields)
    {
        require(frozen_inputs_ && !finished_ && parent > 0 && nodes.size() == 4
                && edges.size() == 12, "invalid frozen input dimensions/state");
        std::set<std::string> names;
        for (const auto& field : fields) {
            const std::string name(field.name);
            require(!name.empty() && name.find_first_not_of("abcdefghijklmnopqrstuvwxyz_") == std::string::npos
                    && names.insert(name).second, "invalid/duplicate input field name");
            require(field.size > 0 && field.data, "empty input field");
            for (std::size_t i = 0; i < field.size; ++i)
                require(std::isfinite(field.data[i]), "nonfinite frozen input");
        }
        out_ << "{\"kind\":\"inputs\",\"parent\":" << parent << ",\"nodes\":";
        array(nodes);
        out_ << ",\"edges\":"; array(edges);
        out_ << ",\"fields\":{";
        bool first = true;
        for (const auto& field : fields) {
            if (!first) out_ << ',';
            first = false;
            out_ << '"' << field.name << "\":[";
            for (std::size_t i = 0; i < field.size; ++i) {
                if (i) out_ << ',';
                out_ << field.data[i];
            }
            out_ << ']';
        }
        out_ << "}}\n";
        ++records_;
    }

    void block(std::uint64_t parent, const std::vector<std::uint64_t>& nodes,
               int components, const std::vector<double>& lhs, const std::vector<double>& rhs)
    {
        require(!finished_, "export already finished");
        const auto n = nodes.size() * std::size_t(components);
        require(parent > 0 && nodes.size() == 4 && (components == 1 || components == 3)
                && lhs.size() == n*n && rhs.size() == n, "invalid Tet4 local block dimensions");
        require(std::set<std::uint64_t>(nodes.begin(), nodes.end()).size() == nodes.size()
                && *std::min_element(nodes.begin(), nodes.end()) > 0, "invalid global node IDs");
        out_ << "{\"kind\":\"block\",\"parent\":" << parent << ",\"components\":" << components
             << ",\"nodes\":";
        array(nodes);
        out_ << ",\"lhs\":"; array(lhs);
        out_ << ",\"rhs\":"; array(rhs);
        out_ << "}\n";
        ++records_;
    }

    void sample(std::uint64_t parent, std::uint64_t left, std::uint64_t right,
                double flux, const std::vector<double>& area)
    {
        require(!finished_ && parent > 0 && left > 0 && right > 0 && left != right
                && area.size() == 3 && std::isfinite(flux), "invalid interior sample");
        out_ << "{\"kind\":\"sample\",\"parent\":" << parent << ",\"left\":" << left
             << ",\"right\":" << right << ",\"flux\":" << flux << ",\"area\":";
        array(area);
        out_ << "}\n";
        ++records_;
    }

    void finish()
    {
        require(!finished_, "export already finished");
        out_ << "{\"kind\":\"end\",\"records\":" << records_ << "}\n";
        out_.flush();
        out_.close();
        finished_ = true;
    }
};
} // namespace mars_reference
