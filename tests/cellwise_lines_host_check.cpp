// Runs cellwise::build_element_lines (backend/distributed/unstructured/solvers/
// mars_cellwise_lines_host.hpp) on E hexahedra read from files and writes its tables, for
// scripts/test_node_lines.py, which emulates the device build of the node blocks from them and
// compares every block entry with the assembled operator.
//
// Usage: cellwise_lines_host_check E corners.bin vid.bin outdir
//   corners.bin: E * 8 * 3 doubles (corner_bits order); vid.bin: E * 8 int64 vertex ids.
// Build from the MARS root:
//   c++ -std=c++20 -O2 -I. tests/cellwise_lines_host_check.cpp -o /tmp/cellwise_lines_host_check
#include <cstdio>
#include <cstdlib>
#include <string>

#include "backend/distributed/unstructured/solvers/mars_cellwise_lines_host.hpp"

namespace {

template <class T>
bool load(const std::string& path, std::vector<T>& v)
{
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    const size_t got = fread(v.data(), sizeof(T), v.size(), f);
    fclose(f);
    return got == v.size();
}

template <class T>
bool save(const std::string& path, const std::vector<T>& v)
{
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return false;
    const size_t put = fwrite(v.data(), sizeof(T), v.size(), f);
    fclose(f);
    return put == v.size();
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc != 5) {
        fprintf(stderr, "usage: %s E corners.bin vid.bin outdir\n", argv[0]);
        return 2;
    }
    const long long E = atoll(argv[1]);
    std::vector<double> corners(E * 24);
    std::vector<long long> vid(E * 8);
    if (!load(argv[2], corners) || !load(argv[3], vid)) {
        fprintf(stderr, "cannot read the mesh files\n");
        return 2;
    }
    const auto L = mars::cellwise::build_element_lines(corners, vid, E);
    const std::string out = std::string(argv[4]) + "/";
    const bool ok = save(out + "start.bin", L.start) && save(out + "elem.bin", L.elem) && save(out + "code.bin", L.code) &&
                    save(out + "face_elem.bin", L.face_elem) && save(out + "face_code.bin", L.face_code) &&
                    save(out + "edge_start.bin", L.edge_start) && save(out + "edge_elem.bin", L.edge_elem) &&
                    save(out + "edge_code.bin", L.edge_code);
    if (!ok) {
        fprintf(stderr, "cannot write the tables to %s\n", argv[4]);
        return 2;
    }
    printf("%lld lines over %lld elements, %zu edge entries\n", L.lines(), E, L.edge_elem.size());
    return 0;
}
