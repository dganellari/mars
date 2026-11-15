#pragma once
#include "mars_nonconformal_scalar_runtime.hpp"
#include <bit>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>

// Input staging is replicated for this small-mesh example. It is not the
// ElementDomain distributed loader; file node IDs remain distinct on both sides.
static Mesh read_mesh(const std::string& path)
{
    if constexpr (std::endian::native != std::endian::little)
        throw std::runtime_error("prepared input requires a little-endian machine");
    std::ifstream in(path, std::ios::binary);
    auto read = [&](void* data, size_t bytes)
    {
        if (!in.read(static_cast<char*>(data), std::streamsize(bytes)))
            throw std::runtime_error("truncated prepared mesh");
    };
    char magic[8];
    uint64_t n, e, t, outer = 0;
    read(magic, 8);
    read(&n, 8);
    read(&e, 8);
    read(&t, 8);
    bool flow = std::memcmp(magic, "MARSNC02", 8) == 0;
    if (flow) read(&outer, 8);
    if ((std::memcmp(magic, "MARSNC01", 8) && !flow) || outer > 4 * e || (flow && !outer) || !e || e > 1000000 || !n ||
        n > 4 * e || !t || t > 4 * e)
        throw std::runtime_error("invalid prepared mesh header");
    if (std::filesystem::file_size(path) != 32 + (flow ? 8 : 0) + 32 * n + 24 * e + 12 * t + 12 * outer)
        throw std::runtime_error("prepared mesh length mismatch");
    static_assert(sizeof(Vec3<double>) == 24 && sizeof(int) == 4 && sizeof(InterfaceFaceTag) == 12);
    Mesh m;
    m.x.resize(n);
    m.cells.resize(4 * e);
    m.k.resize(e);
    m.boundary.resize(n);
    m.side.resize(n);
    m.tags.resize(t);
    read(m.x.data(), 24 * n);
    read(m.cells.data(), 16 * e);
    read(m.k.data(), 8 * e);
    read(m.boundary.data(), 4 * n);
    read(m.side.data(), 4 * n);
    read(m.tags.data(), 12 * t);
    m.outer.resize(outer);
    if (outer) read(m.outer.data(), 12 * outer);
    for (size_t i = 0; i < n; ++i)
        if (!std::isfinite(m.x[i].x) || !std::isfinite(m.x[i].y) || !std::isfinite(m.x[i].z) || m.boundary[i] < 0 ||
            m.boundary[i] > 2 || m.side[i] < 0 || m.side[i] > 1)
            throw std::runtime_error("invalid prepared node");
    for (int id : m.cells)
        if (id < 0 || uint64_t(id) >= n) throw std::runtime_error("invalid prepared connectivity");
    for (double k : m.k)
        if (!std::isfinite(k) || k <= 0) throw std::runtime_error("invalid diffusivity");
    for (auto tag : m.tags)
        if (tag.element < 0 || uint64_t(tag.element) >= e || tag.opposite < 0 || tag.opposite >= 4 || tag.side < 0 ||
            tag.side > 1)
            throw std::runtime_error("invalid prepared face tag");
    for (auto tag : m.outer)
        if (tag.element < 0 || uint64_t(tag.element) >= e || tag.opposite < 0 || tag.opposite >= 4 || tag.side < 0 ||
            tag.side > 2)
            throw std::runtime_error("invalid exterior flow face");
    m.dofs.resize(n);
    std::iota(m.dofs.begin(), m.dofs.end(), 0);
    m.elements.resize(e);
    std::iota(m.elements.begin(), m.elements.end(), 0);
    m.owners.resize(e);
    for (size_t i = 0; i < e; ++i)
        m.owners[i] = int(i % rankCount);
    return m;
}

// Each rank reads only its incident tetrahedra, boundary faces and node coordinates.
static Mesh read_mesh_partition(const std::string& directory)
{
    if constexpr (std::endian::native != std::endian::little)
        throw std::runtime_error("prepared partition requires a little-endian machine");
    auto path = std::filesystem::path(directory) / ("rank" + std::to_string(rankId) + ".mncp");
    std::ifstream in(path, std::ios::binary);
    auto read = [&](void* data, size_t bytes)
    {
        if (bytes && !in.read(static_cast<char*>(data), std::streamsize(bytes)))
            throw std::runtime_error("truncated mesh partition");
    };
    char magic[8];
    uint64_t h[9];
    read(magic, 8);
    read(h, sizeof(h));
    auto [ranks, rank, globalNodes, n, e, t, f, p, globalElements] =
        std::array<uint64_t, 9>{h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], h[8]};
    if (std::memcmp(magic, "MARSNP01", 8) || ranks != uint64_t(rankCount) || rank != uint64_t(rankId) ||
        !globalElements || globalElements > 1000000 || !globalNodes || globalNodes > 4 * globalElements ||
        n > globalNodes || e > globalElements || t > 4 * e || f > 4 * e || p >= ranks)
        throw std::runtime_error("invalid mesh partition header or MPI size");
    if (std::filesystem::file_size(path) != 80 + 40 * n + 36 * e + 12 * (t + f) + 4 * p)
        throw std::runtime_error("mesh partition length mismatch");
    Mesh m;
    m.partitioned = true;
    m.globalNodes = int(globalNodes);
    m.dofs.resize(n);
    m.elements.resize(e);
    m.x.resize(n);
    m.cells.resize(4 * e);
    m.k.resize(e);
    m.boundary.resize(n);
    m.side.resize(n);
    m.owners.resize(e);
    m.tags.resize(t);
    m.outer.resize(f);
    m.candidatePeers.resize(p);
    read(m.dofs.data(), 8 * n);
    read(m.elements.data(), 8 * e);
    read(m.x.data(), 24 * n);
    read(m.cells.data(), 16 * e);
    read(m.k.data(), 8 * e);
    read(m.boundary.data(), 4 * n);
    read(m.side.data(), 4 * n);
    read(m.owners.data(), 4 * e);
    read(m.tags.data(), 12 * t);
    read(m.outer.data(), 12 * f);
    read(m.candidatePeers.data(), 4 * p);
    auto unique_sorted = [](const auto& a)
    { return std::is_sorted(a.begin(), a.end()) && std::adjacent_find(a.begin(), a.end()) == a.end(); };
    if (!unique_sorted(m.dofs) || !unique_sorted(m.elements) || !unique_sorted(m.candidatePeers))
        throw std::runtime_error("partition identities must be sorted and unique");
    for (size_t i = 0; i < n; ++i)
        if (m.dofs[i] < 0 || uint64_t(m.dofs[i]) >= globalNodes || !std::isfinite(m.x[i].x) ||
            !std::isfinite(m.x[i].y) || !std::isfinite(m.x[i].z) || m.boundary[i] < 0 || m.boundary[i] > 2 ||
            m.side[i] < 0 || m.side[i] > 1)
            throw std::runtime_error("invalid partition node");
    for (int id : m.cells)
        if (id < 0 || uint64_t(id) >= n) throw std::runtime_error("invalid partition connectivity");
    for (size_t i = 0; i < e; ++i)
    {
        int firstOwner = rankCount;
        bool incident  = false;
        for (int k = 0; k < 4; ++k)
        {
            int r      = owner(m.dofs[m.cells[4 * i + k]], m.globalNodes);
            firstOwner = std::min(firstOwner, r);
            incident   = incident || r == rankId;
        }
        if (m.elements[i] >= globalElements || m.owners[i] != firstOwner || !incident || !std::isfinite(m.k[i]) ||
            m.k[i] <= 0)
            throw std::runtime_error("invalid partition element");
    }
    for (auto tag : m.tags)
        if (tag.element < 0 || uint64_t(tag.element) >= e || tag.opposite < 0 || tag.opposite >= 4 || tag.side < 0 ||
            tag.side > 1)
            throw std::runtime_error("invalid partition interface");
    for (auto tag : m.outer)
        if (tag.element < 0 || uint64_t(tag.element) >= e || tag.opposite < 0 || tag.opposite >= 4 || tag.side < 0 ||
            tag.side > 2)
            throw std::runtime_error("invalid partition boundary");
    for (int r : m.candidatePeers)
        if (r < 0 || r >= rankCount || r == rankId) throw std::runtime_error("invalid partition peer");
    for (int64_t id = int64_t(globalNodes) * rankId / rankCount; id < int64_t(globalNodes) * (rankId + 1) / rankCount;
         ++id)
        if (!std::binary_search(m.dofs.begin(), m.dofs.end(), id))
            throw std::runtime_error("missing owned partition node");
    return m;
}
