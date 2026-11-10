#pragma once
#include "mars_gpu_runtime.hpp"
#include <algorithm>
#include <cstdint>
#include <tuple>
#include <type_traits>
namespace mars
{
struct TetTag
{
};
struct HexTag
{
};
template<class Tag, class Real, class Key, class Accelerator>
struct ElementDomain
{
    static constexpr int npe = std::is_same_v<Tag, TetTag> ? 4 : 8;
    struct Domain
    {
        size_t begin = 0, end = 0;
        size_t startIndex() const { return begin; }
        size_t endIndex() const { return end; }
    } domain;
    std::vector<uint8_t> own;
    std::vector<Key> keys;
    std::vector<Real> x, y, z;
    std::tuple<std::vector<Key>,
               std::vector<Key>,
               std::vector<Key>,
               std::vector<Key>,
               std::vector<Key>,
               std::vector<Key>,
               std::vector<Key>,
               std::vector<Key>>
        conn;
    ElementDomain(int rank, int ranks)
    {
        const int active = std::min(ranks, 2), nodes = 2 * npe;
        for (int i = 0; i < nodes; ++i)
        {
            const int id = nodes - 1 - i;
            keys.push_back(Key(100 + 3 * id));
            own.push_back(id % active == rank ? 1 : 0);
            x.push_back(Real(id) / 3);
            y.push_back(Real(id) / 7);
            z.push_back(Real(-id) / 11);
        }
        for (int e = 0; e < 2; ++e)
        {
            if (e % active != rank) continue;
            int k = 0;
            std::apply([&](auto&... v) { (v.push_back(Key(nodes - 1 - (e * npe + (k++ % npe)))), ...); }, conn);
            ++domain.end;
        }
    }
    size_t getNodeCount() const { return keys.size(); }
    const auto& getNodeOwnershipMap() const { return own; }
    const auto& getLocalToGlobalSfcMap() const { return keys; }
    const auto& getNodeX() const { return x; }
    const auto& getNodeY() const { return y; }
    const auto& getNodeZ() const { return z; }
    const auto& getDomain() const { return domain; }
    const auto& getElementToNodeConnectivity() const { return conn; }
};
} // namespace mars
