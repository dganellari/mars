#ifndef MARS_DOMAIN_GHOSTING_HPP
#define MARS_DOMAIN_GHOSTING_HPP

// Builds a GhostRegistry over an ElementDomain's nodes -- the MARS equivalent of
// `bulkData.create_ghosting(name)`.
//
// A free function taking the domain, NOT a method on ElementDomain: domain.hpp is a lower layer
// than fem/, so having it include the registry would invert that. Call it as
//
//     auto ghosting = createNodeGhosting(domain, "fluid-solid-interface", comm, participates);
//
// which reads the same as STK's call and keeps the layering intact.
//
// Everything the registry needs is already in the node halo topology:
//   - the OWNER RANK of every node. The domain itself only stores an owned/not-owned flag, but the
//     halo's recv list per peer says exactly which foreign rank each ghost came from, so the owner
//     map falls out of it.
//   - the peer list, which is by construction a superset of the ranks that can own a shared node --
//     a node's owner is always a mesh neighbour.
//
// Setup-time only. The registry's build is host-side (see its header for why), so this copies the
// halo lists and the SFC keys down once. Nothing here runs per step; the per-step exchange is
// GhostRegistry::forward/reverseAdd, which stays on the device.

#include <mpi.h>

#include <cstdint>
#include <string>
#include <vector>

#include "mars_ghost_registry.hpp"

namespace mars
{
namespace fem
{

// `participates` selects which nodes take part. Pass an empty vector to ghost every SHARED node
// (anything this rank sends or receives in the node halo) -- that is the "same as the node halo,
// but as a named registry" case. Pass your own mask to ghost a subset, e.g. only the nodes on one
// interface side, which is the whole point of having named ghostings.
template<typename Domain, typename KeyType>
GhostRegistry<KeyType> createNodeGhosting(const Domain& domain, const std::string& name,
                                          MPI_Comm comm,
                                          const std::vector<uint8_t>& participates = {})
{
    int myRank = 0;
    MPI_Comm_rank(comm, &myRank);

    const auto&  topo      = domain.getNodeHaloTopology();
    const size_t nodeCount = domain.getNodeCount();

    std::vector<int> h_recv(topo.recvNodeIds_.size());
    if (!h_recv.empty())
        cudaMemcpy(h_recv.data(), topo.recvNodeIds_.data(), h_recv.size() * sizeof(int),
                   cudaMemcpyDeviceToHost);
    std::vector<int> h_send(topo.sendNodeIds_.size());
    if (!h_send.empty())
        cudaMemcpy(h_send.data(), topo.sendNodeIds_.data(), h_send.size() * sizeof(int),
                   cudaMemcpyDeviceToHost);

    // Owner rank per node: mine unless it arrived from a peer, in which case that peer owns it.
    std::vector<int> owner(nodeCount, myRank);
    for (size_t p = 0; p < topo.peers_.size(); ++p)
        for (int k = topo.recvOffsets_[p]; k < topo.recvOffsets_[p + 1]; ++k)
        {
            const int local = h_recv[k];
            if (local >= 0 && size_t(local) < nodeCount) owner[local] = topo.peers_[p];
        }

    const auto& d_sfc = domain.getLocalToGlobalSfcMap();
    std::vector<KeyType> key(nodeCount);
    if (nodeCount > 0)
        cudaMemcpy(key.data(), d_sfc.data(), nodeCount * sizeof(KeyType), cudaMemcpyDeviceToHost);

    // Default mask: every node this rank sends or receives. A node must be marked on the OWNER as
    // well as on every rank that ghosts it, or the two sides disagree and isConsistent() says so --
    // taking it from both halo lists satisfies that by construction.
    std::vector<uint8_t> part;
    if (participates.empty())
    {
        part.assign(nodeCount, 0);
        for (int local : h_recv)
            if (local >= 0 && size_t(local) < nodeCount) part[local] = 1;
        for (int local : h_send)
            if (local >= 0 && size_t(local) < nodeCount) part[local] = 1;
    }
    else { part = participates; }

    GhostRegistry<KeyType> reg(name);
    reg.build(nodeCount, owner, key, part, myRank, topo.peers_, comm);
    return reg;
}

}  // namespace fem
}  // namespace mars

#endif  // MARS_DOMAIN_GHOSTING_HPP
