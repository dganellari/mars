#include "tet_gpu_checks.hpp"

namespace amr_test
{
std::unique_ptr<Domain> make_cube(int rank, int ranks)
{
    constexpr int n = 4;
    Domain::HostCoordsTuple coords;
    Domain::HostConnectivityTuple conn;
    auto index = [](int x, int y, int z) { return KeyType(x + (n + 1) * (y + (n + 1) * z)); };
    for (int z = 0; z <= n; ++z)
        for (int y = 0; y <= n; ++y)
            for (int x = 0; x <= n; ++x)
            {
                std::get<0>(coords).push_back(RealType(x) / n);
                std::get<1>(coords).push_back(RealType(y) / n);
                std::get<2>(coords).push_back(RealType(z) / n);
            }
    // The central tet uses globally odd-parity vertices, so adjacent cubes share face diagonals.
    constexpr int pattern[5][4] = {{1, 2, 4, 7}, {0, 1, 2, 4}, {3, 1, 2, 7}, {5, 1, 4, 7}, {6, 2, 4, 7}};
    int element                 = 0;
    for (int z = 0; z < n; ++z)
        for (int y = 0; y < n; ++y)
            for (int x = 0; x < n; ++x)
                for (const auto& cell : pattern)
                {
                    if (element++ % ranks != rank) continue;
                    std::array<KeyType, 4> nodes;
                    std::array<Point, 4> points;
                    for (int c = 0; c < 4; ++c)
                    {
                        int corner = cell[c] ^ ((x + y + z) & 1);
                        int ix = x + (corner & 1), iy = y + ((corner >> 1) & 1), iz = z + ((corner >> 2) & 1);
                        nodes[c]  = index(ix, iy, iz);
                        points[c] = {RealType(ix) / n, RealType(iy) / n, RealType(iz) / n};
                    }
                    if (determinant(points) < 0) std::swap(nodes[0], nodes[1]);
                    std::get<0>(conn).push_back(nodes[0]);
                    std::get<1>(conn).push_back(nodes[1]);
                    std::get<2>(conn).push_back(nodes[2]);
                    std::get<3>(conn).push_back(nodes[3]);
                }
    return std::make_unique<Domain>(coords, conn, rank, ranks, 8, false, 1);
}

std::vector<Tet> snapshot(const std::array<const KeyType*, 4>& conn,
                          size_t start,
                          size_t count,
                          const RealType* x,
                          const RealType* y,
                          const RealType* z,
                          size_t nodes,
                          const KeyType* keys,
                          int rank)
{
    auto hx = read_device(x, nodes), hy = read_device(y, nodes), hz = read_device(z, nodes);
    auto hk = keys ? read_device(keys, nodes) : std::vector<KeyType>(nodes, 0);
    std::array<std::vector<KeyType>, 4> hc;
    for (int c = 0; c < 4; ++c)
        hc[c] = read_device(conn[c] + start, count);
    std::vector<Tet> result(count);
    for (size_t e = 0; e < count; ++e)
    {
        result[e].rank = rank;
        for (int c = 0; c < 4; ++c)
        {
            auto node = hc[c][e];
            require(node < nodes, "Out-of-range tet connectivity");
            result[e].keys[c]   = hk[node];
            result[e].points[c] = {hx[node], hy[node], hz[node]};
        }
    }
    return result;
}

std::vector<Tet> snapshot(Domain& domain, int rank)
{
    const auto& conn = domain.getElementToNodeConnectivity();
    return snapshot(
        {std::get<0>(conn).data(), std::get<1>(conn).data(), std::get<2>(conn).data(), std::get<3>(conn).data()},
        domain.startIndex(), domain.localElementCount(), domain.getNodeX().data(), domain.getNodeY().data(),
        domain.getNodeZ().data(), domain.getNodeCount(), domain.getLocalToGlobalSfcMap().data(), rank);
}

std::array<Values, 2> initial_fields(Domain& domain)
{
    auto x = read_device(domain.getNodeX().data(), domain.getNodeCount());
    auto y = read_device(domain.getNodeY().data(), domain.getNodeCount());
    auto z = read_device(domain.getNodeZ().data(), domain.getNodeCount());
    std::vector<RealType> constant(x.size(), 2.5), linear(x.size());
    for (size_t i = 0; i < x.size(); ++i)
        linear[i] = affine({x[i], y[i], z[i]});
    return {upload(constant), upload(linear)};
}

void check_fields(const RealType* x,
                  const RealType* y,
                  const RealType* z,
                  size_t nodes,
                  const std::array<Values, 2>& fields,
                  RealType tolerance,
                  const std::vector<uint8_t>* current)
{
    require(fields[0].size() == nodes && fields[1].size() == nodes, "Incorrect transferred field size");
    auto hx = read_device(x, nodes), hy = read_device(y, nodes), hz = read_device(z, nodes);
    auto constant = read_device(fields[0].data(), nodes), linear = read_device(fields[1].data(), nodes);
    int invalid = 0;
    for (size_t i = 0; i < nodes; ++i)
        if (!current || (*current)[i])
            invalid |= !std::isfinite(constant[i]) || !std::isfinite(linear[i]) ||
                       std::abs(constant[i] - 2.5) > 1.e-13 ||
                       std::abs(linear[i] - affine({hx[i], hy[i], hz[i]})) > tolerance;
    MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    require(!invalid, "Constant or affine field transfer failed on owned/ghost nodes");
}

std::vector<uint8_t> current_nodes(Domain& domain)
{
    auto current         = read_device(domain.getNodeOwnershipMap().data(), domain.getNodeCount());
    const auto& topology = domain.getNodeHaloTopology();
    if (domain.numRanks() == 1)
        require(topology.peers_.empty() && topology.sendOffsets_ == std::vector<int>{0} &&
                    topology.recvOffsets_ == std::vector<int>{0},
                "Invalid single-rank halo topology");
    const auto& recv_ids = topology.recvNodeIds_;
    for (int node : read_device(recv_ids.data(), recv_ids.size()))
        current.at(node) = 1;
    return current;
}

void check_sparse_device_exchange(int rank, int ranks)
{
    auto value = [](int round, int source, int target, int item)
    {
        return (KeyType(0x5a) << 56) | (KeyType(round) << 48) | (KeyType(source) << 40) | (KeyType(target) << 32) |
               KeyType(item);
    };
    for (int round = 0; round < 3; ++round)
    {
        std::vector<int> peers, counts;
        std::vector<KeyType> packets;
        for (int peer = 0; peer < ranks; ++peer)
            if (peer != rank)
            {
                peers.push_back(peer);
                counts.push_back(11 + 7 * rank + 3 * peer + round * 251);
                for (int i = 0; i < counts.back(); ++i)
                    packets.push_back(value(round, rank, peer, i));
            }
        auto device_packets = upload(packets);
        cstone::DeviceVector<KeyType> device_received;
        auto messages = mars::sparseExchange<KeyType>(peers, counts, device_packets.data(), device_received,
                                                      0x3a10 + round, MPI_COMM_WORLD);
        auto received = read_device(device_received.data(), device_received.size());
        int invalid   = messages.size() != size_t(ranks - 1);
        for (const auto& message : messages)
        {
            int expected_count = 11 + 7 * message.source + 3 * rank + round * 251;
            invalid += message.count != expected_count || message.offset + message.count > received.size();
            for (int i = 0; i < message.count && message.offset + i < received.size(); ++i)
                invalid += received[message.offset + i] != value(round, message.source, rank, i);
        }
        if (invalid)
            std::cerr << "rank " << rank << ": device sparse exchange round=" << round << " invalid=" << invalid
                      << '\n';
        MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        require(!invalid, "Existing device sparse exchange changed records");
    }
    if (rank == 0) std::cout << "Existing device sparse exchange: three variable-size rounds PASS\n";
}

void check_domain_layout(Domain& domain, int rank, int ranks)
{
    // Inspect the domain before its lazy node halo is built, so a failure identifies the broken boundary.
    const auto& partition = domain.getDomain();
    int uncovered = 0;
    for (size_t i = 0; i < domain.getElementCount(); ++i)
    {
        int suppliers = i >= partition.startIndex() && i < partition.endIndex();
        for (const auto& range : partition.incomingHaloIndices())
            suppliers += i >= size_t(range.start()) && i < size_t(range.end());
        uncovered += suppliers != 1;
    }
    MPI_Allreduce(MPI_IN_PLACE, &uncovered, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    require(!uncovered, "Cornerstone element layout has missing or overlapping halo suppliers");
    check_mesh(gather(snapshot(domain, rank)));
    auto keys = read_device(domain.getLocalToGlobalSfcMap().data(), domain.getNodeCount());
    std::vector<KeyType> bounds(ranks + 1);
    for (int r = 0; r <= ranks; ++r)
        bounds[r] = domain.getDomain().assignment()[r];
    mars::SfcNodeOwner<KeyType, RealType> host_owner{domain.getBoundingBox(), domain.getDomain().box(), bounds.data(),
                                                     ranks};
    auto device_owner = domain.sfcNodeOwner();
    auto device_keys  = domain.getLocalToGlobalSfcMap().data();
    cstone::DeviceVector<int> device_owners(keys.size());
    auto output = device_owners.data();
    thrust::for_each_n(thrust::device, thrust::counting_iterator<size_t>(0), keys.size(),
                       [=] __device__(size_t i) { output[i] = ranks == 1 ? 0 : device_owner(device_keys[i]); });
    auto owners = read_device(device_owners.data(), keys.size());
    struct Identity
    {
        KeyType key;
        int owner;
        int present_on;
    };
    std::vector<Identity> local;
    int invalid = 0;
    for (size_t i = 0; i < keys.size(); ++i)
    {
        invalid += owners[i] != (ranks == 1 ? 0 : host_owner(keys[i]));
        local.push_back({keys[i], owners[i], rank});
    }
    MPI_Allreduce(MPI_IN_PLACE, &invalid, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    require(!invalid, "CPU and device SFC owners disagree before the node halo");
    std::map<KeyType, int> owner_by_key;
    std::map<KeyType, std::set<int>> present;
    for (const auto& identity : gather(local))
    {
        auto [it, inserted] = owner_by_key.emplace(identity.key, identity.owner);
        require(inserted || it->second == identity.owner, "Ranks disagree on the SFC owner before the node halo");
        present[identity.key].insert(identity.present_on);
    }
    for (const auto& [key, owner] : owner_by_key)
        require(present.at(key).count(owner), "SFC node owner lacks the node before the node halo");
    if (rank == 0)
        std::cout << "Domain layout: positive parent cells, CPU/device owner parity and owner presence PASS\n";
}

void poison_ghost_fields(Domain& domain, std::array<Values, 2>& fields)
{
    auto owned = read_device(domain.getNodeOwnershipMap().data(), domain.getNodeCount());
    for (auto& field : fields)
    {
        auto values = read_device(field.data(), field.size());
        for (size_t i = 0; i < values.size(); ++i)
            if (!owned[i]) values[i] = -12345;
        field = upload(values);
    }
}

std::vector<Node> check_identity(Domain& domain, int rank)
{
    size_t count = domain.getNodeCount();
    auto keys    = read_device(domain.getLocalToGlobalSfcMap().data(), count);
    auto x = read_device(domain.getNodeX().data(), count), y = read_device(domain.getNodeY().data(), count);
    auto z         = read_device(domain.getNodeZ().data(), count);
    auto ownership = read_device(domain.getNodeOwnershipMap().data(), count);
    std::vector<Node> local(count);
    for (size_t i = 0; i < count; ++i)
        local[i] = {keys[i], {x[i], y[i], z[i]}, ownership[i] == 1 ? rank : -1};
    std::map<KeyType, Point> points;
    std::map<Point, KeyType> identities;
    std::map<KeyType, int> owners;
    auto global = gather(local);
    for (const auto& node : global)
    {
        auto [it, inserted] = points.emplace(node.key, node.point);
        require(inserted || it->second == node.point, "Node identity has inconsistent coordinates across ranks");
        auto [coordinate, new_point] = identities.emplace(node.point, node.key);
        require(new_point || coordinate->second == node.key, "One coordinate has inconsistent node identities");
        if (node.owner >= 0) require(owners.emplace(node.key, node.owner).second, "Node has multiple owners");
    }
    require(owners.size() == points.size(), "A node has no owner");
    return global;
}

RealType affine_quantization_bound(const Domain& domain)
{
    const auto& box                = domain.getBoundingBox();
    const auto bits                = box.getBoxDimBits(cstone::maxTreeLevel<KeyType>{});
    RealType lengths[3]            = {box.lx(), box.ly(), box.lz()};
    constexpr RealType gradient[3] = {2, 3, 0.5};
    RealType bound                 = 0;
    for (int d = 0; d < 3; ++d)
        bound += gradient[d] * lengths[d] / RealType((uint64_t(1) << bits[d]) - 1);
    return bound;
}

void select_device(int rank)
{
    MPI_Comm local;
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &local);
    int local_rank = 0, devices = 0;
    MPI_Comm_rank(local, &local_rank);
    MPI_Comm_free(&local);
    cudaCheckCall(cudaGetDeviceCount(&devices));
    require(devices > 0, "No GPU is visible to this process");
    int device = devices == 1 ? 0 : local_rank % devices;
    cudaCheckCall(cudaSetDevice(device));
    cudaDeviceProp properties;
    cudaCheckCall(cudaGetDeviceProperties(&properties, device));
    char processor[MPI_MAX_PROCESSOR_NAME];
    int processor_length;
    MPI_Get_processor_name(processor, &processor_length);
    std::cout << "rank=" << rank << " host=" << std::string(processor, processor_length)
              << " pci=" << properties.pciDomainID << ':' << properties.pciBusID << ':' << properties.pciDeviceID
              << " local_rank=" << local_rank << " visible_devices=" << devices << " selected_device=" << device
              << " gpu=" << properties.name << "\n";
}

} // namespace amr_test
