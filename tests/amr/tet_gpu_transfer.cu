#include "tet_gpu_checks.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_field_migration.hpp"

namespace amr_test
{
void check_forced_migration(Domain& domain, int rank, int ranks)
{
    auto nodes = check_identity(domain, rank);
    std::map<KeyType, Point> unique;
    std::map<KeyType, int> owners;
    size_t remote = 0;
    for (const auto& node : nodes)
        if (node.owner >= 0)
        {
            unique.emplace(node.key, node.point);
            owners.emplace(node.key, node.owner);
            if (node.owner != 0) ++remote;
        }
    std::vector<KeyType> keys;
    std::array<std::vector<RealType>, 2> values;
    if (rank == 0)
        for (const auto& [key, point] : unique)
        {
            keys.push_back(key);
            values[0].push_back(2.5);
            values[1].push_back(affine(point));
        }
    auto device_keys = upload(keys);
    std::array<Values, 2> migrated;
    for (int f = 0; f < 2; ++f)
    {
        auto device_values = upload(values[f]);
        mars::amr::migrate_amr_field(domain, device_keys, device_values, migrated[f]);
    }
    require(ranks == 1 || remote > 0, "Forced migration fixture has no remote owners");
    auto current = current_nodes(domain);
    check_fields(domain.getNodeX().data(), domain.getNodeY().data(), domain.getNodeZ().data(), domain.getNodeCount(),
                 migrated, 1.e-12, &current);
    if (rank == 0) std::cout << "Forced field migration: " << remote << " remote owner records PASS\n";

    keys.clear();
    values[0].clear();
    values[1].clear();
    std::vector<uint8_t> valid;
    for (const auto& [key, point] : unique)
    {
        bool owned = owners.at(key) == rank;
        keys.push_back(key);
        valid.push_back(owned);
        values[0].push_back(owned ? 2.5 : -12345);
        values[1].push_back(owned ? affine(point) : -12345);
    }
    device_keys       = upload(keys);
    auto device_valid = upload(valid);
    for (int f = 0; f < 2; ++f)
    {
        auto device_values = upload(values[f]);
        mars::amr::migrate_amr_field(domain, device_keys, device_values, migrated[f], device_valid.data());
    }
    check_fields(domain.getNodeX().data(), domain.getNodeY().data(), domain.getNodeZ().data(), domain.getNodeCount(),
                 migrated, 1.e-12, &current);
    if (rank == 0) std::cout << "Authoritative field migration: poisoned duplicate records excluded PASS\n";
}

} // namespace amr_test
