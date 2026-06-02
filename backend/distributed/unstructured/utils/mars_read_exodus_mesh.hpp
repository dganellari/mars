#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <tuple>
#include <algorithm>
#include <cstring>
#include <map>
#include <array>
#include <unordered_set>
#include "mars.hpp"

#ifdef MARS_HAVE_NETCDF
#include <netcdf.h>
#endif

namespace mars
{

#ifdef MARS_HAVE_NETCDF

// Helper macro for netCDF error checking
#define NC_CHECK(call) \
    do { \
        int status = (call); \
        if (status != NC_NOERR) { \
            throw std::runtime_error(std::string("NetCDF error: ") + nc_strerror(status)); \
        } \
    } while(0)

/**
 * @brief Read Exodus II mesh format and return data compatible with mars::ElementDomain
 *
 * This reads ExodusII (.exo) files directly using netCDF without requiring external tools.
 * Returns coordinates and connectivity in the same SoA format expected by ElementDomain.
 *
 * @tparam N Number of nodes per element (4 for tetrahedra, 8 for hexahedra)
 * @tparam RealType Type for coordinate data (float or double)
 * @tparam KeyType Type for connectivity indices
 * @param meshFile Path to Exodus mesh file (.exo)
 * @param rank Current MPI rank
 * @param numRanks Total MPI ranks
 * @return Tuple containing: node count, element count, x_coords, y_coords, z_coords,
 *         connectivity tuple, local_to_global_map, boundary_nodes
 */
template<int N, typename RealType = double, typename KeyType = unsigned>
inline auto readExodusMeshWithElementPartitioning(const std::string& meshFile, int rank, int numRanks)
{
    // Helper to create N-tuple of vectors
    auto createNVectors = [](size_t size) {
        if constexpr (N == 4) {
            return std::tuple<std::vector<KeyType>, std::vector<KeyType>,
                              std::vector<KeyType>, std::vector<KeyType>>(
                std::vector<KeyType>(size), std::vector<KeyType>(size),
                std::vector<KeyType>(size), std::vector<KeyType>(size));
        } else if constexpr (N == 8) {
            return std::tuple<std::vector<KeyType>, std::vector<KeyType>,
                              std::vector<KeyType>, std::vector<KeyType>,
                              std::vector<KeyType>, std::vector<KeyType>,
                              std::vector<KeyType>, std::vector<KeyType>>(
                std::vector<KeyType>(size), std::vector<KeyType>(size),
                std::vector<KeyType>(size), std::vector<KeyType>(size),
                std::vector<KeyType>(size), std::vector<KeyType>(size),
                std::vector<KeyType>(size), std::vector<KeyType>(size));
        }
    };

    // Open netCDF file
    int ncid;
    NC_CHECK(nc_open(meshFile.c_str(), NC_NOWRITE, &ncid));

    // Read dimension for number of nodes
    int num_nodes_dim;
    size_t total_nodes;
    NC_CHECK(nc_inq_dimid(ncid, "num_nodes", &num_nodes_dim));
    NC_CHECK(nc_inq_dimlen(ncid, num_nodes_dim, &total_nodes));

    // Read dimension for number of elements (check multiple possible names)
    int num_elem_dim;
    size_t total_elements = 0;

    // Try different dimension names for total elements
    if (nc_inq_dimid(ncid, "num_elem", &num_elem_dim) == NC_NOERR) {
        NC_CHECK(nc_inq_dimlen(ncid, num_elem_dim, &total_elements));
    } else {
        // Count elements from all element blocks
        int num_el_blk_dim;
        size_t num_el_blk;
        NC_CHECK(nc_inq_dimid(ncid, "num_el_blk", &num_el_blk_dim));
        NC_CHECK(nc_inq_dimlen(ncid, num_el_blk_dim, &num_el_blk));

        for (size_t blk = 1; blk <= num_el_blk; ++blk) {
            char dim_name[64];
            snprintf(dim_name, sizeof(dim_name), "num_el_in_blk%zu", blk);
            int dim_id;
            size_t num_el_in_blk;
            NC_CHECK(nc_inq_dimid(ncid, dim_name, &dim_id));
            NC_CHECK(nc_inq_dimlen(ncid, dim_id, &num_el_in_blk));
            total_elements += num_el_in_blk;
        }
    }

    if (rank == 0) {
        std::cout << "Exodus mesh: " << total_nodes << " nodes, " << total_elements << " elements" << std::endl;
    }

    // Read all coordinates (we need global coordinates for proper partitioning)
    std::vector<double> all_x(total_nodes), all_y(total_nodes), all_z(total_nodes);

    int coord_x_id, coord_y_id, coord_z_id;
    NC_CHECK(nc_inq_varid(ncid, "coordx", &coord_x_id));
    NC_CHECK(nc_inq_varid(ncid, "coordy", &coord_y_id));
    NC_CHECK(nc_inq_varid(ncid, "coordz", &coord_z_id));

    NC_CHECK(nc_get_var_double(ncid, coord_x_id, all_x.data()));
    NC_CHECK(nc_get_var_double(ncid, coord_y_id, all_y.data()));
    NC_CHECK(nc_get_var_double(ncid, coord_z_id, all_z.data()));

    // Read all element connectivity from element blocks
    std::vector<int> all_connectivity(total_elements * N);

    int num_el_blk_dim;
    size_t num_el_blk;
    NC_CHECK(nc_inq_dimid(ncid, "num_el_blk", &num_el_blk_dim));
    NC_CHECK(nc_inq_dimlen(ncid, num_el_blk_dim, &num_el_blk));

    size_t elem_offset = 0;
    for (size_t blk = 1; blk <= num_el_blk; ++blk) {
        char dim_name[64];
        snprintf(dim_name, sizeof(dim_name), "num_el_in_blk%zu", blk);
        int dim_id;
        size_t num_el_in_blk;
        NC_CHECK(nc_inq_dimid(ncid, dim_name, &dim_id));
        NC_CHECK(nc_inq_dimlen(ncid, dim_id, &num_el_in_blk));

        // Read connectivity for this block
        char var_name[64];
        snprintf(var_name, sizeof(var_name), "connect%zu", blk);
        int conn_id;
        NC_CHECK(nc_inq_varid(ncid, var_name, &conn_id));

        std::vector<int> block_conn(num_el_in_blk * N);
        NC_CHECK(nc_get_var_int(ncid, conn_id, block_conn.data()));

        // Copy to global connectivity array
        for (size_t e = 0; e < num_el_in_blk; ++e) {
            for (int n = 0; n < N; ++n) {
                all_connectivity[(elem_offset + e) * N + n] = block_conn[e * N + n];
            }
        }
        elem_offset += num_el_in_blk;
    }

    // Close netCDF file
    NC_CHECK(nc_close(ncid));

    // Partition elements across ranks
    size_t elemPerRank = total_elements / numRanks;
    size_t elemStartIdx = rank * elemPerRank;
    size_t elemEndIdx = (rank == numRanks - 1) ? total_elements : elemStartIdx + elemPerRank;
    size_t elementCount = elemEndIdx - elemStartIdx;

    // Collect unique nodes needed by this rank's elements
    std::vector<KeyType> uniqueNodes;
    {
        std::vector<bool> nodeUsed(total_nodes, false);
        for (size_t e = elemStartIdx; e < elemEndIdx; ++e) {
            for (int n = 0; n < N; ++n) {
                int globalNodeIdx = all_connectivity[e * N + n] - 1; // Exodus is 1-based
                nodeUsed[globalNodeIdx] = true;
            }
        }
        for (size_t i = 0; i < total_nodes; ++i) {
            if (nodeUsed[i]) {
                uniqueNodes.push_back(static_cast<KeyType>(i));
            }
        }
    }

    size_t nodeCount = uniqueNodes.size();

    // Create global-to-local mapping
    std::vector<KeyType> globalToLocal(total_nodes, static_cast<KeyType>(-1));
    for (size_t i = 0; i < nodeCount; ++i) {
        globalToLocal[uniqueNodes[i]] = static_cast<KeyType>(i);
    }

    // Extract coordinates for local nodes
    std::vector<RealType> x_coords(nodeCount), y_coords(nodeCount), z_coords(nodeCount);
    for (size_t i = 0; i < nodeCount; ++i) {
        KeyType globalIdx = uniqueNodes[i];
        x_coords[i] = static_cast<RealType>(all_x[globalIdx]);
        y_coords[i] = static_cast<RealType>(all_y[globalIdx]);
        z_coords[i] = static_cast<RealType>(all_z[globalIdx]);
    }

    // Create local connectivity using local node indices
    auto connectivity = createNVectors(elementCount);

    for (size_t e = 0; e < elementCount; ++e) {
        size_t globalElemIdx = elemStartIdx + e;
        if constexpr (N == 4) {
            std::get<0>(connectivity)[e] = globalToLocal[all_connectivity[globalElemIdx * N + 0] - 1];
            std::get<1>(connectivity)[e] = globalToLocal[all_connectivity[globalElemIdx * N + 1] - 1];
            std::get<2>(connectivity)[e] = globalToLocal[all_connectivity[globalElemIdx * N + 2] - 1];
            std::get<3>(connectivity)[e] = globalToLocal[all_connectivity[globalElemIdx * N + 3] - 1];
        } else if constexpr (N == 8) {
            std::get<0>(connectivity)[e] = globalToLocal[all_connectivity[globalElemIdx * N + 0] - 1];
            std::get<1>(connectivity)[e] = globalToLocal[all_connectivity[globalElemIdx * N + 1] - 1];
            std::get<2>(connectivity)[e] = globalToLocal[all_connectivity[globalElemIdx * N + 2] - 1];
            std::get<3>(connectivity)[e] = globalToLocal[all_connectivity[globalElemIdx * N + 3] - 1];
            std::get<4>(connectivity)[e] = globalToLocal[all_connectivity[globalElemIdx * N + 4] - 1];
            std::get<5>(connectivity)[e] = globalToLocal[all_connectivity[globalElemIdx * N + 5] - 1];
            std::get<6>(connectivity)[e] = globalToLocal[all_connectivity[globalElemIdx * N + 6] - 1];
            std::get<7>(connectivity)[e] = globalToLocal[all_connectivity[globalElemIdx * N + 7] - 1];
        }
    }

    // Local to global map (uniqueNodes already contains this)
    std::vector<KeyType> localToGlobal = uniqueNodes;

    // Boundary nodes (empty for now - could be extracted from side sets)
    std::vector<uint8_t> boundaryNodes(nodeCount, 0);

    if (rank == 0) {
        std::cout << "Rank 0: " << elementCount << " elements, " << nodeCount << " nodes (from Exodus)" << std::endl;
    }

    return std::make_tuple(nodeCount, elementCount,
                           std::move(x_coords), std::move(y_coords), std::move(z_coords),
                           std::move(connectivity), std::move(localToGlobal), std::move(boundaryNodes));
}

// Read all named side-sets from an Exodus mesh and return a map from side-set
// name -> sorted unique list of GLOBAL (1-based-removed -> 0-based) node IDs
// that belong to any face in that side-set.
//
// Every rank reads the SAME global side-set node lists (small enough to
// duplicate -- 7 side-sets * ~80k faces * 4 nodes = ~2.2M global IDs total,
// ~18 MB on each rank). The driver then walks its own local-to-global node
// map post-cstone-partition to mark which of its locally-owned nodes belong
// to which side-set. This avoids re-engineering cstone's partitioner for
// side-set carry-through.
//
// Hex face -> local-node convention follows Exodus II (1-based sides 1..6,
// 1-based local nodes); we use the 0-based form internally.
struct ExodusSideSets
{
    std::map<std::string, std::vector<uint64_t>> nodesByName;   // side-set name -> sorted unique GLOBAL 0-based node IDs
    // side-set name -> per-node coordinate (x,y,z), aligned 1:1 with
    // nodesByName[name]. The driver resolves an Exodus boundary node to its
    // RUNTIME node id by recomputing the SFC key from this coordinate -- because
    // after mesh ingest MARS/cstone index nodes by their SFC-sorted local id,
    // NOT by the Exodus node id. Resolving through the Exodus id (the old g2l)
    // was wrong on >1 rank (the two index spaces differ in order and size).
    std::map<std::string, std::vector<std::array<double, 3>>> nodeCoordsByName;
    // side-set name -> flat list of face triangle coordinates, 3 (x,y,z) per
    // face (size = 3*num_faces), in the same order as trianglesByName. Lets the
    // inlet-normal calc use coordinates directly instead of an index lookup.
    std::map<std::string, std::vector<std::array<double, 3>>> triangleCoordsByName;
    // side-set name -> flat list of face triangles, 3 GLOBAL 0-based node IDs
    // per face (size = 3*num_faces). Used to compute the inlet face normal so
    // an inlet velocity can be applied along the surface normal, not a global
    // axis. Tet faces are triangles; only the tet reader fills this.
    std::map<std::string, std::vector<uint64_t>> trianglesByName;
};

inline ExodusSideSets readExodusSideSetsHex8(const std::string& meshFile, int rank)
{
    ExodusSideSets out;

    // Hex face -> 4-node local indices (0-based). Exodus side k uses Exodus
    // 1-based local nodes; subtract 1 here.
    //   side 1 (y- "south"):   1,2,6,5  ->  0,1,5,4
    //   side 2 (x+ "east"):    2,3,7,6  ->  1,2,6,5
    //   side 3 (y+ "north"):   3,4,8,7  ->  2,3,7,6
    //   side 4 (x- "west"):    4,1,5,8  ->  3,0,4,7
    //   side 5 (z- "bottom"):  1,4,3,2  ->  0,3,2,1
    //   side 6 (z+ "top"):     5,6,7,8  ->  4,5,6,7
    static const int hexFaceNodes[6][4] = {
        {0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {3,0,4,7}, {0,3,2,1}, {4,5,6,7}
    };

    int ncid;
    NC_CHECK(nc_open(meshFile.c_str(), NC_NOWRITE, &ncid));

    // Number of side-sets (may be absent for meshes without boundary tags).
    int num_ss_dim;
    size_t num_side_sets = 0;
    if (nc_inq_dimid(ncid, "num_side_sets", &num_ss_dim) == NC_NOERR)
    {
        NC_CHECK(nc_inq_dimlen(ncid, num_ss_dim, &num_side_sets));
    }
    if (num_side_sets == 0)
    {
        nc_close(ncid);
        if (rank == 0) std::cout << "Exodus: 0 named side-sets in mesh\n";
        return out;
    }

    // Side-set names from ss_names(num_side_sets, len_name).
    int ss_names_id;
    NC_CHECK(nc_inq_varid(ncid, "ss_names", &ss_names_id));
    int len_name_dim;
    size_t len_name = 0;
    NC_CHECK(nc_inq_dimid(ncid, "len_name", &len_name_dim));
    NC_CHECK(nc_inq_dimlen(ncid, len_name_dim, &len_name));
    std::vector<char> raw_names(num_side_sets * len_name, 0);
    NC_CHECK(nc_get_var_text(ncid, ss_names_id, raw_names.data()));
    std::vector<std::string> ssNames(num_side_sets);
    for (size_t k = 0; k < num_side_sets; ++k)
    {
        const char* p = raw_names.data() + k * len_name;
        // length up to first NUL or len_name
        size_t n = strnlen(p, len_name);
        // Exodus pads names to len_name. strnlen stops at the first NUL but
        // KEEPS trailing spaces -- some writers space-pad, which would make the
        // stored key "name   " and silently fail an exact match against a
        // user-passed "name". Trim trailing whitespace (space/tab/CR/NL/NUL).
        while (n > 0 && (p[n-1] == ' ' || p[n-1] == '\t' ||
                         p[n-1] == '\r' || p[n-1] == '\n' || p[n-1] == '\0'))
            --n;
        ssNames[k] = std::string(p, n);
    }

    // Side-set IDs from ss_prop1 (e.g. [1,2,...]).
    int ss_prop1_id;
    NC_CHECK(nc_inq_varid(ncid, "ss_prop1", &ss_prop1_id));
    std::vector<long long> ssIds(num_side_sets);
    NC_CHECK(nc_get_var_longlong(ncid, ss_prop1_id, ssIds.data()));

    // Element-block connectivity: needed to expand (elem, side) -> 4 face nodes.
    // We use the same single-block path the main reader uses (num_el_blk==1
    // here). For multi-block meshes we would need elem-id -> block lookup; the
    // wing mesh is single-block so single-block is enough.
    int num_el_blk_dim;
    size_t num_el_blk = 0;
    NC_CHECK(nc_inq_dimid(ncid, "num_el_blk", &num_el_blk_dim));
    NC_CHECK(nc_inq_dimlen(ncid, num_el_blk_dim, &num_el_blk));
    if (num_el_blk != 1)
    {
        nc_close(ncid);
        throw std::runtime_error("readExodusSideSetsHex8: only single-element-block meshes supported (num_el_blk = "
                                 + std::to_string(num_el_blk) + ")");
    }
    int conn_id;
    NC_CHECK(nc_inq_varid(ncid, "connect1", &conn_id));
    int num_el_blk1_dim;
    size_t num_el_in_blk1 = 0;
    NC_CHECK(nc_inq_dimid(ncid, "num_el_in_blk1", &num_el_blk1_dim));
    NC_CHECK(nc_inq_dimlen(ncid, num_el_blk1_dim, &num_el_in_blk1));
    std::vector<int> connBlock1(num_el_in_blk1 * 8);
    NC_CHECK(nc_get_var_int(ncid, conn_id, connBlock1.data()));

    // For each side-set: read elem_ss<id> + side_ss<id>, expand to 4 face-nodes
    // per face, dedup, sort, store in out.nodesByName[name].
    for (size_t k = 0; k < num_side_sets; ++k)
    {
        std::string ssName = ssNames[k];

        // Exodus indexes the per-side-set variables (elem_ss<n>, side_ss<n>,
        // num_side_ss<n>) by SEQUENCE number n = 1..num_side_sets, NOT by the
        // ss_prop1 "ID" value. Earlier code keyed off the ID and only worked
        // when IDs happened to equal 1..N; meshes whose IDs differ (e.g.
        // 100,200) threw "Invalid dimension ID or name". Use the sequence
        // index; fall back to the ID form for any writer that used it.
        long long seq   = static_cast<long long>(k + 1);
        long long ssId  = ssIds[k];

        char elem_var[64], side_var[64], dim_var[64];
        int dim_id;
        size_t num_faces;

        snprintf(dim_var, sizeof(dim_var), "num_side_ss%lld", seq);
        if (nc_inq_dimid(ncid, dim_var, &dim_id) != NC_NOERR)
        {
            snprintf(dim_var, sizeof(dim_var), "num_side_ss%lld", ssId);
            NC_CHECK(nc_inq_dimid(ncid, dim_var, &dim_id));
            seq = ssId;   // this writer keyed by ID; use it for elem/side too
        }
        NC_CHECK(nc_inq_dimlen(ncid, dim_id, &num_faces));

        snprintf(elem_var, sizeof(elem_var), "elem_ss%lld", seq);
        snprintf(side_var, sizeof(side_var), "side_ss%lld", seq);

        int elem_id, side_id;
        NC_CHECK(nc_inq_varid(ncid, elem_var, &elem_id));
        NC_CHECK(nc_inq_varid(ncid, side_var, &side_id));
        std::vector<long long> elemList(num_faces), sideList(num_faces);
        NC_CHECK(nc_get_var_longlong(ncid, elem_id, elemList.data()));
        NC_CHECK(nc_get_var_longlong(ncid, side_id, sideList.data()));

        std::unordered_set<uint64_t> nodeSet;
        nodeSet.reserve(num_faces * 4);
        for (size_t f = 0; f < num_faces; ++f)
        {
            // 1-based element ID -> 0-based block index (single block here).
            long long elem1 = elemList[f];
            int       side1 = static_cast<int>(sideList[f]);
            if (elem1 < 1 || (uint64_t)elem1 > num_el_in_blk1 || side1 < 1 || side1 > 6) continue;
            size_t elem0 = static_cast<size_t>(elem1 - 1);
            const int* hex = &connBlock1[elem0 * 8];
            const int* faceLocal = hexFaceNodes[side1 - 1];
            for (int j = 0; j < 4; ++j)
            {
                // 1-based global node ID -> 0-based.
                int node1 = hex[faceLocal[j]];
                if (node1 < 1) continue;
                nodeSet.insert(static_cast<uint64_t>(node1 - 1));
            }
        }

        std::vector<uint64_t> nodeVec(nodeSet.begin(), nodeSet.end());
        std::sort(nodeVec.begin(), nodeVec.end());
        if (rank == 0)
        {
            std::cout << "Exodus side-set [" << ssName << "] id=" << ssId
                      << ": " << num_faces << " faces, " << nodeVec.size() << " unique nodes\n";
        }
        out.nodesByName[ssName] = std::move(nodeVec);
    }

    NC_CHECK(nc_close(ncid));
    return out;
}

// Tet4 counterpart of readExodusSideSetsHex8. Same Exodus side-set machinery
// (read names, ids, per-set (elem, side) pairs, expand to face nodes, dedup);
// the only element-specific parts are the connectivity stride (4 nodes/elem),
// the face->local-node table (4 triangular faces), and 3 nodes per face.
//
// Exodus II tet4 face convention (1-based sides 1..4, 1-based local nodes):
//   side 1: 1,2,4   side 2: 2,3,4   side 3: 1,4,3   side 4: 1,3,2
// (each side is the triangle NOT containing the node opposite it). Converted
// to 0-based below.
inline ExodusSideSets readExodusSideSetsTet4(const std::string& meshFile, int rank)
{
    ExodusSideSets out;

    static const int tetFaceNodes[4][3] = {
        {0,1,3}, {1,2,3}, {0,3,2}, {0,2,1}
    };

    int ncid;
    NC_CHECK(nc_open(meshFile.c_str(), NC_NOWRITE, &ncid));

    int num_ss_dim;
    size_t num_side_sets = 0;
    if (nc_inq_dimid(ncid, "num_side_sets", &num_ss_dim) == NC_NOERR)
    {
        NC_CHECK(nc_inq_dimlen(ncid, num_ss_dim, &num_side_sets));
    }
    if (num_side_sets == 0)
    {
        nc_close(ncid);
        if (rank == 0) std::cout << "Exodus: 0 named side-sets in mesh\n";
        return out;
    }

    int ss_names_id;
    NC_CHECK(nc_inq_varid(ncid, "ss_names", &ss_names_id));
    int len_name_dim;
    size_t len_name = 0;
    NC_CHECK(nc_inq_dimid(ncid, "len_name", &len_name_dim));
    NC_CHECK(nc_inq_dimlen(ncid, len_name_dim, &len_name));
    std::vector<char> raw_names(num_side_sets * len_name, 0);
    NC_CHECK(nc_get_var_text(ncid, ss_names_id, raw_names.data()));
    std::vector<std::string> ssNames(num_side_sets);
    for (size_t k = 0; k < num_side_sets; ++k)
    {
        const char* p = raw_names.data() + k * len_name;
        size_t n = strnlen(p, len_name);
        ssNames[k] = std::string(p, n);
    }

    int ss_prop1_id;
    NC_CHECK(nc_inq_varid(ncid, "ss_prop1", &ss_prop1_id));
    std::vector<long long> ssIds(num_side_sets);
    NC_CHECK(nc_get_var_longlong(ncid, ss_prop1_id, ssIds.data()));

    int num_el_blk_dim;
    size_t num_el_blk = 0;
    NC_CHECK(nc_inq_dimid(ncid, "num_el_blk", &num_el_blk_dim));
    NC_CHECK(nc_inq_dimlen(ncid, num_el_blk_dim, &num_el_blk));
    if (num_el_blk != 1)
    {
        nc_close(ncid);
        throw std::runtime_error("readExodusSideSetsTet4: only single-element-block meshes supported (num_el_blk = "
                                 + std::to_string(num_el_blk) + ")");
    }
    int conn_id;
    NC_CHECK(nc_inq_varid(ncid, "connect1", &conn_id));
    int num_el_blk1_dim;
    size_t num_el_in_blk1 = 0;
    NC_CHECK(nc_inq_dimid(ncid, "num_el_in_blk1", &num_el_blk1_dim));
    NC_CHECK(nc_inq_dimlen(ncid, num_el_blk1_dim, &num_el_in_blk1));
    std::vector<int> connBlock1(num_el_in_blk1 * 4);
    NC_CHECK(nc_get_var_int(ncid, conn_id, connBlock1.data()));

    // Global node coordinates, indexed by 0-based global node id (same as the
    // main mesh reader). Needed so each side-set node can be resolved to its
    // runtime SFC local id via its coordinate (the Exodus id is not the runtime
    // id). Read once for the whole file.
    size_t total_nodes = 0;
    {
        int num_nodes_dim;
        NC_CHECK(nc_inq_dimid(ncid, "num_nodes", &num_nodes_dim));
        NC_CHECK(nc_inq_dimlen(ncid, num_nodes_dim, &total_nodes));
    }
    std::vector<double> all_x(total_nodes), all_y(total_nodes), all_z(total_nodes);
    {
        int cx, cy, cz;
        NC_CHECK(nc_inq_varid(ncid, "coordx", &cx));
        NC_CHECK(nc_inq_varid(ncid, "coordy", &cy));
        NC_CHECK(nc_inq_varid(ncid, "coordz", &cz));
        NC_CHECK(nc_get_var_double(ncid, cx, all_x.data()));
        NC_CHECK(nc_get_var_double(ncid, cy, all_y.data()));
        NC_CHECK(nc_get_var_double(ncid, cz, all_z.data()));
    }
    auto coordOf = [&](uint64_t g0) -> std::array<double, 3> {
        return { all_x[g0], all_y[g0], all_z[g0] };
    };

    for (size_t k = 0; k < num_side_sets; ++k)
    {
        std::string ssName = ssNames[k];

        // Exodus indexes the per-side-set variables (elem_ss<n>, side_ss<n>,
        // num_side_ss<n>) by SEQUENCE number n = 1..num_side_sets, NOT by the
        // ss_prop1 "ID" value. Earlier code keyed off the ID and only worked
        // when IDs happened to equal 1..N; meshes whose IDs differ (e.g.
        // 100,200) threw "Invalid dimension ID or name". Use the sequence
        // index; fall back to the ID form for any writer that used it.
        long long seq   = static_cast<long long>(k + 1);
        long long ssId  = ssIds[k];

        char elem_var[64], side_var[64], dim_var[64];
        int dim_id;
        size_t num_faces;

        snprintf(dim_var, sizeof(dim_var), "num_side_ss%lld", seq);
        if (nc_inq_dimid(ncid, dim_var, &dim_id) != NC_NOERR)
        {
            snprintf(dim_var, sizeof(dim_var), "num_side_ss%lld", ssId);
            NC_CHECK(nc_inq_dimid(ncid, dim_var, &dim_id));
            seq = ssId;   // this writer keyed by ID; use it for elem/side too
        }
        NC_CHECK(nc_inq_dimlen(ncid, dim_id, &num_faces));

        snprintf(elem_var, sizeof(elem_var), "elem_ss%lld", seq);
        snprintf(side_var, sizeof(side_var), "side_ss%lld", seq);

        int elem_id, side_id;
        NC_CHECK(nc_inq_varid(ncid, elem_var, &elem_id));
        NC_CHECK(nc_inq_varid(ncid, side_var, &side_id));
        std::vector<long long> elemList(num_faces), sideList(num_faces);
        NC_CHECK(nc_get_var_longlong(ncid, elem_id, elemList.data()));
        NC_CHECK(nc_get_var_longlong(ncid, side_id, sideList.data()));

        std::unordered_set<uint64_t> nodeSet;
        nodeSet.reserve(num_faces * 3);
        std::vector<uint64_t> tris;       // 3 global node IDs per face
        tris.reserve(num_faces * 3);
        for (size_t f = 0; f < num_faces; ++f)
        {
            long long elem1 = elemList[f];
            int       side1 = static_cast<int>(sideList[f]);
            if (elem1 < 1 || (uint64_t)elem1 > num_el_in_blk1 || side1 < 1 || side1 > 4) continue;
            size_t elem0 = static_cast<size_t>(elem1 - 1);
            const int* tet = &connBlock1[elem0 * 4];
            const int* faceLocal = tetFaceNodes[side1 - 1];
            bool faceOk = true;
            uint64_t triNodes[3];
            for (int j = 0; j < 3; ++j)
            {
                int node1 = tet[faceLocal[j]];
                if (node1 < 1) { faceOk = false; continue; }
                uint64_t g0 = static_cast<uint64_t>(node1 - 1);
                nodeSet.insert(g0);
                triNodes[j] = g0;
            }
            // Keep the face triangle in the Exodus winding (tetFaceNodes is
            // outward-oriented), so the cross product gives the outward normal.
            if (faceOk)
                for (int j = 0; j < 3; ++j) tris.push_back(triNodes[j]);
        }

        std::vector<uint64_t> nodeVec(nodeSet.begin(), nodeSet.end());
        std::sort(nodeVec.begin(), nodeVec.end());
        // Per-node coords aligned with nodeVec, and per-triangle-node coords
        // aligned with tris, so the driver can resolve via coordinate->SFC key.
        std::vector<std::array<double, 3>> nodeCoords;
        nodeCoords.reserve(nodeVec.size());
        for (uint64_t g0 : nodeVec) nodeCoords.push_back(coordOf(g0));
        std::vector<std::array<double, 3>> triCoords;
        triCoords.reserve(tris.size());
        for (uint64_t g0 : tris) triCoords.push_back(coordOf(g0));
        if (rank == 0)
        {
            std::cout << "Exodus side-set [" << ssName << "] id=" << ssId
                      << ": " << num_faces << " faces, " << nodeVec.size() << " unique nodes\n";
        }
        out.nodesByName[ssName]         = std::move(nodeVec);
        out.nodeCoordsByName[ssName]    = std::move(nodeCoords);
        out.trianglesByName[ssName]     = std::move(tris);
        out.triangleCoordsByName[ssName] = std::move(triCoords);
    }

    NC_CHECK(nc_close(ncid));
    return out;
}

#undef NC_CHECK

#else // !MARS_HAVE_NETCDF

// Stub implementation when netCDF is not available
template<int N, typename RealType = double, typename KeyType = unsigned>
inline auto readExodusMeshWithElementPartitioning(const std::string& meshFile, int rank, int numRanks)
{
    // Helper to create N-tuple of vectors (must match the real implementation's return type)
    auto createNVectors = [](size_t size) {
        if constexpr (N == 4) {
            return std::tuple<std::vector<KeyType>, std::vector<KeyType>,
                              std::vector<KeyType>, std::vector<KeyType>>(
                std::vector<KeyType>(size), std::vector<KeyType>(size),
                std::vector<KeyType>(size), std::vector<KeyType>(size));
        } else if constexpr (N == 8) {
            return std::tuple<std::vector<KeyType>, std::vector<KeyType>,
                              std::vector<KeyType>, std::vector<KeyType>,
                              std::vector<KeyType>, std::vector<KeyType>,
                              std::vector<KeyType>, std::vector<KeyType>>(
                std::vector<KeyType>(size), std::vector<KeyType>(size),
                std::vector<KeyType>(size), std::vector<KeyType>(size),
                std::vector<KeyType>(size), std::vector<KeyType>(size),
                std::vector<KeyType>(size), std::vector<KeyType>(size));
        }
    };

    throw std::runtime_error("Exodus mesh reading requires netCDF. Please rebuild with netCDF support.");

    // Unreachable, but needed for return type deduction
    return std::make_tuple(size_t(0), size_t(0),
                           std::vector<RealType>(), std::vector<RealType>(), std::vector<RealType>(),
                           createNVectors(0), std::vector<KeyType>(), std::vector<uint8_t>());
}

struct ExodusSideSets
{
    std::map<std::string, std::vector<uint64_t>> nodesByName;
    std::map<std::string, std::vector<uint64_t>> trianglesByName;
};

inline ExodusSideSets readExodusSideSetsHex8(const std::string&, int)
{
    throw std::runtime_error("readExodusSideSetsHex8 requires netCDF.");
}

inline ExodusSideSets readExodusSideSetsTet4(const std::string&, int)
{
    throw std::runtime_error("readExodusSideSetsTet4 requires netCDF.");
}

#endif // MARS_HAVE_NETCDF

} // namespace mars
