#pragma once

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <mpi.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#include <cuda_runtime.h>
#include "backend/distributed/unstructured/domain.hpp"

namespace mars
{
namespace fem
{

// Field kind for multi-field VTU output. Free enum (not nested in the
// templated writer class) so callers can refer to it without spelling out
// VTUParallelWriter<KeyType, RealType>::FieldKind.
enum class VtuFieldKind { PointScalar, PointVector3, CellScalar };

// Parallel VTU writer for hex8 unstructured meshes. Emits:
//   - <prefix>_<step>_r<rank>.vtu (one per rank, owned elements + their nodes)
//   - <prefix>_<step>.pvtu        (master, rank 0 only)
//   - <prefix>.pvd                (time-collection, rank 0 only, appended)
//
// The per-rank file holds OWNED elements and the nodes they touch (including
// ghost nodes that those elements need for vertex coords + scalar values).
// ParaView's pvtu reader merges them; duplicate ghost nodes are stitched by
// position. We don't try to dedupe across ranks - simpler and ParaView handles it.
//
// Hex8 only for now. VTK cell type 12.
template<typename KeyType, typename RealType>
class VTUParallelWriter
{
public:
    using DomainT = ElementDomain<HexTag, RealType, KeyType, cstone::GpuTag>;

    explicit VTUParallelWriter(const std::string& prefix) : prefix_(prefix) {}

    // Write one timestep / AMR-level frame.
    //   step      = frame index (e.g. AMR level)
    //   timeValue = value in the .pvd timeline (e.g. step as double, or wall time)
    //   solution  = device per-node scalar (size = nodeCount, includes ghost slots)
    //   field     = name of the scalar in the output ("u", "pressure", etc.)
    void writeFrame(int step,
                    double timeValue,
                    const DomainT& domain,
                    const cstone::DeviceVector<RealType>& solution,
                    const std::string& field = "u") const
    {
        int rank = 0, numRanks = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

        writeRankPiece(rank, step, domain, solution, field);

        if (rank == 0)
        {
            writeMaster(step, numRanks, field);
            appendToPvd(step, timeValue);
        }
    }

    // Multi-field writeFrame: emits any combination of per-node scalars,
    // per-node 3-vectors, and per-element scalars into one .vtu file. Use this
    // for richer ParaView output (velocity vector, vorticity, refinement
    // level, etc.) without writing 5 separate PVD timelines.
    //
    // Each field is described by a FieldDesc struct.
    using FieldKind = VtuFieldKind;  // alias so existing callers can spell FieldKind
    struct FieldDesc
    {
        std::string name;
        VtuFieldKind kind;
        const cstone::DeviceVector<RealType>* data;
        const cstone::DeviceVector<RealType>* data_y;  // PointVector3 only
        const cstone::DeviceVector<RealType>* data_z;  // PointVector3 only
        using Kind = VtuFieldKind;                     // FD::Kind alias
    };

    void writeMultiFieldFrame(int step,
                              double timeValue,
                              const DomainT& domain,
                              const std::vector<FieldDesc>& fields) const
    {
        int rank = 0, numRanks = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

        writeRankPieceMulti(rank, step, domain, fields);

        if (rank == 0)
        {
            writeMasterMulti(step, numRanks, fields);
            appendToPvd(step, timeValue);
        }
    }

private:
    std::string prefix_;

    static std::string rankPieceName(const std::string& prefix, int step, int rank)
    {
        std::ostringstream oss;
        oss << prefix << "_step" << std::setw(4) << std::setfill('0') << step
            << "_r" << std::setw(4) << std::setfill('0') << rank << ".vtu";
        return oss.str();
    }

    static std::string masterName(const std::string& prefix, int step)
    {
        std::ostringstream oss;
        oss << prefix << "_step" << std::setw(4) << std::setfill('0') << step << ".pvtu";
        return oss.str();
    }

    static std::string pvdName(const std::string& prefix) { return prefix + ".pvd"; }

    void writeRankPiece(int rank,
                        int step,
                        const DomainT& domain,
                        const cstone::DeviceVector<RealType>& solution,
                        const std::string& field) const
    {
        // Pull mesh + solution to host. cstone partitions by element index; owned
        // range is [startIndex, endIndex). We write only owned elements; the
        // nodes they reference are emitted via direct local index (including
        // ghosts). Across ranks the same node appears in multiple pieces;
        // ParaView pvtu merges by coordinate.
        const auto& cstoneDom = domain.getDomain();
        const size_t elemStart = static_cast<size_t>(cstoneDom.startIndex());
        const size_t elemEnd   = static_cast<size_t>(cstoneDom.endIndex());
        const size_t numOwnedElem = elemEnd - elemStart;
        const size_t nodeCount    = domain.getNodeCount();

        // Coords (full local node array; we'll renumber to a compact set below).
        // cstone::DeviceVector::begin()/end() aren't const, so go through
        // rawPtr + cudaMemcpy which works on const refs.
        const auto& d_x = domain.getNodeX();
        const auto& d_y = domain.getNodeY();
        const auto& d_z = domain.getNodeZ();
        std::vector<RealType> hx(nodeCount), hy(nodeCount), hz(nodeCount);
        cudaMemcpy(hx.data(), thrust::raw_pointer_cast(d_x.data()), nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hy.data(), thrust::raw_pointer_cast(d_y.data()), nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hz.data(), thrust::raw_pointer_cast(d_z.data()), nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);

        // Solution
        std::vector<RealType> hu(nodeCount);
        cudaMemcpy(hu.data(), thrust::raw_pointer_cast(solution.data()), nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);

        // Connectivity (8 per hex, full element range; we slice owned range)
        const auto& d_conn = domain.getElementToNodeConnectivity();
        size_t fullElemCount = std::get<0>(d_conn).size();
        std::vector<KeyType> hc[8];
        for (int k = 0; k < 8; ++k) hc[k].resize(fullElemCount);
        cudaMemcpy(hc[0].data(), thrust::raw_pointer_cast(std::get<0>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[1].data(), thrust::raw_pointer_cast(std::get<1>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[2].data(), thrust::raw_pointer_cast(std::get<2>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[3].data(), thrust::raw_pointer_cast(std::get<3>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[4].data(), thrust::raw_pointer_cast(std::get<4>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[5].data(), thrust::raw_pointer_cast(std::get<5>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[6].data(), thrust::raw_pointer_cast(std::get<6>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[7].data(), thrust::raw_pointer_cast(std::get<7>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);

        // Mark which local node IDs are referenced by owned elements; compact-renumber.
        std::vector<int> remap(nodeCount, -1);
        std::vector<int> compact;
        compact.reserve(numOwnedElem * 8);  // upper bound; many duplicates
        for (size_t e = elemStart; e < elemEnd; ++e)
        {
            for (int k = 0; k < 8; ++k)
            {
                int nid = static_cast<int>(hc[k][e]);
                if (remap[nid] < 0)
                {
                    remap[nid] = static_cast<int>(compact.size());
                    compact.push_back(nid);
                }
            }
        }
        const size_t numWrittenNodes = compact.size();

        // Open file
        const std::string fname = rankPieceName(prefix_, step, rank);
        std::ofstream f(fname);
        if (!f) throw std::runtime_error("Cannot open " + fname);

        f << "<?xml version=\"1.0\"?>\n";
        f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        f << "<UnstructuredGrid>\n";
        f << "<Piece NumberOfPoints=\"" << numWrittenNodes
          << "\" NumberOfCells=\"" << numOwnedElem << "\">\n";

        // Points
        f << "<Points>\n<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (size_t i = 0; i < numWrittenNodes; ++i)
        {
            int n = compact[i];
            f << hx[n] << " " << hy[n] << " " << hz[n] << "\n";
        }
        f << "</DataArray>\n</Points>\n";

        // Cells: connectivity in VTK hex ordering. Mesh is already VTK-ordered
        // (i0..i7 corners), so we emit as-is, remapped to compact indices.
        f << "<Cells>\n";
        f << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (size_t e = elemStart; e < elemEnd; ++e)
        {
            for (int k = 0; k < 8; ++k)
            {
                f << remap[static_cast<int>(hc[k][e])] << " ";
            }
            f << "\n";
        }
        f << "</DataArray>\n";

        f << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        for (size_t i = 1; i <= numOwnedElem; ++i) f << (8 * i) << "\n";
        f << "</DataArray>\n";

        f << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (size_t i = 0; i < numOwnedElem; ++i) f << "12\n";  // VTK_HEXAHEDRON
        f << "</DataArray>\n";
        f << "</Cells>\n";

        // Point data: scalar solution. Emitted as the only data array so
        // ParaView auto-selects it; an additional cell array (e.g. rank id)
        // gets picked first by some ParaView versions and hides the solution.
        f << "<PointData Scalars=\"" << field << "\">\n";
        f << "<DataArray type=\"Float32\" Name=\"" << field << "\" format=\"ascii\">\n";
        for (size_t i = 0; i < numWrittenNodes; ++i) f << hu[compact[i]] << "\n";
        f << "</DataArray>\n";
        f << "</PointData>\n";

        f << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    }

    // Multi-field per-rank piece. Mirrors writeRankPiece but emits multiple
    // PointData and CellData arrays in one .vtu file.
    void writeRankPieceMulti(int rank,
                             int step,
                             const DomainT& domain,
                             const std::vector<FieldDesc>& fields) const
    {
        const auto& cstoneDom = domain.getDomain();
        const size_t elemStart    = static_cast<size_t>(cstoneDom.startIndex());
        const size_t elemEnd      = static_cast<size_t>(cstoneDom.endIndex());
        const size_t numOwnedElem = elemEnd - elemStart;
        const size_t nodeCount    = domain.getNodeCount();

        const auto& d_x = domain.getNodeX();
        const auto& d_y = domain.getNodeY();
        const auto& d_z = domain.getNodeZ();
        std::vector<RealType> hx(nodeCount), hy(nodeCount), hz(nodeCount);
        cudaMemcpy(hx.data(), thrust::raw_pointer_cast(d_x.data()), nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hy.data(), thrust::raw_pointer_cast(d_y.data()), nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hz.data(), thrust::raw_pointer_cast(d_z.data()), nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);

        // Pull each field to host.
        // PointScalar / PointVector3 fields are sized nodeCount; CellScalar
        // fields are sized full element count and we slice owned range.
        struct HostField {
            const FieldDesc* desc;
            std::vector<RealType> x, y, z;  // y, z used only for Vector3
        };
        std::vector<HostField> hfields;
        hfields.reserve(fields.size());
        for (const auto& fd : fields)
        {
            HostField hf{ &fd, {}, {}, {} };
            using K = typename FieldDesc::Kind;
            if (fd.kind == K::PointScalar)
            {
                hf.x.resize(nodeCount);
                cudaMemcpy(hf.x.data(), thrust::raw_pointer_cast(fd.data->data()),
                           nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
            }
            else if (fd.kind == K::PointVector3)
            {
                hf.x.resize(nodeCount); hf.y.resize(nodeCount); hf.z.resize(nodeCount);
                cudaMemcpy(hf.x.data(), thrust::raw_pointer_cast(fd.data->data()),
                           nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
                cudaMemcpy(hf.y.data(), thrust::raw_pointer_cast(fd.data_y->data()),
                           nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
                cudaMemcpy(hf.z.data(), thrust::raw_pointer_cast(fd.data_z->data()),
                           nodeCount * sizeof(RealType), cudaMemcpyDeviceToHost);
            }
            else if (fd.kind == K::CellScalar)
            {
                size_t sz = fd.data->size();
                hf.x.resize(sz);
                cudaMemcpy(hf.x.data(), thrust::raw_pointer_cast(fd.data->data()),
                           sz * sizeof(RealType), cudaMemcpyDeviceToHost);
            }
            hfields.push_back(std::move(hf));
        }

        const auto& d_conn = domain.getElementToNodeConnectivity();
        size_t fullElemCount = std::get<0>(d_conn).size();
        std::vector<KeyType> hc[8];
        for (int k = 0; k < 8; ++k) hc[k].resize(fullElemCount);
        cudaMemcpy(hc[0].data(), thrust::raw_pointer_cast(std::get<0>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[1].data(), thrust::raw_pointer_cast(std::get<1>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[2].data(), thrust::raw_pointer_cast(std::get<2>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[3].data(), thrust::raw_pointer_cast(std::get<3>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[4].data(), thrust::raw_pointer_cast(std::get<4>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[5].data(), thrust::raw_pointer_cast(std::get<5>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[6].data(), thrust::raw_pointer_cast(std::get<6>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);
        cudaMemcpy(hc[7].data(), thrust::raw_pointer_cast(std::get<7>(d_conn).data()), fullElemCount * sizeof(KeyType), cudaMemcpyDeviceToHost);

        std::vector<int> remap(nodeCount, -1);
        std::vector<int> compact;
        compact.reserve(numOwnedElem * 8);
        for (size_t e = elemStart; e < elemEnd; ++e)
            for (int k = 0; k < 8; ++k) {
                int nid = static_cast<int>(hc[k][e]);
                if (remap[nid] < 0) { remap[nid] = static_cast<int>(compact.size()); compact.push_back(nid); }
            }
        const size_t numWrittenNodes = compact.size();

        const std::string fname = rankPieceName(prefix_, step, rank);
        std::ofstream f(fname);
        if (!f) throw std::runtime_error("Cannot open " + fname);

        f << "<?xml version=\"1.0\"?>\n";
        f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        f << "<UnstructuredGrid>\n";
        f << "<Piece NumberOfPoints=\"" << numWrittenNodes
          << "\" NumberOfCells=\"" << numOwnedElem << "\">\n";

        // Points
        f << "<Points>\n<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (size_t i = 0; i < numWrittenNodes; ++i) {
            int n = compact[i];
            f << hx[n] << " " << hy[n] << " " << hz[n] << "\n";
        }
        f << "</DataArray>\n</Points>\n";

        // Cells
        f << "<Cells>\n";
        f << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (size_t e = elemStart; e < elemEnd; ++e) {
            for (int k = 0; k < 8; ++k) f << remap[static_cast<int>(hc[k][e])] << " ";
            f << "\n";
        }
        f << "</DataArray>\n";
        f << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        for (size_t i = 1; i <= numOwnedElem; ++i) f << (8 * i) << "\n";
        f << "</DataArray>\n";
        f << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (size_t i = 0; i < numOwnedElem; ++i) f << "12\n";
        f << "</DataArray>\n";
        f << "</Cells>\n";

        // PointData: scalars + vectors
        std::string scalarsAttr, vectorsAttr;
        for (const auto& hf : hfields) {
            using K = typename FieldDesc::Kind;
            if (hf.desc->kind == K::PointScalar)
                scalarsAttr = scalarsAttr.empty() ? hf.desc->name : scalarsAttr;
            else if (hf.desc->kind == K::PointVector3)
                vectorsAttr = vectorsAttr.empty() ? hf.desc->name : vectorsAttr;
        }
        f << "<PointData";
        if (!scalarsAttr.empty()) f << " Scalars=\"" << scalarsAttr << "\"";
        if (!vectorsAttr.empty()) f << " Vectors=\"" << vectorsAttr << "\"";
        f << ">\n";
        for (const auto& hf : hfields) {
            using K = typename FieldDesc::Kind;
            if (hf.desc->kind == K::PointScalar) {
                f << "<DataArray type=\"Float64\" Name=\"" << hf.desc->name << "\" format=\"ascii\">\n";
                for (size_t i = 0; i < numWrittenNodes; ++i) f << hf.x[compact[i]] << "\n";
                f << "</DataArray>\n";
            } else if (hf.desc->kind == K::PointVector3) {
                f << "<DataArray type=\"Float32\" Name=\"" << hf.desc->name
                  << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
                for (size_t i = 0; i < numWrittenNodes; ++i) {
                    int n = compact[i];
                    f << hf.x[n] << " " << hf.y[n] << " " << hf.z[n] << "\n";
                }
                f << "</DataArray>\n";
            }
        }
        f << "</PointData>\n";

        // CellData: scalars (per-owned-element)
        std::string cellScalarsAttr;
        for (const auto& hf : hfields) {
            using K = typename FieldDesc::Kind;
            if (hf.desc->kind == K::CellScalar) {
                cellScalarsAttr = cellScalarsAttr.empty() ? hf.desc->name : cellScalarsAttr;
                break;
            }
        }
        f << "<CellData";
        if (!cellScalarsAttr.empty()) f << " Scalars=\"" << cellScalarsAttr << "\"";
        f << ">\n";
        for (const auto& hf : hfields) {
            using K = typename FieldDesc::Kind;
            if (hf.desc->kind == K::CellScalar) {
                f << "<DataArray type=\"Float64\" Name=\"" << hf.desc->name << "\" format=\"ascii\">\n";
                for (size_t e = elemStart; e < elemEnd; ++e) f << hf.x[e] << "\n";
                f << "</DataArray>\n";
            }
        }
        f << "</CellData>\n";

        f << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    }

    // .pvtu master for multi-field. Lists every field array.
    void writeMasterMulti(int step, int numRanks, const std::vector<FieldDesc>& fields) const
    {
        const std::string fname = masterName(prefix_, step);
        std::ofstream f(fname);
        if (!f) throw std::runtime_error("Cannot open " + fname);

        f << "<?xml version=\"1.0\"?>\n";
        f << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        f << "<PUnstructuredGrid GhostLevel=\"0\">\n";
        f << "<PPoints>\n<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n</PPoints>\n";

        std::string scalarsAttr, vectorsAttr, cellScalarsAttr;
        for (const auto& fd : fields) {
            using K = typename FieldDesc::Kind;
            if (fd.kind == K::PointScalar  && scalarsAttr.empty())     scalarsAttr     = fd.name;
            if (fd.kind == K::PointVector3 && vectorsAttr.empty())     vectorsAttr     = fd.name;
            if (fd.kind == K::CellScalar   && cellScalarsAttr.empty()) cellScalarsAttr = fd.name;
        }
        f << "<PPointData";
        if (!scalarsAttr.empty()) f << " Scalars=\"" << scalarsAttr << "\"";
        if (!vectorsAttr.empty()) f << " Vectors=\"" << vectorsAttr << "\"";
        f << ">\n";
        for (const auto& fd : fields) {
            using K = typename FieldDesc::Kind;
            if (fd.kind == K::PointScalar)
                f << "<PDataArray type=\"Float32\" Name=\"" << fd.name << "\"/>\n";
            else if (fd.kind == K::PointVector3)
                f << "<PDataArray type=\"Float32\" Name=\"" << fd.name
                  << "\" NumberOfComponents=\"3\"/>\n";
        }
        f << "</PPointData>\n";

        f << "<PCellData";
        if (!cellScalarsAttr.empty()) f << " Scalars=\"" << cellScalarsAttr << "\"";
        f << ">\n";
        for (const auto& fd : fields) {
            using K = typename FieldDesc::Kind;
            if (fd.kind == K::CellScalar)
                f << "<PDataArray type=\"Float32\" Name=\"" << fd.name << "\"/>\n";
        }
        f << "</PCellData>\n";

        for (int r = 0; r < numRanks; ++r) {
            std::string piece = rankPieceName(prefix_, step, r);
            auto slash = piece.find_last_of("/\\");
            std::string basename = (slash == std::string::npos) ? piece : piece.substr(slash + 1);
            f << "<Piece Source=\"" << basename << "\"/>\n";
        }
        f << "</PUnstructuredGrid>\n</VTKFile>\n";
    }

    // Rank 0 writes the .pvtu master that references each rank's .vtu piece.
    // Path in the Source attribute is the basename only; ParaView resolves
    // relative to the .pvtu location.
    void writeMaster(int step, int numRanks, const std::string& field) const
    {
        const std::string fname = masterName(prefix_, step);
        std::ofstream f(fname);
        if (!f) throw std::runtime_error("Cannot open " + fname);

        f << "<?xml version=\"1.0\"?>\n";
        f << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        f << "<PUnstructuredGrid GhostLevel=\"0\">\n";

        f << "<PPoints>\n";
        f << "<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
        f << "</PPoints>\n";

        f << "<PPointData Scalars=\"" << field << "\">\n";
        f << "<PDataArray type=\"Float32\" Name=\"" << field << "\"/>\n";
        f << "</PPointData>\n";

        for (int r = 0; r < numRanks; ++r)
        {
            std::string piece = rankPieceName(prefix_, step, r);
            // Strip directory if present so ParaView resolves relative to the pvtu
            auto slash = piece.find_last_of("/\\");
            std::string basename = (slash == std::string::npos) ? piece : piece.substr(slash + 1);
            f << "<Piece Source=\"" << basename << "\"/>\n";
        }

        f << "</PUnstructuredGrid>\n</VTKFile>\n";
    }

    // Append (or recreate) the .pvd time collection. We rewrite each call
    // because ParaView's pvd format is a single closed XML document; appending
    // requires rewriting. Step count is small for AMR (≤10) so cost is trivial.
    void appendToPvd(int step, double timeValue) const
    {
        // Read existing entries if present
        std::vector<std::pair<int, double>> entries;
        {
            std::ifstream in(pvdName(prefix_));
            if (in)
            {
                std::string line;
                while (std::getline(in, line))
                {
                    // `step="` is a SUBSTRING of `timestep="`, so we must
                    // require a leading space (or other token-break) for the
                    // standalone step= attribute or every line collapses to
                    // step=0 and the PVD points all timesteps at frame 0.
                    auto p = line.find("timestep=\"");
                    auto s = line.find(" step=\"");
                    if (p != std::string::npos && s != std::string::npos)
                    {
                        double t = std::stod(line.substr(p + 10));
                        int st   = std::stoi(line.substr(s + 7));
                        entries.emplace_back(st, t);
                    }
                }
            }
        }
        // Drop any prior entry at the same step (in case we're overwriting)
        entries.erase(
            std::remove_if(entries.begin(), entries.end(),
                           [step](auto& e) { return e.first == step; }),
            entries.end());
        entries.emplace_back(step, timeValue);
        std::sort(entries.begin(), entries.end(),
                  [](auto& a, auto& b) { return a.first < b.first; });

        std::ofstream f(pvdName(prefix_), std::ios::trunc);
        if (!f) throw std::runtime_error("Cannot open " + pvdName(prefix_));
        f << "<?xml version=\"1.0\"?>\n";
        f << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        f << "<Collection>\n";
        for (auto& [st, t] : entries)
        {
            std::string master = masterName(prefix_, st);
            auto slash = master.find_last_of("/\\");
            std::string basename = (slash == std::string::npos) ? master : master.substr(slash + 1);
            f << "<DataSet timestep=\"" << t << "\" step=\"" << st
              << "\" file=\"" << basename << "\"/>\n";
        }
        f << "</Collection>\n</VTKFile>\n";
    }
};

}  // namespace fem
}  // namespace mars
