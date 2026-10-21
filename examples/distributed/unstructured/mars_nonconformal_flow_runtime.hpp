#pragma once
#include "mars_nonconformal_mesh_input.hpp"
#include "backend/distributed/unstructured/fem/mars_nonconformal_flow.hpp"
#include "backend/distributed/unstructured/fem/segregated/mars_segregated_geometry.hpp"
#include "backend/distributed/unstructured/fem/segregated/mars_segregated_tet_interior.hpp"
#include "backend/distributed/unstructured/fem/segregated/mars_segregated_boundary.hpp"
#include "backend/distributed/unstructured/fem/segregated/mars_segregated_node.hpp"

namespace ncflow
{
inline void finish()
{
#if defined(__CUDACC__) || defined(__HIPCC__)
    bool ok = cudaGetLastError() == cudaSuccess;
    ok      = (cudaStreamSynchronize(0) == cudaSuccess) && ok;
    check(ok, "flow stream completion");
#endif
}

using namespace mars::segregated;
using PatchRow = nonconformal_search_detail::PatchRow<double>;
struct Node
{
    InterfaceFlowState<double> flow;
    double increment{}, volume{}, divergence{};
};
struct Face
{
    int element, ordinal, kind;
};
struct Cell
{
    int node[4];
    double xyz[12];
    TetGeometry<double> geometry;
};
struct View
{
    const Cell* cells;
    const Face* faces;
    const PatchRow* patches;
    const int *cellOffsets, *cellIds, *faceOffsets, *faceIds, *patchOffsets, *patchIds;
    const int *slot, *boundary;
    const int64_t *rows, *columns;
    const int* offsets;
    const Vec3<double>* coordinates;
    Node* state;
    double *interiorMass, *boundaryMass, *interfaceMass;
    double density, viscosity, alpha, dt, speed;
    bool plug;
    const int* rowLocal;
    const int64_t *globalDofs, *stateIds;
    int ownedCount, stateCount;
    MARS_FS_FN int findSlot(int64_t id) const
    {
        int pos = nonconformal_detail::columnSlot(stateIds, 0, ownedCount, id);
        return pos >= 0 ? pos : nonconformal_detail::columnSlot(stateIds, ownedCount, stateCount, id);
    }
    MARS_FS_FN double prescribed(int node, int c) const
    {
        if (c != 0) return 0;
        if (plug) return speed;
        if (boundary[node] == 2) return 0;
        const auto x = coordinates[node];
        return 16 * speed * x.y * (1 - x.y) * x.z * (1 - x.z);
    }
    MARS_FS_FN TetInteriorInput interior(int e, int stage) const
    {
        TetInteriorInput in{};
        in.stage         = stage;
        const auto& cell = cells[e];
        for (int k = 0; k < 4; ++k)
        {
            const auto& s   = state[slot[cell.node[k]]].flow;
            in.density[k]   = density;
            in.viscosity[k] = viscosity;
            in.pressure[k]  = s.pressure;
            for (int j = 0; j < 3; ++j)
            {
                in.coordinates[3 * k + j]       = cell.xyz[3 * k + j];
                in.velocity[3 * k + j]          = nonconformal_flow_detail::component(s.velocity, j);
                in.pressure_gradient[3 * k + j] = nonconformal_flow_detail::component(s.pressureGradient, j);
                in.influence_lhs[3 * k + j]     = in.influence_rhs[3 * k + j] =
                    nonconformal_flow_detail::component(s.influence, j);
            }
        }
        for (int q = 0; q < 6; ++q)
        {
            in.edges[2 * q]     = tet_edge_node(q, 0);
            in.edges[2 * q + 1] = tet_edge_node(q, 1);
            tet_sample_shape(q, false, in.velocity_shape + 4 * q);
            tet_sample_shape(q, false, in.coordinate_shape + 4 * q);
            in.stored_flux[q] = interiorMass[6 * e + q];
            for (int j = 0; j < 3; ++j)
                in.area[3 * q + j] = cell.geometry.area[3 * q + j];
            for (int k = 0; k < 12; ++k)
                in.shape_gradient[12 * q + k] = cell.geometry.gradient[k];
        }
        return in;
    }
    MARS_FS_FN BoundaryInput boundaryInput(int f, bool momentum) const
    {
        const auto face  = faces[f];
        const auto& cell = cells[face.element];
        BoundaryInput in{};
        in.stage = momentum ? 4 : (face.kind == 1 ? 1 : 0);
        for (int k = 0; k < 4; ++k)
        {
            const auto& s       = state[slot[cell.node[k]]].flow;
            in.pressure[k]      = s.pressure;
            in.bc_multiplier[k] = 1;
            for (int j = 0; j < 3; ++j)
            {
                in.velocity[3 * k + j]          = nonconformal_flow_detail::component(s.velocity, j);
                in.pressure_gradient[3 * k + j] = nonconformal_flow_detail::component(s.pressureGradient, j);
            }
        }
        double area[3];
        tet_boundary_area(cell.geometry, face.ordinal, area);
        for (int q = 0; q < 3; ++q)
        {
            int k = tet_face_node(face.ordinal, q), node = cell.node[k];
            in.face_nodes[q] = in.nearest[q] = k;
            in.opposing[q]                   = tet_opposite_node(face.ordinal);
            in.viscosity[q]                  = viscosity;
            in.density[q]                    = density;
            in.stored_flux[q]                = boundaryMass[3 * f + q];
            tri_sample_shape(q, false, in.shape + 3 * q);
            if (face.kind == 1)
            {
                in.pressure[k]      = 0;
                in.bc_multiplier[k] = 0;
            }
            for (int j = 0; j < 3; ++j)
            {
                in.area[3 * q + j]              = area[j];
                in.boundary_velocity[3 * q + j] = prescribed(node, j);
                in.influence_lhs[3 * q + j]     = in.influence_rhs[3 * q + j] =
                    nonconformal_flow_detail::component(state[slot[node]].flow.influence, j);
            }
            for (int j = 0; j < 12; ++j)
                in.gradient[12 * q + j] = cell.geometry.gradient[j];
        }
        return in;
    }
    MARS_FS_FN InterfaceFlowFlux<double> interface(int q) const
    {
        InterfaceFlowState<double> s[8];
        for (int k = 0; k < 8; ++k)
            s[k] = state[findSlot(patches[q].patch.column[k])].flow;
        return nonconformalFlowFlux(patches[q].patch, s, density, viscosity, viscosity, interfaceMass[q]);
    }
};
struct Gradient
{
    View v;
    bool correction;
    Vec3<double>* output;
    MARS_FS_FN void operator()(int r) const
    {
        Vec3<double> result{};
        double volume = 0;
        for (int t = v.cellOffsets[r]; t < v.cellOffsets[r + 1]; ++t)
        {
            const auto& cell = v.cells[v.cellIds[t]];
            double weight    = cell.geometry.volume * .25;
            volume += weight;
            for (int k = 0; k < 4; ++k)
            {
                const auto& s = v.state[v.slot[cell.node[k]]];
                double value  = correction ? s.increment : s.flow.pressure;
                result.x += weight * cell.geometry.gradient[3 * k] * value;
                result.y += weight * cell.geometry.gradient[3 * k + 1] * value;
                result.z += weight * cell.geometry.gradient[3 * k + 2] * value;
            }
        }
        for (int t = v.patchOffsets[r]; t < v.patchOffsets[r + 1]; ++t)
        {
            const auto& entry = v.patches[v.patchIds[t]];
            const auto& p     = entry.patch;
            double jump       = 0;
            for (int k = 0; k < 8; ++k)
            {
                const auto& s = v.state[v.findSlot(p.column[k])];
                double value  = correction ? s.increment : s.flow.pressure;
                jump += (k < 4 ? -1. : 1.) * p.shape[k] * value;
            }
            // Common trace minus the one-sided trace, with the side's outward normal.
            result = nonconformal_detail::add(result, nonconformal_detail::scale(p.normal, .5 * p.area * jump));
        }
        // The outlet has prescribed physical pressure zero, including its increment.
        // Replace the one-sided boundary trace in the Gauss gradient by that value.
        for (int t = v.faceOffsets[r]; t < v.faceOffsets[r + 1]; ++t)
        {
            int f           = v.faceIds[t];
            const auto face = v.faces[f];
            if (face.kind != 1) continue;
            const auto& cell = v.cells[face.element];
            int q            = 0;
            while (cell.node[tet_face_node(face.ordinal, q)] != v.rowLocal[r])
                ++q;
            double shape[3], area[3], value = 0;
            tri_sample_shape(q, false, shape);
            tet_boundary_area(cell.geometry, face.ordinal, area);
            for (int k = 0; k < 3; ++k)
            {
                const auto& s = v.state[v.slot[cell.node[tet_face_node(face.ordinal, k)]]];
                value += shape[k] * (correction ? s.increment : s.flow.pressure);
            }
            result.x -= area[0] * value;
            result.y -= area[1] * value;
            result.z -= area[2] * value;
        }
        output[r]         = nonconformal_detail::scale(result, 1. / volume);
        v.state[r].volume = volume;
    }
};
struct StoreGradient
{
    Node* state;
    const Vec3<double>* gradient;
    MARS_FS_FN void operator()(int r) const { state[r].flow.pressureGradient = gradient[r]; }
};
struct Assemble
{
    View v;
    bool momentum;
    double *values, *rhs;
    int* valid;
    Vec3<double>* influence;
    MARS_FS_FN void add(int r, int64_t column, int i, int j, double value) const
    {
        int pos = nonconformal_detail::columnSlot(v.columns, v.offsets[r], v.offsets[r + 1], column);
        if (pos < 0)
        {
            valid[r] = 0;
            return;
        }
        if (momentum)
        {
            int len = v.offsets[r + 1] - v.offsets[r];
            values[9 * v.offsets[r] + 3 * len * i + 3 * (pos - v.offsets[r]) + j] += value;
        }
        else
            values[pos] += value;
    }
    MARS_FS_FN void operator()(int r) const
    {
        int node = v.rowLocal[r], width = momentum ? 3 : 1;
        valid[r] = 1;
        for (int j = width * width * v.offsets[r]; j < width * width * v.offsets[r + 1]; ++j)
            values[j] = 0;
        for (int i = 0; i < width; ++i)
            rhs[width * r + i] = 0;
        if (momentum && v.boundary[node])
        {
            for (int i = 0; i < 3; ++i)
            {
                add(r, v.rows[r], i, i, 1);
                rhs[3 * r + i] =
                    v.prescribed(node, i) - nonconformal_flow_detail::component(v.state[r].flow.velocity, i);
            }
            influence[r] = {v.dt * v.alpha / v.density, v.dt * v.alpha / v.density, v.dt * v.alpha / v.density};
            return;
        }
        for (int t = v.cellOffsets[r]; t < v.cellOffsets[r + 1]; ++t)
        {
            int e            = v.cellIds[t];
            const auto& cell = v.cells[e];
            int local        = 0;
            while (cell.node[local] != node)
                ++local;
            auto in = v.interior(e, momentum ? 1 : 0);
            TetInteriorOutput out;
            tet_interior(in, out);
            for (int i = 0; i < width; ++i)
            {
                rhs[width * r + i] += out.rhs[width * local + i];
                for (int k = 0; k < 4; ++k)
                    for (int j = 0; j < width; ++j)
                        add(r, v.globalDofs[cell.node[k]], i, j,
                            out.lhs[(width * local + i) * 4 * width + width * k + j]);
            }
        }
        for (int t = v.faceOffsets[r]; t < v.faceOffsets[r + 1]; ++t)
        {
            int f = v.faceIds[t];
            if (momentum && v.faces[f].kind != 1) continue;
            const auto& cell = v.cells[v.faces[f].element];
            int local        = 0;
            while (cell.node[local] != node)
                ++local;
            auto in = v.boundaryInput(f, momentum);
            BoundaryOutput out;
            boundary_block(in, out);
            for (int i = 0; i < width; ++i)
            {
                rhs[width * r + i] += out.rhs[width * local + i];
                for (int k = 0; k < 4; ++k)
                    for (int j = 0; j < width; ++j)
                        add(r, v.globalDofs[cell.node[k]], i, j,
                            out.lhs[(width * local + i) * 4 * width + width * k + j]);
            }
        }
        for (int t = v.patchOffsets[r]; t < v.patchOffsets[r + 1]; ++t)
        {
            int q         = v.patchIds[t];
            const auto& p = v.patches[q];
            auto out      = v.interface(q);
            if (!out.valid)
            {
                valid[r] = 0;
                continue;
            }
            double sign = p.side == 0 ? 1. : -1.;
            for (int i = 0; i < width; ++i)
            {
                rhs[width * r + i] -= sign * (momentum ? out.momentum[i] : out.mass);
                for (int k = 0; k < 8; ++k)
                    for (int j = 0; j < width; ++j)
                        add(r, p.patch.column[k], i, j,
                            sign * (momentum ? out.momentumCoefficient[k][i][j] : out.pressureCoefficient[k]));
            }
        }
        if (momentum)
        {
            auto& s          = v.state[r];
            double transient = v.density * s.volume / v.dt;
            for (int i = 0; i < 3; ++i)
            {
                rhs[3 * r + i] -= s.volume * nonconformal_flow_detail::component(s.flow.pressureGradient, i);
                add(r, v.rows[r], i, i, transient);
                int pos      = nonconformal_detail::columnSlot(v.columns, v.offsets[r], v.offsets[r + 1], v.rows[r]);
                int len      = v.offsets[r + 1] - v.offsets[r];
                int diagonal = 9 * v.offsets[r] + 3 * len * i + 3 * (pos - v.offsets[r]) + i;
                values[diagonal] /= v.alpha;
                double d = s.volume / values[diagonal];
                if (!(d > 0) || !std::isfinite(d)) valid[r] = 0;
                if (i == 0) influence[r].x = d;
                if (i == 1) influence[r].y = d;
                if (i == 2) influence[r].z = d;
            }
        }
    }
};
struct Jacobi
{
    InterfaceCsrView<double> matrix;
    double* diagonal;
    int* valid;
    MARS_FS_FN void operator()(int r) const
    {
        int pos =
            nonconformal_detail::columnSlot(matrix.columns, matrix.offsets[r], matrix.offsets[r + 1], matrix.rowIds[r]);
        valid[r]    = pos >= 0 && std::isfinite(matrix.values[pos]) && matrix.values[pos] != 0;
        diagonal[r] = valid[r] ? 1. / matrix.values[pos] : 0;
    }
};
struct StoreInfluence
{
    Node* state;
    const Vec3<double>* influence;
    MARS_FS_FN void operator()(int r) const { state[r].flow.influence = influence[r]; }
};
struct AddVelocity
{
    Node* state;
    const double* increment;
    MARS_FS_FN void operator()(int r) const
    {
        state[r].flow.velocity.x += increment[3 * r];
        state[r].flow.velocity.y += increment[3 * r + 1];
        state[r].flow.velocity.z += increment[3 * r + 2];
    }
};
struct SetIncrement
{
    Node* state;
    const double* increment;
    MARS_FS_FN void operator()(int r) const { state[r].increment = increment[r]; }
};
struct Correct
{
    View v;
    const Vec3<double>* gradient;
    double alphaPressure;
    MARS_FS_FN void operator()(int r) const
    {
        auto& s = v.state[r];
        s.flow.pressure += alphaPressure * s.increment;
        if (!v.boundary[v.rowLocal[r]])
        {
            s.flow.velocity.x -= s.flow.influence.x * gradient[r].x;
            s.flow.velocity.y -= s.flow.influence.y * gradient[r].y;
            s.flow.velocity.z -= s.flow.influence.z * gradient[r].z;
        }
    }
};
// Flux histories are evaluated before velocity correction, using the complete
// pressure increment. Pressure underrelaxation must not weaken mass correction.
struct RefreshMass
{
    View v;
    MARS_FS_FN void operator()(int r) const
    {
        double divergence = 0;
        int node          = v.rowLocal[r];
        for (int t = v.cellOffsets[r]; t < v.cellOffsets[r + 1]; ++t)
        {
            int e   = v.cellIds[t];
            auto in = v.interior(e, 0);
            TetInteriorOutput out;
            tet_interior(in, out);
            const auto& cell = v.cells[e];
            for (int q = 0; q < 6; ++q)
            {
                int a = tet_edge_node(q, 0), b = tet_edge_node(q, 1);
                if (cell.node[a] != node && cell.node[b] != node) continue;
                double mass = out.flux[q];
                for (int k = 0; k < 4; ++k)
                    for (int j = 0; j < 3; ++j)
                    {
                        double d = 0;
                        for (int l = 0; l < 4; ++l)
                            d += in.velocity_shape[4 * q + l] * in.influence_lhs[3 * l + j];
                        mass -= v.density * d * cell.geometry.gradient[3 * k + j] * cell.geometry.area[3 * q + j] *
                                v.state[v.slot[cell.node[k]]].increment;
                    }
                divergence += (cell.node[a] == node ? 1. : -1.) * mass;
            }
        }
        for (int t = v.faceOffsets[r]; t < v.faceOffsets[r + 1]; ++t)
        {
            int f   = v.faceIds[t];
            auto in = v.boundaryInput(f, false);
            BoundaryOutput out;
            boundary_block(in, out);
            const auto& cell = v.cells[v.faces[f].element];
            int local        = 0;
            while (cell.node[local] != node)
                ++local;
            for (int q = 0; q < 3; ++q)
                if (tet_face_node(v.faces[f].ordinal, q) == local)
                {
                    double mass = out.flux[q];
                    for (int k = 0; k < 4; ++k)
                        mass += out.lhs[4 * local + k] * v.state[v.slot[cell.node[k]]].increment;
                    divergence += mass;
                }
        }
        for (int t = v.patchOffsets[r]; t < v.patchOffsets[r + 1]; ++t)
        {
            int q         = v.patchIds[t];
            auto out      = v.interface(q);
            const auto& p = v.patches[q];
            double mass   = out.mass;
            for (int k = 0; k < 8; ++k)
                mass += out.pressureCoefficient[k] * v.state[v.findSlot(p.patch.column[k])].increment;
            divergence += (p.side == 0 ? 1. : -1.) * mass;
        }
        v.state[r].divergence = divergence;
    }
};
} // namespace ncflow
namespace ncflow
{
struct CellMass
{
    View v;
    const int* active;
    MARS_FS_FN void operator()(int index) const
    {
        int e   = active[index];
        auto in = v.interior(e, 0);
        TetInteriorOutput out;
        tet_interior(in, out);
        const auto& cell = v.cells[e];
        for (int q = 0; q < 6; ++q)
        {
            double mass = out.flux[q];
            for (int k = 0; k < 4; ++k)
                for (int j = 0; j < 3; ++j)
                {
                    double d = 0;
                    for (int l = 0; l < 4; ++l)
                        d += in.velocity_shape[4 * q + l] * in.influence_lhs[3 * l + j];
                    mass -= v.density * d * cell.geometry.gradient[3 * k + j] * cell.geometry.area[3 * q + j] *
                            v.state[v.slot[cell.node[k]]].increment;
                }
            v.interiorMass[6 * e + q] = mass;
        }
    }
};
struct BoundaryMass
{
    View v;
    const int* active;
    MARS_FS_FN void operator()(int index) const
    {
        int f   = active[index];
        auto in = v.boundaryInput(f, false);
        BoundaryOutput out;
        boundary_block(in, out);
        const auto& cell = v.cells[v.faces[f].element];
        for (int q = 0; q < 3; ++q)
        {
            int local   = tet_face_node(v.faces[f].ordinal, q);
            double mass = out.flux[q];
            for (int k = 0; k < 4; ++k)
                mass += out.lhs[4 * local + k] * v.state[v.slot[cell.node[k]]].increment;
            v.boundaryMass[3 * f + q] = mass;
        }
    }
};
struct PatchMass
{
    View v;
    MARS_FS_FN void operator()(int q) const
    {
        auto out      = v.interface(q);
        const auto& p = v.patches[q];
        double mass   = out.mass;
        for (int k = 0; k < 8; ++k)
            mass += out.pressureCoefficient[k] * v.state[v.findSlot(p.patch.column[k])].increment;
        v.interfaceMass[q] = mass;
    }
};
struct CellGeometry
{
    Cell* cells;
    int* valid;
    MARS_FS_FN void operator()(int e) const
    {
        auto& cell = cells[e];
        InterfaceTetFace<double> face{};
        face.opposite    = 3;
        face.diffusivity = 1;
        for (int k = 0; k < 4; ++k)
        {
            face.x[k]   = {cell.xyz[3 * k], cell.xyz[3 * k + 1], cell.xyz[3 * k + 2]};
            face.dof[k] = cell.node[k];
        }
        nonconformal_detail::Geometry<double> g;
        valid[e] = nonconformal_detail::geometry(face, g) == InterfaceStatus::valid;
        if (!valid[e]) return;
        cell.geometry.volume =
            std::abs(dot(sub(face.x[1], face.x[0]),
                         nonconformal_detail::cross(sub(face.x[2], face.x[0]), sub(face.x[3], face.x[0])))) /
            6;
        for (int k = 0; k < 4; ++k)
            for (int j = 0; j < 3; ++j)
                cell.geometry.gradient[3 * k + j] = nonconformal_flow_detail::component(g.gradient[k], j);
        for (int q = 0; q < 6; ++q)
            for (int j = 0; j < 3; ++j)
                cell.geometry.area[3 * q + j] = cell.geometry.volume * .25 *
                                                (cell.geometry.gradient[3 * tet_edge_node(q, 1) + j] -
                                                 cell.geometry.gradient[3 * tet_edge_node(q, 0) + j]);
    }
};
struct PatchIdentity
{
    const PatchRow* patches;
    int64_t* row;
    MARS_FS_FN void operator()(int q) const { row[q] = patches[q].patch.row[patches[q].side]; }
};
inline void validFlags(const Buffer<int>& flags, const char* message)
{
#if defined(__CUDACC__) || defined(__HIPCC__)
    check(thrust::count(flags.begin(), flags.end(), 0) == 0, message);
#else
    check(std::find(flags.begin(), flags.end(), 0) == flags.end(), message);
#endif
}
inline void incidence(const std::vector<std::vector<int>>& lists, Buffer<int>& offsets, Buffer<int>& ids)
{
    std::vector<int> off{0}, all;
    for (const auto& list : lists)
    {
        all.insert(all.end(), list.begin(), list.end());
        off.push_back(int(all.size()));
    }
    offsets = Buffer<int>(off.begin(), off.end());
    ids     = Buffer<int>(all.begin(), all.end());
}
struct BlockGraphRows
{
    const int64_t *rows, *columns;
    const int* offsets;
    int64_t *blockRows, *blockColumns;
    int* blockOffsets;
    int count;
    MARS_FS_FN void operator()(int r) const
    {
        if (r == count)
        {
            blockOffsets[3 * r] = 9 * offsets[r];
            return;
        }
        int length = offsets[r + 1] - offsets[r];
        for (int a = 0; a < 3; ++a)
        {
            blockRows[3 * r + a]    = 3 * rows[r] + a;
            int start               = 9 * offsets[r] + 3 * a * length;
            blockOffsets[3 * r + a] = start;
            for (int j = 0; j < length; ++j)
                for (int b = 0; b < 3; ++b)
                    blockColumns[start + 3 * j + b] = 3 * columns[offsets[r] + j] + b;
        }
    }
};
inline Graph blockGraph(const Graph& node)
{
    check(node.columns.size() <= size_t(INT_MAX / 9) && node.rowIds.size() < size_t(INT_MAX / 3),
          "block graph exceeds local index capacity");
    Graph result;
    result.rowIds.resize(3 * node.rowIds.size());
    result.columns.resize(9 * node.columns.size());
    result.offsets.resize(3 * node.rowIds.size() + 1);
    result.values.resize(result.columns.size());
    each<device>(int(node.rowIds.size()) + 1,
                 BlockGraphRows{pointer(node.rowIds), pointer(node.columns), pointer(node.offsets),
                                pointer(result.rowIds), pointer(result.columns), pointer(result.offsets),
                                int(node.rowIds.size())});
    finish();
    return result;
}

struct FlowReport
{
    double sums[11]{}, maxima[4]{};
};
struct ReportSum
{
    MARS_FS_FN FlowReport operator()(FlowReport a, const FlowReport& b) const
    {
        for (int i = 0; i < 11; ++i)
            a.sums[i] += b.sums[i];
        for (int i = 0; i < 4; ++i)
            if (b.maxima[i] > a.maxima[i]) a.maxima[i] = b.maxima[i];
        return a;
    }
};
struct Measure
{
    View v;
    MARS_FS_FN FlowReport operator()(int r) const
    {
        FlowReport out;
        int node      = v.rowLocal[r];
        const auto& s = v.state[r];
        out.sums[8]   = std::abs(s.divergence);
        out.maxima[0] = std::abs(s.divergence);
        out.maxima[1] = std::max(std::abs(s.flow.velocity.x - v.speed),
                                 std::max(std::abs(s.flow.velocity.y), std::abs(s.flow.velocity.z)));
        out.maxima[2] = std::abs(s.flow.pressure);
        bool finite   = nonconformal_flow_detail::finite(s.flow.velocity) &&
                      nonconformal_flow_detail::finite(s.flow.pressureGradient) && std::isfinite(s.flow.pressure) &&
                      std::isfinite(s.divergence);
        for (int t = v.faceOffsets[r]; t < v.faceOffsets[r + 1]; ++t)
        {
            int f            = v.faceIds[t];
            const auto& face = v.faces[f];
            const auto& cell = v.cells[face.element];
            for (int q = 0; q < 3; ++q)
                if (cell.node[tet_face_node(face.ordinal, q)] == node) out.sums[face.kind] += v.boundaryMass[3 * f + q];
        }
        for (int t = v.patchOffsets[r]; t < v.patchOffsets[r + 1]; ++t)
        {
            int q             = v.patchIds[t];
            const auto& entry = v.patches[q];
            double sign       = entry.side == 0 ? 1. : -1.;
            out.sums[3 + entry.side] += sign * v.interfaceMass[q];
            auto flux = v.interface(q);
            finite    = finite && flux.valid;
            for (int j = 0; j < 3; ++j)
                out.sums[5 + j] += sign * (flux.momentum[j] + flux.pressureForce[j]);
            Vec3<double> jump{};
            for (int k = 0; k < 8; ++k)
            {
                auto u = v.state[v.findSlot(entry.patch.column[k])].flow.velocity;
                jump   = nonconformal_detail::add(
                    jump, nonconformal_detail::scale(u, (k < 4 ? 1. : -1.) * entry.patch.shape[k]));
            }
            out.maxima[3] = std::max(out.maxima[3], std::sqrt(dot(jump, jump)));
        }
        double rawDivergence = 0;
        for (int t = v.cellOffsets[r]; t < v.cellOffsets[r + 1]; ++t)
        {
            const auto& cell = v.cells[v.cellIds[t]];
            for (int k = 0; k < 4; ++k)
                for (int j = 0; j < 3; ++j)
                    rawDivergence +=
                        cell.geometry.volume * .25 * cell.geometry.gradient[3 * k + j] *
                        nonconformal_flow_detail::component(v.state[v.slot[cell.node[k]]].flow.velocity, j);
        }
        out.sums[10] = std::abs(v.density * rawDivergence);
        out.sums[9]  = finite ? 0 : 1;
        return out;
    }
};
class Channel
{
public:
    Channel(const Mesh& m, bool plug, double speed, double viscosity, double density = 1.)
        : plan(MPI_COMM_WORLD)
        , halo(MPI_COMM_WORLD)
    {
        auto own = ownedRows(m);
        n        = int(own.size());
        auto ex  = extract(m);
        check(plan.build(ex.faces, ex.ids, ex.owners, own, peers(m), 1e-11, 1e-10, true), "flow interface coverage");
        std::vector<std::set<int64_t>> neighbors(n);
        std::vector<std::vector<int>> cellsByRow(n), facesByRow(n), patchesByRow(n);
        std::vector<Cell> hc(m.elements.size());
        std::set<int> activeCells;
        int first = int(int64_t(global_node_count(m)) * rankId / rankCount), last = first + n;
        std::vector<int> ownLocal;
        for (auto id : own)
        {
            auto it = std::lower_bound(m.dofs.begin(), m.dofs.end(), id);
            if (it == m.dofs.end() || *it != id) throw std::runtime_error("missing owned input node");
            ownLocal.push_back(int(it - m.dofs.begin()));
        }
        rowLocal   = Buffer<int>(ownLocal.begin(), ownLocal.end());
        globalDofs = Buffer<int64_t>(m.dofs.begin(), m.dofs.end());
        for (int e = 0; e < int(hc.size()); ++e)
        {
            for (int k = 0; k < 4; ++k)
            {
                int id               = m.cells[4 * e + k];
                hc[e].node[k]        = id;
                hc[e].xyz[3 * k]     = m.x[id].x;
                hc[e].xyz[3 * k + 1] = m.x[id].y;
                hc[e].xyz[3 * k + 2] = m.x[id].z;
                if (m.dofs[id] >= first && m.dofs[id] < last)
                {
                    cellsByRow[m.dofs[id] - first].push_back(e);
                    activeCells.insert(e);
                    for (int j = 0; j < 4; ++j)
                        neighbors[m.dofs[id] - first].insert(m.dofs[m.cells[4 * e + j]]);
                }
            }
        }
        std::vector<int64_t> cols;
        std::vector<int> offsets{0};
        for (auto& row : neighbors)
        {
            cols.insert(cols.end(), row.begin(), row.end());
            offsets.push_back(int(cols.size()));
        }
        Graph base;
        base.rowIds  = Buffer<int64_t>(own.begin(), own.end());
        base.columns = Buffer<int64_t>(cols.begin(), cols.end());
        base.offsets = Buffer<int>(offsets.begin(), offsets.end());
        base.values.resize(cols.size());
#if defined(__CUDACC__) || defined(__HIPCC__)
        check(
            buildNonconformalScalarGraph(base.view(), int(plan.stencils().size()), pointer(plan.stencils()), pressure),
            "flow graph");
#else
        check(buildNonconformalScalarGraphHost(base.view(), int(plan.stencils().size()), pointer(plan.stencils()),
                                               pressure),
              "flow graph");
#endif
        momentum     = blockGraph(pressure);
        auto allCols = host<int64_t>(pressure.columns);
        std::set<int64_t> ghosts;
        for (auto c : allCols)
            if (c < first || c >= last) ghosts.insert(c);
        std::vector<uint64_t> ids(own.begin(), own.end());
        ids.insert(ids.end(), ghosts.begin(), ghosts.end());
        std::vector<int> owners, slots(m.x.size(), -1);
        for (int i = 0; i < int(ids.size()); ++i)
        {
            owners.push_back(owner(ids[i], global_node_count(m)));
            auto it = std::lower_bound(m.dofs.begin(), m.dofs.end(), ids[i]);
            if (it != m.dofs.end() && *it == ids[i]) slots[it - m.dofs.begin()] = i;
        }
        check(halo.buildFromHost(ids, owners, peers(m)), "flow state halo");
        state.resize(ids.size());
        stateIds = Buffer<int64_t>(ids.begin(), ids.end());
        slot     = Buffer<int>(slots.begin(), slots.end());
        Buffer<int64_t> patchRows(plan.patches().size());
        each<device>(int(patchRows.size()), PatchIdentity{pointer(plan.patches()), pointer(patchRows)});
        finish();
        auto hp = host<int64_t>(patchRows);
        for (int q = 0; q < int(hp.size()); ++q)
            patchesByRow.at(hp[q] - first).push_back(q);
        int localOuter = int(m.outer.size()), totalOuter = 0;
        MPI_Allreduce(&localOuter, &totalOuter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        check(totalOuter > 0, "flow requires input prepared with --flow-boundaries");
        std::vector<Face> hf;
        std::vector<int> codes(m.x.size(), 0);
        std::set<int> activeFaces;
        for (auto tag : m.outer)
        {
            int f = 0;
            while (tet_opposite_node(f) != tag.opposite)
                ++f;
            int e = tag.element, kind = tag.side, index = int(hf.size());
            hf.push_back({e, f, kind});
            for (int j = 0; j < 3; ++j)
            {
                int id = hc[e].node[tet_face_node(f, j)];
                if (kind != 1) codes[id] = std::max(codes[id], kind == 0 ? 1 : 2);
                if (m.dofs[id] >= first && m.dofs[id] < last)
                {
                    facesByRow[m.dofs[id] - first].push_back(index);
                    activeFaces.insert(index);
                }
            }
        }
        cells       = Buffer<Cell>(hc.begin(), hc.end());
        faces       = Buffer<Face>(hf.begin(), hf.end());
        coordinates = Buffer<Vec3<double>>(m.x.begin(), m.x.end());
        boundary    = Buffer<int>(codes.begin(), codes.end());
        activeCell  = Buffer<int>(activeCells.begin(), activeCells.end());
        activeFace  = Buffer<int>(activeFaces.begin(), activeFaces.end());
        incidence(cellsByRow, cellOffsets, cellIds);
        incidence(facesByRow, faceOffsets, faceIds);
        incidence(patchesByRow, patchOffsets, patchIds);
        Buffer<int> valid(hc.size());
        each<device>(int(hc.size()), CellGeometry{pointer(cells), pointer(valid)});
        finish();
        validFlags(valid, "flow volume geometry");
        interiorMass.resize(6 * hc.size());
        boundaryMass.resize(3 * hf.size());
        interfaceMass.resize(plan.patches().size());
        gradient.resize(n);
        influence.resize(n);
        view = {pointer(cells),
                pointer(faces),
                pointer(plan.patches()),
                pointer(cellOffsets),
                pointer(cellIds),
                pointer(faceOffsets),
                pointer(faceIds),
                pointer(patchOffsets),
                pointer(patchIds),
                pointer(slot),
                pointer(boundary),
                pointer(pressure.rowIds),
                pointer(pressure.columns),
                pointer(pressure.offsets),
                pointer(coordinates),
                pointer(state),
                pointer(interiorMass),
                pointer(boundaryMass),
                pointer(interfaceMass),
                density,
                viscosity,
                .7,
                1.,
                speed,
                plug,
                pointer(rowLocal),
                pointer(globalDofs),
                pointer(stateIds),
                n,
                int(state.size())};
    }
    void exchange()
    {
#if defined(__CUDACC__) || defined(__HIPCC__)
        check(halo.forward(pointer(state), state.size()), "flow state exchange");
#else
        check(halo.forwardHost(state), "flow state exchange");
#endif
    }
    double assemble(bool mu, ChannelOps& ops)
    {
        Buffer<int> valid(n);
        auto& graph = mu ? momentum : pressure;
        each<device>(n, Assemble{view, mu, pointer(graph.values), ops.rhs(), pointer(valid), pointer(influence)});
        finish();
        validFlags(valid, "flow coefficients");
        if (mu) each<device>(n, StoreInfluence{pointer(state), pointer(influence)});
        valid.resize(ops.n);
        each<device>(ops.n, Jacobi{graph.view(), pointer(ops.diagonal), pointer(valid)});
        finish();
        validFlags(valid, "flow diagonal");
        return ops.norm(ops.rhs());
    }
    bool solve(const Mesh& m, int budget, double tolerance)
    {
        ChannelOps u(momentum, m, 100, 3, false), p(pressure, m, 100, 1, false);
        double initialU = -1, initialP = -1;
        for (int iteration = 0; iteration < budget; ++iteration)
        {
            exchange();
            each<device>(n, Gradient{view, false, pointer(gradient)});
            each<device>(n, StoreGradient{pointer(state), pointer(gradient)});
            exchange();
            double momentumNorm = assemble(true, u);
            exchange();
            if (initialU < 0) initialU = std::max(momentumNorm, 1e-20);
            u.zero(u.solution());
            auto ur = outlet_fgmres(u, 100, 3000, 1e-10, std::numeric_limits<double>::epsilon());
            if (!ur.converged)
            {
                if (rankId == 0) std::printf("momentum linear failure %g\n", ur.relative_residual);
                return false;
            }
            each<device>(n, AddVelocity{pointer(state), u.solution()});
            exchange();
            double pressureNorm = assemble(false, p);
            if (initialP < 0) initialP = std::max(pressureNorm, 1e-20);
            p.zero(p.solution());
            auto pr = outlet_fgmres(p, 100, 3000, 1e-10, std::numeric_limits<double>::epsilon());
            if (!pr.converged)
            {
                if (rankId == 0) std::printf("pressure linear failure %g\n", pr.relative_residual);
                return false;
            }
            each<device>(n, SetIncrement{pointer(state), p.solution()});
            exchange();
            each<device>(n, Gradient{view, true, pointer(gradient)});
            each<device>(n, RefreshMass{view});
            each<device>(int(activeCell.size()), CellMass{view, pointer(activeCell)});
            each<device>(int(activeFace.size()), BoundaryMass{view, pointer(activeFace)});
            each<device>(int(plan.patches().size()), PatchMass{view});
            each<device>(n, Correct{view, pointer(gradient), .3});
            finish();
            if (rankId == 0 && (iteration % 10 == 0 || iteration + 1 == budget))
            {
                std::printf("iteration %d momentum %.4e continuity %.4e linear %d/%d\n", iteration + 1,
                            momentumNorm / initialU, pressureNorm / initialP, ur.iterations, pr.iterations);
                std::fflush(stdout);
            }
            if (momentumNorm <= std::max(tolerance * initialU, 1e-12 * view.viscosity * view.speed) &&
                pressureNorm <= std::max(tolerance * initialP, 1e-12 * view.density * view.speed))
            {
                if (rankId == 0)
                    std::printf("converged iteration %d momentum %.6e continuity %.6e\n", iteration + 1,
                                momentumNorm / initialU, pressureNorm / initialP);
                return verify();
            }
        }
        return false;
    }
    bool verify()
    {
        exchange();
        FlowReport local;
#if defined(__CUDACC__) || defined(__HIPCC__)
        local = thrust::transform_reduce(thrust::device, thrust::counting_iterator<int>(0),
                                         thrust::counting_iterator<int>(n), Measure{view}, FlowReport{}, ReportSum{});
#else
        for (int i = 0; i < n; ++i)
            local = ReportSum{}(local, Measure{view}(i));
#endif
        FlowReport global;
        MPI_Allreduce(local.sums, global.sums, 11, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(local.maxima, global.maxima, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        double scale            = std::max(std::abs(global.sums[0]), 1e-30);
        double boundaryBalance  = std::abs(global.sums[0] + global.sums[1] + global.sums[2]) / scale;
        double interfaceBalance = std::abs(global.sums[3] + global.sums[4]) / scale;
        double momentumBalance  = 0;
        for (int j = 0; j < 3; ++j)
            momentumBalance += std::abs(global.sums[5 + j]);
        bool ok = global.sums[9] == 0 && global.sums[0] < 0 && global.sums[1] > 0 && boundaryBalance < 1e-8 &&
                  interfaceBalance < 1e-10 && global.sums[8] / scale < 1e-8 &&
                  momentumBalance / (scale * view.speed) < 1e-8;
        if (view.plug)
            ok = ok && global.maxima[1] / view.speed < 1e-6 &&
                 global.maxima[2] / (view.density * view.speed * view.speed) < 1e-5;
        if (rankId == 0)
            std::printf("%s mass inlet=%.10e outlet=%.10e boundary_balance=%.3e interface_balance=%.3e "
                        "continuity_l1=%.3e momentum_balance=%.3e velocity_jump=%.3e pressure_max=%.3e\n",
                        ok ? "PASS" : "FAIL", global.sums[0], global.sums[1], boundaryBalance, interfaceBalance,
                        global.sums[8] / scale, momentumBalance, global.maxima[3], global.maxima[2]);
        if (rankId == 0)
            std::printf("raw_velocity_divergence_l1/inlet_mass=%.6e (separate from corrected mass flux)\n",
                        global.sums[10] / scale);
        if (rankId == 0 && view.plug)
            std::printf("analytic plug velocity_error=%.6e pressure_error=%.6e\n", global.maxima[1], global.maxima[2]);
        return ok;
    }
    void write(const Mesh& m, const std::string& prefix)
    {
        if (prefix.empty()) return;
        std::vector<Node> states(n);
#if defined(__CUDACC__) || defined(__HIPCC__)
        thrust::copy(state.begin(), state.begin() + n, states.begin());
#else
        std::copy(state.begin(), state.begin() + n, states.begin());
#endif
        auto rows        = host<int64_t>(pressure.rowIds);
        std::string path = prefix + ".rank" + std::to_string(rankId) + ".csv";
        if (std::filesystem::exists(path)) throw std::runtime_error("output already exists");
        std::ofstream out(path);
        out.precision(17);
        out << "node_id,x,y,z,velocity_x,velocity_y,velocity_z,pressure\n";
        for (int r = 0; r < n; ++r)
        {
            int id        = int(rows[r]);
            const auto& x = m.x[std::lower_bound(m.dofs.begin(), m.dofs.end(), id) - m.dofs.begin()];
            const auto& s = states[r].flow;
            out << id + 1 << ',' << x.x << ',' << x.y << ',' << x.z << ',' << s.velocity.x << ',' << s.velocity.y << ','
                << s.velocity.z << ',' << s.pressure << '\n';
        }
        out.close();
        if (!out) throw std::runtime_error("flow output failed");
    }
    int n;
    NonconformalInterfacePlan<double, device> plan;
    NonconformalDonorExchange<double, Node> halo;
    Graph pressure, momentum;
    Buffer<Cell> cells;
    Buffer<Face> faces;
    Buffer<Node> state;
    Buffer<Vec3<double>> coordinates, gradient, influence;
    Buffer<int64_t> globalDofs, stateIds;
    Buffer<int> rowLocal;
    Buffer<int> slot, boundary, cellOffsets, cellIds, faceOffsets, faceIds, patchOffsets, patchIds, activeCell,
        activeFace;
    Buffer<double> interiorMass, boundaryMass, interfaceMass;
    View view;
};
} // namespace ncflow
