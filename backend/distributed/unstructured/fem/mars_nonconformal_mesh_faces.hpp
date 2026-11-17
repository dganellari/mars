#pragma once
#include "mars_nonconformal_scalar.hpp"
#include <cstddef>

namespace mars::fem
{
// Side-set entries name a local tetrahedron and its opposite corner, avoiding
// coordinate matching of coincident nodes on the two independent mesh blocks.
struct InterfaceFaceTag
{
    int element, opposite, side;
};
template<class T>
struct InterfaceTetMeshView
{
    int nodes, elements;
    const Vec3<T>* coordinates;
    const int* connectivity; // four local node indices per tetrahedron
    const int64_t* globalDofs;
    const uint64_t* globalElements;
    const int* elementOwners;
    const T* diffusivity;
};
template<class T>
struct ExtractInterfaceFace
{
    InterfaceTetMeshView<T> mesh;
    const InterfaceFaceTag* tags;
    InterfaceTetFace<T>* faces;
    uint64_t* ids;
    int* owners;
    int* valid;
    MARS_FS_FN void operator()(int i) const
    {
        valid[i] = 0;
        auto tag = tags[i];
        if (tag.element < 0 || tag.element >= mesh.elements || tag.opposite < 0 || tag.opposite >= 4 || tag.side < 0 ||
            tag.side > 1)
            return;
        uint64_t element = mesh.globalElements[tag.element];
        if (element > (std::numeric_limits<uint64_t>::max() - 3) / 4) return;
        InterfaceTetFace<T> f{};
        f.id          = 4 * element + tag.opposite;
        f.side        = tag.side;
        f.opposite    = tag.opposite;
        f.diffusivity = mesh.diffusivity[tag.element];
        for (int j = 0; j < 4; ++j)
        {
            int node = mesh.connectivity[4 * std::size_t(tag.element) + j];
            if (node < 0 || node >= mesh.nodes) return;
            f.x[j]   = mesh.coordinates[node];
            f.dof[j] = mesh.globalDofs[node];
        }
        nonconformal_detail::Geometry<T> g;
        if (nonconformal_detail::geometry(f, g) != InterfaceStatus::valid) return;
        faces[i]  = f;
        ids[i]    = f.id;
        owners[i] = mesh.elementOwners[tag.element];
        valid[i]  = 1;
    }
};
} // namespace mars::fem
