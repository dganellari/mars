
#include "mars_vtk_IO.hpp"

namespace mars {
    namespace vtk {

        using Hex8Mesh = ::mars::DistributedMesh<::mars::ElementType::Hex8>;
        using Hex8DofHandler = ::mars::DofHandler<Hex8Mesh, 1>;
        using Hex8FEDofMap = ::mars::FEDofMap<Hex8DofHandler>;

        template class IO<Hex8FEDofMap>;

    }  // namespace vtk

}  // namespace mars