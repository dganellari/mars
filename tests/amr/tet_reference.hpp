#pragma once

#include "tet_checks.hpp"
#include "backend/distributed/unstructured/amr/mars_amr_lepp.hpp"

namespace amr_test
{
using HostResult = mars::amr::TetLeppResult<mars::amr::LeppHostExecution, KeyType, RealType>;
HostResult refine_reference(mars::amr::TetLeppView<KeyType, RealType> view,
                            const cstone::Box<RealType>& box,
                            const std::vector<mars::amr::LeppTaggedFace>& faces = {});
} // namespace amr_test
