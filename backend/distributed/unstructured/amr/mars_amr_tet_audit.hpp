#pragma once

#include "mars_amr_tet_audit_math.hpp"
#include <cstone/cuda/device_vector.h>
#include <cstone/sfc/box.hpp>
#include <array>
#include <cstdlib>
#include <string>

namespace mars::amr
{
// Explicitly enabled only for the public tet diagnostic; normal domain rebuilds do no audit work.
inline bool tet_audit_enabled()
{
    const char* value = std::getenv("MARS_AMR_AUDIT_TET");
    return value && std::string(value) == "1";
}

template<class KeyType>
struct TetAuditSnapshot
{
    TetAuditTotals totals;
    cstone::DeviceVector<char> records;
    std::string stage;
};

// With coordinates, conn contains node indices. Otherwise it contains node keys.
// keyed_geometry makes volume comparable across redistribution, independently of coordinate decoding.
template<class KeyType, class RealType>
TetAuditSnapshot<KeyType> audit_tets(const char* stage,
                                     std::array<const KeyType*, 4> conn,
                                     size_t begin,
                                     size_t end,
                                     const cstone::Box<RealType>& box,
                                     const RealType* x        = nullptr,
                                     const RealType* y        = nullptr,
                                     const RealType* z        = nullptr,
                                     size_t num_nodes         = 0,
                                     bool keyed_geometry      = true,
                                     const KeyType* node_keys = nullptr);

template<class KeyType>
void audit_require(bool condition,
                   const char* reason,
                   const TetAuditSnapshot<KeyType>& current,
                   const TetAuditSnapshot<KeyType>* previous = nullptr);

template<class KeyType>
void audit_compare(const TetAuditSnapshot<KeyType>& before, const TetAuditSnapshot<KeyType>& after);

template<class KeyType>
void audit_red_children(const TetAuditSnapshot<KeyType>& parents,
                        const TetAuditSnapshot<KeyType>& children,
                        const uint8_t* marks,
                        const uint64_t* marked_prefix);
} // namespace mars::amr
