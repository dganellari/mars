#ifndef MARS_DISTRIBUTED_DOF_HPP_
#define MARS_DISTRIBUTED_DOF_HPP_

#include "mars_globals.hpp"

namespace mars {
    struct Dof {
        Integer gid;
        Integer owner_proc;
        bool valid = true;

        MARS_INLINE_FUNCTION
        Dof(Integer g, Integer p) : gid(g), owner_proc(p) {}

        MARS_INLINE_FUNCTION
        Dof() = default;

        MARS_INLINE_FUNCTION
        void set_gid(const Integer g) { gid = g; }

        MARS_INLINE_FUNCTION
        Integer get_gid() const { return gid; }

        MARS_INLINE_FUNCTION
        void set_proc(const Integer p) { owner_proc = p; }

        MARS_INLINE_FUNCTION
        Integer get_proc() const { return owner_proc; }

        MARS_INLINE_FUNCTION
        void set_invalid() { valid = false; }

        MARS_INLINE_FUNCTION
        bool is_valid() const { return valid; }
    };
}  // namespace mars

#endif
