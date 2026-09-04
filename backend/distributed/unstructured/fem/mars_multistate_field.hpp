#pragma once

// Generic previous-timestep / previous-iteration field history, generalizing the hand-rolled
// pattern already proven in mars_ns_solver.hpp's BDF2/EXT2 velocity + advection-accumulator
// history (NSStepper's d_u_nm1/d_v_nm1/d_w_nm1 and d_advU_n/d_advU_nm1 etc., shuffled each step
// via explicit thrust::copy at mars_ns_solver.hpp ~9742-9777). Closes the "no FieldState" gap
// flagged against STK's stk::mesh::FieldState/field_of_state -- see project memory
// project_stk_mars_gap_eval / internal-notes/stk-mars.md. Any new caller that needs N historical
// copies of a field declares one MultiStateField instead of hand-rolling N separate buffers plus a
// manual copy shuffle. The existing NS solver is NOT migrated onto this (its hand-rolled pattern
// is validated in production); this is for new code.
//
// Implemented as index rotation over NStates fixed buffers, not a copy shuffle: advance() just
// rotates which buffer is "current", so it costs O(1) regardless of field size, unlike the
// thrust::copy pattern it generalizes (which moves the whole field every step).

#include <array>
#include <cstddef>
#include "domain.hpp"

namespace mars
{
namespace fem
{

// RealType/AcceleratorTag pick the storage type via ElementDomain's VectorSelector (DeviceVector on
// GPU, std::vector on CPU), so this works with either backend. NStates is compile-time: 2 covers
// BDF2/EXT2 (current + n-1); pass more for higher-order multistep schemes.
template<typename RealType, typename AcceleratorTag, size_t NStates>
class MultiStateField
{
public:
    static_assert(NStates >= 1, "MultiStateField needs at least one state");

    using Buffer = typename VectorSelector<RealType, AcceleratorTag>::type;

    void resize(size_t n)
    {
        for (auto& buf : buffers_) buf.resize(n);
    }

    size_t size() const { return buffers_[ringIndex(0)].size(); }

    // state(0) = current ("now"), state(1) = previous, state(2) = previous-previous, ...
    Buffer&       state(size_t stepsBack)       { return buffers_[ringIndex(stepsBack)]; }
    const Buffer& state(size_t stepsBack) const { return buffers_[ringIndex(stepsBack)]; }

    Buffer&       current()       { return state(0); }
    const Buffer& current() const { return state(0); }

    // Rotate the ring by one step: every state(k) becomes state(k+1), and the buffer that WAS the
    // oldest state (NStates-1) becomes the new state(0).
    //
    // WARNING -- this differs from the copy-shuffle it replaces: right after advance(), current()
    // holds STALE data (the values from NStates steps ago), NOT the previous step's values. The
    // caller must FULLY overwrite current() before reading it. A caller that needs the previous
    // value while computing the new one reads state(1), or calls advance() only after the new
    // value is complete.
    void advance() { head_ = (head_ + NStates - 1) % NStates; }

private:
    size_t ringIndex(size_t stepsBack) const { return (head_ + stepsBack) % NStates; }

    std::array<Buffer, NStates> buffers_;
    size_t head_ = 0;   // index of the state(0) ("current") buffer
};

} // namespace fem
} // namespace mars
