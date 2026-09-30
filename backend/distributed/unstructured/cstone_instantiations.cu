#ifdef MARS_ENABLE_CUDA
// Cornerstone instantiates these for other key/real combinations only; MARS also uses unsigned keys with double
#include <cstone/traversal/collisions_gpu.cu>

namespace cstone
{
MARK_MACS_GPU(unsigned, double);
FIND_HALOS_GPU(uint32_t, double);
} // namespace cstone
#endif
