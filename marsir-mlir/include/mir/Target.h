#ifndef MIR_TARGET_H
#define MIR_TARGET_H

#include "mlir/IR/Operation.h"
#include "llvm/ADT/StringRef.h"

#include <cstdint>

namespace mir {

// The GPU properties the MARSIR passes depend on. They come from the module's
// dlti.target_system_spec, device "GPU" (set by --mir-set-target); a module without one
// gets the NVIDIA values every committed kernel was generated with, so passes that read
// the target change nothing for those kernels.
struct GpuTarget {
  llvm::StringRef vendor = "nvidia";   // "nvidia" | "amd"
  int64_t waveSize = 32;               // lanes that execute together (warp / wavefront)
  int64_t maxSharedBytes = 49152;      // static shared memory (LDS) per block
  bool isAmd() const { return vendor == "amd"; }
};

inline constexpr const char *kVendorKey = "mir.vendor";
inline constexpr const char *kWaveSizeKey = "mir.wave_size";
inline constexpr const char *kMaxSharedBytesKey = "mir.max_shared_bytes";

// The target of the closest enclosing module (any op inside it may be passed).
GpuTarget getGpuTarget(mlir::Operation *op);

// Registers --mir-set-target=target=nvidia|amd.
void registerSetTargetPass();

}  // namespace mir

#endif  // MIR_TARGET_H
