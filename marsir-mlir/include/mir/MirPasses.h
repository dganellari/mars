#ifndef MIR_MIRPASSES_H
#define MIR_MIRPASSES_H

namespace mir {
// Registers the mir lowering passes (convert-mir-to-linalg) with the pass
// registry so mir-opt can run them by name.
void registerMirPasses();
// Registers --mir-hoist-transfer-pairs (accumulator registerization).
void registerHoistTransferPairsPass();
// Registers --mir-lower-copies (memref.copy -> scf.for loops, GPU-legal).
void registerLowerCopiesPass();
// Registers --mir-warp-distribute (lane distribution of warp regions).
void registerWarpDistributePass();
}  // namespace mir

#endif  // MIR_MIRPASSES_H
