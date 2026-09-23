#ifndef MIR_MIRPASSES_H
#define MIR_MIRPASSES_H

namespace mir {
// Set by --mir-chain-contracts on a C-fragment transfer: lane L touches only
// its own entries [i, 2k..2k+1] (i = L/4, k = L%4) of the window its indices
// select. --mir-warp-barriers relies on it to skip barriers between accesses
// that provably touch the same entries in the same lanes.
inline constexpr const char *kLaneOwnedAttr = "mir.lane_owned";

// Registers the mir lowering passes (convert-mir-to-linalg) with the pass
// registry so mir-opt can run them by name.
void registerMirPasses();
// Registers --mir-hoist-transfer-pairs (accumulator registerization).
void registerHoistTransferPairsPass();
// Registers --mir-lower-copies (memref.copy -> scf.for loops, GPU-legal).
void registerLowerCopiesPass();
// Registers --mir-warp-distribute (lane distribution of warp regions).
void registerWarpDistributePass();
// Registers --mir-warp-wrap (wrap per-element bodies in warp regions).
void registerWarpWrapPass();
// Registers --mir-chain-contracts (register-resident mma+shuffle chaining).
void registerChainContractsPass();

// --mir-batch-elements
void registerBatchElementsPass();

// --mir-gpu-wrap
void registerGpuWrapPass();

// --mir-forward-transfers
void registerForwardTransfersPass();
void registerEmulateWarpPass();
void registerWorkgroupBuffersPass();
void registerDistributeFillsPass();
void registerWarpBarriersPass();
void registerHoistInvariantReadsPass();
void registerUnrollLoopsPass();
void registerForwardOwnedPass();
}  // namespace mir

#endif  // MIR_MIRPASSES_H
