// The FP64 tensor-core (DMMA) target for the FEM sum-factorization contraction.
// m8n8k4 f64 mma.sync is exactly the instruction MARS hand-writes in
// mars_ho_laplacian_dmma.hpp, and the 8x8 tile matches the p=7 reference operator.
//
// Proof that MLIR 19 lowers it to the NVVM tensor-core intrinsic:
//   mlir-opt %s --convert-nvgpu-to-nvvm | grep nvvm.mma.sync
//
// Warp-distributed f64 m8n8k4 fragments: A 8x4=32 elts (1/thread), B 4x8=32
// (1/thread), C/D 8x8=64 (2/thread, shaped 1x2).
func.func @dmma_f64(%a: vector<1x1xf64>, %b: vector<1x1xf64>, %c: vector<1x2xf64>)
    -> vector<1x2xf64> {
  %0 = nvgpu.mma.sync(%a, %b, %c) {mmaShape = [8, 8, 4]}
       : (vector<1x1xf64>, vector<1x1xf64>, vector<1x2xf64>) -> vector<1x2xf64>
  return %0 : vector<1x2xf64>
}
