// One direction's face contract chain in TARGET warp form (what the emitter
// must produce): dt = iface . D  (axis-0, 2 k-slabs, standard arm)
//                dts = dt * g     (flux, elementwise, fuses in-register)
//                tmp = dts . W    (axis-1 => transposed-B arm, 2 slabs)
//                intf = tmp . W   (axis-0 => standard arm, 2 slabs)
// staged through scratch memrefs S1/S2 between contract stages (real design =
// shared memory + barrier; here plain memrefs suffice for the distribution +
// PTX-gen proof). Expect 6 mma.sync.
#a  = affine_map<(m, n, k) -> (m, k)>
#b  = affine_map<(m, n, k) -> (k, n)>
#bt = affine_map<(m, n, k) -> (n, k)>
#c  = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @one_dir(%If: memref<8x8xf64>, %D: memref<8x8xf64>, %G: memref<8x8xf64>,
                    %W: memref<8x8xf64>, %S1: memref<8x8xf64>, %S2: memref<8x8xf64>,
                    %Out: memref<8x8xf64>) kernel {
    %c0 = arith.constant 0 : index
    %c4 = arith.constant 4 : index
    %z = arith.constant 0.0 : f64
    %zc = arith.constant dense<0.0> : vector<8x8xf64>
    %lane = gpu.thread_id x
    // Stage 1: dt = D @ If  (axis-0: out[m,n]=sum_k D[m,k]*If[k,n]); scale by g.
    vector.warp_execute_on_lane_0(%lane)[32] {
      %d0 = vector.transfer_read %D[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %d1 = vector.transfer_read %D[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %i0 = vector.transfer_read %If[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<4x8xf64>
      %i1 = vector.transfer_read %If[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<4x8xf64>
      %r0 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %d0, %i0, %zc : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %r1 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %d1, %i1, %r0 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %g  = vector.transfer_read %G[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
      %rs = arith.mulf %r1, %g : vector<8x8xf64>
      vector.transfer_write %rs, %S1[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64>
    }
    gpu.barrier
    // Stage 2: tmp = dts @ W^T  (axis-1: out[m,n]=sum_k dts[m,k]*W[n,k]) -> S2.
    vector.warp_execute_on_lane_0(%lane)[32] {
      %a0 = vector.transfer_read %S1[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %a1 = vector.transfer_read %S1[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %w0 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %w1 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %s0 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a0, %w0, %zc : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      %s1 = vector.contract {indexing_maps=[#a,#bt,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1, %w1, %s0 : vector<8x4xf64>, vector<8x4xf64> into vector<8x8xf64>
      vector.transfer_write %s1, %S2[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64>
    }
    gpu.barrier
    // Stage 3: intf = W @ tmp  (axis-0, standard) -> Out.
    vector.warp_execute_on_lane_0(%lane)[32] {
      %w0 = vector.transfer_read %W[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %w1 = vector.transfer_read %W[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %t0 = vector.transfer_read %S2[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<4x8xf64>
      %t1 = vector.transfer_read %S2[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<4x8xf64>
      %o0 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %w0, %t0, %zc : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      %o1 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %w1, %t1, %o0 : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %o1, %Out[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64>
    }
    gpu.return
  }
}
