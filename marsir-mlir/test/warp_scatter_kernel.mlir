// GPU gate for the +/- plane SCATTER (the census's hardest blocker), direction
// 0, two faces with the loop-carried plane overlap: plane l+1 is written by
// face l's += then face l+1's -=. Each face:
//   R_compute: intf = A[:,4l:4l+4] @ B[4l:4l+4,:]  -> scratch Sf   (contract used ONCE)
//   R_scatter: y[l] -= Sf ; y[l+1] += Sf           (Sf re-read TWICE, distinct SSA)
// Staging intf through Sf is REQUIRED: a contract result feeding two in-region
// elementwise ops crashes upstream WarpOpElementwise (MLIR 19). Barriers order
// the smem/global RMW. Real design uses shared memory for Sf and Y.
// Launch: grid 1, block (32,1,1). Result: y0=-intf0, y1=intf0-intf1, y2=intf1.
#a = affine_map<(m, n, k) -> (m, k)>
#b = affine_map<(m, n, k) -> (k, n)>
#c = affine_map<(m, n, k) -> (m, n)>
gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {
  gpu.func @scatter(%A: memref<8x8xf64>, %B: memref<8x8xf64>,
                    %Sf: memref<8x8xf64>, %Y: memref<8x8x8xf64>) kernel {
    %c0 = arith.constant 0 : index
    %c4 = arith.constant 4 : index
    %z  = arith.constant 0.0 : f64
    %zc = arith.constant dense<0.0> : vector<8x8xf64>
    %lane = gpu.thread_id x
    %y0 = memref.subview %Y[0,0,0][1,8,8][1,1,1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8,1], offset: 0>>
    %y1 = memref.subview %Y[1,0,0][1,8,8][1,1,1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8,1], offset: 64>>
    %y2 = memref.subview %Y[2,0,0][1,8,8][1,1,1] : memref<8x8x8xf64> to memref<8x8xf64, strided<[8,1], offset: 128>>

    // face 0: intf0 = A[:,0:4] @ B[0:4,:] -> Sf
    vector.warp_execute_on_lane_0(%lane)[32] {
      %a0 = vector.transfer_read %A[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b0 = vector.transfer_read %B[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<4x8xf64>
      %f0 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a0, %b0, %zc : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %f0, %Sf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa = vector.transfer_read %Sf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
      %r0 = vector.transfer_read %y0[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8,1], offset: 0>>, vector<8x8xf64>
      %m0 = arith.subf %r0, %fa : vector<8x8xf64>
      vector.transfer_write %m0, %y0[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8,1], offset: 0>>
      %fb = vector.transfer_read %Sf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
      %r1 = vector.transfer_read %y1[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8,1], offset: 64>>, vector<8x8xf64>
      %p1 = arith.addf %r1, %fb : vector<8x8xf64>
      vector.transfer_write %p1, %y1[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8,1], offset: 64>>
    }
    gpu.barrier

    // face 1: intf1 = A[:,4:8] @ B[4:8,:] -> Sf
    vector.warp_execute_on_lane_0(%lane)[32] {
      %a1 = vector.transfer_read %A[%c0, %c4], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x4xf64>
      %b1 = vector.transfer_read %B[%c4, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<4x8xf64>
      %f1 = vector.contract {indexing_maps=[#a,#b,#c], iterator_types=["parallel","parallel","reduction"], kind=#vector.kind<add>} %a1, %b1, %zc : vector<8x4xf64>, vector<4x8xf64> into vector<8x8xf64>
      vector.transfer_write %f1, %Sf[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64>
    }
    gpu.barrier
    vector.warp_execute_on_lane_0(%lane)[32] {
      %fa = vector.transfer_read %Sf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
      %r1 = vector.transfer_read %y1[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8,1], offset: 64>>, vector<8x8xf64>
      %m1 = arith.subf %r1, %fa : vector<8x8xf64>
      vector.transfer_write %m1, %y1[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8,1], offset: 64>>
      %fb = vector.transfer_read %Sf[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64>, vector<8x8xf64>
      %r2 = vector.transfer_read %y2[%c0, %c0], %z {in_bounds=[true,true]} : memref<8x8xf64, strided<[8,1], offset: 128>>, vector<8x8xf64>
      %p2 = arith.addf %r2, %fb : vector<8x8xf64>
      vector.transfer_write %p2, %y2[%c0, %c0] {in_bounds=[true,true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8,1], offset: 128>>
    }
    gpu.return
  }
}
