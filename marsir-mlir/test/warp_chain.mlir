func.func @batched_flux(%A: memref<?x8x8xf64>, %B: memref<?x8x8xf64>, %O: memref<?x8x8xf64>) {
  %c0 = arith.constant 0 : index
  %c1 = arith.constant 1 : index
  %c8 = arith.constant 8 : index
  %E = memref.dim %A, %c0 : memref<?x8x8xf64>
  %Bk = arith.divui %E, %c8 : index
  scf.parallel (%b) = (%c0) to (%Bk) step (%c1) {
    %be = arith.muli %b, %c8 : index
    scf.parallel (%w) = (%c0) to (%c8) step (%c1) {
      %e = arith.addi %be, %w : index
      %as = memref.subview %A[%e, 0, 0] [1, 8, 8] [1, 1, 1] : memref<?x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %bs = memref.subview %B[%e, 0, 0] [1, 8, 8] [1, 1, 1] : memref<?x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %os = memref.subview %O[%e, 0, 0] [1, 8, 8] [1, 1, 1] : memref<?x8x8xf64> to memref<8x8xf64, strided<[8, 1], offset: ?>>
      %z = arith.constant 0.0 : f64
      %va = vector.transfer_read %as[%c0, %c0], %z {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %vb = vector.transfer_read %bs[%c0, %c0], %z {in_bounds = [true, true]} : memref<8x8xf64, strided<[8, 1], offset: ?>>, vector<8x8xf64>
      %m = arith.mulf %va, %vb : vector<8x8xf64>
      %s = arith.addf %m, %va : vector<8x8xf64>
      vector.transfer_write %s, %os[%c0, %c0] {in_bounds = [true, true]} : vector<8x8xf64>, memref<8x8xf64, strided<[8, 1], offset: ?>>
      scf.reduce
    } {mapping = [#gpu.loop_dim_map<processor = thread_y, map = (d0) -> (d0), bound = (d0) -> (d0)>]}
    scf.reduce
  } {mapping = [#gpu.loop_dim_map<processor = block_x, map = (d0) -> (d0), bound = (d0) -> (d0)>]}
  return
}
