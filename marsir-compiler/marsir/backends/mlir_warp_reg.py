# Register-resident emission for the HO operator: the sum-factorization face
# chain emitted as DIRECT per-lane nvgpu.mma.sync + gpu.shuffle -- NO
# warp_execute_on_lane_0 region, NO shared memory, NO barriers. Intermediates
# (dt1, dt2, flux, tmp, intf) stay as per-lane C fragments (vector<1x2>) in
# registers; only the inputs (read from memory) and intf (written) touch memory.
# This empties the L1/shared-instruction pipe that bounds the warp-region path.
#
# Proven foundation: the C-fragment -> A-fragment gpu.shuffle relayout is
# GH200-bit-exact (test/warp_chain_reg.mlir, run_chain_reg). This module builds
# the full face chain on it.
#
# m8n8k4 f64 fragment conventions (row.col), lane L, i=L/4, k=L%4:
#   A-frag       lane L = A[i, k]                (read mem [i, 4s+k])
#   B-frag std   lane L = B[k, i]                (read mem [4s+k, i])
#   B-frag transp (stored [n,k]) = read mem [i, 4s+k]  (same index as A)
#   C-frag       lane L = C[i,2k], C[i,2k+1]     (vector<1x2>, read [i,2k])
# Relayouts (register-only, gpu.shuffle idx, width 32):
#   C->A slab s: src=(L&28)|(2s)|((L&3)>>1),  elem=L&1
#   C->B std slab s: src=4*(4s+L%4)+((L/4)>>1)=16s+4*(L&3)+(L>>3), elem=(L>>2)&1

FRAG1 = "vector<1x1xf64>"
FRAG2 = "vector<1x2xf64>"
MR8 = "memref<8x8xf64>"


class Ctr(object):
    def __init__(self):
        self.i = 0

    def fresh(self, tag):
        self.i += 1
        return "%%%s%d" % (tag, self.i)


def _read1(L, ctr, J, mem, mem_ty, row, col):
    r = ctr.fresh("f")
    L.append("%s%s = vector.transfer_read %s[%s, %s], %%z {in_bounds=[true,true]}"
             " : %s, %s" % (J, r, mem, row, col, mem_ty, FRAG1))
    return r


def _read2(L, ctr, J, mem, mem_ty, row, col):
    r = ctr.fresh("cf")
    L.append("%s%s = vector.transfer_read %s[%s, %s], %%z {in_bounds=[true,true]}"
             " : %s, %s" % (J, r, mem, row, col, mem_ty, FRAG2))
    return r


def _mma(L, ctr, J, a, b, acc):
    r = ctr.fresh("m")
    L.append("%s%s = nvgpu.mma.sync(%s, %s, %s) {mmaShape=[8,8,4]} : "
             "(%s, %s, %s) -> %s" % (J, r, a, b, acc, FRAG1, FRAG1, FRAG2, FRAG2))
    return r


def _contract2(L, ctr, J, n, a_reader, b_reader):
    """Two k-slabs chained: acc = slab0 + slab1, starting from a zero C-frag.
    a_reader/b_reader(slab) return the A/B fragment SSA for that slab."""
    acc = "%zc2"
    for s in range(n // 4):
        a = a_reader(s)
        b = b_reader(s)
        acc = _mma(L, ctr, J, a, b, acc)
    return acc


def _relayout_c2a(L, ctr, J, cfrag, slab):
    """C-fragment -> A-fragment for k-slab `slab`: A[i,k]=C[i,4s+k].
    src=(L&28)|(2s)|((L&3)>>1), elem=L&1. Register-only gpu.shuffle."""
    lo = ctr.fresh("lo"); hi = ctr.fresh("hi")
    L.append("%s%s = vector.extract %s[0, 0] : f64 from %s" % (J, lo, cfrag, FRAG2))
    L.append("%s%s = vector.extract %s[0, 1] : f64 from %s" % (J, hi, cfrag, FRAG2))
    base = ctr.fresh("sb")
    # src = (li & 28) | (2*slab) | ((li & 3) >> 1)
    m28 = ctr.fresh("t"); L.append("%s%s = arith.andi %%li, %%j28 : i32" % (J, m28))
    m3 = ctr.fresh("t"); L.append("%s%s = arith.andi %%li, %%j3 : i32" % (J, m3))
    sh = ctr.fresh("t"); L.append("%s%s = arith.shrui %s, %%j1 : i32" % (J, sh, m3))
    o1 = ctr.fresh("t"); L.append("%s%s = arith.ori %s, %s : i32" % (J, o1, m28, sh))
    if slab:
        src = ctr.fresh("t"); L.append("%s%s = arith.ori %s, %%j2 : i32" % (J, src, o1))
    else:
        src = o1
    s0 = ctr.fresh("s"); L.append("%s%s, %%unused%d = gpu.shuffle idx %s, %s, %%j32 : f64"
                                  % (J, s0, ctr.i, lo, src)); ctr.i += 1
    s1 = ctr.fresh("s"); L.append("%s%s, %%unused%d = gpu.shuffle idx %s, %s, %%j32 : f64"
                                  % (J, s1, ctr.i, hi, src)); ctr.i += 1
    par = ctr.fresh("t"); L.append("%s%s = arith.andi %%li, %%j1 : i32" % (J, par))
    odd = ctr.fresh("t"); L.append("%s%s = arith.cmpi eq, %s, %%j1 : i32" % (J, odd, par))
    sel = ctr.fresh("sv"); L.append("%s%s = arith.select %s, %s, %s : f64" % (J, sel, odd, s1, s0))
    a = ctr.fresh("a"); L.append("%s%s = vector.insert %s, %%zc1 [0, 0] : f64 into %s"
                                 % (J, a, sel, FRAG1))
    return a


def _relayout_c2b(L, ctr, J, cfrag, slab):
    """C-fragment -> B-fragment (standard) for k-slab `slab`: B[k,n]=C[4s+k, n]
    i.e. the mma consumes tmp[4s+L%4, L/4]. src=16s+4*(L&3)+(L>>3), elem=(L>>2)&1."""
    lo = ctr.fresh("lo"); hi = ctr.fresh("hi")
    L.append("%s%s = vector.extract %s[0, 0] : f64 from %s" % (J, lo, cfrag, FRAG2))
    L.append("%s%s = vector.extract %s[0, 1] : f64 from %s" % (J, hi, cfrag, FRAG2))
    m3 = ctr.fresh("t"); L.append("%s%s = arith.andi %%li, %%j3 : i32" % (J, m3))
    q4 = ctr.fresh("t"); L.append("%s%s = arith.muli %s, %%j4 : i32" % (J, q4, m3))
    hi3 = ctr.fresh("t"); L.append("%s%s = arith.shrui %%li, %%j3 : i32" % (J, hi3))
    o1 = ctr.fresh("t"); L.append("%s%s = arith.addi %s, %s : i32" % (J, o1, q4, hi3))
    if slab:
        src = ctr.fresh("t"); L.append("%s%s = arith.addi %s, %%j16 : i32" % (J, src, o1))
    else:
        src = o1
    s0 = ctr.fresh("s"); L.append("%s%s, %%unused%d = gpu.shuffle idx %s, %s, %%j32 : f64"
                                  % (J, s0, ctr.i, lo, src)); ctr.i += 1
    s1 = ctr.fresh("s"); L.append("%s%s, %%unused%d = gpu.shuffle idx %s, %s, %%j32 : f64"
                                  % (J, s1, ctr.i, hi, src)); ctr.i += 1
    # element = (L>>2) & 1
    sh2 = ctr.fresh("t"); L.append("%s%s = arith.shrui %%li, %%j2 : i32" % (J, sh2))
    par = ctr.fresh("t"); L.append("%s%s = arith.andi %s, %%j1 : i32" % (J, par, sh2))
    odd = ctr.fresh("t"); L.append("%s%s = arith.cmpi eq, %s, %%j1 : i32" % (J, odd, par))
    sel = ctr.fresh("sv"); L.append("%s%s = arith.select %s, %s, %s : f64" % (J, sel, odd, s1, s0))
    b = ctr.fresh("b"); L.append("%s%s = vector.insert %s, %%zc1 [0, 0] : f64 into %s"
                                 % (J, b, sel, FRAG1))
    return b


def emit_face_reg(L, ctr, J, n, mems, write=True):
    """Register-resident face chain. If write, stores intf (C-frag) to
    mems['intf']; ALWAYS returns the intf C-frag SSA so a register-resident
    scatter can consume it without a memory round-trip.
       dt2 = interp @ Dm^T (transp) ; dt1 = Dm @ interp (std) ; no relayout
       flux = g2*deriv + g0*dt2 + g1*dt1  (per-lane C-frag)
       tmp  = flux @ W^T (transp)  -- flux C->A relayout
       intf = W @ tmp  (std)       -- tmp  C->B relayout
    mems: {interp,deriv,Dm,W,g0,g1,g2[,intf]} nxn memref SSA names (+ '_ty')."""
    def ty(kk):
        return mems.get(kk + "_ty", MR8)
    Im, Dm, Wm = mems["interp"], mems["Dm"], mems["W"]
    # dt1 = Dm @ interp (standard): A=Dm[i,4s+k], B=interp[4s+k,i]
    dt1 = _contract2(L, ctr, J, n,
        lambda s: _read1(L, ctr, J, Dm, ty("Dm"), "%i", "%k" if s == 0 else "%k4"),
        lambda s: _read1(L, ctr, J, Im, ty("interp"), "%k" if s == 0 else "%k4", "%i"))
    # dt2 = interp @ Dm^T (transposed): A=interp[i,4s+k], B(transp)=Dm[i,4s+k]
    dt2 = _contract2(L, ctr, J, n,
        lambda s: _read1(L, ctr, J, Im, ty("interp"), "%i", "%k" if s == 0 else "%k4"),
        lambda s: _read1(L, ctr, J, Dm, ty("Dm"), "%i", "%k" if s == 0 else "%k4"))
    # flux = g2*deriv + g0*dt2 + g1*dt1  (C-frag pointwise)
    dv = _read2(L, ctr, J, mems["deriv"], ty("deriv"), "%i", "%tc")
    g2c = _read2(L, ctr, J, mems["g2"], ty("g2"), "%i", "%tc")
    g0c = _read2(L, ctr, J, mems["g0"], ty("g0"), "%i", "%tc")
    g1c = _read2(L, ctr, J, mems["g1"], ty("g1"), "%i", "%tc")
    t = ctr.fresh("x"); L.append("%s%s = arith.mulf %s, %s : %s" % (J, t, g2c, dv, FRAG2))
    p0 = ctr.fresh("x"); L.append("%s%s = arith.mulf %s, %s : %s" % (J, p0, g0c, dt2, FRAG2))
    p1 = ctr.fresh("x"); L.append("%s%s = arith.mulf %s, %s : %s" % (J, p1, g1c, dt1, FRAG2))
    a01 = ctr.fresh("x"); L.append("%s%s = arith.addf %s, %s : %s" % (J, a01, t, p0, FRAG2))
    flux = ctr.fresh("fl"); L.append("%s%s = arith.addf %s, %s : %s" % (J, flux, a01, p1, FRAG2))
    # tmp = flux @ W^T (transposed): A=relayout_C2A(flux,s), B(transp)=W[i,4s+k]
    tmp = _contract2(L, ctr, J, n,
        lambda s: _relayout_c2a(L, ctr, J, flux, s),
        lambda s: _read1(L, ctr, J, Wm, ty("W"), "%i", "%k" if s == 0 else "%k4"))
    # intf = W @ tmp (standard): A=W[i,4s+k], B=relayout_C2B(tmp,s)
    intf = _contract2(L, ctr, J, n,
        lambda s: _read1(L, ctr, J, Wm, ty("W"), "%i", "%k" if s == 0 else "%k4"),
        lambda s: _relayout_c2b(L, ctr, J, tmp, s))
    if write:
        L.append("%svector.transfer_write %s, %s[%%i, %%tc] {in_bounds=[true,true]}"
                 " : %s, %s" % (J, intf, mems["intf"], FRAG2, ty("intf")))
    return intf


def emit_scatter_reg(L, ctr, J, n, intf, ym, yty, yp, ypty):
    """Register-resident +/- plane scatter: y[l] -= intf ; y[l+1] += intf, per
    lane on C-frags (no smem, no warp region). ym/yp are the plane l / l+1
    memref subviews (their types yty/ypty). Each lane owns disjoint (i,2k) cells
    across all planes, so no cross-lane hazard within a direction; sequential
    program order covers the l/l+1 loop-carried dep. A barrier is only needed
    BETWEEN directions (different lane->cell mapping), handled by the caller."""
    r0 = ctr.fresh("y")
    L.append("%s%s = vector.transfer_read %s[%%i, %%tc], %%z {in_bounds=[true,true]}"
             " : %s, %s" % (J, r0, ym, yty, FRAG2))
    m0 = ctr.fresh("y")
    L.append("%s%s = arith.subf %s, %s : %s" % (J, m0, r0, intf, FRAG2))
    L.append("%svector.transfer_write %s, %s[%%i, %%tc] {in_bounds=[true,true]} : "
             "%s, %s" % (J, m0, ym, FRAG2, yty))
    r1 = ctr.fresh("y")
    L.append("%s%s = vector.transfer_read %s[%%i, %%tc], %%z {in_bounds=[true,true]}"
             " : %s, %s" % (J, r1, yp, ypty, FRAG2))
    p1 = ctr.fresh("y")
    L.append("%s%s = arith.addf %s, %s : %s" % (J, p1, r1, intf, FRAG2))
    L.append("%svector.transfer_write %s, %s[%%i, %%tc] {in_bounds=[true,true]} : "
             "%s, %s" % (J, p1, yp, FRAG2, ypty))



if __name__ == "__main__":
    import sys
    n = 8
    ctr = Ctr()
    args = ", ".join("%%%s: %s" % (a, MR8) for a in
                     ["interp", "deriv", "Dm", "W", "g0", "g1", "g2", "intf"])
    L = ['gpu.module @mir_kernels [#nvvm.target<chip = "sm_90", O = 3>] {',
         "  gpu.func @face_reg(%s) kernel {" % args,
         "    %z = arith.constant 0.0 : f64",
         "    %zc1 = arith.constant dense<0.0> : " + FRAG1,
         "    %zc2 = arith.constant dense<0.0> : " + FRAG2,
         "    %c2i = arith.constant 2 : index",
         "    %c4i = arith.constant 4 : index",
         "    %lane = gpu.thread_id x",
         "    %i = arith.divui %lane, %c4i : index",
         "    %k = arith.remui %lane, %c4i : index",
         "    %k4 = arith.addi %k, %c4i : index",
         "    %tc = arith.muli %k, %c2i : index",
         "    %li = arith.index_cast %lane : index to i32",
         "    %j1 = arith.constant 1 : i32",
         "    %j2 = arith.constant 2 : i32",
         "    %j3 = arith.constant 3 : i32",
         "    %j4 = arith.constant 4 : i32",
         "    %j16 = arith.constant 16 : i32",
         "    %j28 = arith.constant 28 : i32",
         "    %j32 = arith.constant 32 : i32"]
    mems = {kk: "%" + kk for kk in
            ["interp", "deriv", "Dm", "W", "g0", "g1", "g2", "intf"]}
    emit_face_reg(L, ctr, "    ", n, mems)
    L += ["    gpu.return", "  }", "}"]
    sys.stdout.write("\n".join(L) + "\n")
