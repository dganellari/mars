# Warp-per-element emission for the HO matrix-free operator: the sum-
# factorization contracts become m8n8k4 f64 tensor-core ops instead of scalar
# thread-per-element loops. Every construct here is a hand-verified atom (see
# marsir-mlir/test/warp_*.mlir, all GH200 bit-exact):
#   - an nxn contract along axis 0 (standard A@B) or axis 1 (transposed B@ .^T)
#     emitted as n/4 chained vector.contract k-slabs inside a warp region;
#   - the pointwise flux, which fuses in-register with its producing contract;
#   - the +/- plane scatter as a per-lane strided-subview RMW.
#
# Load-bearing constraints found while proving the atoms (violate -> upstream
# MLIR 19 crashes or silent non-distribution):
#   * a CONTRACT RESULT must feed AT MOST ONE in-region elementwise op; the
#     scatter's two planes each RE-READ intf from scratch (distinct SSA).
#   * B/scratch memrefs read inside a warp region must be DEFINED OUTSIDE it
#     (moveScalarUniformCode hoists the subviews; we emit them outside anyway).
#   * scratch + y live in workgroup (shared) memory; barriers order the stages.
#
# Contract axis semantics (matches mir.contract / the Knaus oracle):
#   axis 0 standard:   out[m,n] = sum_k M[m,k] * X[k,n]   (M @ X)
#   axis 1 transposed: out[m,n] = sum_k X[m,k] * M[n,k]   (X @ M^T)
# For the emitter M is the small operator matrix (Dm/W), X the data.
#
# Operator matrices (Knaus Alg-2 sum-factorization symbols -- same names as the
# hand CUDA kernel, mlir_ir.py, exec_gate, and buildHoCvfemOperators, so the
# vocabulary matches across all of MARS and the paper):
#   Btil  (B-tilde)  trial-basis interpolation at the quadrature points
#   Dtil  (D-tilde)  basis derivative at the quadrature points
#   D  /  Dm         1D differentiation matrix (tangential derivative)
#   W                quadrature / integration weight matrix

WS = "#gpu.address_space<workgroup>"


def _mr(n, space=None):
    return "memref<%dx%dxf64%s>" % (n, n, ", " + space if space else "")


class Ctr(object):
    def __init__(self):
        self.i = 0

    def fresh(self, tag):
        self.i += 1
        return "%%%s%d" % (tag, self.i)


def index_consts(indent, upto, step=4):
    """Declare %c0,%c<step>,...  index constants once at the top of a kernel.
    The emit_* helpers reference index offsets by the deterministic name %c<N>
    (N = the integer value), so the caller just declares the covering set here;
    canonicalize DCEs any that go unused. `upto` is exclusive."""
    return ["%s%%c%d = arith.constant %d : index" % (indent, v, v)
            for v in range(0, upto, step)]


def emit_contract(L, ctr, indent, a_mem, a_ty, b_mem, b_ty, out_mem, out_ty,
                  n, transposed, scale_mem=None, scale_ty=None):
    """Emit one nxn contract as n/4 chained k-slab vector.contract ops inside a
    fresh warp region, writing the nxn result to out_mem. If scale_mem is given,
    the result is elementwise-multiplied by it before the write (flux fuses in
    register -- ONE elementwise use of the contract result, which is legal).
    n must be a multiple of 4 (p=3 -> n=4, p=7 -> n=8)."""
    assert n % 4 == 0, "n must be a multiple of 4 (pad smaller orders)"
    slabs = n // 4
    I = indent
    c = lambda name: "%%c%s" % name
    va = "vector<%dx4xf64>" % n
    vb_std = "vector<4x%dxf64>" % n
    vb_t = "vector<%dx4xf64>" % n
    vfull = "vector<%dx%dxf64>" % (n, n)

    std_maps = ('#a = affine_map<(m, n, k) -> (m, k)>\n'
                '#b = affine_map<(m, n, k) -> (k, n)>\n'
                '#c = affine_map<(m, n, k) -> (m, n)>')
    # emitted once at module scope by the caller; here we just reference #a/#b/#bt/#c

    L.append("%svector.warp_execute_on_lane_0(%%lane)[32] {" % I)
    J = I + "  "
    acc = ctr.fresh("acc")
    L.append("%s%s = arith.constant dense<0.0> : %s" % (J, acc, vfull))
    for s in range(slabs):
        koff = 4 * s
        # A slab: M[:, koff:koff+4]  (n x 4)  -- the (m,k) operand
        aslab = ctr.fresh("a")
        L.append("%s%s = vector.transfer_read %s[%%c0, %%c%d], %%z "
                 "{in_bounds=[true,true]} : %s, %s"
                 % (J, aslab, a_mem, koff, a_ty, va))
        # B slab: standard X[koff:koff+4, :] (4 x n); transposed M[:, koff:koff+4] (n x 4)
        bslab = ctr.fresh("b")
        if transposed:
            L.append("%s%s = vector.transfer_read %s[%%c0, %%c%d], %%z "
                     "{in_bounds=[true,true]} : %s, %s"
                     % (J, bslab, b_mem, koff, b_ty, vb_t))
            bmap, bty = "#bt", vb_t
        else:
            L.append("%s%s = vector.transfer_read %s[%%c%d, %%c0], %%z "
                     "{in_bounds=[true,true]} : %s, %s"
                     % (J, bslab, b_mem, koff, b_ty, vb_std))
            bmap, bty = "#b", vb_std
        res = ctr.fresh("r")
        L.append("%s%s = vector.contract {indexing_maps=[#a,%s,#c], "
                 "iterator_types=[\"parallel\",\"parallel\",\"reduction\"], "
                 "kind=#vector.kind<add>} %s, %s, %s : %s, %s into %s"
                 % (J, res, bmap, aslab, bslab, acc, va, bty, vfull))
        acc = res
    out_val = acc
    if scale_mem is not None:
        g = ctr.fresh("g")
        L.append("%s%s = vector.transfer_read %s[%%c0, %%c0], %%z "
                 "{in_bounds=[true,true]} : %s, %s"
                 % (J, g, scale_mem, scale_ty, vfull))
        sc = ctr.fresh("sc")
        L.append("%s%s = arith.mulf %s, %s : %s" % (J, sc, out_val, g, vfull))
        out_val = sc
    L.append("%svector.transfer_write %s, %s[%%c0, %%c0] "
             "{in_bounds=[true,true]} : %s, %s"
             % (J, out_val, out_mem, vfull, out_ty))
    L.append("%s}" % I)


def emit_matmul(L, ctr, indent, a_mem, a_ty, b_mem, b_ty, out_mem, out_ty,
                m, k, ncols, transposed, scale_mem=None, scale_ty=None):
    """Rectangular m x k @ k x ncols -> m x ncols, column-tiled to m8n8k4:
    ncols/8 output col-tiles (each a fresh warp region), k/4 chained slabs. Used
    for the B-sweep (Btil/Dtil P->8 padded @ U reshaped n x n^2). m must be 8,
    k a multiple of 4, ncols a multiple of 8. transposed selects the B arm:
      std   B is (k x ncols),   read [koff, coloff]  as 4x8
      transp B is (ncols x k),  read [coloff, koff]  as 8x4  (X @ M^T)."""
    assert m == 8 and k % 4 == 0 and ncols % 8 == 0
    slabs, ntiles = k // 4, ncols // 8
    I = indent
    va = "vector<8x4xf64>"
    vb_std = "vector<4x8xf64>"
    vb_t = "vector<8x4xf64>"
    v8 = "vector<8x8xf64>"
    for ct in range(ntiles):
        coloff = 8 * ct
        L.append("%svector.warp_execute_on_lane_0(%%lane)[32] {" % I)
        J = I + "  "
        acc = ctr.fresh("acc")
        L.append("%s%s = arith.constant dense<0.0> : %s" % (J, acc, v8))
        for s in range(slabs):
            koff = 4 * s
            aslab = ctr.fresh("a")
            L.append("%s%s = vector.transfer_read %s[%%c0, %%c%d], %%z "
                     "{in_bounds=[true,true]} : %s, %s"
                     % (J, aslab, a_mem, koff, a_ty, va))
            bslab = ctr.fresh("b")
            if transposed:
                L.append("%s%s = vector.transfer_read %s[%%c%d, %%c%d], %%z "
                         "{in_bounds=[true,true]} : %s, %s"
                         % (J, bslab, b_mem, coloff, koff, b_ty, vb_t))
                bmap, bty = "#bt", vb_t
            else:
                L.append("%s%s = vector.transfer_read %s[%%c%d, %%c%d], %%z "
                         "{in_bounds=[true,true]} : %s, %s"
                         % (J, bslab, b_mem, koff, coloff, b_ty, vb_std))
                bmap, bty = "#b", vb_std
            res = ctr.fresh("r")
            L.append("%s%s = vector.contract {indexing_maps=[#a,%s,#c], "
                     "iterator_types=[\"parallel\",\"parallel\",\"reduction\"], "
                     "kind=#vector.kind<add>} %s, %s, %s : %s, %s into %s"
                     % (J, res, bmap, aslab, bslab, acc, va, bty, v8))
            acc = res
        out_val = acc
        if scale_mem is not None:
            g = ctr.fresh("g")
            L.append("%s%s = vector.transfer_read %s[%%c0, %%c%d], %%z "
                     "{in_bounds=[true,true]} : %s, %s"
                     % (J, g, scale_mem, coloff, scale_ty, v8))
            sc = ctr.fresh("sc")
            L.append("%s%s = arith.mulf %s, %s : %s" % (J, sc, out_val, g, v8))
            out_val = sc
        L.append("%svector.transfer_write %s, %s[%%c0, %%c%d] "
                 "{in_bounds=[true,true]} : %s, %s"
                 % (J, out_val, out_mem, coloff, v8, out_ty))
        L.append("%s}" % I)


def emit_flux3(L, ctr, indent, n, deriv_mem, deriv_ty, g2_mem, g2_ty,
               dt2g_mem, dt2g_ty, dt1g_mem, dt1g_ty, out_mem, out_ty):
    """flux = g2*deriv + dt2g + dt1g, where dt2g = g0*dt2 and dt1g = g1*dt1 were
    already scaled in their producing contracts. All operands are scratch/input
    READS (no contract result), so several elementwise ops in one region are
    safe (a contract result feeding 2 elementwise ops is what crashes upstream).
    """
    v = "vector<%dx%dxf64>" % (n, n)
    I, J = indent, indent + "  "
    L.append("%svector.warp_execute_on_lane_0(%%lane)[32] {" % I)
    d = ctr.fresh("d")
    L.append("%s%s = vector.transfer_read %s[%%c0, %%c0], %%z "
             "{in_bounds=[true,true]} : %s, %s" % (J, d, deriv_mem, deriv_ty, v))
    g2 = ctr.fresh("g2")
    L.append("%s%s = vector.transfer_read %s[%%c0, %%c0], %%z "
             "{in_bounds=[true,true]} : %s, %s" % (J, g2, g2_mem, g2_ty, v))
    t = ctr.fresh("t")
    L.append("%s%s = arith.mulf %s, %s : %s" % (J, t, d, g2, v))
    a = ctr.fresh("a")
    L.append("%s%s = vector.transfer_read %s[%%c0, %%c0], %%z "
             "{in_bounds=[true,true]} : %s, %s" % (J, a, dt2g_mem, dt2g_ty, v))
    t2 = ctr.fresh("t2")
    L.append("%s%s = arith.addf %s, %s : %s" % (J, t2, t, a, v))
    b = ctr.fresh("b")
    L.append("%s%s = vector.transfer_read %s[%%c0, %%c0], %%z "
             "{in_bounds=[true,true]} : %s, %s" % (J, b, dt1g_mem, dt1g_ty, v))
    fl = ctr.fresh("fl")
    L.append("%s%s = arith.addf %s, %s : %s" % (J, fl, t2, b, v))
    L.append("%svector.transfer_write %s, %s[%%c0, %%c0] "
             "{in_bounds=[true,true]} : %s, %s" % (J, fl, out_mem, v, out_ty))
    L.append("%s}" % I)


def emit_face(L, ctr, indent, n, mems, scratch):
    """One face of the Knaus apply: the four contracts + the 3-term flux, all
    staged through shared memory, writing intf to mems['intf']. Matches the
    oracle:
       dt2  = interp @ Dm^T   (transposed) ; scaled by g0 -> dt2g
       dt1  = Dm @ interp     (standard)   ; scaled by g1 -> dt1g
       flux = g2*deriv + dt2g + dt1g
       tmp  = flux @ W^T      (transposed)
       intf = W @ tmp         (standard)
    `mems`   : dict of input/output nxn memrefs {interp,deriv,Dm,W,g0,g1,g2,intf}
    `scratch`: dict of workgroup nxn scratch memrefs {dt1g,dt2g,flux,tmp} + their
               type strings under key+'_ty' (workgroup addr space)."""
    mr = _mr(n)
    def ty(k):
        return scratch.get(k + "_ty", mr)
    # dt1 = Dm @ interp (standard), fused *g1 -> dt1g
    emit_contract(L, ctr, indent, mems["Dm"], mr, mems["interp"], mr,
                  scratch["dt1g"], ty("dt1g"), n, transposed=False,
                  scale_mem=mems["g1"], scale_ty=mr)
    # dt2 = interp @ Dm^T (transposed), fused *g0 -> dt2g
    emit_contract(L, ctr, indent, mems["interp"], mr, mems["Dm"], mr,
                  scratch["dt2g"], ty("dt2g"), n, transposed=True,
                  scale_mem=mems["g0"], scale_ty=mr)
    L.append("%sgpu.barrier" % indent)
    emit_flux3(L, ctr, indent, n, mems["deriv"], mr, mems["g2"], mr,
               scratch["dt2g"], ty("dt2g"), scratch["dt1g"], ty("dt1g"),
               scratch["flux"], ty("flux"))
    L.append("%sgpu.barrier" % indent)
    # tmp = flux @ W^T (transposed)
    emit_contract(L, ctr, indent, scratch["flux"], ty("flux"), mems["W"], mr,
                  scratch["tmp"], ty("tmp"), n, transposed=True)
    L.append("%sgpu.barrier" % indent)
    # intf = W @ tmp (standard)
    emit_contract(L, ctr, indent, mems["W"], mr, scratch["tmp"], ty("tmp"),
                  mems["intf"], mr, n, transposed=False)


def module_maps():
    """The three affine maps + the transposed-B map, emitted once at top."""
    return ('#a = affine_map<(m, n, k) -> (m, k)>\n'
            '#b = affine_map<(m, n, k) -> (k, n)>\n'
            '#bt = affine_map<(m, n, k) -> (n, k)>\n'
            '#c = affine_map<(m, n, k) -> (m, n)>')


if __name__ == "__main__":
    # Self-test: emit a kernel doing tmp = M @ X (std) then out = X2 @ M^T
    # (transposed) staged through shared memory, and print it. A driver script
    # checks it distributes to the right mma count.
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else "contract"
    n = 8
    mrg = _mr(n)
    ctr = Ctr()
    if which == "matmul":
        # B-sweep shape: interp = Btil(8x8) @ U(8x64) -> 8x64, column-tiled.
        m8 = "memref<8x8xf64>"
        m64 = "memref<8x64xf64>"
        L = [module_maps(),
             "gpu.module @mir_kernels [#nvvm.target<chip = \"sm_90\", O = 3>] {",
             "  gpu.func @bsweep(%%Bt: %s, %%U: %s, %%Interp: %s) kernel {"
             % (m8, m64, m64),
             "    %z = arith.constant 0.0 : f64",
             "    %lane = gpu.thread_id x"]
        L += index_consts("    ", 64)   # covers koffs {0,4} + coloffs {0,8..56}
        emit_matmul(L, ctr, "    ", "%Bt", m8, "%U", m64, "%Interp", m64,
                    m=8, k=8, ncols=64, transposed=False)
        L += ["    gpu.return", "  }", "}"]
    elif which == "face":
        # One face of the Knaus apply: interp,deriv,Dm,W,g0,g1,g2 in -> intf out.
        # dt1g/dt2g/flux/tmp are workgroup scratch.
        wsty = _mr(n, WS)
        args = ", ".join("%%%s: %s" % (a, mrg)
                         for a in ["interp", "deriv", "Dm", "W",
                                   "g0", "g1", "g2", "intf"])
        wg = ", ".join("%%%s: %s" % (s, wsty)
                       for s in ["dt1g", "dt2g", "flux", "tmp"])
        L = [module_maps(),
             "gpu.module @mir_kernels [#nvvm.target<chip = \"sm_90\", O = 3>] {",
             "  gpu.func @face(%s) workgroup(%s) kernel {" % (args, wg),
             "    %c0 = arith.constant 0 : index",
             "    %c4 = arith.constant 4 : index",
             "    %z = arith.constant 0.0 : f64",
             "    %lane = gpu.thread_id x"]
        mems = {k: "%" + k for k in
                ["interp", "deriv", "Dm", "W", "g0", "g1", "g2", "intf"]}
        scratch = {}
        for s in ["dt1g", "dt2g", "flux", "tmp"]:
            scratch[s] = "%" + s
            scratch[s + "_ty"] = wsty
        emit_face(L, ctr, "    ", n, mems, scratch)
        L += ["    gpu.return", "  }", "}"]
    else:
        L = [module_maps(),
             "gpu.module @mir_kernels [#nvvm.target<chip = \"sm_90\", O = 3>] {",
             "  gpu.func @two(%%M: %s, %%X: %s, %%Out: %s)" % (mrg, mrg, mrg),
             "      workgroup(%%S1: %s) kernel {" % _mr(n, WS),
             "    %c0 = arith.constant 0 : index",
             "    %c4 = arith.constant 4 : index",
             "    %z = arith.constant 0.0 : f64",
             "    %lane = gpu.thread_id x"]
        s1ty = _mr(n, WS)
        emit_contract(L, ctr, "    ", "%M", mrg, "%X", mrg, "%S1", s1ty, n,
                      transposed=False)
        L.append("    gpu.barrier")
        emit_contract(L, ctr, "    ", "%S1", s1ty, "%M", mrg, "%Out", mrg, n,
                      transposed=True)
        L += ["    gpu.return", "  }", "}"]
    sys.stdout.write("\n".join(L) + "\n")
