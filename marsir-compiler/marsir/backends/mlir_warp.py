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


def _chained_contract_body(L, ctr, J, n, a_mem, a_ty, b_mem, b_ty, transposed):
    """Emit the n/4 chained k-slab vector.contracts INSIDE an already-open warp
    region and return the accumulator SSA (a register vector<8x8>, not written
    to smem). Same math as emit_contract's inner loop -- factored so several
    contracts can share one region and stay register-resident."""
    va, vb_std, vb_t = "vector<8x4xf64>", "vector<4x8xf64>", "vector<8x4xf64>"
    v8 = "vector<%dx%dxf64>" % (n, n)
    acc = ctr.fresh("acc")
    L.append("%s%s = arith.constant dense<0.0> : %s" % (J, acc, v8))
    for s in range(n // 4):
        koff = 4 * s
        a = ctr.fresh("a")
        L.append("%s%s = vector.transfer_read %s[%%c0, %%c%d], %%z "
                 "{in_bounds=[true,true]} : %s, %s" % (J, a, a_mem, koff, a_ty, va))
        b = ctr.fresh("b")
        if transposed:
            L.append("%s%s = vector.transfer_read %s[%%c0, %%c%d], %%z "
                     "{in_bounds=[true,true]} : %s, %s" % (J, b, b_mem, koff, b_ty, vb_t))
            bmap, bty = "#bt", vb_t
        else:
            L.append("%s%s = vector.transfer_read %s[%%c%d, %%c0], %%z "
                     "{in_bounds=[true,true]} : %s, %s" % (J, b, b_mem, koff, b_ty, vb_std))
            bmap, bty = "#b", vb_std
        r = ctr.fresh("r")
        L.append("%s%s = vector.contract {indexing_maps=[#a,%s,#c], "
                 "iterator_types=[\"parallel\",\"parallel\",\"reduction\"], "
                 "kind=#vector.kind<add>} %s, %s, %s : %s, %s into %s"
                 % (J, r, bmap, a, b, acc, va, bty, v8))
        acc = r
    return acc


def emit_fused_dt_flux(L, ctr, indent, n, interp_mem, interp_ty, Dm_mem, Dm_ty,
                       deriv_mem, deriv_ty, g0_mem, g1_mem, g2_mem, g_ty,
                       flux_mem, flux_ty):
    """WIN2+: dt1, dt2 and the 3-term flux fused into ONE warp region so the
    dt1g/dt2g intermediates stay REGISTER-resident (never touch smem). Each
    contract result feeds exactly one mulf (legal); the flux is a relayout-free
    C-fragment junction. Drops 2 of 7 smem buffers + 1 barrier per face.
       dt1 = Dm @ interp (std) *g1 ; dt2 = interp @ Dm^T (transp) *g0
       flux = g2*deriv + dt2g + dt1g   -> flux_mem"""
    v8 = "vector<%dx%dxf64>" % (n, n)
    I, J = indent, indent + "  "
    L.append("%svector.warp_execute_on_lane_0(%%lane)[32] {" % I)
    acc1 = _chained_contract_body(L, ctr, J, n, Dm_mem, Dm_ty, interp_mem,
                                  interp_ty, transposed=False)
    g1v = ctr.fresh("g")
    L.append("%s%s = vector.transfer_read %s[%%c0, %%c0], %%z "
             "{in_bounds=[true,true]} : %s, %s" % (J, g1v, g1_mem, g_ty, v8))
    dt1g = ctr.fresh("dt")
    L.append("%s%s = arith.mulf %s, %s : %s" % (J, dt1g, acc1, g1v, v8))
    acc2 = _chained_contract_body(L, ctr, J, n, interp_mem, interp_ty, Dm_mem,
                                  Dm_ty, transposed=True)
    g0v = ctr.fresh("g")
    L.append("%s%s = vector.transfer_read %s[%%c0, %%c0], %%z "
             "{in_bounds=[true,true]} : %s, %s" % (J, g0v, g0_mem, g_ty, v8))
    dt2g = ctr.fresh("dt")
    L.append("%s%s = arith.mulf %s, %s : %s" % (J, dt2g, acc2, g0v, v8))
    dv = ctr.fresh("dv")
    L.append("%s%s = vector.transfer_read %s[%%c0, %%c0], %%z "
             "{in_bounds=[true,true]} : %s, %s" % (J, dv, deriv_mem, deriv_ty, v8))
    g2v = ctr.fresh("g")
    L.append("%s%s = vector.transfer_read %s[%%c0, %%c0], %%z "
             "{in_bounds=[true,true]} : %s, %s" % (J, g2v, g2_mem, g_ty, v8))
    t = ctr.fresh("t")
    L.append("%s%s = arith.mulf %s, %s : %s" % (J, t, dv, g2v, v8))
    t2 = ctr.fresh("t2")
    L.append("%s%s = arith.addf %s, %s : %s" % (J, t2, t, dt2g, v8))
    fl = ctr.fresh("fl")
    L.append("%s%s = arith.addf %s, %s : %s" % (J, fl, t2, dt1g, v8))
    L.append("%svector.transfer_write %s, %s[%%c0, %%c0] "
             "{in_bounds=[true,true]} : %s, %s" % (J, fl, flux_mem, v8, flux_ty))
    L.append("%s}" % I)


def emit_face(L, ctr, indent, n, mems, scratch, tys=None):
    """One face of the Knaus apply: the four contracts + the 3-term flux, all
    staged through shared memory, writing intf to mems['intf']. Matches the
    oracle:
       dt2  = interp @ Dm^T   (transposed) ; scaled by g0 -> dt2g
       dt1  = Dm @ interp     (standard)   ; scaled by g1 -> dt1g
       flux = g2*deriv + dt2g + dt1g
       tmp  = flux @ W^T      (transposed)
       intf = W @ tmp         (standard)
    `mems`   : dict of input/output nxn memref SSA names
               {interp,deriv,Dm,W,g0,g1,g2,intf}
    `tys`    : dict of the ACTUAL memref TYPE strings for each mems key (plane
               subviews carry a strided + workgroup layout, not plain nxn); any
               key absent defaults to a plain nxn memref.
    `scratch`: dict of workgroup nxn scratch memref SSA names {dt1g,dt2g,flux,
               tmp} + their type strings under key+'_ty'."""
    mr = _mr(n)
    tys = tys or {}
    def mt(k):
        return tys.get(k, mr)
    def st(k):
        return scratch.get(k + "_ty", mr)
    # WIN2+: dt1 + dt2 + flux fused into one region (dt1g/dt2g register-resident)
    emit_fused_dt_flux(L, ctr, indent, n, mems["interp"], mt("interp"),
                       mems["Dm"], mt("Dm"), mems["deriv"], mt("deriv"),
                       mems["g0"], mems["g1"], mems["g2"], mt("g0"),
                       scratch["flux"], st("flux"))
    L.append("%sgpu.barrier" % indent)
    # tmp = flux @ W^T (transposed)
    emit_contract(L, ctr, indent, scratch["flux"], st("flux"), mems["W"],
                  mt("W"), scratch["tmp"], st("tmp"), n, transposed=True)
    L.append("%sgpu.barrier" % indent)
    # intf = W @ tmp (standard)
    emit_contract(L, ctr, indent, mems["W"], mt("W"), scratch["tmp"], st("tmp"),
                  mems["intf"], mt("intf"), n, transposed=False)


def _plane(space):
    return (", " + space) if space else ""


def reinterpret_2d_to_3d(L, ctr, indent, src, n, space=None):
    """View a (P-padded 8) x n^2 buffer as (8, n, n): row l of the matmul output
    is face l flattened (s*n+r), so the 3-D view's [l, s, r] is contiguous with
    the 2-D [l, s*n+r]. Returns the new 3-D SSA name + its type."""
    dst = ctr.fresh("v3")
    ty = "memref<8x%dx%dxf64%s>" % (n, n, _plane(space))
    L.append("%s%s = memref.reinterpret_cast %s to offset: [0], "
             "sizes: [8, %d, %d], strides: [%d, %d, 1] : memref<8x%dxf64%s> to %s"
             % (indent, dst, src, n, n, n * n, n, n * n, _plane(space), ty))
    return dst, ty


def subview_plane(L, ctr, indent, src, n, l, axis, space=None, src_ty=None,
                  dyn=False):
    """Rank-reduced n x n plane l of an 8xnxn buffer along `axis` (0/1/2). The
    +/- plane scatter slices y on a different axis per direction; the B-sweep
    face slice is always axis 0. When the source is a per-element subview of a
    batched buffer its OFFSET is dynamic (pass dyn=True + src_ty); the plane's
    offset is then also dynamic. Strides are unchanged either way."""
    off = [0, 0, 0]; off[axis] = l
    sz = [8, n, n]; sz[axis] = 1
    strides = [n * n, n, 1]
    keep = [strides[i] for i in range(3) if i != axis]
    r = ctr.fresh("sv")
    boff = "?" if dyn else str(l * strides[axis])
    lay = "strided<[%d, %d], offset: %s>" % (keep[0], keep[1], boff)
    dst_ty = "memref<%dx%dxf64, %s%s>" % (n, n, lay, _plane(space))
    if src_ty is None:
        src_ty = "memref<8x%dx%dxf64%s>" % (n, n, _plane(space))
    L.append("%s%s = memref.subview %s[%d, %d, %d] [%d, %d, %d] [1, 1, 1] : "
             "%s to %s"
             % (indent, r, src, off[0], off[1], off[2], sz[0], sz[1], sz[2],
                src_ty, dst_ty))
    return r, dst_ty


def subview_elem(L, ctr, indent, src, e, n, space=None):
    """Per-element rank-4 -> rank-3 subview src[e,:,:,:] of a batched Ex8xnxn
    buffer; the offset is dynamic (e is a runtime block id)."""
    r = ctr.fresh("el")
    # rank-reduced 8xnxn result -> 3 strides (the leading E dim is dropped)
    strd = "strided<[%d, %d, 1], offset: ?>" % (n * n, n)
    src_ty = "memref<?x8x%dx%dxf64%s>" % (n, n, _plane(space))
    dst_ty = "memref<8x%dx%dxf64, %s%s>" % (n, n, strd, _plane(space))
    L.append("%s%s = memref.subview %s[%s, 0, 0, 0] [1, 8, %d, %d] "
             "[1, 1, 1, 1] : %s to %s"
             % (indent, r, src, e, n, n, src_ty, dst_ty))
    return r, dst_ty


def subview_elem_2d(L, ctr, indent, src, e, n, space=None):
    """Per-element rank-3 -> rank-2 subview src[e,:,:] of a batched Exnxn buffer
    (the per-element metric planes). Dynamic offset."""
    r = ctr.fresh("el")
    strd = "strided<[%d, 1], offset: ?>" % n
    src_ty = "memref<?x%dx%dxf64%s>" % (n, n, _plane(space))
    dst_ty = "memref<%dx%dxf64, %s%s>" % (n, n, strd, _plane(space))
    L.append("%s%s = memref.subview %s[%s, 0, 0] [1, %d, %d] [1, 1, 1] : %s to %s"
             % (indent, r, src, e, n, n, src_ty, dst_ty))
    return r, dst_ty


def emit_scatter(L, ctr, indent, n, intf_mem, intf_ty, y3, l, axis,
                 y3_ty=None, dyn=False):
    """y[l] -= intf ; y[l+1] += intf. intf is read from scratch ONCE and reused
    for both planes -- the contract-result-double-use crash is specific to a
    CONTRACT result feeding two elementwise ops; a transfer_read result feeds
    both subf and addf fine (WIN1: one fewer smem load per face). y3 plane l is
    sliced along `axis` (0/1/2 per direction); batched y3 -> pass type + dyn."""
    v = "vector<%dx%dxf64>" % (n, n)
    ym, yty = subview_plane(L, ctr, indent, y3, n, l, axis, src_ty=y3_ty, dyn=dyn)
    yp, ypty = subview_plane(L, ctr, indent, y3, n, l + 1, axis, src_ty=y3_ty,
                             dyn=dyn)
    I, J = indent, indent + "  "
    L.append("%svector.warp_execute_on_lane_0(%%lane)[32] {" % I)
    fa = ctr.fresh("fa")
    L.append("%s%s = vector.transfer_read %s[%%c0, %%c0], %%z "
             "{in_bounds=[true,true]} : %s, %s" % (J, fa, intf_mem, intf_ty, v))
    r0 = ctr.fresh("r"); L.append("%s%s = vector.transfer_read %s[%%c0, %%c0], "
             "%%z {in_bounds=[true,true]} : %s, %s" % (J, r0, ym, yty, v))
    m0 = ctr.fresh("m"); L.append("%s%s = arith.subf %s, %s : %s" % (J, m0, r0, fa, v))
    L.append("%svector.transfer_write %s, %s[%%c0, %%c0] {in_bounds=[true,true]} "
             ": %s, %s" % (J, m0, ym, v, yty))
    r1 = ctr.fresh("r"); L.append("%s%s = vector.transfer_read %s[%%c0, %%c0], "
             "%%z {in_bounds=[true,true]} : %s, %s" % (J, r1, yp, ypty, v))
    p1 = ctr.fresh("p"); L.append("%s%s = arith.addf %s, %s : %s" % (J, p1, r1, fa, v))
    L.append("%svector.transfer_write %s, %s[%%c0, %%c0] {in_bounds=[true,true]} "
             ": %s, %s" % (J, p1, yp, v, ypty))
    L.append("%s}" % I)


def emit_bsweep(L, ctr, indent, n, u_mem, u_ty, op_mem, out_mem, out_ty,
                pax, transposed, dyn=False):
    """B-sweep out = op @ (u viewed for this direction), reading u PLANES
    directly -- no transpose buffer. Col-tile ct = interp_all[:, ct*8:ct*8+8]
    comes from plane ct of u along axis `pax`; the direction's moveaxis is
    absorbed into (pax, arm):
       dir0: pax=1 (u[:,ct,:]), standard    out[l,r]=sum_q op[l,q]*u[q,ct,r]
       dir1: pax=0 (u[ct,:,:]), standard    out[l,r]=sum_q op[l,q]*u[ct,q,r]
       dir2: pax=0 (u[ct,:,:]), transposed  out[l,r]=sum_q op[l,q]*u[ct,r,q]
    op is Btil or Dtil (nxn). u is 8xnxn. out is 8xn^2 (col-tiles of 8)."""
    mr = _mr(n)
    va = "vector<8x4xf64>"
    vb_std = "vector<4x8xf64>"
    vb_t = "vector<8x4xf64>"
    v8 = "vector<8x8xf64>"
    slabs = n // 4
    for ct in range(n):
        plane, pty = subview_plane(L, ctr, indent, u_mem, n, ct, pax,
                                   src_ty=u_ty, dyn=dyn)
        L.append("%svector.warp_execute_on_lane_0(%%lane)[32] {" % indent)
        J = indent + "  "
        acc = ctr.fresh("acc")
        L.append("%s%s = arith.constant dense<0.0> : %s" % (J, acc, v8))
        for s in range(slabs):
            koff = 4 * s
            a = ctr.fresh("a")
            L.append("%s%s = vector.transfer_read %s[%%c0, %%c%d], %%z "
                     "{in_bounds=[true,true]} : %s, %s" % (J, a, op_mem, koff, mr, va))
            b = ctr.fresh("b")
            if transposed:
                L.append("%s%s = vector.transfer_read %s[%%c0, %%c%d], %%z "
                         "{in_bounds=[true,true]} : %s, %s" % (J, b, plane, koff, pty, vb_t))
                bmap, bty = "#bt", vb_t
            else:
                L.append("%s%s = vector.transfer_read %s[%%c%d, %%c0], %%z "
                         "{in_bounds=[true,true]} : %s, %s" % (J, b, plane, koff, pty, vb_std))
                bmap, bty = "#b", vb_std
            r = ctr.fresh("r")
            L.append("%s%s = vector.contract {indexing_maps=[#a,%s,#c], "
                     "iterator_types=[\"parallel\",\"parallel\",\"reduction\"], "
                     "kind=#vector.kind<add>} %s, %s, %s : %s, %s into %s"
                     % (J, r, bmap, a, b, acc, va, bty, v8))
            acc = r
        L.append("%svector.transfer_write %s, %s[%%c0, %%c%d] "
                 "{in_bounds=[true,true]} : %s, %s"
                 % (J, acc, out_mem, 8 * ct, v8, out_ty))
        L.append("%s}" % indent)


def emit_direction_u(L, ctr, indent, n, P, mems, scratch, g_keys,
                     scatter_axis, pax, transposed, dyn=False, u_ty=None,
                     g_ty=None, y_ty=None, metric_per_face=True):
    """Direction with the transpose INTERNALIZED: B-sweeps read u planes via
    emit_bsweep (single u input, no host-presented u_flat). For a batched kernel
    u/metric/y are per-element subviews with dynamic offsets -- pass dyn=True and
    their strided types (u_ty, g_ty, y_ty)."""
    emit_bsweep(L, ctr, indent, n, mems["u"], u_ty, mems["Btil"],
                scratch["interp_all"], scratch["interp_all_ty"], pax, transposed,
                dyn=dyn)
    emit_bsweep(L, ctr, indent, n, mems["u"], u_ty, mems["Dtil"],
                scratch["deriv_all"], scratch["deriv_all_ty"], pax, transposed,
                dyn=dyn)
    L.append("%sgpu.barrier" % indent)
    ip3, _ = reinterpret_2d_to_3d(L, ctr, indent, scratch["interp_all"], n, WS)
    dv3, _ = reinterpret_2d_to_3d(L, ctr, indent, scratch["deriv_all"], n, WS)
    g0k, g1k, g2k = g_keys
    for l in range(P):
        # interp/deriv planes come from workgroup scratch (static, not dyn).
        ifc, ifc_ty = subview_plane(L, ctr, indent, ip3, n, l, 0, WS)
        dfc, dfc_ty = subview_plane(L, ctr, indent, dv3, n, l, 0, WS)
        if metric_per_face:
            g0, g0_ty = subview_plane(L, ctr, indent, mems[g0k], n, l, 0,
                                      src_ty=g_ty, dyn=dyn)
            g1, g1_ty = subview_plane(L, ctr, indent, mems[g1k], n, l, 0,
                                      src_ty=g_ty, dyn=dyn)
            g2, g2_ty = subview_plane(L, ctr, indent, mems[g2k], n, l, 0,
                                      src_ty=g_ty, dyn=dyn)
        else:
            # per-element metric (nxn), shared across the P faces -> read once
            # per direction; ~8x less metric traffic (the affine/cheap-metric
            # path that lifts arithmetic intensity above the FP64 ridge).
            g0, g0_ty = mems[g0k], (g_ty or _mr(n))
            g1, g1_ty = mems[g1k], (g_ty or _mr(n))
            g2, g2_ty = mems[g2k], (g_ty or _mr(n))
        fmems = {"interp": ifc, "deriv": dfc, "Dm": mems["Dm"], "W": mems["W"],
                 "g0": g0, "g1": g1, "g2": g2, "intf": scratch["intf"]}
        ftys = {"interp": ifc_ty, "deriv": dfc_ty, "g0": g0_ty, "g1": g1_ty,
                "g2": g2_ty, "intf": scratch["intf_ty"]}
        emit_face(L, ctr, indent, n, fmems, scratch, tys=ftys)
        L.append("%sgpu.barrier" % indent)
        emit_scatter(L, ctr, indent, n, scratch["intf"], scratch["intf_ty"],
                     mems["y3"], l, scatter_axis, y3_ty=y_ty, dyn=dyn)
        L.append("%sgpu.barrier" % indent)


def emit_direction(L, ctr, indent, n, P, mems, scratch, u2_key, g_keys,
                   scatter_axis):
    """One direction of the Knaus apply. B-sweep interp/deriv = Btil/Dtil @
    u_flat (mems[u2_key], already in this direction's layout), then per face l
    apply emit_face and scatter intf into y along `scatter_axis` (0/1/2 for the
    three directions -- moveaxis(y,d,0)). The face slice of the B-sweep output is
    ALWAYS axis 0 (the face dimension). g_keys = (g0,g1,g2) mems keys for this
    direction's per-face metric. interp_all/deriv_all + per-face buffers are
    reused across directions (barriers order them)."""
    m64 = "memref<8x%dxf64>" % (n * n)
    wg64 = "memref<8x%dxf64, %s>" % (n * n, WS)
    emit_matmul(L, ctr, indent, mems["Btil"], _mr(n), mems[u2_key], m64,
                scratch["interp_all"], wg64, m=8, k=n, ncols=n * n,
                transposed=False)
    emit_matmul(L, ctr, indent, mems["Dtil"], _mr(n), mems[u2_key], m64,
                scratch["deriv_all"], wg64, m=8, k=n, ncols=n * n,
                transposed=False)
    L.append("%sgpu.barrier" % indent)
    ip3, _ = reinterpret_2d_to_3d(L, ctr, indent, scratch["interp_all"], n, WS)
    dv3, _ = reinterpret_2d_to_3d(L, ctr, indent, scratch["deriv_all"], n, WS)
    g0k, g1k, g2k = g_keys
    for l in range(P):
        ifc, ifc_ty = subview_plane(L, ctr, indent, ip3, n, l, 0, WS)
        dfc, dfc_ty = subview_plane(L, ctr, indent, dv3, n, l, 0, WS)
        g0, g0_ty = subview_plane(L, ctr, indent, mems[g0k], n, l, 0)
        g1, g1_ty = subview_plane(L, ctr, indent, mems[g1k], n, l, 0)
        g2, g2_ty = subview_plane(L, ctr, indent, mems[g2k], n, l, 0)
        fmems = {"interp": ifc, "deriv": dfc, "Dm": mems["Dm"], "W": mems["W"],
                 "g0": g0, "g1": g1, "g2": g2, "intf": scratch["intf"]}
        ftys = {"interp": ifc_ty, "deriv": dfc_ty, "g0": g0_ty, "g1": g1_ty,
                "g2": g2_ty, "intf": scratch["intf_ty"]}
        emit_face(L, ctr, indent, n, fmems, scratch, tys=ftys)
        L.append("%sgpu.barrier" % indent)
        emit_scatter(L, ctr, indent, n, scratch["intf"], scratch["intf_ty"],
                     mems["y3"], l, scatter_axis)
        L.append("%sgpu.barrier" % indent)


def emit_direction0(L, ctr, indent, n, P, mems, scratch):
    """Direction 0 (no input transpose, scatter axis 0). Thin wrapper kept for
    the single-direction gate."""
    emit_direction(L, ctr, indent, n, P, mems, scratch, "u2",
                   ("g0all", "g1all", "g2all"), scatter_axis=0)


def build_dir0_kernel(n, P):
    """Full single-direction (dir 0, no input transpose) operator kernel string:
    B-sweep interp/deriv = Btil/Dtil @ u_flat, then the P faces (emit_face +
    scatter into y). Inputs u2 (8xn^2), Btil/Dtil/Dm/W (nxn), g0all/g1all/g2all
    (8xnxn per-face metric components), out y3 (8xnxn). interp_all/deriv_all +
    per-face buffers are workgroup (shared) scratch."""
    ctr = Ctr()
    mrg = _mr(n)
    m64 = "memref<8x%dxf64>" % (n * n)
    wg64 = "memref<8x%dxf64, %s>" % (n * n, WS)
    wsty = _mr(n, WS)
    gall = "memref<8x%dx%dxf64>" % (n, n)
    inargs = [("u2", m64), ("Btil", mrg), ("Dtil", mrg), ("Dm", mrg),
              ("W", mrg), ("g0all", gall), ("g1all", gall), ("g2all", gall),
              ("y3", gall)]
    argstr = ", ".join("%%%s: %s" % (a, t) for a, t in inargs)
    scratch_names = ["interp_all", "deriv_all", "flux", "tmp", "intf"]
    wg = ", ".join("%%%s: %s" % (s, wg64 if s in ("interp_all", "deriv_all")
                                 else wsty) for s in scratch_names)
    L = [module_maps(),
         "gpu.module @mir_kernels [#nvvm.target<chip = \"sm_90\", O = 3>] {",
         "  gpu.func @dir0(%s) workgroup(%s) kernel {" % (argstr, wg),
         "    %z = arith.constant 0.0 : f64",
         "    %lane = gpu.thread_id x"]
    L += index_consts("    ", n * n)
    mems = {k: "%" + k for k in
            ["u2", "Btil", "Dtil", "Dm", "W", "g0all", "g1all", "g2all", "y3"]}
    scratch = {}
    for s in scratch_names:
        scratch[s] = "%" + s
        scratch[s + "_ty"] = wg64 if s in ("interp_all", "deriv_all") else wsty
    emit_direction0(L, ctr, "    ", n, P, mems, scratch)
    L += ["    gpu.return", "  }", "}"]
    return "\n".join(L) + "\n"


def build_full_kernel(n, P):
    """Full 3-direction operator (correctness form): the input is presented in
    the three moveaxis layouts u0/u1/u2 (each 8xn^2) so the transpose is out of
    the kernel for now; each direction runs emit_direction and scatters into the
    SAME y along its own axis (0/1/2), accumulating. Per-direction metric is 9
    arrays g{d}{c}all (8xnxn). Internalizing the transpose is the next step."""
    ctr = Ctr()
    mrg = _mr(n)
    m64 = "memref<8x%dxf64>" % (n * n)
    wg64 = "memref<8x%dxf64, %s>" % (n * n, WS)
    wsty = _mr(n, WS)
    gall = "memref<8x%dx%dxf64>" % (n, n)
    uflat = ["u0", "u1", "u2"]
    gnames = ["g%d%dall" % (d, c) for d in range(3) for c in range(3)]
    inargs = [(u, m64) for u in uflat]
    inargs += [("Btil", mrg), ("Dtil", mrg), ("Dm", mrg), ("W", mrg)]
    inargs += [(g, gall) for g in gnames]
    inargs += [("y3", gall)]
    argstr = ", ".join("%%%s: %s" % (a, t) for a, t in inargs)
    scratch_names = ["interp_all", "deriv_all", "flux", "tmp", "intf"]
    wg = ", ".join("%%%s: %s" % (s, wg64 if s in ("interp_all", "deriv_all")
                                 else wsty) for s in scratch_names)
    L = [module_maps(),
         "gpu.module @mir_kernels [#nvvm.target<chip = \"sm_90\", O = 3>] {",
         "  gpu.func @full(%s) workgroup(%s) kernel {" % (argstr, wg),
         "    %z = arith.constant 0.0 : f64",
         "    %lane = gpu.thread_id x"]
    L += index_consts("    ", n * n)
    mems = {k: "%" + k for k, _ in inargs}
    scratch = {}
    for s in scratch_names:
        scratch[s] = "%" + s
        scratch[s + "_ty"] = wg64 if s in ("interp_all", "deriv_all") else wsty
    for d in range(3):
        gk = ("g%d0all" % d, "g%d1all" % d, "g%d2all" % d)
        emit_direction(L, ctr, "    ", n, P, mems, scratch, uflat[d], gk,
                       scatter_axis=d)
    L += ["    gpu.return", "  }", "}"]
    return "\n".join(L) + "\n"


def build_full_single(n, P):
    """Full 3-direction operator with the transpose INTERNALIZED: a single u
    (8xnxn) input; each direction reads u planes via emit_bsweep (pax/arm absorb
    the moveaxis) and scatters into y along its axis. This is the real
    single-input operator -- the head-to-head measurement target."""
    ctr = Ctr()
    mrg = _mr(n)
    wg64 = "memref<8x%dxf64, %s>" % (n * n, WS)
    wsty = _mr(n, WS)
    gall = "memref<8x%dx%dxf64>" % (n, n)
    u3 = "memref<8x%dx%dxf64>" % (n, n)
    gnames = ["g%d%dall" % (d, c) for d in range(3) for c in range(3)]
    inargs = [("u", u3), ("Btil", mrg), ("Dtil", mrg), ("Dm", mrg), ("W", mrg)]
    inargs += [(g, gall) for g in gnames]
    inargs += [("y3", gall)]
    argstr = ", ".join("%%%s: %s" % (a, t) for a, t in inargs)
    scratch_names = ["interp_all", "deriv_all", "flux", "tmp", "intf"]
    wg = ", ".join("%%%s: %s" % (s, wg64 if s in ("interp_all", "deriv_all")
                                 else wsty) for s in scratch_names)
    L = [module_maps(),
         "gpu.module @mir_kernels [#nvvm.target<chip = \"sm_90\", O = 3>] {",
         "  gpu.func @full(%s) workgroup(%s) kernel {" % (argstr, wg),
         "    %z = arith.constant 0.0 : f64",
         "    %lane = gpu.thread_id x"]
    L += index_consts("    ", n * n)
    mems = {k: "%" + k for k, _ in inargs}
    scratch = {}
    for s in scratch_names:
        scratch[s] = "%" + s
        scratch[s + "_ty"] = wg64 if s in ("interp_all", "deriv_all") else wsty
    # (pax, transposed, scatter_axis) per direction
    dirs = [(1, False, 0), (0, False, 1), (0, True, 2)]
    for d, (pax, transp, sax) in enumerate(dirs):
        gk = ("g%d0all" % d, "g%d1all" % d, "g%d2all" % d)
        emit_direction_u(L, ctr, "    ", n, P, mems, scratch, gk, sax, pax, transp)
    L += ["    gpu.return", "  }", "}"]
    return "\n".join(L) + "\n"


def build_full_batched(n, P):
    """Batched operator: one WARP per element over E elements (grid = E blocks,
    block = 32). e = gpu.block_id x selects the element; U/metric/Y are batched
    (Ex8xnxn) and per-element subviews (dynamic offset) feed the same body. The
    operator matrices Btil/Dtil/Dm/W are shared. This is the throughput-
    measurement kernel (launch grid=(E,1,1), block=(32,1,1))."""
    ctr = Ctr()
    mrg = _mr(n)
    wg64 = "memref<8x%dxf64, %s>" % (n * n, WS)
    wsty = _mr(n, WS)
    batched = "memref<?x8x%dx%dxf64>" % (n, n)
    gnames = ["g%d%dall" % (d, c) for d in range(3) for c in range(3)]
    inargs = [("U", batched), ("Btil", mrg), ("Dtil", mrg), ("Dm", mrg),
              ("W", mrg)]
    inargs += [(g, batched) for g in gnames]
    inargs += [("Y", batched)]
    argstr = ", ".join("%%%s: %s" % (a, t) for a, t in inargs)
    scratch_names = ["interp_all", "deriv_all", "flux", "tmp", "intf"]
    wg = ", ".join("%%%s: %s" % (s, wg64 if s in ("interp_all", "deriv_all")
                                 else wsty) for s in scratch_names)
    L = [module_maps(),
         "gpu.module @mir_kernels [#nvvm.target<chip = \"sm_90\", O = 3>] {",
         "  gpu.func @full_batched(%s) workgroup(%s) kernel {" % (argstr, wg),
         "    %z = arith.constant 0.0 : f64",
         "    %lane = gpu.thread_id x",
         "    %e = gpu.block_id x"]
    L += index_consts("    ", n * n)
    ue, elem_ty = subview_elem(L, ctr, "    ", "%U", "%e", n)
    ye, _ = subview_elem(L, ctr, "    ", "%Y", "%e", n)
    mems = {"u": ue, "Btil": "%Btil", "Dtil": "%Dtil", "Dm": "%Dm", "W": "%W",
            "y3": ye}
    for g in gnames:
        ge, _ = subview_elem(L, ctr, "    ", "%" + g, "%e", n)
        mems[g] = ge
    scratch = {}
    for s in scratch_names:
        scratch[s] = "%" + s
        scratch[s + "_ty"] = wg64 if s in ("interp_all", "deriv_all") else wsty
    dirs = [(1, False, 0), (0, False, 1), (0, True, 2)]
    for d, (pax, transp, sax) in enumerate(dirs):
        gk = ("g%d0all" % d, "g%d1all" % d, "g%d2all" % d)
        emit_direction_u(L, ctr, "    ", n, P, mems, scratch, gk, sax, pax,
                         transp, dyn=True, u_ty=elem_ty, g_ty=elem_ty,
                         y_ty=elem_ty)
    L += ["    gpu.return", "  }", "}"]
    return "\n".join(L) + "\n"


def build_full_batched_affine(n, P):
    """Batched operator with a CHEAP per-element metric: g{d}{c} are Exnxn (one
    nxn per element, shared across the P faces) instead of Ex8xnxn per-face.
    Cuts metric traffic ~8x -> arithmetic intensity above the FP64 ridge
    (compute-bound regime). Same warp-per-element structure otherwise; this is
    the Track-A speedup AND the regime where tensor cores could matter."""
    ctr = Ctr()
    mrg = _mr(n)
    wg64 = "memref<8x%dxf64, %s>" % (n * n, WS)
    wsty = _mr(n, WS)
    batched = "memref<?x8x%dx%dxf64>" % (n, n)
    metric = "memref<?x%dx%dxf64>" % (n, n)          # Exnxn per-element metric
    gnames = ["g%d%dall" % (d, c) for d in range(3) for c in range(3)]
    inargs = [("U", batched), ("Btil", mrg), ("Dtil", mrg), ("Dm", mrg),
              ("W", mrg)]
    inargs += [(g, metric) for g in gnames]
    inargs += [("Y", batched)]
    argstr = ", ".join("%%%s: %s" % (a, t) for a, t in inargs)
    scratch_names = ["interp_all", "deriv_all", "flux", "tmp", "intf"]
    wg = ", ".join("%%%s: %s" % (s, wg64 if s in ("interp_all", "deriv_all")
                                 else wsty) for s in scratch_names)
    L = [module_maps(),
         "gpu.module @mir_kernels [#nvvm.target<chip = \"sm_90\", O = 3>] {",
         "  gpu.func @full_batched_affine(%s) workgroup(%s) kernel {"
         % (argstr, wg),
         "    %z = arith.constant 0.0 : f64",
         "    %lane = gpu.thread_id x",
         "    %e = gpu.block_id x"]
    L += index_consts("    ", n * n)
    ue, elem_ty = subview_elem(L, ctr, "    ", "%U", "%e", n)
    ye, _ = subview_elem(L, ctr, "    ", "%Y", "%e", n)
    mems = {"u": ue, "Btil": "%Btil", "Dtil": "%Dtil", "Dm": "%Dm", "W": "%W",
            "y3": ye}
    g_ty = None
    for g in gnames:
        ge, g_ty = subview_elem_2d(L, ctr, "    ", "%" + g, "%e", n)
        mems[g] = ge
    scratch = {}
    for s in scratch_names:
        scratch[s] = "%" + s
        scratch[s + "_ty"] = wg64 if s in ("interp_all", "deriv_all") else wsty
    dirs = [(1, False, 0), (0, False, 1), (0, True, 2)]
    for d, (pax, transp, sax) in enumerate(dirs):
        gk = ("g%d0all" % d, "g%d1all" % d, "g%d2all" % d)
        emit_direction_u(L, ctr, "    ", n, P, mems, scratch, gk, sax, pax,
                         transp, dyn=True, u_ty=elem_ty, g_ty=g_ty,
                         y_ty=elem_ty, metric_per_face=False)
    L += ["    gpu.return", "  }", "}"]
    return "\n".join(L) + "\n"


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
    if which == "dir0":
        sys.stdout.write(build_dir0_kernel(n, 7))
        sys.exit(0)
    if which == "full":
        sys.stdout.write(build_full_kernel(n, 7))
        sys.exit(0)
    if which == "fulls":
        sys.stdout.write(build_full_single(n, 7))
        sys.exit(0)
    if which == "batched":
        sys.stdout.write(build_full_batched(n, 7))
        sys.exit(0)
    if which == "batched_affine":
        sys.stdout.write(build_full_batched_affine(n, 7))
        sys.exit(0)
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
                       for s in ["flux", "tmp"])
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
        for s in ["flux", "tmp"]:
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
