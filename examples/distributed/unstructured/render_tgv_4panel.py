#!/usr/bin/env pvbatch
# Render a 4-panel ParaView animation of the TGV simulation.
#
# Inputs (all derived from --pvd PREFIX):
#   PREFIX.pvd  -- multi-field PVD timeline written by VTUParallelWriter::
#                  writeMultiFieldFrame. Must contain:
#                    - PointData scalar  "u"          (x-velocity)
#                    - PointData scalar  "p"          (pressure)
#                    - PointData scalar  "omega"      (vorticity magnitude)
#                    - PointData vector  "velocity"   (3-component velocity)
#                    - CellData   scalar "level"      (refinement level)
#
# Output:
#   <PREFIX>_4panel/frame_NNNN.png    -- per-timestep frames
#   <PREFIX>_4panel.mp4               -- assembled animation (if ffmpeg in PATH)
#
# Usage on Mac:
#   /Applications/ParaView-5.13.2.app/Contents/bin/pvbatch \
#       render_tgv_4panel.py --pvd /tmp/tgv/run --width 1920 --height 1080 \
#       --fps 30
#
# Usage on Alps (offscreen pvbatch from a Spack module env):
#   pvbatch render_tgv_4panel.py --pvd /tmp/tgv/run
#
# Layout (clockwise from top-left, MFlex/Mavriplis-style):
#   TL: mesh outline + refinement level color (CellData)
#   TR: refinement level isosurface (volume cells colored by level)
#   BL: vorticity magnitude (volume rendering or isosurface)
#   BR: u-velocity slice + surface pressure or pressure isosurface

import argparse, os, sys, glob, math

ap = argparse.ArgumentParser()
ap.add_argument("--pvd",     required=True, help="prefix of the .pvd file (without extension)")
ap.add_argument("--width",   type=int, default=1920)
ap.add_argument("--height",  type=int, default=1080)
ap.add_argument("--fps",     type=int, default=30)
ap.add_argument("--frames",  type=int, default=0, help="0 = all timesteps in PVD")
ap.add_argument("--no-mp4",  action="store_true", help="skip ffmpeg assembly")
ap.add_argument("--vort-range", default="0,30",
                help="vorticity magnitude color range LOW,HI")
ap.add_argument("--p-range",    default="-5,5",
                help="pressure color range LOW,HI")
ap.add_argument("--u-range",    default="-1,1",
                help="u-velocity color range LOW,HI")
ap.add_argument("--vort-iso",   type=float, default=5.0,
                help="vorticity-magnitude isosurface value for the BL panel")
args = ap.parse_args()

pvd_path = args.pvd if args.pvd.endswith(".pvd") else (args.pvd + ".pvd")
if not os.path.exists(pvd_path):
    print(f"PVD not found: {pvd_path}", file=sys.stderr)
    sys.exit(1)

prefix    = pvd_path[:-4]
out_dir   = prefix + "_4panel"
os.makedirs(out_dir, exist_ok=True)

vlo, vhi = [float(x) for x in args.vort_range.split(",")]
plo, phi = [float(x) for x in args.p_range.split(",")]
ulo, uhi = [float(x) for x in args.u_range.split(",")]

# --- ParaView setup ----------------------------------------------------------
from paraview.simple import *

paraview.simple._DisableFirstRenderCameraReset()

reader = PVDReader(FileName=pvd_path)
reader.UpdatePipeline()
all_times = list(reader.TimestepValues or [])
if args.frames > 0:
    all_times = all_times[:args.frames]
print(f"PVD: {pvd_path}  timesteps: {len(all_times)}")

# Bounds
reader.UpdatePipeline(time=all_times[0])
bounds = reader.GetDataInformation().GetBounds()
cx = 0.5 * (bounds[0] + bounds[1])
cy = 0.5 * (bounds[2] + bounds[3])
cz = 0.5 * (bounds[4] + bounds[5])
L  = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
print(f"bounds: x=[{bounds[0]:.3f},{bounds[1]:.3f}]  y=[{bounds[2]:.3f},{bounds[3]:.3f}]  z=[{bounds[4]:.3f},{bounds[5]:.3f}]")

# --- LUTs ---------------------------------------------------------------------
def lut(name, low, high, preset="Cool to Warm"):
    L = GetColorTransferFunction(name)
    L.RescaleTransferFunction(low, high)
    L.ApplyPreset(preset, True)
    return L

# Vorticity uses "Cool to Warm (Extended)" -- the ONERA TGV video uses a
# rainbow blue-cyan-green-yellow LUT; the extended Cool-to-Warm preset matches
# that look while staying perceptually monotonic. Bright across the entire
# range, so the iso surface stays visible as the field decays.
vortLUT = lut("omega",    vlo, vhi, "Blue Orange (divergent)")
pLUT    = lut("p",        plo, phi, "Blue Orange (divergent)")
uLUT    = lut("u",        ulo, uhi, "Blue Orange (divergent)")
# |velocity| LUT. "Cool to Warm" (non-extended) maps low to a *medium* blue
# (#3b4cc0), not dark navy / black. Important because the TGV initial field
# has |velocity| = 0 at the cube corners and 8 isolated interior points;
# Turbo and the *Extended* variant of Cool-to-Warm both bottom out near
# black there, which made the AMR-refined corners look bad. The plain
# "Cool to Warm" preset stays in a saturated mid-range across all values.
vmagLUT = lut("vmag",     0.0, 1.0, "Cool to Warm")
# Discrete 3-bin LUT for refinement level (0, 1, 2). Set the colors here;
# `InterpretValuesAsCategories` is enabled on TR's display side below.
levelLUT= lut("level",    0,   2,   "Viridis (matplotlib)")

# --- Filters reused across panels --------------------------------------------
# Detect whether CellData "level" exists (AMR run) or not (uniform mesh).
# The mars_tgv writer only emits 'level' AFTER the first adapt step, so frame 0
# of an AMR run won't have it. Probe a frame near the middle of the PVD to
# get a reliable answer; if found, treat the entire run as AMR.
_has_level = False
_probe_t   = all_times[min(len(all_times) - 1, len(all_times) // 2)]
reader.UpdatePipeline(time=_probe_t)
_pinfo = reader.GetDataInformation().GetCellDataInformation()
for _ai in range(_pinfo.GetNumberOfArrays()):
    if _pinfo.GetArrayInformation(_ai).GetName() == "level":
        _has_level = True
        break
print(f"AMR level CellData {'present' if _has_level else 'NOT present (uniform mesh)'} (probed t={_probe_t})")
# Restore the reader to the first frame for downstream filter setup.
reader.UpdatePipeline(time=all_times[0])

# Compute |velocity| as a PointData scalar so we can color the iso surface by
# velocity magnitude (ONERA aesthetic: vortex ropes colored by |u|, not |omega|).
vmagCalc = Calculator(Input=reader)
vmagCalc.ResultArrayName = "vmag"
vmagCalc.Function        = "sqrt(u*u + v*v + w*w)"

# Isosurface of vorticity magnitude. For periodic TGV at Re=20, |omega| ranges
# 0..~15 at t=0 and decays. The isosurface is built on the omega field, but
# COLORED by the velocity-magnitude field that the Calculator computed (this
# matches ONERA's render where iso tubes are |u|-shaded, not |omega|-shaded).
vortIso = Contour(Input=vmagCalc)
vortIso.ContourBy   = ["POINTS", "omega"]
vortIso.Isosurfaces = [args.vort_iso]

# Slice through middle of domain on z-mid plane for the u-velocity panel
uSlice = Slice(Input=reader)
uSlice.SliceType        = "Plane"
uSlice.SliceType.Origin = [cx, cy, cz]
uSlice.SliceType.Normal = [0, 0, 1]

# Pressure isosurface (low-p lobes are vortex cores). Pressure range for the
# TGV setup we're rendering is ~|p|<=0.5, so use ±0.2 as default contour values.
pIso = Contour(Input=reader)
pIso.ContourBy   = ["POINTS", "p"]
pIso.Isosurfaces = [-0.2, 0.2]

# --- Render views (4 viewports composed into one layout) ---------------------
def make_view(w, h, bg=(0.0, 0.0, 0.0)):
    rv = CreateView("RenderView")
    rv.ViewSize             = [w, h]
    rv.Background           = list(bg)
    # FXAA only -- MSAA breaks offscreen pvbatch on macOS Metal (returns
    # all-black frames). FXAA is shader-based and works reliably.
    try:
        rv.UseFXAA            = 1
    except Exception:
        pass
    # ParaView 5.13+ replaced UseGradientBackground with BackgroundColorMode.
    try:
        rv.BackgroundColorMode = "Single Color"
    except Exception:
        pass
    rv.OrientationAxesVisibility = 0
    return rv

half_w = args.width  // 2
half_h = args.height // 2

# Per-panel padding so the cube renders don't touch the panel borders. The
# simulation panels (TL/TR/BL) render at a tighter inner box; the AMR panel
# (BR) renders at a smaller box still so it reads as a SUMMARY view, not a
# competing main panel.
panel_pad   = 24   # padding around all panels
amr_shrink  = 96   # extra inset just for the BR (AMR) panel

# Standard isometric-ish camera shared across all four views.
cam_focal = [cx, cy, cz]
# Camera at ~2.2 L: pulled back slightly from 2.0 L because the bottom
# corner of the cube was clipping out of the view at the previous distance.
cam_pos   = [cx + 2.25 * L, cy - 2.00 * L, cz + 1.74 * L]
cam_up    = [0, 0, 1]
def aim(view):
    view.CameraFocalPoint = cam_focal
    view.CameraPosition   = cam_pos
    view.CameraViewUp     = cam_up
    view.CameraParallelProjection = 0
    view.Update()

# AMR: configure the level LUT as a categorical 3-color (coarsened/base/refined).
# cube16 starts at cstone tree level 4; coarsen->3, refine->5. Annotations are
# read by the color legend; IndexedColors picks the discrete fill per category.
if _has_level:
    levelLUT.InterpretValuesAsCategories = 1
    levelLUT.AnnotationsInitialized      = 1
    levelLUT.Annotations  = ["3", "coarsened",
                             "4", "base",
                             "5", "refined"]
    # Yellow refinement (#facc15). Strongest contrast against the dominant
    # blue palette (TL Cool-to-Warm blue end, TR Blue-Orange omega blue cloud,
    # BL Cool-to-Warm Extended iso). Cyan blended into the blue cloud and
    # was hard to read; yellow is the complementary high-contrast pick.
    levelLUT.IndexedColors = [0.21, 0.49, 0.72,   # level 3 -> tab:blue
                              0.85, 0.85, 0.85,   # level 4 -> light gray
                              0.98, 0.80, 0.08]   # level 5 -> bright yellow

# Top-Left: 3D mesh + |velocity| surface. Using |u| instead of u so all six
# cube faces show signal. The x-component u is identically zero on the x=0
# and x=1 faces (because sin(k*0) = sin(k*L) = 0 for one full period),
# which left two faces unlit when we colored by raw u.
view_TL = make_view(half_w - 2 * panel_pad, half_h - 2 * panel_pad)
dispTL = Show(vmagCalc, view_TL)
dispTL.Representation = "Surface With Edges"
# Light gray edges: visible against both the dark-blue low-velocity faces
# and the warm high-velocity ones. Pure black edges merged into the dark
# blue corners of the cube and the AMR-refined corners looked black on black.
dispTL.EdgeColor      = [0.75, 0.75, 0.75]
ColorBy(dispTL, ("POINTS", "vmag"))
dispTL.RescaleTransferFunctionToDataRange(True, False)
aim(view_TL)

# Top-Right: vorticity volume render (the yellow panel) for ALL runs.
view_TR = make_view(half_w - 2 * panel_pad, half_h - 2 * panel_pad)
dispTR = Show(reader, view_TR)
dispTR.Representation = "Volume"
ColorBy(dispTR, ("POINTS", "omega"))
dispTR.RescaleTransferFunctionToDataRange(True, False)
aim(view_TR)

# Bottom-Left: vorticity isosurface colored by |velocity|. Phong-shaded for the
# ONERA "twisted ropes" look -- specular highlight from a single key light.
view_BL = make_view(half_w - 2 * panel_pad, half_h - 2 * panel_pad)
dispBL = Show(vortIso, view_BL)
dispBL.Representation = "Surface"
ColorBy(dispBL, ("POINTS", "vmag"))
# Specular makes the iso surface look metallic / waxed, picks up highlights
# along the curvature -- closest match to ONERA's rope appearance.
try:
    dispBL.Specular        = 0.6
    dispBL.SpecularPower   = 60.0
    dispBL.SpecularColor   = [1.0, 1.0, 1.0]
    dispBL.Ambient         = 0.15
    dispBL.Diffuse         = 0.8
except Exception:
    pass
aim(view_BL)

# Bottom-Right: refined cells only (AMR runs) or u-slice + p-iso fallback.
# For AMR the refined-cells-only Threshold view is the "AMR-only" panel the
# user asked for: it floats just the refined cells in empty space so the
# spatial pattern of where AMR cares is immediately readable.
view_BR = make_view(half_w - 2 * panel_pad - amr_shrink,
                    half_h - 2 * panel_pad - amr_shrink)
if _has_level:
    amrRefined = Threshold(Input=reader)
    amrRefined.Scalars         = ["CELLS", "level"]
    amrRefined.ThresholdMethod = "Between"
    amrRefined.LowerThreshold  = 5
    amrRefined.UpperThreshold  = 100
    dispBR = Show(amrRefined, view_BR)
    dispBR.Representation = "Surface With Edges"
    dispBR.EdgeColor      = [0.2, 0.2, 0.2]
    ColorBy(dispBR, ("CELLS", "level"))
    # Boost ambient lighting so the refined-cells color reads as bright orange
    # not brown. Phong shading on small cubes dims them substantially otherwise.
    try:
        dispBR.Ambient  = 0.55
        dispBR.Diffuse  = 0.55
        dispBR.Specular = 0.05
    except Exception:
        pass
else:
    dispBR_slice = Show(uSlice, view_BR)
    dispBR_slice.Representation = "Surface"
    ColorBy(dispBR_slice, ("POINTS", "u"))
    dispBR_iso = Show(pIso, view_BR)
    dispBR_iso.Representation = "Surface"
    ColorBy(dispBR_iso, ("POINTS", "p"))
aim(view_BR)

# Bottom AMR strip removed per user request. The previous AM (full mesh) view
# duplicated TL once TL switched to coloring by level, and AS was just a KE
# plot that's now a smaller inset in the top header. Layout is plain 2x2.
view_AMR_3D    = None
view_AMR_slice = None

# --- Compose layout: 2x2 grid -------------------------------------------------
# pvbatch can't easily save a layout (multiview-per-step is tricky in headless
# mode); the conventional approach is render each view to its own image then
# composite with PIL. We do that.
try:
    from PIL import Image, ImageDraw, ImageFont
    have_pil = True
except ImportError:
    have_pil = False
    print("WARNING: PIL not available; per-panel PNGs will be written without compositing")

# Load a reasonably-sized font for the overlay. Falls back to PIL default if
# Helvetica isn't found (e.g. on Alps/Linux pvbatch).
def _load_font():
    if not have_pil:
        return None
    for path in (
        "/System/Library/Fonts/Helvetica.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf",
    ):
        try:
            return ImageFont.truetype(path, 22)
        except Exception:
            continue
    return ImageFont.load_default()
_font = _load_font()
_font_big   = None
_font_small = None
if have_pil:
    try:
        _font_big = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 40)
    except Exception:
        _font_big = _font
    try:
        _font_small = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 14)
    except Exception:
        _font_small = _font

# Pull dataset stats out of the current reader state. KE is a lumped sum that
# matches what the solver writes (0.5 * |u|^2 over cells), computed locally
# in pvbatch from the per-node u,v,w arrays.
def _frame_stats():
    info = reader.GetDataInformation()
    n_cells = info.GetNumberOfCells()
    n_pts   = info.GetNumberOfPoints()
    # Approximate KE via vtkArray range/mean: cheap, no full reduction needed.
    # We use the L2-norm range of |u|^2 + |v|^2 + |w|^2 to get a representative
    # magnitude that tracks KE decay even if not exact.
    pdi = info.GetPointDataInformation()
    KE_proxy = 0.0
    for i in range(pdi.GetNumberOfArrays()):
        a = pdi.GetArrayInformation(i)
        if a.GetName() == "velocity":
            for c in range(3):
                r = a.GetComponentRange(c)
                KE_proxy += 0.5 * max(abs(r[0]), abs(r[1])) ** 2
            break
    return n_cells, n_pts, KE_proxy

def _frame_has_level():
    info = reader.GetDataInformation().GetCellDataInformation()
    for i in range(info.GetNumberOfArrays()):
        if info.GetArrayInformation(i).GetName() == "level":
            return True
    return False

# Pre-compute KE(t) at every PVD timestep via PointDataToCellData +
# IntegrateVariables. KE = 0.5 * sum_cells( |u|^2 * V_cell ).
# This is the same lumped-mass quantity the solver reports in its step log,
# so the plotted curve matches what the user saw in the run output.
print("Computing KE(t) for all timesteps...")
keVsT = []
_vcalc = Calculator(Input=reader)
_vcalc.ResultArrayName = "vsq"
_vcalc.Function        = "u*u + v*v + w*w"
_p2c    = PointDatatoCellData(Input=_vcalc)
_p2c.PointDataArraytoprocess = ["vsq"]
_integ  = IntegrateVariables(Input=_p2c)
for t in all_times:
    _integ.UpdatePipeline(time=t)
    info = _integ.GetCellDataInformation()  # IntegrateVariables emits per-cell totals
    cdi  = servermanager.Fetch(_integ)  # pull the (single-cell) integrated result host-side
    cd   = cdi.GetCellData()
    a    = cd.GetArray("vsq")
    ke   = 0.5 * a.GetTuple1(0) if a is not None else 0.0
    keVsT.append(ke)
print(f"KE range: {min(keVsT):.4e} -> {max(keVsT):.4e}")

# Output canvas: 2x2 (4 panels) on top + bottom strip (2 AMR panels) when AMR is on.
# Top half is half_w x half_h per cell -> args.width x args.height.
# Bottom strip is the same width with half_h height -> adds half_h to total height.
# Plain 2x2 layout. Bottom AMR strip was dropped per user request; the AMR
# story now lives in the TL (mesh + level) and BR (refined cells only) panels.
canvas_h = args.height
canvas_w = args.width

def render_step(step_idx, t):
    reader.UpdatePipeline(time=t)
    # Hybrid iso threshold: use the user value when it's a sensible fraction of
    # the current |omega|max (so early-time renders are stable and the user can
    # tune); fall back to 15% of running max when the user value is too high
    # for the current data (late-time, after decay). This keeps the BL panel
    # visually full at all times while still showing decay -- you see the
    # PATTERN of vortex sheets change even when their amplitude shrinks.
    om_info = reader.GetDataInformation().GetPointDataInformation()
    om_max = 1.0
    for i in range(om_info.GetNumberOfArrays()):
        a = om_info.GetArrayInformation(i)
        if a.GetName() == "omega":
            om_max = max(1e-6, a.GetComponentRange(0)[1])
            break
    iso_val = min(args.vort_iso, 0.15 * om_max)
    vortIso.Isosurfaces = [iso_val]
    all_views = [view_TL, view_TR, view_BL, view_BR]
    if view_AMR_3D    is not None: all_views.append(view_AMR_3D)
    if view_AMR_slice is not None: all_views.append(view_AMR_slice)
    for v in all_views:
        v.ViewTime = t
        v.Update()

    # LUTs are FROZEN to step-0 ranges. As the field decays, colors fade
    # against the (unchanging) colorbar -- that fading IS the decay signal.
    # Per-frame rescaling kept colors saturated and hid the energy loss.

    paths = {}
    panel_list = [("TL", view_TL), ("TR", view_TR), ("BL", view_BL), ("BR", view_BR)]
    if view_AMR_3D    is not None: panel_list.append(("AM", view_AMR_3D))
    if view_AMR_slice is not None: panel_list.append(("AS", view_AMR_slice))
    for name, v in panel_list:
        p = os.path.join(out_dir, f"_{name}_{step_idx:04d}.png")
        SaveScreenshot(p, v)
        paths[name] = p

    if have_pil:
        imgs = {k: Image.open(p) for k, p in paths.items()}
        canvas = Image.new("RGB", (canvas_w, canvas_h), (0, 0, 0))
        # Each panel sits inside its slot with `panel_pad` margin. The BR slot
        # is the AMR view: drawn smaller so it reads as a summary, not a main
        # panel of the simulation.
        canvas.paste(imgs["TL"], (panel_pad,           panel_pad))
        canvas.paste(imgs["TR"], (half_w + panel_pad,  panel_pad))
        canvas.paste(imgs["BL"], (panel_pad,           half_h + panel_pad))
        amr_off = amr_shrink // 2
        canvas.paste(imgs["BR"], (half_w + panel_pad + amr_off,
                                  half_h + panel_pad + amr_off))

        draw = ImageDraw.Draw(canvas)

        # KE-vs-time inset: centered at the canvas mid-cross (the boundary
        # intersection of the 4 panels). That spot is empty in all views
        # because each panel renders its own background there, and the inset
        # straddling all 4 panels reads as "global summary".
        if len(keVsT) > 0:
            # Bigger inset (was 320x180); ticks and title use _font_small so
            # the text doesn't look bloated at the larger plot size. Inset
            # is centered horizontally; vertically offset DOWN into the
            # bottom-row panels so the top-row labels (which now sit at
            # y=half_h-30) are not occluded.
            iw, ih = 440, 270
            px = (canvas_w - iw) // 2
            py = (canvas_h - ih) // 2
            pad_l, pad_r, pad_t, pad_b = 70, 18, 32, 38
            plot_x0 = px + pad_l
            plot_y0 = py + pad_t
            plot_w  = iw - pad_l - pad_r
            plot_h  = ih - pad_t - pad_b
            # White panel, thin black axes.
            draw.rectangle([(px, py), (px + iw, py + ih)],
                           fill=(245, 245, 245), outline=(60, 60, 60), width=1)
            draw.rectangle([(plot_x0, plot_y0),
                            (plot_x0 + plot_w, plot_y0 + plot_h)],
                           outline=(60, 60, 60), width=1)
            import math
            tmin, tmax = float(all_times[0]), float(all_times[-1])
            ke_pos = [max(v, 1e-12) for v in keVsT]
            ymin = math.log10(min(ke_pos))
            ymax = math.log10(max(ke_pos))
            yspan = max(1e-6, ymax - ymin)
            ymin -= 0.05 * yspan
            ymax += 0.05 * yspan
            def _xy(ti, kei):
                u = (ti - tmin) / max(1e-12, tmax - tmin)
                logk = math.log10(max(kei, 1e-12))
                v = (logk - ymin) / (ymax - ymin)
                return (plot_x0 + u * plot_w, plot_y0 + (1.0 - v) * plot_h)
            pts = [_xy(all_times[k], keVsT[k]) for k in range(len(keVsT))]
            # Dark-blue line on white (clean look).
            draw.line(pts, fill=(40, 80, 200), width=2)
            cur = _xy(t, keVsT[step_idx])
            # Cursor: orange dot on white outline -- not the saturated red the
            # user pushed back on.
            draw.ellipse([(cur[0] - 5, cur[1] - 5), (cur[0] + 5, cur[1] + 5)],
                         fill=(245, 150, 30), outline=(255, 255, 255), width=1)
            # 3 x-ticks, 3 y-ticks (sparse so the inset stays readable).
            # 5 x-ticks and 5 y-ticks for the bigger inset, small font.
            for tk in range(5):
                tt = tmin + tk * (tmax - tmin) / 4
                tx = plot_x0 + tk * plot_w / 4
                draw.line([(tx, plot_y0 + plot_h), (tx, plot_y0 + plot_h + 3)],
                          fill=(60, 60, 60), width=1)
                draw.text((tx - 12, plot_y0 + plot_h + 6), f"{tt:.2f}",
                          fill=(40, 40, 40), font=_font_small)
            for tk in range(5):
                ly  = ymin + tk * (ymax - ymin) / 4
                ty  = plot_y0 + (1.0 - tk / 4) * plot_h
                draw.line([(plot_x0 - 3, ty), (plot_x0, ty)],
                          fill=(60, 60, 60), width=1)
                draw.text((plot_x0 - 62, ty - 7), f"1e{ly:+.1f}",
                          fill=(40, 40, 40), font=_font_small)
            # Title (small).
            draw.text((px + 10, py + 6), "KE vs t (log)",
                      fill=(20, 20, 20), font=_font_small)

        # Per-frame overlay: time, n_cells, n_points + panel labels.
        n_cells, n_pts, _ = _frame_stats()
        KE_true = keVsT[step_idx] if step_idx < len(keVsT) else 0.0
        header = f"t = {t:.4f}     elements = {n_cells:>7}     nodes = {n_pts:>7}     KE = {KE_true:.3e}"
        draw.rectangle([(0, 0), (canvas_w, 56)], fill=(0, 0, 0, 220))
        draw.text((14, 8), header, fill=(255, 255, 255), font=_font_big)

        labels = {
            "TL": "mesh + |velocity|",
            "TR": "vorticity volume",
            "BL": "vorticity isosurface",
            "BR": ("AMR refined cells (level 5)" if _has_level
                   else "u-slice + pressure isosurfaces"),
        }
        # Labels at the BOTTOM-LEFT corner of each panel (per user request).
        # TL bottom-left, BL/BR bottom-left of their respective bottom-row
        # panels. TR's natural bottom-left (half_w+10, half_h-30) lands
        # inside the centered KE inset, so TR's label is shifted to the
        # bottom-RIGHT of its panel where there's empty space.
        positions = {
            "TL": (10,                  half_h - 30),
            "TR": (canvas_w - 250,      half_h - 30),
            "BL": (10,                  canvas_h - 30),
            "BR": (half_w + 10,         canvas_h - 30),
        }
        # Labels for the four real panels (no longer a bottom AMR strip).
        label_slots = [("TL",), ("TR",), ("BL",), ("BR",)]
        for k in label_slots:
            kname = k[0]
            if kname not in labels: continue
            x, y = positions[kname]
            label = labels[kname]
            draw.text((x+1, y+1), label, fill=(0, 0, 0), font=_font)
            draw.text((x,   y),   label, fill=(255, 255, 255), font=_font)

        out = os.path.join(out_dir, f"frame_{step_idx:04d}.png")
        canvas.save(out)
        for p in paths.values():
            try: os.remove(p)
            except OSError: pass

print(f"rendering {len(all_times)} frames at {args.width}x{args.height}")
for i, t in enumerate(all_times):
    render_step(i, t)
    print(f"  frame {i+1}/{len(all_times)}  t={t:.4f}")

print(f"frames in {out_dir}/")

# --- Assemble mp4 -------------------------------------------------------------
if not args.no_mp4 and have_pil:
    import shutil, subprocess
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg:
        mp4 = prefix + "_4panel.mp4"
        cmd = [ffmpeg, "-y", "-framerate", str(args.fps),
               "-i", os.path.join(out_dir, "frame_%04d.png"),
               "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
               mp4]
        print("ffmpeg:", " ".join(cmd))
        subprocess.run(cmd, check=False)
        print(f"video: {mp4}")
    else:
        print("ffmpeg not found in PATH; skipping video assembly")
