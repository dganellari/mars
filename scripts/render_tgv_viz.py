#!/usr/bin/env pvbatch
# Render a TGV multi-rank video from MARS VTU/PVD output.
#
# Run on the cluster (headless ParaView, e.g. via pvbatch):
#   pvbatch scripts/render_tgv_viz.py --pvd tgv.pvd --frames-dir tgv_frames --field umag
# Then ffmpeg locally:
#   ffmpeg -framerate 12 -i tgv_frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p -crf 18 tgv.mp4
#
# MARS TGV writer (mars_vtu_parallel_writer.hpp) writes point arrays:
# u, v, w, p, omega, velocity.
# Field can be: umag (computed from u,v,w), p, omega, or any name above.

import argparse
import math
import os
import sys

if "--gpu-filters" not in sys.argv:
    os.environ.setdefault("VTK_DISABLE_VISKORES_OVERRIDES", "1")
    os.environ.setdefault("VTKM_DEVICE", "Serial")
    os.environ.setdefault("PARAVIEW_USE_ACCELERATED_FILTERS", "0")

from paraview.simple import *  # noqa: F401,F403,E402


def main():
    ap = argparse.ArgumentParser(description="Render MARS TGV multi-rank video")
    ap.add_argument("--pvd", required=True, help="path to .pvd from --vtu-output")
    ap.add_argument("--frames-dir", default="tgv_frames", help="PNG frames dir")
    ap.add_argument("--field", default="umag",
                    help="field to color by: umag (default), p, omega, u/v/w")
    ap.add_argument("--res", default="1280x720", help="WxH render resolution")
    ap.add_argument("--orbit-deg", type=float, default=45.0,
                    help="total azimuth orbit over the run (0 = fixed cam)")
    ap.add_argument("--preset", default="Viridis (matplotlib)",
                    help="initial color preset; falls back if unknown")
    ap.add_argument("--zoom", type=float, default=1.3,
                    help="camera Dolly factor (>1 zooms in)")
    ap.add_argument("--vmax", type=float, default=0.0,
                    help="clamp color range max (0=auto). For TGV V0=1, try 1.5; "
                         "auto picks the last frame's max which gets dominated by "
                         "the seam blow-up so early frames look black.")
    ap.add_argument("--vmin", type=float, default=0.0,
                    help="clamp color range min (default 0)")
    ap.add_argument("--gpu-filters", action="store_true")
    ap.add_argument("--show-outline", action="store_true", default=True,
                    help="show bounding-box outline")
    args = ap.parse_args()

    if not os.path.exists(args.pvd):
        sys.exit("ERROR: pvd not found: %s" % args.pvd)
    W, H = (int(x) for x in args.res.lower().split("x"))

    reader = PVDReader(FileName=args.pvd)
    reader.UpdatePipeline()
    tvals = list(reader.TimestepValues) if reader.TimestepValues else [0.0]

    # Discover available arrays for safety
    pd = reader.PointData
    arr_names = [pd.GetArray(i).GetName() for i in range(pd.GetNumberOfArrays())]
    print("[render] Point arrays found:", arr_names)

    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [W, H]
    view.Background = [0.10, 0.12, 0.16]
    view.Background2 = [0.22, 0.26, 0.32]
    view.UseColorPaletteForBackground = 0
    try:
        view.UseGradientBackground = 1
    except Exception:
        pass
    view.OrientationAxesVisibility = 1

    bounds = reader.GetDataInformation().GetBounds()
    cx = 0.5 * (bounds[0] + bounds[1])
    cy = 0.5 * (bounds[2] + bounds[3])
    cz = 0.5 * (bounds[4] + bounds[5])
    diag = math.sqrt((bounds[1] - bounds[0])**2 +
                     (bounds[3] - bounds[2])**2 +
                     (bounds[5] - bounds[4])**2)

    if args.show_outline:
        try:
            outline = Outline(Input=reader)
            odisp = Show(outline, view)
            odisp.DiffuseColor = [0.7, 0.75, 0.82]
            odisp.AmbientColor = [0.7, 0.75, 0.82]
            try:
                odisp.LineWidth = 2
            except Exception:
                pass
        except Exception as e:
            print("[render] outline skipped:", e)

    # Build the color source: either an existing scalar/vector array, or compute
    # |u| from u,v,w via a Calculator. The display pipeline then reads "color_field".
    color_field = args.field
    source = reader

    if args.field == "umag":
        # compute |u| from u,v,w (works whether 'velocity' vector exists or not)
        if all(n in arr_names for n in ("u", "v", "w")):
            mag = Calculator(Input=reader)
            mag.ResultArrayName = "umag"
            mag.Function = "sqrt(u*u + v*v + w*w)"
            source = mag
        elif "velocity" in arr_names:
            mag = Calculator(Input=reader)
            mag.ResultArrayName = "umag"
            mag.Function = "mag(velocity)"
            source = mag
        else:
            print("[render] cannot compute umag; falling back to first array")
            color_field = arr_names[0] if arr_names else None
    elif args.field not in arr_names:
        print("[render] WARNING: field '%s' not found; trying it anyway" % args.field)

    if color_field is None:
        sys.exit("ERROR: no fields available to render")

    # Robust color range: evaluate on the LAST (developed) frame only
    print("[render] computing %s range on developed frame..." % color_field)
    try:
        reader.UpdatePipeline(tvals[-1])
        source.UpdatePipeline(tvals[-1])
    except Exception as e:
        print("[render] UpdatePipeline failed:", e)
    arr = source.PointData.GetArray(color_field) if hasattr(source, "PointData") else None
    if arr is None:
        # Refresh
        try:
            source.UpdatePipeline(tvals[-1])
            arr = source.PointData.GetArray(color_field)
        except Exception:
            pass
    if arr:
        rng = arr.GetRange()
        smin, smax = rng[0], (rng[1] if rng[1] > rng[0] else rng[0] + 1.0)
    else:
        smin, smax = 0.0, 1.0
    print("[render] %s data range: [%g, %g]" % (color_field, smin, smax))
    # User override (--vmin/--vmax): clamp so blow-up frames don't crush the color
    # map for the early clean steps. For TGV V0=1, --vmax=1.5 keeps the color band
    # readable across the whole sequence.
    if args.vmax > 0:
        smax = args.vmax
    if args.vmin != 0.0 or args.vmax > 0:
        smin = args.vmin
    print("[render] %s color range: [%g, %g]" % (color_field, smin, smax))

    display = Show(source, view)
    display.Representation = "Surface"
    ColorBy(display, ("POINTS", color_field))

    lut = GetColorTransferFunction(color_field)
    # Try a list of common preset names; ParaView versions disagree on naming
    preset_candidates = [args.preset, "Viridis (matplotlib)", "viridis", "Viridis",
                         "Cool to Warm (Extended)", "Cool to Warm", "Turbo",
                         "Inferno (matplotlib)", "Plasma (matplotlib)",
                         "Rainbow Uniform", "Black-Body Radiation"]
    for p in preset_candidates:
        try:
            lut.ApplyPreset(p, True)
            print("[render] color preset: %s" % p)
            break
        except Exception:
            continue
    try:
        lut.RescaleTransferFunction(smin, smax)
    except Exception:
        pass

    display.SetScalarBarVisibility(view, True)
    bar = GetScalarBar(lut, view)
    bar.Title = color_field
    bar.ComponentTitle = ""
    bar.TitleColor = [1, 1, 1]
    bar.LabelColor = [0.9, 0.9, 0.9]
    bar.TitleFontSize = 16
    bar.LabelFontSize = 13
    for attr, val in (("WindowLocation", "Lower Right Corner"),
                      ("ScalarBarLength", 0.33)):
        try:
            setattr(bar, attr, val)
        except Exception:
            pass

    # Title + time annotation
    try:
        title = Text()
        title.Text = "MARS TGV  --  4-rank periodic Taylor-Green Vortex"
        tdisp = Show(title, view)
        tdisp.Color = [1, 1, 1]
        tdisp.FontSize = 18
        try:
            tdisp.WindowLocation = "Upper Left Corner"
        except Exception:
            pass
    except Exception as e:
        print("[render] title skipped:", e)

    try:
        tann = AnnotateTimeFilter(Input=reader)
        for fmt in ("t = {time:.4f}", "t = %f"):
            try:
                tann.Format = fmt
                break
            except Exception:
                continue
        adisp = Show(tann, view)
        adisp.Color = [0.85, 0.9, 1.0]
        adisp.FontSize = 16
        try:
            adisp.WindowLocation = "Lower Left Corner"
        except Exception:
            pass
    except Exception as e:
        print("[render] time annotation skipped:", e)

    ResetCamera(view)
    cam = GetActiveCamera()
    cam.Elevation(20)
    cam.Azimuth(30)
    view.CameraFocalPoint = [cx, cy, cz]
    cam.Dolly(args.zoom)
    view.CameraParallelProjection = 0
    # ResetCameraClippingRange() was removed in ParaView 6.x; ResetCamera()
    # above already recomputed the clipping range.

    os.makedirs(args.frames_dir, exist_ok=True)
    total = len(tvals)
    az_per_frame = args.orbit_deg / max(1, total - 1)
    print("[render] %d frames, %dx%d, %s range [%g,%g]"
          % (total, W, H, color_field, smin, smax))

    written = 0
    for i, t in enumerate(tvals):
        view.ViewTime = t
        try:
            reader.UpdatePipeline(t)
            if source is not reader:
                source.UpdatePipeline(t)
        except Exception as e:
            print("  [warn] UpdatePipeline t=%g failed (%s)" % (t, e))
        cam.Azimuth(az_per_frame)
        fn = os.path.join(args.frames_dir, "frame_%04d.png" % i)
        try:
            Render(view)
            SaveScreenshot(fn, view, ImageResolution=[W, H])
            written += 1
        except Exception as e:
            print("  [warn] frame %d failed: %s" % (i, e))
        if i % 5 == 0:
            print("  frame %d/%d  t=%g" % (i + 1, total, t))

    print("[render] wrote %d/%d PNG frames to %s" % (written, total, args.frames_dir))
    print("[render] now run locally:")
    print("  ffmpeg -framerate 12 -i %s/frame_%%04d.png \\" % args.frames_dir)
    print("    -c:v libx264 -pix_fmt yuv420p -crf 18 tgv.mp4")


if __name__ == "__main__":
    main()
