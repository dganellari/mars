#!/usr/bin/env pvbatch
# Render a GREAT pump NS-solver flow video from the MARS pump VTU/PVD output.
#
# A pump video should look alive and read clearly. This pipeline builds:
#   - a clean translucent geometry SHELL for spatial context (Fresnel-ish faint),
#   - inlet->outlet STREAMLINES as tubes colored by |velocity| (the flow path +
#     recirculation), with the seed riding the high-speed inlet region,
#   - the .pvd TIME animation so the jet is seen DEVELOPING from rest,
#   - a slow CINEMATIC camera orbit, ambient-occlusion depth, a dark studio
#     background, a clean color bar + a time/flow-through annotation,
#   - high-res frames encoded to mp4.
#
# Run on the cluster (headless ParaView):
#   pvbatch scripts/render_pump_flow.py --pvd pump.pvd --out pump_flow.mp4
# Frames-only (encode yourself for full codec control):
#   pvbatch scripts/render_pump_flow.py --pvd pump.pvd --frames-dir frames
#   ffmpeg -framerate 30 -i frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p -crf 18 pump_flow.mp4
#
# Pump writer (mars_vtu_parallel_writer.hpp): <prefix>.pvd -> per-step .pvtu;
# point scalars u,v,w,p and point vector "velocity".

import argparse
import math
import os
import sys

from paraview.simple import *  # noqa: F401,F403


def main():
    ap = argparse.ArgumentParser(description="Render a great MARS pump flow video")
    ap.add_argument("--pvd", required=True, help="path to <prefix>.pvd from the pump run")
    ap.add_argument("--out", default="pump_flow.mp4", help="output mp4")
    ap.add_argument("--frames-dir", default="pump_frames", help="PNG frames dir (always written)")
    ap.add_argument("--res", default="1920x1080", help="WxH render resolution")
    ap.add_argument("--field", default="velocity", help="vector field (default velocity)")
    ap.add_argument("--n-streamlines", type=int, default=900, help="streamline seed count")
    ap.add_argument("--orbit-deg", type=float, default=120.0, help="total azimuth sweep over the run")
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--hold-frames", type=int, default=20,
                    help="extra orbit-only frames at the end on the final (developed) state")
    ap.add_argument("--preset", default="Viridis (matplotlib)", help="color preset")
    ap.add_argument("--streamlines", action="store_true",
                    help="use StreamTracer tubes (prettier but SLOW/can hang on big tet meshes; "
                         "default is a fast volume clip colored by |velocity|)")
    args = ap.parse_args()

    if not os.path.exists(args.pvd):
        sys.exit("ERROR: pvd not found: %s" % args.pvd)
    W, H = (int(x) for x in args.res.lower().split("x"))

    reader = PVDReader(FileName=args.pvd)
    reader.UpdatePipeline()
    tvals = list(reader.TimestepValues) if reader.TimestepValues else [0.0]

    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [W, H]
    view.Background = [0.03, 0.035, 0.05]      # near-black studio
    view.Background2 = [0.10, 0.12, 0.16]      # subtle gradient top
    view.UseColorPaletteForBackground = 0
    try:
        view.UseGradientBackground = 1
    except Exception:
        pass
    view.OrientationAxesVisibility = 0
    # depth + soft shadows for a 3D feel (ParaView 5.9+; guarded)
    for attr, val in (("UseAmbientOcclusion", 1), ("EnableRayTracing", 0)):
        try:
            setattr(view, attr, val)
        except Exception:
            pass

    bounds = reader.GetDataInformation().GetBounds()
    cx = 0.5 * (bounds[0] + bounds[1]); cy = 0.5 * (bounds[2] + bounds[3]); cz = 0.5 * (bounds[4] + bounds[5])
    diag = math.sqrt((bounds[1]-bounds[0])**2 + (bounds[3]-bounds[2])**2 + (bounds[5]-bounds[4])**2)

    # ---- translucent geometry shell ----
    surf = Show(reader, view)
    surf.Representation = "Surface"
    surf.Opacity = 0.10
    surf.AmbientColor = [0.55, 0.62, 0.72]
    surf.DiffuseColor = [0.55, 0.62, 0.72]
    surf.Specular = 0.3
    ColorBy(surf, None)

    # ---- |velocity| ----
    mag = Calculator(Input=reader)
    mag.ResultArrayName = "speed"
    mag.Function = "mag(%s)" % args.field
    mag.UpdatePipeline()

    # robust color range: use the LAST frame (developed jet) so colors don't wash
    # out early. Rescale to ~98th percentile feel via the data range of frame -1.
    reader.UpdatePipeline(tvals[-1])
    mag.UpdatePipeline(tvals[-1])
    srng = mag.PointData.GetArray("speed").GetRange() if mag.PointData.GetArray("speed") else (0.0, 1.0)
    smin, smax = srng[0], (srng[1] if srng[1] > srng[0] else srng[0] + 1.0)
    lut = GetColorTransferFunction("speed")
    # Preset names vary across ParaView builds; try the requested one then a list
    # of common synonyms, and just keep the default if none apply.
    _preset_candidates = [args.preset, "Viridis (matplotlib)", "viridis", "Viridis",
                          "Cool to Warm (Extended)", "Cool to Warm", "Turbo",
                          "Inferno (matplotlib)", "Plasma (matplotlib)", "Rainbow Uniform"]
    for _p in _preset_candidates:
        try:
            lut.ApplyPreset(_p, True)
            print("[render] color preset: %s" % _p)
            break
        except Exception:
            continue
    lut.RescaleTransferFunction(smin, smax)
    try:
        lut.NanColor = [0.2, 0.2, 0.2]
    except Exception:
        pass

    # ---- flow visualization ----
    # DEFAULT (fast, robust): a CLIP through the volume colored by |velocity| --
    # shows the inlet->outlet jet without the GPU StreamTracer (which hangs on
    # large tet meshes in headless ParaView). Streamlines are opt-in via
    # --streamlines (slow but prettier).
    if args.streamlines:
        stream = StreamTracer(Input=mag, SeedType="Point Cloud")
        stream.Vectors = ["POINTS", args.field]
        stream.IntegrationDirection = "BOTH"
        try:
            stream.IntegratorType = "Runge-Kutta 4-5"
        except Exception:
            pass
        stream.MaximumStreamlineLength = diag * 5.0
        stream.SeedType.Center = [cx, cy, cz]
        stream.SeedType.Radius = diag * 0.48
        stream.SeedType.NumberOfPoints = args.n_streamlines
        stream.UpdatePipeline(tvals[-1])
        tubes = Tube(Input=stream)
        tubes.Radius = diag * 0.0016
        tubes.Capping = 1
        tubes.UpdatePipeline(tvals[-1])
        flow = tubes
    else:
        # Clip the volume in half along the dominant flow axis so the interior
        # jet is visible; color the cut volume by |velocity|.
        clip = Clip(Input=mag)
        try:
            clip.ClipType = "Plane"
            clip.ClipType.Origin = [cx, cy, cz]
            # cut perpendicular to the largest extent (the flow usually runs along it)
            ext = [bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4]]
            axis = ext.index(max(ext))
            nrm = [0, 0, 0]; nrm[axis] = 1
            clip.ClipType.Normal = nrm
        except Exception:
            pass
        clip.UpdatePipeline(tvals[-1])
        flow = clip

    sdisp = Show(flow, view)
    sdisp.Representation = "Surface"
    ColorBy(sdisp, ("POINTS", "speed"))
    sdisp.Specular = 0.4
    sdisp.SetScalarBarVisibility(view, True)

    bar = GetScalarBar(lut, view)
    bar.Title = "|velocity|"
    bar.ComponentTitle = ""
    bar.TitleColor = [1, 1, 1]; bar.LabelColor = [0.9, 0.9, 0.9]
    bar.TitleFontSize = 16; bar.LabelFontSize = 13
    for attr, val in (("WindowLocation", "Lower Right Corner"), ("ScalarBarLength", 0.33)):
        try:
            setattr(bar, attr, val)
        except Exception:
            pass

    # ---- title + time annotation (all optional/cosmetic; never fail the render) ----
    try:
        title = Text()
        title.Text = "MARS pump  --  incompressible NS (CVFEM/DDT)"
        tdisp = Show(title, view)
        tdisp.Color = [1, 1, 1]; tdisp.FontSize = 18
        try:
            tdisp.WindowLocation = "Upper Left Corner"
        except Exception:
            pass
    except Exception as e:  # noqa: BLE001
        print("[render] title annotation skipped: %s" % e)

    try:
        tann = AnnotateTimeFilter(Input=reader)
        for fmt in ("t = {time:.4f}", "t = %f"):
            try:
                tann.Format = fmt
                break
            except Exception:
                continue
        adisp = Show(tann, view)
        adisp.Color = [0.85, 0.9, 1.0]; adisp.FontSize = 16
        try:
            adisp.WindowLocation = "Lower Left Corner"
        except Exception:
            pass
    except Exception as e:  # noqa: BLE001
        print("[render] time annotation skipped: %s" % e)

    # ---- camera ----
    ResetCamera(view)
    cam = GetActiveCamera()
    cam.Elevation(16); cam.Azimuth(25)
    view.CameraFocalPoint = [cx, cy, cz]
    # pull back slightly for headroom
    cam.Dolly(0.85)
    view.CameraParallelProjection = 0

    os.makedirs(args.frames_dir, exist_ok=True)
    total = len(tvals) + max(0, args.hold_frames)
    az_per_frame = args.orbit_deg / max(1, total - 1)
    print("[render] %d frames (%d time + %d hold), %dx%d, seeds=%d, speed[%.3g,%.3g]"
          % (total, len(tvals), args.hold_frames, W, H, args.n_streamlines, smin, smax))

    fidx = 0
    for t in tvals:
        view.ViewTime = t
        reader.UpdatePipeline(t)
        cam.Azimuth(az_per_frame)
        Render(view)
        SaveScreenshot(os.path.join(args.frames_dir, "frame_%04d.png" % fidx), view, ImageResolution=[W, H])
        if fidx % 10 == 0:
            print("  frame %d/%d  t=%.4g" % (fidx + 1, total, t))
        fidx += 1
    # hold on the developed state, keep orbiting
    for _ in range(max(0, args.hold_frames)):
        cam.Azimuth(az_per_frame)
        Render(view)
        SaveScreenshot(os.path.join(args.frames_dir, "frame_%04d.png" % fidx), view, ImageResolution=[W, H])
        fidx += 1

    if args.out:
        try:
            SaveAnimation(args.out, view, ImageResolution=[W, H], FrameRate=args.fps,
                          FrameWindow=[0, total - 1])
            print("[render] wrote %s" % args.out)
        except Exception as e:  # noqa: BLE001
            print("[render] SaveAnimation unavailable (%s). Encode the frames:" % e)
            print("  ffmpeg -framerate %d -i %s/frame_%%04d.png -c:v libx264 -pix_fmt yuv420p -crf 18 %s"
                  % (args.fps, args.frames_dir, args.out))


if __name__ == "__main__":
    main()
