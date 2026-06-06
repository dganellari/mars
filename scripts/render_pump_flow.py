#!/usr/bin/env pvbatch
# Render a pump NS-solver flow video from the MARS pump VTU/PVD output.
#
# Best-for-a-pump visualization: the inlet->outlet JET is the thing to show, so
# we combine (a) a translucent geometry surface for context, (b) inlet-seeded
# STREAMLINES traced through the volume colored by velocity magnitude (the flow
# path + recirculation), and (c) the time animation over the .pvd so the jet is
# seen DEVELOPING. A slow camera orbit gives depth.
#
# Run on the cluster (headless ParaView):
#   pvbatch scripts/render_pump_flow.py --pvd pump.pvd --out pump_flow.mp4
# or render PNG frames + encode yourself:
#   pvbatch scripts/render_pump_flow.py --pvd pump.pvd --frames-dir frames
#
# The pump writer (mars_vtu_parallel_writer.hpp) emits <prefix>.pvd referencing
# per-step .pvtu; fields are point scalars u,v,w,p and point vector "velocity".

import argparse
import os
import sys

from paraview.simple import *  # noqa: F401,F403


def main():
    ap = argparse.ArgumentParser(description="Render MARS pump flow video")
    ap.add_argument("--pvd", required=True, help="path to <prefix>.pvd from the pump run")
    ap.add_argument("--out", default="pump_flow.mp4", help="output mp4 (if encoder available)")
    ap.add_argument("--frames-dir", default=None,
                    help="if set, write PNG frames here instead of/in addition to mp4")
    ap.add_argument("--res", default="1600x900", help="WxH render resolution")
    ap.add_argument("--field", default="velocity", help="vector field to seed/color (default velocity)")
    ap.add_argument("--n-streamlines", type=int, default=400, help="streamline seed count")
    ap.add_argument("--orbit-deg", type=float, default=90.0,
                    help="total camera azimuth sweep over the animation (deg)")
    ap.add_argument("--fps", type=int, default=15)
    ap.add_argument("--no-surface", action="store_true", help="hide the translucent geometry surface")
    args = ap.parse_args()

    if not os.path.exists(args.pvd):
        sys.exit("ERROR: pvd not found: %s" % args.pvd)
    W, H = (int(x) for x in args.res.lower().split("x"))

    # ---- load the time series ----
    reader = PVDReader(FileName=args.pvd)
    reader.UpdatePipeline()
    tvals = list(reader.TimestepValues) if reader.TimestepValues else [0.0]

    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [W, H]
    view.Background = [0.08, 0.09, 0.11]   # dark slate, makes the jet pop
    view.OrientationAxesVisibility = 0

    # ---- (a) translucent geometry surface for spatial context ----
    if not args.no_surface:
        surf = Show(reader, view)
        surf.Representation = "Surface"
        surf.Opacity = 0.12
        surf.AmbientColor = [0.6, 0.65, 0.7]
        surf.DiffuseColor = [0.6, 0.65, 0.7]
        ColorBy(surf, None)

    # ---- velocity magnitude (Calculator) for coloring ----
    mag = Calculator(Input=reader)
    mag.ResultArrayName = "speed"
    mag.Function = "mag(%s)" % args.field
    mag.UpdatePipeline()

    speed_range = mag.PointData.GetArray("speed").GetRange() if mag.PointData.GetArray("speed") else (0.0, 1.0)
    lut = GetColorTransferFunction("speed")
    lut.ApplyPreset("Cool to Warm (Extended)", True)
    lut.RescaleTransferFunction(speed_range[0], speed_range[1])

    # ---- (b) inlet-seeded streamlines through the volume ----
    # Seed from a point cloud near the high-speed (inlet) region. Using the full
    # bounds with a dense seed + both-direction integration captures the jet and
    # any recirculation without needing the exact inlet coords.
    bounds = reader.GetDataInformation().GetBounds()
    cx = 0.5 * (bounds[0] + bounds[1])
    cy = 0.5 * (bounds[2] + bounds[3])
    cz = 0.5 * (bounds[4] + bounds[5])
    diag = ((bounds[1]-bounds[0])**2 + (bounds[3]-bounds[2])**2 + (bounds[5]-bounds[4])**2) ** 0.5

    stream = StreamTracer(Input=mag, SeedType="Point Cloud")
    stream.Vectors = ["POINTS", args.field]
    stream.IntegrationDirection = "BOTH"
    stream.MaximumStreamlineLength = diag * 4.0
    stream.SeedType.Center = [cx, cy, cz]
    stream.SeedType.Radius = diag * 0.45
    stream.SeedType.NumberOfPoints = args.n_streamlines
    stream.UpdatePipeline()

    # tubes for visibility
    tubes = Tube(Input=stream)
    tubes.Radius = diag * 0.0015
    tubes.UpdatePipeline()

    sdisp = Show(tubes, view)
    sdisp.Representation = "Surface"
    ColorBy(sdisp, ("POINTS", "speed"))
    sdisp.RescaleTransferFunctionToDataRange(False, True)
    sdisp.SetScalarBarVisibility(view, True)
    bar = GetScalarBar(lut, view)
    bar.Title = "|velocity|"
    bar.ComponentTitle = ""
    bar.TitleColor = [1, 1, 1]
    bar.LabelColor = [1, 1, 1]

    # ---- camera ----
    ResetCamera(view)
    cam = GetActiveCamera()
    cam.Azimuth(20); cam.Elevation(18)
    view.CameraFocalPoint = [cx, cy, cz]

    # ---- animate over time + orbit ----
    scene = GetAnimationScene()
    scene.UpdateAnimationUsingDataTimeSteps()
    nframes = len(tvals)
    az_per_frame = args.orbit_deg / max(1, nframes - 1)

    frames_dir = args.frames_dir
    if frames_dir:
        os.makedirs(frames_dir, exist_ok=True)

    print("[render] %d frames, res %dx%d, field=%s, seeds=%d" % (nframes, W, H, args.field, args.n_streamlines))
    for i, t in enumerate(tvals):
        scene.AnimationTime = t
        view.ViewTime = t
        reader.UpdatePipeline(t)
        cam.Azimuth(az_per_frame)
        Render(view)
        if frames_dir:
            fn = os.path.join(frames_dir, "frame_%04d.png" % i)
            SaveScreenshot(fn, view, ImageResolution=[W, H])
        print("  frame %d/%d  t=%.4g" % (i + 1, nframes, t))

    # ---- encode mp4 (ParaView's built-in encoder; falls back to frames if unavailable) ----
    if args.out:
        try:
            SaveAnimation(args.out, view, ImageResolution=[W, H], FrameRate=args.fps,
                          FrameWindow=[0, nframes - 1])
            print("[render] wrote %s" % args.out)
        except Exception as e:  # noqa: BLE001
            print("[render] SaveAnimation failed (%s); use the PNG frames + ffmpeg:" % e)
            print("  ffmpeg -framerate %d -i %s/frame_%%04d.png -pix_fmt yuv420p %s"
                  % (args.fps, frames_dir or "frames", args.out))


if __name__ == "__main__":
    main()
