#!/usr/bin/env pvbatch
# Render the Poiseuille channel animation from mars_poiseuille_flow VTU output.
#
# Input (derived from --pvd PREFIX): PREFIX_<field>.pvd. Default field is u
# (the validated streamwise velocity). NOTE: umag = sqrt(u^2+v^2+w^2) is
# contaminated by the incremental-pressure artifact leaking into v/w -- it
# renders as saturated garbage even when the u physics is perfect. The channel is one
# element thick in z, so a single top-down (x-y) view shows everything: the
# uniform plug inflow on the left bending into the parabolic profile
# downstream (fast red core, slow blue walls).
#
# Output:
#   <PREFIX>_movie/frame_NNNN.png   -- per-timestep frames
#   <PREFIX>_movie.mp4              -- assembled animation (if ffmpeg in PATH)
#
# Usage on Mac:
#   /Applications/ParaView-5.13.2.app/Contents/bin/pvbatch \
#       scripts/render_poiseuille.py --pvd poiseuille
# Usage on Alps (offscreen pvbatch):
#   pvbatch scripts/render_poiseuille.py --pvd poiseuille
#
# For a smooth video run the driver with --vtu-every=10 (121 frames over the
# 1200-step run) rather than the default 50.

import argparse, os, sys, shutil, subprocess

ap = argparse.ArgumentParser()
ap.add_argument("--pvd",    required=True, help="driver --vtu-output prefix (reads PREFIX_umag.pvd)")
ap.add_argument("--width",  type=int, default=1920)
ap.add_argument("--height", type=int, default=540)
ap.add_argument("--fps",    type=int, default=12)
ap.add_argument("--field", default="u",
                help="which PVD timeline to render: u (default; the validated streamwise "
                     "component) or umag (CAUTION: contains v/w garbage from the "
                     "incremental-pressure artifact -- looks wrong even when u is correct)")
ap.add_argument("--u-range", default="0,1.5",
                help="color range LOW,HI (default 0,1.5 = up to U_max)")
ap.add_argument("--no-mp4", action="store_true", help="skip ffmpeg assembly")
args = ap.parse_args()

pvd_path = args.pvd + "_" + args.field + ".pvd"
if not os.path.exists(pvd_path):
    print(f"PVD not found: {pvd_path}", file=sys.stderr)
    sys.exit(1)

out_dir = args.pvd + "_movie"
os.makedirs(out_dir, exist_ok=True)
ulo, uhi = [float(x) for x in args.u_range.split(",")]

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

reader = PVDReader(FileName=pvd_path)
reader.UpdatePipeline()
times = reader.TimestepValues or [0.0]

view = CreateView("RenderView")
view.ViewSize = [args.width, args.height]
view.Background = [1, 1, 1]
view.OrientationAxesVisibility = 0
view.UseFXAA = 1

disp = Show(reader, view)
disp.Representation = "Surface"
ColorBy(disp, ("POINTS", args.field))
lut = GetColorTransferFunction(args.field)
lut.ApplyPreset("Turbo", True)
lut.RescaleTransferFunction(ulo, uhi)
bar = GetScalarBar(lut, view)
bar.Title = args.field
bar.ComponentTitle = ""
bar.TitleColor = [0, 0, 0]
bar.LabelColor = [0, 0, 0]
bar.WindowLocation = "Any Location"
bar.Position = [0.92, 0.25]
bar.ScalarBarLength = 0.5

annot = AnnotateTimeFilter(Input=reader)
annot.Format = "t = {time:.2f}"
adisp = Show(annot, view)
adisp.Color = [0, 0, 0]
adisp.WindowLocation = "Any Location"
adisp.Position = [0.02, 0.88]

# Top-down camera fitted to the 11:1 channel with parallel projection.
b = reader.GetDataInformation().GetBounds()
cx, cy = 0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3])
view.InteractionMode = "2D"
view.CameraPosition = [cx, cy, 10]
view.CameraFocalPoint = [cx, cy, 0]
view.CameraViewUp = [0, 1, 0]
view.CameraParallelProjection = 1
aspect = float(args.width) / float(args.height)
half_x = 0.5 * (b[1] - b[0]) * 1.04
half_y = 0.5 * (b[3] - b[2]) * 1.30
view.CameraParallelScale = max(half_y, half_x / aspect)

for i, t in enumerate(times):
    view.ViewTime = t
    Render(view)
    SaveScreenshot(os.path.join(out_dir, f"frame_{i:04d}.png"), view,
                   ImageResolution=[args.width, args.height])
    print(f"frame {i + 1}/{len(times)}  t={t:.2f}")

if not args.no_mp4:
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg:
        mp4 = args.pvd + "_movie.mp4"
        subprocess.run([ffmpeg, "-y", "-framerate", str(args.fps),
                        "-i", os.path.join(out_dir, "frame_%04d.png"),
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
                        mp4], check=True)
        print(f"wrote {mp4}")
    else:
        print("ffmpeg not in PATH; frames left in " + out_dir)
