#!/usr/bin/env pvbatch
# Render the Poiseuille channel validation animation from mars_poiseuille_flow
# VTU output. Pump-renderer style (render_pump_4panel.py): per-view PNGs
# composited with PIL into a labeled multi-panel frame + header, then ffmpeg.
#
# Default layout (--layout panels):
#   TOP    : u colormap of the whole channel (top-down x-y view)
#   BOT-L  : velocity profile u(y) at the probe station vs the EXACT parabola
#            u = G/(2mu) y(h-y) -- the validation figure, animated
#   BOT-R  : centerline u(x) -- entrance development toward U_max
# --layout simple gives the old single-panel field view.
#
# Input: PREFIX_<field>.pvd (default field u -- the validated component;
# umag/p carry the incremental-pressure artifact and render as garbage).
#
# Usage:
#   pvbatch scripts/render_poiseuille.py --pvd poiseuille_final
#   /Applications/ParaView-6.1.1.app/Contents/bin/pvbatch \
#       scripts/render_poiseuille.py --pvd /tmp/poise_render/poiseuille_final

import argparse, os, sys, shutil, subprocess

ap = argparse.ArgumentParser()
ap.add_argument("--pvd",    required=True, help="driver --vtu-output prefix")
ap.add_argument("--field",  default="u",
                help="PVD timeline to render (default u; umag/p are artifact-contaminated)")
ap.add_argument("--layout", default="panels", choices=["panels", "simple"])
ap.add_argument("--width",  type=int, default=1920)
ap.add_argument("--height", type=int, default=1080)
ap.add_argument("--fps",    type=int, default=12)
ap.add_argument("--frames", type=int, default=0, help="0 = all timesteps")
ap.add_argument("--u-range", default="0,1.5", help="color/axis range LOW,HI")
ap.add_argument("--umax",   type=float, default=1.5, help="exact centerline U_max")
ap.add_argument("--profile-x", type=float, default=9.0, help="profile station x")
ap.add_argument("--no-mp4", action="store_true")
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
if args.frames > 0:
    times = times[:args.frames]

b = reader.GetDataInformation().GetBounds()
cx, cy, cz = 0.5*(b[0]+b[1]), 0.5*(b[2]+b[3]), 0.5*(b[4]+b[5])
h = b[3] - b[2]
y0 = b[2]

# --- TOP: field view ----------------------------------------------------------
def make_render_view(w, hgt):
    rv = CreateView("RenderView")
    rv.ViewSize = [w, hgt]
    rv.Background = [1, 1, 1]
    try: rv.UseColorPaletteForBackground = 0
    except Exception: pass
    try: rv.BackgroundColorMode = "Single Color"
    except Exception: pass
    rv.OrientationAxesVisibility = 0
    try: rv.UseFXAA = 1
    except Exception: pass
    return rv

panel_pad = 10
top_h    = int(args.height * 0.40)
bot_h    = args.height - top_h - 56      # 56 px header strip
half_w   = args.width // 2

simple_mode = (args.layout == "simple")
top_w  = args.width if simple_mode else args.width - 2*panel_pad
view_T = make_render_view(top_w, (args.height if simple_mode else top_h) - 2*panel_pad)

disp = Show(reader, view_T)
disp.Representation = "Surface"
ColorBy(disp, ("POINTS", args.field))
lut = GetColorTransferFunction(args.field)
lut.ApplyPreset("Turbo", True)
lut.RescaleTransferFunction(ulo, uhi)
bar = GetScalarBar(lut, view_T)
bar.Title = args.field
bar.ComponentTitle = ""
bar.TitleColor = [0, 0, 0]
bar.LabelColor = [0, 0, 0]
bar.WindowLocation = "Any Location"
bar.Position = [0.93, 0.20]
bar.ScalarBarLength = 0.6

view_T.InteractionMode = "2D"
view_T.CameraPosition = [cx, cy, 10]
view_T.CameraFocalPoint = [cx, cy, 0]
view_T.CameraViewUp = [0, 1, 0]
view_T.CameraParallelProjection = 1
aspect = float(view_T.ViewSize[0]) / float(view_T.ViewSize[1])
view_T.CameraParallelScale = max(0.5*(b[3]-b[2])*1.12, 0.5*(b[1]-b[0])*1.02 / aspect)

# --- Charts (panels layout only) ----------------------------------------------
if not simple_mode:
    # Profile at the probe station: vertical line through the channel.
    profLine = PlotOverLine(Input=reader)
    for attr, val in (("Point1", [args.profile_x, b[2], cz]),
                      ("Point2", [args.profile_x, b[3], cz])):
        try: setattr(profLine, attr, val)
        except Exception: setattr(profLine.Source, attr, val)
    try: profLine.Resolution = 200
    except Exception: pass

    # Exact parabola on the same line: u = Umax * 4 (y-y0)(h-(y-y0)) / h^2.
    exactCalc = Calculator(Input=profLine)
    exactCalc.ResultArrayName = "u_exact"
    exactCalc.Function = (f"{args.umax}*4*(coordsY-{y0})*({h}-(coordsY-{y0}))/({h}*{h})")

    view_BL = CreateView("XYChartView")
    view_BL.ViewSize = [half_w - 2*panel_pad, bot_h - 2*panel_pad]
    view_BL.ChartTitle = "profile at x=%g vs exact parabola" % args.profile_x
    view_BL.BottomAxisTitle = "u"
    view_BL.LeftAxisTitle = "y"
    dBL_exact = Show(exactCalc, view_BL)
    dBL_exact.UseIndexForXAxis = 0
    dBL_exact.XArrayName = "u_exact"
    dBL_exact.SeriesVisibility = ["Points_Y", "1"]
    dBL_exact.SeriesColor = ["Points_Y", "0.12", "0.38", "0.85"]
    dBL_exact.SeriesLineThickness = ["Points_Y", "4"]
    dBL_exact.SeriesLabel = ["Points_Y", "exact"]
    dBL_sol = Show(profLine, view_BL)
    dBL_sol.UseIndexForXAxis = 0
    dBL_sol.XArrayName = args.field
    dBL_sol.SeriesVisibility = ["Points_Y", "1"]
    dBL_sol.SeriesColor = ["Points_Y", "0.85", "0.10", "0.10"]
    dBL_sol.SeriesLineThickness = ["Points_Y", "2"]
    dBL_sol.SeriesLineStyle = ["Points_Y", "0"]
    dBL_sol.SeriesMarkerStyle = ["Points_Y", "4"]
    dBL_sol.SeriesMarkerSize = ["Points_Y", "7"]
    dBL_sol.SeriesLabel = ["Points_Y", "MARS (solved)"]
    view_BL.BottomAxisUseCustomRange = 1
    view_BL.BottomAxisRangeMinimum = ulo - 0.05
    view_BL.BottomAxisRangeMaximum = uhi * 1.08
    view_BL.LeftAxisUseCustomRange = 1
    view_BL.LeftAxisRangeMinimum = b[2] - 0.03*h
    view_BL.LeftAxisRangeMaximum = b[3] + 0.03*h

    # Centerline development u(x) at mid-height.
    centerLine = PlotOverLine(Input=reader)
    ymid = 0.5 * (b[2] + b[3])
    for attr, val in (("Point1", [b[0], ymid, cz]),
                      ("Point2", [b[1], ymid, cz])):
        try: setattr(centerLine, attr, val)
        except Exception: setattr(centerLine.Source, attr, val)
    try: centerLine.Resolution = 300
    except Exception: pass

    view_BR = CreateView("XYChartView")
    view_BR.ViewSize = [half_w - 2*panel_pad, bot_h - 2*panel_pad]
    view_BR.ChartTitle = "centerline u(x): development toward U_max=%g" % args.umax
    view_BR.BottomAxisTitle = "x"
    view_BR.LeftAxisTitle = "u centerline"
    dBR = Show(centerLine, view_BR)
    dBR.UseIndexForXAxis = 0
    dBR.XArrayName = "Points_X"
    dBR.SeriesVisibility = [args.field, "1"]
    dBR.SeriesColor = [args.field, "0.85", "0.10", "0.10"]
    dBR.SeriesLineThickness = [args.field, "3"]
    dBR.SeriesLabel = [args.field, "u(x, y=H/2)"]
    view_BR.LeftAxisUseCustomRange = 1
    view_BR.LeftAxisRangeMinimum = 0
    view_BR.LeftAxisRangeMaximum = uhi * 1.08
    view_BR.BottomAxisUseCustomRange = 1
    view_BR.BottomAxisRangeMinimum = b[0]
    view_BR.BottomAxisRangeMaximum = b[1]

# --- Compositing ----------------------------------------------------------------
try:
    from PIL import Image, ImageDraw, ImageFont
    have_pil = True
except ImportError:
    have_pil = False
    if not simple_mode:
        print("WARNING: PIL unavailable -- falling back to simple layout")
        simple_mode = True

def _font(size):
    for path in ("/System/Library/Fonts/Helvetica.ttc",
                 "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
                 "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf"):
        try: return ImageFont.truetype(path, size)
        except Exception: continue
    return ImageFont.load_default()

font_big = _font(34) if have_pil else None

for i, t in enumerate(times):
    if simple_mode:
        view_T.ViewTime = t
        Render(view_T)
        SaveScreenshot(os.path.join(out_dir, f"frame_{i:04d}.png"), view_T,
                       ImageResolution=list(view_T.ViewSize))
    else:
        views = [("T", view_T), ("BL", view_BL), ("BR", view_BR)]
        paths = {}
        for name, v in views:
            v.ViewTime = t
            try: v.Update()
            except Exception: pass
            p = os.path.join(out_dir, f"_{name}_{i:04d}.png")
            SaveScreenshot(p, v, ImageResolution=list(v.ViewSize))
            paths[name] = p
        canvas = Image.new("RGB", (args.width, args.height), (255, 255, 255))
        imgs = {k: Image.open(p) for k, p in paths.items()}
        canvas.paste(imgs["T"],  (panel_pad, 56 + panel_pad))
        canvas.paste(imgs["BL"], (panel_pad, 56 + top_h + panel_pad))
        canvas.paste(imgs["BR"], (half_w + panel_pad, 56 + top_h + panel_pad))
        draw = ImageDraw.Draw(canvas)
        draw.rectangle([(0, 0), (args.width, 56)], fill=(25, 25, 25))
        draw.text((16, 9), f"MARS Poiseuille validation    t = {t:6.2f}",
                  fill=(255, 255, 255), font=font_big)
        canvas.save(os.path.join(out_dir, f"frame_{i:04d}.png"))
        for p in paths.values():
            try: os.remove(p)
            except OSError: pass
    print(f"frame {i+1}/{len(times)}  t={t:.2f}")

if not args.no_mp4:
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg:
        mp4 = args.pvd + "_movie.mp4"
        subprocess.run([ffmpeg, "-y", "-loglevel", "error",
                        "-framerate", str(args.fps),
                        "-i", os.path.join(out_dir, "frame_%04d.png"),
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
                        mp4], check=True)
        print(f"wrote {mp4}")
    else:
        print("ffmpeg not in PATH; frames left in " + out_dir)
