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
                    help="use StreamTracer tubes (the real flow-path video; run under srun on a GPU "
                         "compute node so it doesn't crash; default is the fast jet-threshold)")
    ap.add_argument("--seed-inlet", default="",
                    help="seed streamlines at 'x,y,z' (the inlet); default = auto peak-speed point")
    ap.add_argument("--seed-radius", type=float, default=0.0,
                    help="streamline seed-cloud radius as a fraction of bbox diagonal (default auto)")
    ap.add_argument("--jet-frac", type=float, default=0.08,
                    help="threshold isovolume keeps speed > jet_frac*max (the visible jet); lower=more fluid shown")
    ap.add_argument("--color-frac", type=float, default=0.6,
                    help="color range clamped to [0, color_frac*max] so the jet stands out (not washed flat)")
    ap.add_argument("--glyphs", action="store_true",
                    help="add velocity-arrow glyphs (can crash on big tet meshes in headless PV; OFF by default)")
    ap.add_argument("--n-glyphs", type=int, default=4000, help="max velocity-arrow glyphs")
    ap.add_argument("--glyph-stride", type=int, default=20, help="glyph subsample stride (Every Nth Point)")
    ap.add_argument("--show-surface", dest="show_surface", action="store_true", default=True,
                    help="show translucent geometry shell (default on)")
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
    # The pump is a small jet into a big mostly-stagnant chamber: a plain clip is
    # ~all one color (most fluid is slow). To SHOW THE FLOW we combine, by default:
    #   (1) a THRESHOLD isovolume of the fast-moving fluid (speed > frac*max) --
    #       this is literally "where the jet is", and it GROWS over time as the jet
    #       penetrates, so the animation shows flow developing;
    #   (2) velocity GLYPHS (arrows) sampled through the volume, colored by speed,
    #       showing direction + magnitude (the actual flow vectors).
    # Color is clamped to a meaningful band (not [0,max]) so slow fluid isn't flat.
    # --streamlines switches to StreamTracer tubes (prettier, but can hang on big
    # tet meshes in headless ParaView -- run it on a GPU compute node).
    flow_props = []  # (proxy, colorByName) to Show

    if args.streamlines:
        # Seed center: by default the high-speed (inlet/jet) region of the LAST
        # frame, so lines trace the actual pumped path; override with --seed-inlet.
        seed_c = [cx, cy, cz]
        seed_r = diag * (args.seed_radius if args.seed_radius > 0 else 0.45)
        if args.seed_inlet:
            try:
                seed_c = [float(x) for x in args.seed_inlet.split(",")]
            except Exception:
                print("[render] bad --seed-inlet '%s', using center" % args.seed_inlet)
        else:
            # find the location of peak speed on the developed frame to seed there
            try:
                reader.UpdatePipeline(tvals[-1]); mag.UpdatePipeline(tvals[-1])
                from paraview import servermanager as _sm
                data = _sm.Fetch(mag)
                pts = data.GetPoints(); arr = data.GetPointData().GetArray("speed")
                if pts and arr:
                    best_i, best_v = 0, -1.0
                    n = arr.GetNumberOfTuples()
                    step = max(1, n // 200000)  # subsample for speed
                    for i in range(0, n, step):
                        v = arr.GetValue(i)
                        if v > best_v:
                            best_v, best_i = v, i
                    seed_c = list(pts.GetPoint(best_i))
                    seed_r = diag * 0.18  # tighter cloud around the jet
                    print("[render] seeding streamlines at peak-speed point %s (|u|=%.3g)" % (seed_c, best_v))
            except Exception as e:  # noqa: BLE001
                print("[render] peak-seed failed (%s), using center cloud" % e)

        stream = StreamTracer(Input=mag, SeedType="Point Cloud")
        stream.Vectors = ["POINTS", args.field]
        stream.IntegrationDirection = "BOTH"
        try:
            stream.IntegratorType = "Runge-Kutta 4-5"
        except Exception:
            pass
        stream.MaximumStreamlineLength = diag * 5.0
        stream.SeedType.Center = seed_c
        stream.SeedType.Radius = seed_r
        stream.SeedType.NumberOfPoints = args.n_streamlines
        stream.UpdatePipeline(tvals[-1])
        tubes = Tube(Input=stream)
        tubes.Radius = diag * 0.0018
        tubes.Capping = 1
        tubes.UpdatePipeline(tvals[-1])
        flow_props.append((tubes, "speed"))
    else:
        # (1) fast-jet isovolume: keep cells where speed exceeds a fraction of the
        # peak. This is the moving jet; it must re-evaluate per timestep so it
        # grows with the flow.
        thr = Threshold(Input=mag)
        thr.Scalars = ["POINTS", "speed"]
        lo = max(args.jet_frac * smax, 1e-12)
        for setter in ("LowerThreshold", "UpperThreshold"):
            try:
                setattr(thr, setter, lo if setter == "LowerThreshold" else smax * 1.05)
            except Exception:
                pass
        # ParaView >= 5.10 uses a single ThresholdMethod/range; guard both APIs
        try:
            thr.ThresholdMethod = "Above Lower Threshold"
        except Exception:
            pass
        thr.UpdatePipeline(tvals[-1])
        flow_props.append((thr, "speed"))

        # (2) OPTIONAL velocity glyphs (--glyphs). The Glyph filter with spatial
        # sampling on a large tet mesh can hard-crash (core dump) in headless
        # ParaView when re-run per frame, so it is OFF by default. The threshold
        # isovolume alone already shows the jet developing.
        if args.glyphs:
            glyph = Glyph(Input=mag, GlyphType="Arrow")
            glyph.OrientationArray = ["POINTS", args.field]
            glyph.ScaleArray = ["POINTS", "speed"]
            try:
                glyph.ScaleFactor = diag * 0.04
            except Exception:
                pass
            for mode in ("Every Nth Point", "Uniform Spatial Distribution"):
                try:
                    glyph.GlyphMode = mode
                    break
                except Exception:
                    continue
            try:
                glyph.MaximumNumberOfSamplePoints = args.n_glyphs
            except Exception:
                pass
            try:
                glyph.Stride = max(1, args.glyph_stride)
            except Exception:
                pass
            glyph.UpdatePipeline(tvals[-1])
            flow_props.append((glyph, "speed"))

    # clamp color to a readable band: [0, jet_color_frac*max] so the jet stands out
    cmax = max(args.color_frac * smax, smin + 1e-9)
    lut.RescaleTransferFunction(smin, cmax)

    sdisp = None
    for proxy, cname in flow_props:
        d = Show(proxy, view)
        d.Representation = "Surface"
        try:
            ColorBy(d, ("POINTS", cname))
        except Exception:
            pass
        d.Specular = 0.3
        sdisp = d
    if sdisp is not None:
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

    def save_frame(idx):
        # Guard each frame: one bad timestep/render must NOT abort the whole run.
        fn = os.path.join(args.frames_dir, "frame_%04d.png" % idx)
        try:
            Render(view)
            SaveScreenshot(fn, view, ImageResolution=[W, H])
            return True
        except Exception as e:  # noqa: BLE001
            print("  [warn] frame %d failed (%s) -- skipping" % (idx, e))
            return False

    fidx = 0
    written = 0
    for t in tvals:
        view.ViewTime = t
        try:
            reader.UpdatePipeline(t)
        except Exception as e:  # noqa: BLE001
            print("  [warn] UpdatePipeline t=%.4g failed (%s)" % (t, e))
        cam.Azimuth(az_per_frame)
        if save_frame(fidx):
            written += 1
        if fidx % 10 == 0:
            print("  frame %d/%d  t=%.4g" % (fidx + 1, total, t))
        fidx += 1
    # hold on the developed state, keep orbiting
    for _ in range(max(0, args.hold_frames)):
        cam.Azimuth(az_per_frame)
        if save_frame(fidx):
            written += 1
        fidx += 1

    print("[render] wrote %d/%d PNG frames to %s" % (written, total, args.frames_dir))

    # Encode FROM the PNG frames with ffmpeg -- do NOT use SaveAnimation (it
    # re-renders every frame and would crash the same way / waste the work).
    if args.out:
        import shutil
        import subprocess
        ff = shutil.which("ffmpeg")
        ffcmd = ("%s -y -framerate %d -i %s/frame_%%04d.png -c:v libx264 -pix_fmt yuv420p -crf 18 %s"
                 % (ff or "ffmpeg", args.fps, args.frames_dir, args.out))
        if ff:
            print("[render] encoding: %s" % ffcmd)
            try:
                subprocess.run(ffcmd.split(), check=True)
                print("[render] wrote %s" % args.out)
            except Exception as e:  # noqa: BLE001
                print("[render] ffmpeg failed (%s). Run manually:\n  %s" % (e, ffcmd))
        else:
            print("[render] ffmpeg not found. Encode the frames yourself:\n  %s" % ffcmd)


if __name__ == "__main__":
    main()
