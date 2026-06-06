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

# The GPU/viskores backend HANGS on StreamTracer/ParticleTracer over large tet
# meshes in headless ParaView (no progress, no crash). Force the CPU/serial VTK
# filter backend BEFORE importing paraview, unless the user opts back in with
# --gpu-filters. This makes the integrators slower-but-actually-finish.
if "--gpu-filters" not in sys.argv:
    os.environ.setdefault("VTK_DISABLE_VISKORES_OVERRIDES", "1")
    os.environ.setdefault("VTKM_DEVICE", "Serial")
    os.environ.setdefault("PARAVIEW_USE_ACCELERATED_FILTERS", "0")

from paraview.simple import *  # noqa: F401,F403,E402


def main():
    ap = argparse.ArgumentParser(description="Render a great MARS pump flow video")
    ap.add_argument("--pvd", required=True, help="path to <prefix>.pvd from the pump run")
    ap.add_argument("--out", default="", help="output mp4 (needs ffmpeg). Leave empty to write frames only "
                    "and print the ffmpeg command to run elsewhere (e.g. your laptop).")
    ap.add_argument("--frames-dir", default="pump_frames", help="PNG frames dir (always written)")
    ap.add_argument("--res", default="1920x1080", help="WxH render resolution")
    ap.add_argument("--field", default="velocity", help="vector field (default velocity)")
    ap.add_argument("--n-streamlines", type=int, default=900, help="streamline/particle seed count")
    ap.add_argument("--gpu-filters", action="store_true",
                    help="allow the GPU/viskores filter backend (default forces CPU/serial because the "
                         "GPU StreamTracer/ParticleTracer HANGS on large tet meshes in headless ParaView)")
    ap.add_argument("--pathlines", action="store_true",
                    help="PARTICLE TRACKS: release particles at the inlet and let the flow CARRY them "
                         "over time (looks like real fluid moving through). Needs the time series.")
    ap.add_argument("--orbit-deg", type=float, default=30.0,
                    help="total azimuth sweep over the run (gentle sway; 0 = fixed camera)")
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
    ap.add_argument("--no-shell", action="store_true",
                    help="skip the translucent geometry surface (heavy to render on big meshes)")
    ap.add_argument("--slice", action="store_true",
                    help="show a cutting-plane cross-section colored by |velocity| (MRI-like; "
                         "integration-free so it never hangs, unlike streamlines)")
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

    # ---- translucent geometry shell (heavy on big meshes -> opt-OUT with --no-shell) ----
    if not args.no_shell:
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

    # robust color range: evaluate ONCE at the last (developed) frame. (Do not
    # also UpdatePipeline with no time -- on big meshes that doubles the cost.)
    print("[render] computing |velocity| range on developed frame (this can take a while on big meshes)...")
    reader.UpdatePipeline(tvals[-1])
    mag.UpdatePipeline(tvals[-1])
    print("[render] |velocity| range computed")
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
    pathline_tracer = None  # set if --pathlines: advanced per-frame in the loop

    if args.pathlines:
        # PARTICLE TRACKS: seed particles at the inlet (peak-speed point) and let
        # ParticleTracer integrate them ACROSS timesteps -- they get carried by the
        # flow, so the animation shows real fluid moving through. Heavier than
        # streamlines; run under srun on a GPU node.
        seed_c = [cx, cy, cz]; seed_r = diag * 0.15
        if args.seed_inlet:
            try:
                seed_c = [float(x) for x in args.seed_inlet.split(",")]
            except Exception:
                pass
        else:
            try:
                reader.UpdatePipeline(tvals[-1]); mag.UpdatePipeline(tvals[-1])
                from paraview import servermanager as _sm
                data = _sm.Fetch(mag); pts = data.GetPoints(); arr = data.GetPointData().GetArray("speed")
                if pts and arr:
                    n = arr.GetNumberOfTuples(); step = max(1, n // 200000)
                    bi, bv = 0, -1.0
                    for i in range(0, n, step):
                        v = arr.GetValue(i)
                        if v > bv:
                            bv, bi = v, i
                    seed_c = list(pts.GetPoint(bi))
                    print("[render] seeding particles at peak-speed point %s (|u|=%.3g)" % (seed_c, bv))
            except Exception as e:  # noqa: BLE001
                print("[render] peak-seed failed (%s)" % e)
        try:
            tracer = ParticleTracer(Input=mag, SeedSource="Point Cloud")
            tracer.SelectInputVectors = ["POINTS", args.field]
            tracer.SeedSource.Center = seed_c
            tracer.SeedSource.Radius = seed_r
            tracer.SeedSource.NumberOfPoints = args.n_streamlines
            # render tracks as fading tails via a Tube on the tracer output
            ptubes = Tube(Input=tracer)
            ptubes.Radius = diag * 0.0016
            pathline_tracer = (tracer, ptubes)
            flow_props.append((ptubes, "speed"))
            print("[render] pathline (ParticleTracer) mode")
        except Exception as e:  # noqa: BLE001
            print("[render] ParticleTracer unavailable (%s); falling back to streamlines" % e)
            args.streamlines = True
            args.pathlines = False

    if args.streamlines and not args.pathlines:
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
    if args.slice and not args.streamlines and not args.pathlines:
        # SLICE: a cutting plane through the pump colored by |velocity| -- an
        # MRI-like cross-section of the flow. Integration-free, so it never hangs;
        # shows the velocity field on the cut and evolves with time.
        sl = Slice(Input=mag)
        try:
            sl.SliceType = "Plane"
            sl.SliceType.Origin = [cx, cy, cz]
            ext = [bounds[1]-bounds[0], bounds[3]-bounds[2], bounds[5]-bounds[4]]
            # cut perpendicular to the SHORTEST extent so the slice spans the long flow axis
            axis = ext.index(min(ext))
            nrm = [0, 0, 0]; nrm[axis] = 1
            sl.SliceType.Normal = nrm
        except Exception:
            pass
        sl.UpdatePipeline(tvals[-1])
        flow_props.append((sl, "speed"))

    if not args.slice and not args.streamlines and not args.pathlines:
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
        print("[render] building jet threshold (speed > %.3g)..." % lo)
        thr.UpdatePipeline(tvals[-1])
        print("[render] threshold built")
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
        # advance the particle tracer in time so particles get carried forward
        if pathline_tracer is not None:
            try:
                pathline_tracer[0].UpdatePipeline(t)
                pathline_tracer[1].UpdatePipeline(t)
            except Exception as e:  # noqa: BLE001
                print("  [warn] pathline update t=%.4g failed (%s)" % (t, e))
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
    import shutil
    out = args.out or "pump_flow.mp4"
    ff = shutil.which("ffmpeg")
    ffcmd = ("ffmpeg -y -framerate %d -i %s/frame_%%04d.png -c:v libx264 -pix_fmt yuv420p -crf 18 %s"
             % (args.fps, args.frames_dir, out))
    if args.out and ff:
        import subprocess
        print("[render] encoding: %s" % ffcmd)
        try:
            subprocess.run(ffcmd.split(), check=True)
            print("[render] wrote %s" % out)
        except Exception as e:  # noqa: BLE001
            print("[render] ffmpeg failed (%s). Run manually:\n  %s" % (e, ffcmd))
    else:
        # frames-only (no ffmpeg on the cluster): print the command to run on a
        # machine that has ffmpeg (e.g. your laptop, after rsync'ing the frames).
        print("[render] frames-only. Encode where ffmpeg exists (e.g. your laptop):")
        print("  %s" % ffcmd)


if __name__ == "__main__":
    main()
