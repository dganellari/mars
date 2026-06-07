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
                         "over time, rendered as little spheres -- looks like real liquid streaming "
                         "through. Needs the time series. Slow on big meshes; let it run.")
    ap.add_argument("--no-outline", action="store_true",
                    help="never draw the bounding-box outline (the white box)")
    ap.add_argument("--particle-size", type=float, default=0.003,
                    help="particle sphere size as a fraction of the bbox diagonal (smaller=more liquid-like)")
    ap.add_argument("--reinject", type=int, default=2,
                    help="re-inject particles every N steps (smaller=denser continuous stream)")
    ap.add_argument("--seed-whole", action="store_true",
                    help="seed particles across the WHOLE domain (not just the jet) so both sides fill; "
                         "use if particles only appear on one side")
    ap.add_argument("--mark-io", action="store_true",
                    help="auto-detect inlet (most inward boundary flow) + outlet (most outward) and mark "
                         "them with a green/red sphere + a legend, so the intended flow path is clear")
    ap.add_argument("--marker-size", type=float, default=0.025,
                    help="inlet/outlet marker sphere size as a fraction of the bbox diagonal (default 0.025)")
    ap.add_argument("--orbit-deg", type=float, default=0.0,
                    help="total azimuth sweep over the run (default 0 = fixed camera; set e.g. 30 for a gentle sway)")
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--hold-frames", type=int, default=20,
                    help="extra orbit-only frames at the end on the final (developed) state")
    ap.add_argument("--preset", default="",
                    help="color preset name; default is a custom light-blue->deep-blue water ramp "
                         "(no red). Pass a name to override, e.g. 'Cool to Warm', 'Viridis (matplotlib)'")
    ap.add_argument("--streamlines", action="store_true",
                    help="use StreamTracer tubes (the real flow-path video; run under srun on a GPU "
                         "compute node so it doesn't crash; default is the fast jet-threshold)")
    ap.add_argument("--seed-inlet", default="",
                    help="seed streamlines at 'x,y,z' (the inlet); default = auto peak-speed point")
    ap.add_argument("--seed-radius", type=float, default=0.0,
                    help="streamline seed-cloud radius as a fraction of bbox diagonal (default auto)")
    ap.add_argument("--jet-frac", type=float, default=0.08,
                    help="threshold isovolume keeps speed > jet_frac*max (the visible jet); lower=more fluid shown")
    ap.add_argument("--color-frac", type=float, default=1.0,
                    help="color range = [0, color_frac*max]. 1.0 uses the full speed range so you see the "
                         "blue->red gradient across the whole flow; lower saturates the fast jet to red sooner")
    ap.add_argument("--shell", action="store_true",
                    help="add the full translucent geometry surface (prettier but HEAVY on big meshes; "
                         "default is just a cheap always-visible bounding outline)")
    ap.add_argument("--shell-opacity", type=float, default=0.08,
                    help="translucent pump-body opacity (0..1; lower = more see-through, default 0.08)")
    ap.add_argument("--wireframe", action="store_true",
                    help="draw the pump body as a wireframe cage instead of a translucent surface "
                         "(does not block the interior at all -- particles fully visible)")
    ap.add_argument("--zoom", type=float, default=2.0,
                    help="camera zoom: Dolly factor >1 moves IN so the pump fills the frame (default 2.0)")
    ap.add_argument("--frame-stride", type=int, default=1,
                    help="render every Nth timestep (e.g. 4 = 4x fewer frames, 4x faster; the video "
                         "is still smooth). Use this if the render hits the job time limit.")
    ap.add_argument("--max-frames", type=int, default=0,
                    help="cap the number of frames rendered (0 = no cap); another way to fit the time budget")
    ap.add_argument("--volume", action="store_true",
                    help="color the whole pump SURFACE/volume by |velocity| (the 'colored tetrahedra' "
                         "look). Integration-free, fast, robust. Use this for the second video "
                         "(--pathlines for the first); run as a separate job.")
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
    _ntotal = len(tvals)
    if args.frame_stride > 1:
        # always keep the last frame (the developed state) even if stride skips it
        last = tvals[-1]
        tvals = tvals[::args.frame_stride]
        if tvals[-1] != last:
            tvals.append(last)
    if args.max_frames > 0 and len(tvals) > args.max_frames:
        # subsample evenly down to max_frames
        step = len(tvals) / float(args.max_frames)
        tvals = [tvals[int(i * step)] for i in range(args.max_frames)]
    if len(tvals) != _ntotal:
        print("[render] rendering %d of %d timesteps (stride=%d, max=%s)"
              % (len(tvals), _ntotal, args.frame_stride, args.max_frames or "none"))

    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = [W, H]
    # lighter slate-blue gradient (pure-black hid the pump when flow was low)
    view.Background = [0.16, 0.18, 0.22]
    view.Background2 = [0.30, 0.34, 0.40]
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

    # ---- always-visible pump context ----
    # A bounding-box OUTLINE shows the domain extent even at t=0. It's only useful
    # WITHOUT the translucent surface; with --shell the body already shows the
    # geometry, so skip the outline (it reads as a distracting white box).
    if not args.shell and not args.no_outline:
        try:
            outline = Outline(Input=reader)
            odisp = Show(outline, view)
            odisp.DiffuseColor = [0.7, 0.75, 0.82]
            odisp.AmbientColor = [0.7, 0.75, 0.82]
            try:
                odisp.LineWidth = 2
            except Exception:
                pass
        except Exception as e:  # noqa: BLE001
            print("[render] outline skipped: %s" % e)

    if args.shell:
        surf = Show(reader, view)
        if args.wireframe:
            # wireframe shows the pump shape WITHOUT blocking the interior at all
            # -> particles fully visible inside; reads as a glass cage.
            surf.Representation = "Wireframe"
            surf.Opacity = max(args.shell_opacity, 0.25)
            surf.LineWidth = 1
        else:
            surf.Representation = "Surface"
            surf.Opacity = args.shell_opacity
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
    if args.preset:
        # explicit preset requested
        for _p in (args.preset, "Cool to Warm", "Viridis (matplotlib)", "viridis"):
            try:
                lut.ApplyPreset(_p, True)
                print("[render] color preset: %s" % _p)
                break
            except Exception:
                continue
        lut.RescaleTransferFunction(smin, smax)
    else:
        # custom WATER ramp: light/baby blue (slow) -> deep blue (fast). No red.
        # RGBPoints = [value, r, g, b, ...] across [smin, smax].
        lo, hi = smin, (smax if smax > smin else smin + 1.0)
        def _mix(a, b, t):
            return a + (b - a) * t
        # baby-blue DOMINANT: most of the range stays light/baby blue; only the
        # very fastest fluid deepens to a medium sky blue (never dark navy).
        stops = [(0.0, 0.86, 0.94, 1.00),   # very light baby blue (slow)
                 (0.7, 0.70, 0.86, 0.98),   # baby blue (most of the field)
                 (1.0, 0.40, 0.68, 0.92)]   # medium sky blue (fastest)
        rgb_points = []
        for frac, r, g, b in stops:
            rgb_points += [_mix(lo, hi, frac), r, g, b]
        try:
            lut.RGBPoints = rgb_points
            lut.ColorSpace = "RGB"
            print("[render] custom water-blue colormap (light->deep blue, no red)")
        except Exception as e:  # noqa: BLE001
            print("[render] custom colormap failed (%s), falling back to Cool to Warm" % e)
            try:
                lut.ApplyPreset("Cool to Warm", True)
            except Exception:
                pass
            lut.RescaleTransferFunction(lo, hi)
    try:
        lut.NanColor = [0.2, 0.2, 0.2]
    except Exception:
        pass

    # ---- auto-detect inlet & outlet from the boundary flow (NDA-safe: all done
    # inside ParaView, coords never leave the render) and mark them ----
    # Inlet  = boundary point with the most INWARD flow (velocity . inward normal).
    # Outlet = boundary point with the most OUTWARD flow.
    if args.mark_io:
        try:
            surf = ExtractSurface(Input=mag)
            # the surface-normals filter has different names across PV builds
            # (PV 6.1 does NOT have GenerateSurfaceNormals); try the known ones.
            import paraview.simple as _ps
            norm = None
            for fname in ("GenerateSurfaceNormals", "SurfaceNormals",
                          "PolyDataNormals", "Normals", "ComputeNormals"):
                ctor = getattr(_ps, fname, None)
                if ctor is None:
                    continue
                try:
                    norm = ctor(Input=surf)
                    print("[render] surface normals via %s" % fname)
                    break
                except Exception:
                    norm = None
            if norm is None:
                # last resort: the surface may already carry Normals, or we
                # compute v.n with a fabricated outward direction (skip normals).
                norm = surf
                print("[render] no surface-normals filter found; using surface as-is")
            try:
                norm.ComputePointNormals = 1
            except Exception:
                pass
            # materialize the normals BEFORE the Calculator parses the expression
            norm.UpdatePipeline(tvals[-1])
            # auto-discover the point arrays present (normals name varies by build)
            pdn = norm.PointData
            pnames = [pdn.GetArray(i).GetName() for i in range(pdn.GetNumberOfArrays())]
            ncand = [n for n in pnames if n.lower() in ("normals", "surface normals")] \
                or [n for n in pnames if "normal" in n.lower()]
            nname = ncand[0] if ncand else "Normals"
            print("[render] surface point arrays: %s -> normals: %s" % (pnames, nname))
            if args.field not in pnames:
                raise RuntimeError("'%s' not a point array on the surface (%s)" % (args.field, pnames))
            # velocity . normal: outward>0 (outlet), inward<0 (inlet)
            dotc = Calculator(Input=norm)
            try:
                dotc.AttributeType = "Point Data"
            except Exception:
                pass
            dotc.ResultArrayName = "vdotn"
            dotc.Function = "dot(%s, %s)" % (args.field, nname)
            dotc.UpdatePipeline(tvals[-1])
            if dotc.PointData.GetArray("vdotn") is None:
                raise RuntimeError("vdotn not computed")
            from paraview import servermanager as _sm
            sd = _sm.Fetch(dotc)
            spts = sd.GetPoints()
            vn = sd.GetPointData().GetArray("vdotn")
            if spts and vn:
                nn = vn.GetNumberOfTuples(); st = max(1, nn // 300000)
                out_i, out_v, in_i, in_v = -1, -1e30, -1, 1e30
                for i in range(0, nn, st):
                    v = vn.GetValue(i)
                    if v > out_v:
                        out_v, out_i = v, i
                    if v < in_v:
                        in_v, in_i = v, i
                io_marks = []
                # colors chosen to stand out against blue water + each other, no plain red
                if in_i >= 0:
                    io_marks.append(("INLET", list(spts.GetPoint(in_i)), [1.0, 0.85, 0.1]))   # gold
                if out_i >= 0:
                    io_marks.append(("OUTLET", list(spts.GetPoint(out_i)), [1.0, 0.2, 0.8]))   # magenta
                print("[render] inlet found=%s (gold marker), outlet found=%s (magenta marker)"
                      % (in_i >= 0, out_i >= 0))
                for label, pos, col in io_marks:
                    ball = Sphere()
                    ball.Center = pos
                    ball.Radius = diag * args.marker_size
                    ball.ThetaResolution = 24; ball.PhiResolution = 24
                    bdisp = Show(ball, view)
                    bdisp.DiffuseColor = col; bdisp.AmbientColor = col
                    bdisp.Specular = 0.5
                    # NO floating 3D text (it rendered as huge letters and the
                    # auto-detected point can land on a spurious mesh corner).
                    # The small corner legend below carries the color key.
                # small fixed corner legend (the color key)
                try:
                    legend = Text()
                    legend.Text = "gold sphere = INLET    magenta sphere = OUTLET"
                    ld = Show(legend, view)
                    ld.Color = [1, 1, 1]; ld.FontSize = 16
                    try:
                        ld.WindowLocation = "Upper Right Corner"
                    except Exception:
                        pass
                except Exception:
                    pass
                print("[render] marked inlet (green sphere) + outlet (red sphere) + 3D labels")
        except Exception as e:  # noqa: BLE001
            print("[render] inlet/outlet marking skipped (%s)" % e)

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
        # Seed a BROAD cloud around the jet so particles enter across the inlet and
        # get carried through the WHOLE pump, not just a tight knot at one point.
        seed_c = [cx, cy, cz]
        seed_r = diag * (args.seed_radius if args.seed_radius > 0 else 0.4)
        if args.seed_whole:
            # cover the whole domain so particles enter everywhere -> both sides
            # fill (the bbox center + half-diagonal radius; no mesh reading).
            seed_c = [cx, cy, cz]
            seed_r = diag * 0.5
            print("[render] seeding particles across the WHOLE domain")
        elif args.seed_inlet:
            try:
                seed_c = [float(x) for x in args.seed_inlet.split(",")]
                seed_r = diag * (args.seed_radius if args.seed_radius > 0 else 0.12)
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
                    # cloud big enough to span the inlet jet, not a single point
                    seed_r = diag * (args.seed_radius if args.seed_radius > 0 else 0.3)
                    print("[render] seeding particles around jet at %s (|u|=%.3g), radius=%.3g"
                          % (seed_c, bv, seed_r))
            except Exception as e:  # noqa: BLE001
                print("[render] peak-seed failed (%s)" % e)
        try:
            # SeedSource is a PROXY OBJECT, not a string. Create a PointSource and
            # assign it (ParaView 5.12+/6.x). This is the GUI's "Point Cloud" seed.
            seed = PointSource()
            seed.Center = seed_c
            seed.Radius = seed_r
            seed.NumberOfPoints = args.n_streamlines

            tracer = ParticleTracer(Input=mag, SeedSource=seed)
            # array used for integration (list form; bare-string fallback below)
            try:
                tracer.SelectInputVectors = ["POINTS", args.field]
            except Exception:
                tracer.SelectInputVectors = args.field
            # re-inject so the inlet keeps "spraying" a continuous stream
            for attr, val in (("ForceReinjectionEveryNSteps", max(1, args.reinject)),
                              ("StaticSeeds", 1), ("ComputeVorticity", 0)):
                try:
                    setattr(tracer, attr, val)
                except Exception:
                    pass

            # Render advected particles as little SPHERES -> liquid droplets.
            blobs = Glyph(Input=tracer, GlyphType="Sphere")
            try:
                blobs.GlyphMode = "All Points"      # one sphere per particle
            except Exception:
                pass
            try:
                blobs.ScaleArray = ["POINTS", "No scale array"]  # uniform size
            except Exception:
                pass
            try:
                blobs.ScaleFactor = diag * args.particle_size
            except Exception:
                pass
            pathline_tracer = (tracer, blobs)
            flow_props.append((blobs, "speed"))
            print("[render] pathline/particle mode (PointSource seed + sphere glyphs)")
        except Exception as e:  # noqa: BLE001
            print("[render] ParticleTracer setup failed (%s); falling back to streamlines" % e)
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
    if args.volume and not args.streamlines and not args.pathlines:
        # VOLUME: the whole pump SURFACE colored by |velocity| -- the "colored
        # tetrahedra" look. Integration-free (fast, robust). Shows the velocity
        # field over the geometry, evolving with time. This is the second video
        # type (run a separate job with --volume; --pathlines for the first).
        vsurf = ExtractSurface(Input=mag)
        vsurf.UpdatePipeline(tvals[-1])
        flow_props.append((vsurf, "speed"))

    if args.slice and not args.streamlines and not args.pathlines and not args.volume:
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

    if not args.slice and not args.streamlines and not args.pathlines and not args.volume:
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
    # Reset to the PUMP DATA bounds only (not the markers/legend), so stray
    # markers don't make ResetCamera zoom out and shrink the pump.
    try:
        ResetCamera(view, bounds[0], bounds[1], bounds[2], bounds[3], bounds[4], bounds[5])
    except Exception:
        ResetCamera(view)
    cam = GetActiveCamera()
    cam.Elevation(16); cam.Azimuth(25)
    view.CameraFocalPoint = [cx, cy, cz]
    # zoom IN so the pump fills the frame (Dolly>1 moves the camera TOWARD the
    # focal point; the earlier 0.85 pulled it back and made the pump look tiny).
    cam.Dolly(args.zoom)
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
