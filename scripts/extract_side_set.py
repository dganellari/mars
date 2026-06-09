#!/usr/bin/env pvbatch
# Extract ONE side set from a MARS pump VTU/PVD into its own .vtu you can open
# directly in ParaView (no Threshold needed). The MARS pump writer tags every
# node with side_set_id (1,2,3,... per Exodus side set; the run prints the
# legend "N = [name]"). This reads that field and writes only the nodes/cells
# of the requested set.
#
# Run on the cluster (headless ParaView):
#   pvbatch scripts/extract_side_set.py --pvd pump_big/pump_big.pvd --id 1 --out inlet.vtu
#   pvbatch scripts/extract_side_set.py --pvd pump_big/pump_big.pvd --id 2 --out outlet.vtu
# Then open inlet.vtu / outlet.vtu directly in the ParaView GUI.
#
# Needs a VTU written by a binary that includes the side_set_id field (commit
# e668200+). If side_set_id is absent the script says so and exits.

import argparse
import sys

from paraview.simple import *  # noqa: F401,F403,E402


def main():
    ap = argparse.ArgumentParser(description="Extract one side set from a MARS VTU by side_set_id")
    ap.add_argument("--pvd", required=True, help="path to <prefix>.pvd (or a single .vtu/.pvtu)")
    ap.add_argument("--id", type=int, required=True,
                    help="side_set_id to extract (the integer from the run's 'N = [name]' legend)")
    ap.add_argument("--out", required=True, help="output .vtu path to write the isolated set")
    ap.add_argument("--field", default="side_set_id", help="tag field name (default side_set_id)")
    ap.add_argument("--last", action="store_true",
                    help="extract from the LAST timestep only (default: all timesteps)")
    args = ap.parse_args()

    reader = OpenDataFile(args.pvd)
    if reader is None:
        print("[extract] could not open %s" % args.pvd)
        sys.exit(1)
    reader.UpdatePipeline()

    # Confirm the tag field is present (point data). If not, the VTU predates the
    # side_set_id feature -> the user needs to re-run with a current binary.
    info = reader.GetPointDataInformation()
    names = [info.GetArray(i).GetName() for i in range(info.GetNumberOfArrays())]
    if args.field not in names:
        print("[extract] field '%s' NOT in this file. Point arrays present: %s" % (args.field, names))
        print("[extract] -> re-run the simulation with a binary that writes %s (commit e668200+)." % args.field)
        sys.exit(2)

    tvals = list(reader.TimestepValues) if hasattr(reader, "TimestepValues") and reader.TimestepValues else []

    # Threshold on the tag. side_set_id is a float (1.0, 2.0, ...); bracket the
    # integer so float compare never misses, and keep a cell if ANY node is on
    # the set (boundary node fields don't fill whole cells).
    thr = Threshold(Input=reader)
    thr.Scalars = ["POINTS", args.field]
    lo, hi = float(args.id) - 0.5, float(args.id) + 0.5
    # ParaView API for the range varies across versions; set whichever exists.
    for attr in ("LowerThreshold", "UpperThreshold"):
        try:
            setattr(thr, attr, lo if attr == "LowerThreshold" else hi)
        except Exception:
            pass
    try:
        thr.ThresholdRange = [lo, hi]   # older API
    except Exception:
        pass
    # keep a cell if ANY point passes (not all)
    for attr, val in (("AllScalars", 0), ("ComponentMode", "Selected")):
        try:
            setattr(thr, attr, val)
        except Exception:
            pass

    t = tvals[-1] if (args.last and tvals) else None
    if t is not None:
        thr.UpdatePipeline(t)
    else:
        thr.UpdatePipeline()

    # Report what survived so the user knows the id was right.
    n = thr.GetDataInformation().GetNumberOfPoints()
    print("[extract] side_set_id==%d -> %d points survived" % (args.id, n))
    if n == 0:
        print("[extract] nothing matched id=%d. Check the run's legend for the right id "
              "(and that this .pvd is from the tagged run)." % args.id)

    writer = XMLUnstructuredGridWriter(Input=thr, FileName=args.out)
    try:
        writer.DataMode = "Appended"
    except Exception:
        pass
    writer.UpdatePipeline()
    print("[extract] wrote %s (open it directly in ParaView)" % args.out)


if __name__ == "__main__":
    main()
