#!/usr/bin/env python3
"""The production native reader on duct_mesh.py files (needs netCDF, so it runs where MARS has it).

  test_duct_exodus.py CHECK_BINARY MPIEXEC NUMPROC_FLAG [launcher args...]

duct_exodus_check must accept each generated mesh on 1 and 2 ranks (arrays equal to the C++
lattice, side sets as inlet/outlet/walls) and fail when told another lattice.
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import duct_mesh as dm  # noqa: E402


def main(argv):
    if len(argv) < 4:
        print("usage: test_duct_exodus.py CHECK MPIEXEC NUMPROC_FLAG [launcher args...]")
        return 2
    check, mpiexec, flag, extra = argv[1], argv[2], argv[3], argv[4:]
    work = tempfile.mkdtemp(prefix="duct-exodus-")
    failures = 0
    try:
        for cells, args in ((4, []), (6, ["--length", "7.5", "--stretch", "2.5"])):
            exo = os.path.join(work, "duct-%d.exo" % cells)
            if dm.main(["--cells", str(cells), "--output", exo] + args) != 0:
                return 1
            for ranks in (1, 2):
                command = [mpiexec, flag, str(ranks)] + extra + [check, exo, "--cells", str(cells)] + args
                good = subprocess.call(command) == 0
                failures += not good
                print("%s: %s" % ("ok" if good else "FAIL", " ".join(command)))
        wrong = [mpiexec, flag, "1"] + extra + [check, os.path.join(work, "duct-4.exo"), "--cells", "4", "--length", "6"]
        rejected = subprocess.call(wrong, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) != 0
        failures += not rejected
        print("%s: another lattice is rejected" % ("ok" if rejected else "FAIL"))
    finally:
        shutil.rmtree(work)
    print("PASS" if not failures else "FAIL: %d" % failures)
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
