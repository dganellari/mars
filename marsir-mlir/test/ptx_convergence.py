#!/usr/bin/env python3
"""Gate: no warp-collective instruction (mma.sync, shfl.sync, bar) may sit where
only some lanes reach it. LLVM does not treat the mma intrinsic as convergent,
so it may sink an mma whose result feeds only a lane-guarded store into that
guard -- undefined behavior for mma.sync.aligned, which every lane must execute.

Lane-dependent predicates are traced from %tid.x through the integer data flow;
a collective between a branch on such a predicate and the branch target is
reported.

Usage: python3 test/ptx_convergence.py kernel.ptx"""
import re
import sys


def main(path):
    lines = open(path).read().splitlines()
    lane = set()        # registers whose value depends on the lane id
    for l in lines:
        m = re.match(r"\s*(?:@!?%p\d+\s+)?([a-z][\w.]*)\s+(%[\w]+),\s*(.*);", l)
        if not m:
            continue
        op, dst, srcs = m.groups()
        if "%tid" in srcs or "%laneid" in srcs or any(
                s in lane for s in re.findall(r"%\w+", srcs)):
            if not op.startswith(("ld.", "st.", "mma", "shfl")):
                lane.add(dst)
    labels = {l.strip()[:-1]: i for i, l in enumerate(lines)
              if re.match(r"\s*\$L__\w+:", l)}
    bad = []
    for i, l in enumerate(lines):
        m = re.match(r"\s*@!?(%p\d+)\s+bra(?:\.uni)?\s+(\$L__\w+);", l)
        if not m or m.group(1) not in lane:
            continue
        end = labels.get(m.group(2), i)
        if end <= i:
            continue    # a backward branch: a loop, not a guard
        for j in range(i + 1, end):
            if re.search(r"\b(mma\.sync|shfl\.sync|bar\.sync|barrier\.sync)\b", lines[j]):
                bad.append((j + 1, lines[j].strip()[:60], i + 1))
    for j, txt, g in bad:
        print(f"  line {j}: {txt}   (guarded by line {g})")
    print(f"{len(bad)} warp-collective instruction(s) under a lane-dependent branch")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
