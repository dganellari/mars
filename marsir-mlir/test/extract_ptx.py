#!/usr/bin/env python3
"""Extract PTX from a serialized gpu.binary (assembly = "..." attr) on stdin."""
import re
import sys

txt = sys.stdin.read()
blobs = re.findall(r'assembly = "(.*?)">', txt, re.S)
if not blobs:
    sys.exit("no PTX assembly found in gpu.binary")
# MLIR escapes non-printables as \XX (two hex digits), plus \\ and \".
ptx = re.sub(r"\\([0-9A-Fa-f]{2})", lambda m: chr(int(m.group(1), 16)), blobs[0])
ptx = ptx.replace('\\"', '"').replace("\\\\", "\\")
out = sys.argv[1] if len(sys.argv) > 1 else "generated/sumfac_sm90.ptx"
open(out, "w").write(ptx)
n = ptx.count("mma.sync.aligned.m8n8k4")
print("wrote %s (%d chars, %d x mma.sync.aligned.m8n8k4 f64)" % (out, len(ptx), n))
min_mma = int(sys.argv[2]) if len(sys.argv) > 2 else 1
assert n >= min_mma, "expected >= %d DMMA instructions in PTX" % min_mma
