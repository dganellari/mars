#!/usr/bin/env python3
"""BP3 in the column form (docs/design/marsir_bp3_columns_2026-11-21.md): the emitted operator
with every x and y contraction done plane by plane, from a one-plane copy, so each 3-D buffer is
touched only at the (x, y) column of the thread that owns it. With --mir-bp-schedule columns and
lane-owned, those buffers stay in registers and shared memory holds only planes.
Usage: bp_columns.py <p> <out.mir> [--whole=interp,middle,proj]   (q = p + 2; --whole keeps those
parts as whole-cube contractions, for bisecting)"""
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
whole = next((a.split("=", 1)[1].split(",") for a in sys.argv if a.startswith("--whole=")), [])
p = int(sys.argv[1])
n, q = p + 1, p + 2
src = subprocess.run([sys.executable, "-c", f"""
import sys; sys.path.insert(0, {os.path.join(HERE, '..', '..', 'marsir-compiler')!r})
from marsir import parse_spec_file, synthesize
from marsir.backends import mlir_ir
ga = synthesize(parse_spec_file({os.path.join(HERE, '..', '..', 'marsir-compiler', 'specs', 'bp3.op')!r}))
sys.stdout.write(mlir_ir.emit_galerkin(ga, p={p}, q={q}, real='f64'))
"""], capture_output=True, text=True, check=True).stdout
flux = re.search(r"%x13:3 = mir\.flux ins\(.*?\) : \(.*?\) -> \(.*?\) \{\n(.*?)\n  \}", src, re.S).group(1)


def T(*s):
    return f"tensor<{'x'.join(map(str, s))}xf64>"


B, Bt, D = T(q, n), T(n, q), T(q, q)
copy = "{\n  ^bb0(%a: f64):\n    mir.yield %a : f64\n  }"
add = "{\n  ^bb0(%a0: f64, %a1: f64):\n    %t = arith.addf %a0, %a1 : f64\n    mir.yield %t : f64\n  }"


def planes(name, cube, cube_t, plane_t, count, body, out_t, result_plane_t):
    """out = for z < count: plane z of cube -> one-plane copy -> body -> plane z of out."""
    pz = f"[0, 0, %z] [{plane_t[0]}, {plane_t[1]}, 1] [1, 1, 1]"
    oz = f"[0, 0, %z] [{result_plane_t[0]}, {result_plane_t[1]}, 1] [1, 1, 1]"
    return f"""  %{name}_e = tensor.empty() : {T(*out_t)}
  %{name} = scf.for %z = %i0 to %i{count} step %i1 iter_args(%A = %{name}_e) -> ({T(*out_t)}) {{
    %pl = tensor.extract_slice %{cube}{pz} : {T(*cube_t)} to {T(*plane_t)}
    %s = mir.flux ins(%pl) : ({T(*plane_t)}) -> {T(*plane_t)} {copy}
{body}
    %A1 = tensor.insert_slice %r into %A{oz} : {T(*result_plane_t)} into {T(*out_t)}
    scf.yield %A1 : {T(*out_t)}
  }} {{mir.planes}}
"""


Q2 = ", ".join([T(q, q)] * 9)
gs = "\n".join(f"    %g{k} = tensor.extract_slice %G[{k}, 0, 0, %z] [1, {q}, {q}, 1] [1, 1, 1, 1] : "
               f"{T(6, q, q, q)} to {T(q, q)}" for k in range(6))
pz = f"[0, 0, %z] [{q}, {q}, 1] [1, 1, 1]"
out = f"""func.func @bp3_apply(%u: {T(n, n, n)}, %B: {B}, %Bt: {Bt}, %D: {D}, %Dt: {D}, %G: {T(6, q, q, q)}) -> {T(n, n, n)} {{
  %i0 = arith.constant 0 : index
  %i1 = arith.constant 1 : index
  %i{n} = arith.constant {n} : index
""" + (f"  %i{q} = arith.constant {q} : index\n") + f"""  %c1 = mir.contract %u, %B {{axis = 0 : i64}} : ({T(n, n, n)}, {B}) -> {T(q, n, n)}
""" + (f"  %c2 = mir.contract %c1, %B {{axis = 1 : i64}} : ({T(q, n, n)}, {B}) -> {T(q, q, n)}\n"
       if "interp" in whole else planes("c2", "c1", (q, n, n), (q, n), n,
             f"    %r = mir.contract %s, %B {{axis = 1 : i64}} : ({T(q, n)}, {B}) -> {T(q, q)}",
             (q, q, n), (q, q))) + f"""  %c3 = mir.contract %c2, %B {{axis = 2 : i64}} : ({T(q, q, n)}, {B}) -> {T(q, q, q)}
  %c6 = mir.contract %c3, %D {{axis = 2 : i64}} : ({T(q, q, q)}, {D}) -> {T(q, q, q)}
  %we = tensor.empty() : {T(q, q, q)}
  %ve = tensor.empty() : {T(q, q, q)}
  %m:2 = scf.for %z = %i0 to %i{q} step %i1 iter_args(%W = %we, %V = %ve) -> ({T(q, q, q)}, {T(q, q, q)}) {{
    %pl = tensor.extract_slice %c3{pz} : {T(q, q, q)} to {T(q, q)}
    %s = mir.flux ins(%pl) : ({T(q, q)}) -> {T(q, q)} {copy}
    %d0 = mir.contract %s, %D {{axis = 0 : i64}} : ({T(q, q)}, {D}) -> {T(q, q)}
    %d1 = mir.contract %s, %D {{axis = 1 : i64}} : ({T(q, q)}, {D}) -> {T(q, q)}
    %d2 = tensor.extract_slice %c6{pz} : {T(q, q, q)} to {T(q, q)}
{gs}
    %y:3 = mir.flux ins(%d0, %d1, %d2, %g0, %g1, %g2, %g3, %g4, %g5) : ({Q2}) -> ({T(q, q)}, {T(q, q)}, {T(q, q)}) {{
{flux}
    }}
    %t0 = mir.contract %y#0, %Dt {{axis = 0 : i64}} : ({T(q, q)}, {D}) -> {T(q, q)}
    %t1 = mir.contract %y#1, %Dt {{axis = 1 : i64}} : ({T(q, q)}, {D}) -> {T(q, q)}
    %ws = mir.flux ins(%t0, %t1) : ({T(q, q)}, {T(q, q)}) -> {T(q, q)} {add}
    %W1 = tensor.insert_slice %ws into %W{pz} : {T(q, q)} into {T(q, q, q)}
    %V1 = tensor.insert_slice %y#2 into %V{pz} : {T(q, q)} into {T(q, q, q)}
    scf.yield %W1, %V1 : {T(q, q, q)}, {T(q, q, q)}
  }} {{mir.planes}}
  %c17 = mir.contract %m#1, %Dt {{axis = 2 : i64}} : ({T(q, q, q)}, {D}) -> {T(q, q, q)}
  %w18 = mir.flux ins(%m#0, %c17) : ({T(q, q, q)}, {T(q, q, q)}) -> {T(q, q, q)} {add}
  %c19 = mir.contract %w18, %Bt {{axis = 2 : i64}} : ({T(q, q, q)}, {Bt}) -> {T(q, q, n)}
""" + (f"  %c20 = mir.contract %c19, %Bt {{axis = 1 : i64}} : ({T(q, q, n)}, {Bt}) -> {T(q, n, n)}\n"
       f"  %c21 = mir.contract %c20, %Bt {{axis = 0 : i64}} : ({T(q, n, n)}, {Bt}) -> {T(n, n, n)}\n"
       if "proj" in whole else planes("c20", "c19", (q, q, n), (q, q), n,
             f"    %r = mir.contract %s, %Bt {{axis = 1 : i64}} : ({T(q, q)}, {Bt}) -> {T(q, n)}",
             (q, n, n), (q, n)) + planes("c21", "c20", (q, n, n), (q, n), n,
             f"    %r = mir.contract %s, %Bt {{axis = 0 : i64}} : ({T(q, n)}, {Bt}) -> {T(n, n)}",
             (n, n, n), (n, n))) + "  return %c21 : " + T(n, n, n) + "\n}\n"
open(sys.argv[2], "w").write(out)
print(f"BP3 p = {p} in the column form -> {sys.argv[2]}")
