# MARSIR z-march schedule for the Galerkin Laplacian (BP3), 21 November 2026

## Why

On MI300A the generated BP3 operator is limited by how many elements are resident per CU,
and shared memory per element sets that number (Beverin 757263: the same p = 3 kernel takes
1.02 / 1.58 / 2.81 ms at 4 / 2 / 1 blocks per CU). Today four q^3 buffers are live at once:
u at the quadrature points (`uq`) and its three derivatives, because each derivative phase
puts its points on different lanes and the pointwise flux must see all three at one point.
Inside the current phase model (values cross phases only through memory) that peak cannot
drop. The z-march changes the loop structure instead.

## The schedule

Interpolation (three B contractions) and the projection back (three Bt contractions) stay
as they are. The middle section, three D contractions of `uq` along axes 0, 1, 2, the flux,
three Dt contractions and the adds, becomes a loop over the planes c of axis 2 that every
lane of the element runs:

```
w = 0                                                    phase over (a, b, k)
for c in 0..q-1:                                         every lane, same trip count
    phase over (a, b):                                   lane-cyclic, as any phase
        d0 = sum_m D[a,m] uq[m,b,c]                      reads other lanes' points of plane c
        d1 = sum_m D[b,m] uq[a,m,c]
        d2 = sum_m D[c,m] uq[a,b,m]                      the lane's own column
        (dv0, dv1, dv2) = flux(G[:,a,b,c], d0, d1, d2)
        plane0[a,b] = dv0;  plane1[a,b] = dv1
        for k: w[a,b,k] += D[c,k] dv2                    the z part, the lane's own column
    barrier
    phase over (a, b):
        w[a,b,c] += sum_m D[m,a] plane0[m,b] + sum_m D[m,b] plane1[a,m]
    barrier                                              the planes are rewritten next
```

The summation order of the Dt contractions changes (the z part accumulates plane by plane),
so the result agrees with today's kernel to rounding, not bit for bit; the gate checks it
against host quadrature as usual.

## Shared memory per element (doubles)

| | today | z-march |
|---|---|---|
| interpolation (c2 + uq) | 225 / 1,377 | same |
| middle section | 4 q^3 = 500 / 2,916 | 2 q^3 + 2 q^2 = 300 / 1,620 |
| projection (w + c19) | 225 / 1,377 | same |
| peak, p = 3 / p = 7 | 4,000 B / 23,328 B | 2,400 B / 12,960 B |

Elements per 64 KB CU: 16 -> 27 at p = 3, 2 -> 5 at p = 7. `uq` and `w` overlap during the
march (d2 at plane c reads `uq` at every plane, and `w` at plane c is final only after the
last plane), which is why the middle section keeps two cubes.

## Compiler changes

1. `--mir-warp-per-element`: an `scf.for` with constant bounds in the element body may hold
   phases. Barriers inside it follow the same rule as between top-level phases; at the end
   of each iteration a barrier goes in when a phase of the next iteration accesses a buffer
   this iteration writes (the planes). Every lane reaches every barrier: the trip count is
   the same for all lanes, and the phases inside are lane-cyclic as before.
2. `--mir-cellwise-schedule=zmarch=true`: recognise the section on the bufferized linalg ops
   (three contractions of one buffer with one matrix along axes 0, 1, 2; a pointwise op
   with three outputs; three contractions of those with the transposed matrix and the
   adds into one result) and emit the loop above, with `plane0`, `plane1` as new
   `{mir.scratch}` buffers. Off by default, so every committed kernel stays byte-identical.
3. CPU warp emulator gate (`bp_gpu_pipeline.py --emulate`), then the operator gate and the
   full solve on MI300A and GH200.

## Open questions

- Lanes per element: the march phases have q^2 points (25 at p = 3, 81 at p = 7). With 16
  lanes at p = 3 some lanes take two columns, as the current phases already do.
- The z-part update writes q values of `w` per point each plane (q^2 extra shared writes
  per plane, lane-local). Keeping a lane's `w` column in registers would remove them, but
  a lane with two columns would index registers by a loop variable; measure first.

## Result (Beverin 757772, 21 November 2026): correct, slower

The prototype (middle section as an `scf.for {mir.planes}` of plane-sized `mir.contract` /
`mir.flux`; compiler support in commits ac40d416d and 249fd8e2e) passes the CPU emulator at
p = 3 and 7, with no race under ThreadSanitizer, and lowers shared memory per element from
4,000 to 3,000 B (p = 3) and from 23,328 to 14,904 B (p = 7). On MI300A, two rounds:

| BP3 | current AMD default | z-march prototype |
|---|---|---|
| p = 3, 262,144 elements | 1.025 / 1.029 ms, 4 blocks per CU | 1.739 / 1.737 ms, 5 blocks per CU |
| p = 7, 32,768 elements | 1.557 / 1.562 ms, 1 block per CU | 2.529 / 2.524 ms, 2 blocks per CU |

Residency rose as planned, and the extra work outweighed it:
- the z part is a contraction into the whole `w` cube on every plane: q^3 read-modify-writes
  per plane, q times the traffic of the one sweep it replaces;
- the plane phases are small (25 points per element at p = 3), so lanes idle, and each plane
  adds barriers;
- the 2-D derivative phases are not register-blocked (a plane has fewer rows than lanes), and
  at p = 7 D and Dt lose the even-odd folding (they are read through row and column views).

Keeping a lane's `w` column in registers across the planes would remove the first cost, but
that is lane-private state carried across phases, which the phase model (values cross phases
only through memory) cannot express. So within the current compiler the z-march does not pay
on MI300A; it is not a default. The compiler support stays (correct, gated, and the unit-
dimension and per-table-folding fixes found on the way matter on their own).
