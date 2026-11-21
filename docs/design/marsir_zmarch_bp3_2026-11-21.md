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

## Next: lane-owned buffers (design, 21 November 2026)

The main cost above is that the z part reads and writes the whole result cube `w` in shared
memory on every plane. Within the plane loop every access to `w` is at a point the accessing
lane owns: the z part writes `w[a,b,k]` for the lane's column (a, b), the plane update reads
and writes `w[a,b,c]`. Such a buffer can live in the lane's registers.

A buffer is lane-owned over a region when every phase in the region that touches it has the
same iteration space (so point q runs on lane q mod L in all of them) and indexes it with that
space's induction variables in the same positions, plus constants for dimensions kept in the
thread. Each lane then holds only its own points: J = ceil(points / L) slots (2 for the BP3
march at p = 3 and 7) times the constant-indexed extent, e.g. 2 x q doubles.

`--mir-warp-per-element` would:
1. find lane-owned scratch buffers (the same "same point" rule FusePhases uses, with constant
   indices allowed in both accesses);
2. map their phases as `for j in 0..J: q = glane + j L; if (q < points) { ... }` and replace
   each access by a private array indexed [j, constants...];
3. fully unroll those lane loops and the loops of phases around them, so every private index
   is a constant and LLVM keeps the array in registers;
4. lower a fill of the buffer to a per-lane initialisation of its private entries.

The projection after the march reads `w` across lanes, so the march ends with one copy of
the lane-owned columns into a shared `w` cube, in the owner's mapping (phase over (a, b),
k in the thread); the schedule's copy and fill lowering must keep that blocking. Shared memory
in the march is then `uq` plus the planes; `w` takes a shared slot only for the projection,
when `uq` is dead.

## Lane-owned result (Beverin 757797): correct, still slower

Implemented as `--mir-warp-per-element=lane-owned` (gate `test/lane_owned_gate.py`), with the
projection back reordered to run along z first so `w` is touched only lane by lane. Shared
memory per element at p = 3: 2,000 B (half of today). The first GPU run (757788) failed: the
private arrays were allocas inside the element loop, i.e. dynamic stack allocations on every
element, which overrun AMDGPU's private segment (the CPU emulator, 3 elements, cannot see it;
`test/amd_target.py` now rejects any alloca outside the entry block). With them before the loop
(all promoted to registers, no scratch), on MI300A, two rounds:

| BP3 | current default | z-march, w shared | z-march, w lane-owned |
|---|---|---|---|
| p = 3 | 1.023 / 1.027 ms (4 blocks per CU) | 1.74 ms | 1.937 / 1.939 ms (8 blocks per CU) |
| p = 7 | 1.554 / 1.559 ms | 2.53 ms | 2.724 / 2.736 ms |

Residency doubled at p = 3 as designed, and the kernel is still 1.9x slower: the march's many
small plane phases and barriers, the unblocked 2-D derivative phases, and the register cost of
the lane-owned columns (132 / 274 VGPRs, guarded unrolled slots) outweigh it. Lane-owned
scratch stays an option; the BP3 schedule stays as it is. The same check found the trilinear BP
kernels (bp1_tri, bp3_tri) with per-lane scalars inside the element loop on AMD, now moved out
(their first MI300A run: Beverin 757799).
