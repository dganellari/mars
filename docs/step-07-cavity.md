# Step 7: the lid-driven cavity

Step 7 leaves the model problems and solves the incompressible Navier–Stokes equations with the
MARS flow solver. The case is the lid-driven cavity: fluid at rest in a closed box, set in motion by
its top wall, which slides with velocity U in x. It is the smallest flow example and the classic
first test of an incompressible solver, because it has no inlet and no outlet.

The program is `examples/distributed/unstructured/mars_lid_driven_cavity.cu`. It uses the same
solver (`fem/mars_navier_stokes.hpp`) as the Poiseuille and Taylor–Green steps that follow, and it
needs MARS built with Hypre (`-DMARS_ENABLE_HYPRE=ON`): the pressure and velocity systems are solved
with Hypre's algebraic multigrid.

## Run it

`--cells=32` makes the program generate a cube of 32³ hexahedra on every rank, so no mesh file is
needed:

```bash
mpirun -np 4 ./build/examples/distributed/unstructured/mars_lid_driven_cavity \
    --cells=32 --nu=0.01 --dt=0.005 --num-steps=400 --report-every=50
```

With a unit lid speed and a unit box, the Reynolds number is Re = U L / ν = 100. Every 50 steps the
program prints a line like

```
Step     50  t=0.2500  KE=...  continuity=...  amg(u,v,w,p)=.../.../.../...
```

| step | t | kinetic energy | continuity |
|---:|---:|---:|---:|
| 50 | 0.25 | 1.526e-2 | 9e-13 |
| 100 | 0.5 | 1.879e-2 | 2e-13 |
| 200 | 1.0 | 2.258e-2 | 7e-14 |
| 400 | 2.0 | 2.590e-2 | 4e-14 |

The kinetic energy still grows at t = 2, more slowly each time: the flow is approaching its steady
state. The continuity stays at round-off level. The multigrid solves take 4 to 5 iterations for each
velocity component and 20 to 22 for the pressure. On NVIDIA and AMD GPUs the kinetic energy agrees to
8 digits. The release check `marsReleaseCavity` runs 100 steps on a 16³ mesh and requires the same
kinetic energy on 1 and N ranks.

- **KE** is the kinetic energy of the flow, ½ Σ V_i |u_i|² over the nodes' control volumes V_i. It
  starts at zero and grows while the lid drags the fluid along, then levels off as the flow
  approaches its steady state.
- **continuity** is the largest net flow out of a control volume, divided by its volume: how far
  the discrete velocity is from divergence-free. An incompressible solver must keep it small at
  every step.
- **amg(u,v,w,p)** are the multigrid iterations of the three velocity solves and the pressure solve.

## What the program sets up

### Boundary conditions

```cpp
auto cavity = [box, eps, U] __device__(RealType x, RealType y, RealType z) {
    NodeCondition<RealType> node;
    bool side   = ... x or y on the box faces ...;
    bool bottom = fabs(z - box.lo[2]) < eps;
    bool lid    = fabs(z - box.hi[2]) < eps && !side;
    node.velocityFixed = side || bottom || lid;
    node.velocity[0]   = lid ? U : 0;
    return node;
};
```

Every boundary node is a wall: the velocity is prescribed there, zero except on the lid. The
function runs on the GPU for every node, which is how all MARS flow examples state their boundary
conditions. Where the lid meets a side wall, the side wall wins, so the corners stay at rest.

With no inlet or outlet, the pressure is fixed only up to a constant, which the solver removes.

### The time loop

```cpp
for (int step = 1; step <= opt.numSteps; ++step)
{
    solver.step();
    ...   // every report-every steps: kineticEnergy(), maxContinuity()
}
```

Each `step()` advances the velocity and pressure by one time step with the control-volume method of
[step 4](step-04-cvfem.md). The [Poiseuille step](poiseuille_tutorial.md) explains the method in
detail.

## Try this

- `--vtu-output=cavity --vtu-every=50` writes velocity and pressure frames; open `cavity.pvd` in
  ParaView to watch the main vortex form.
- `--nu=0.001` gives Re = 1000, a stronger and less viscous flow; use a finer mesh (`--cells=64`) and
  a smaller time step, and compare the frames with the Re = 100 run.

Next: [step 8](poiseuille_tutorial.md) validates the solver against an exact solution.
