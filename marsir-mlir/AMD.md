# MARSIR on AMD (MI300A, gfx942)

Status, 20 November 2026. The generated BP3 operator runs correctly on MI300A, built from
the same MARSIR passes as the NVIDIA kernels. It is not yet tuned for AMD.

## Where things stand

| Kernel | GH200 (NVIDIA) | MI300A (AMD) |
|---|---|---|
| BP3 p=3, 262,144 elements | 0.661 ms, 25.4 GDoF/s | 1.638 ms, 10.2 GDoF/s |
| BP3 p=7, 32,768 elements | 0.578 ms, 29.0 GDoF/s | 1.986 ms, 8.4 GDoF/s |

Same gate (`mars_marsir_galerkin_gate`), same source, median of 20 warm applies; both pass the
check against host quadrature. GDoF/s counts element-local values. Jobs: Daint 5410396,
Beverin 757090.

Both kernels are limited by memory traffic: per element they read the geometry tensor G
(6 values per quadrature point) and u, and write y, 7,024 B at p=3 and 43,184 B at p=7. That
gives 2.8 / 2.4 TB/s on GH200 and 1.1 / 0.7 TB/s on MI300A, so MI300A reaches a much smaller
share of its memory bandwidth.

Occupancy is not the cause. MI300A holds 2 blocks per CU against 6 to 8 per SM on GH200, but
an AMD ablation at p=3 (Beverin 757093, lanes per element 16 / 32 / 64, each run twice)
raised it to 4 and 16 blocks per CU and was slower: 1.65 / 2.38 / 1.83 ms.

How to build an AMD kernel:

```bash
python3 test/amd_target.py bp_gpu_pipeline bp3 3     # -> generated/bp3_p3_gfx942.s
llvm-mc -triple=amdgcn-amd-amdhsa -mcpu=gfx942 -filetype=obj generated/bp3_p3_gfx942.s -o bp3_p3.o
# on the AMD machine, with ROCm:
clang++ -target amdgcn-amd-amdhsa -mcpu=gfx942 bp3_p3.o -o bp3_p3_gfx942.hsaco
```

MARS loads the `.hsaco` through `marsir/mars_marsir_driver.hpp`, which maps the CUDA driver
calls the loaders use to HIP's module API.

## Four rough edges, and what we did about them

Upstream MLIR (LLVM 19) can target AMD, but in four places the AMD path is less smooth than
the NVIDIA one. Each workaround lives in `test/amd_target.py`.

### 1. MLIR's packaging step wants ROCm's device libraries

MLIR's last step, `gpu-module-to-binary`, turns the kernel into a GPU binary. For AMD it
first links ROCm's "device libraries" (`ockl.bc`, `ocml.bc`): small precompiled helper
functions for things like math and kernel queries. It requires them even when the kernel
uses none. A Mac has no ROCm, so the step stops.

Workaround: skip that step. Translate the kernel to plain LLVM IR and compile it with `llc`
for gfx942.

Proper fix: build where ROCm is installed and point MLIR at its libraries, or remove the
only reason we needed them (see 2).

### 2. Asking for the grid size calls a ROCm library function

Our kernels loop over elements in steps of the grid size (the total number of thread blocks).
MLIR lowers "how many blocks are there?" (`gpu.grid_dim`) to a call to a ROCm library
function, `__ockl_get_num_groups`. On NVIDIA the same question is a single hardware register
read. Without the library the call is left unresolved, so linking fails: the "undefined
symbol" error in the second Beverin job.

Workaround: we define that function ourselves inside the kernel. It reads two numbers from
the *dispatch packet*, a small record the GPU receives with every kernel launch that holds
the grid size and the block size, and divides them.

Proper fix: lower `gpu.grid_dim` straight to that read in MLIR, with no library call.

### 3. Scratch arrays land in the wrong address space ("alloca address space")

Some kernels keep small per-thread scratch arrays, for example in the p=7 tensor-core
kernel. In compiler terms such an array is an `alloca`: memory on the thread's own stack.

GPUs have several kinds of memory, each with its own numbered *address space*:

| Address space (AMD) | Memory | Who sees it |
|---|---|---|
| 1 | global (HBM) | every thread |
| 3 | shared / LDS | threads of one block |
| 4 | constant | every thread, read-only |
| 5 | private | one thread (its stack or registers) |
| 0 | "generic" | a pointer that can point into any of the above |

On AMD, a thread's own stack must be in address space 5; LLVM's AMD data layout says so
(`A5`). MLIR's memory lowering creates the scratch arrays in the generic space 0 instead.
NVIDIA's back end accepts that and sorts it out itself; AMD's instruction selector stops
with `Cannot select: FrameIndex` ("I don't know where this stack slot lives").

Workaround: after translating to LLVM IR we rewrite each such array: allocate it in space 5
and convert the pointer to a generic one for its uses. This works, but the conversion stops
LLVM from turning the small array into registers, so the p=7 tensor-core kernel keeps
2.5 KB of slow scratch memory per thread.

Proper fix: create these arrays in the private space from the start in MLIR
(`#gpu.address_space<private>`), so LLVM sees private memory and can promote it to registers,
as it does on NVIDIA.

### 4. Our newer LLVM and Beverin's ROCm disagree on a directive

The kernel's assembly starts with `.amdhsa_code_object_version 5`, which names the format of
the binary. LLVM 19 writes it; the older assembler in Beverin's ROCm does not recognize it, so
assembling on Beverin fails. That was the first Beverin job.

Workaround: assemble on the laptop with LLVM 19, and only link on Beverin. Both code object
versions 5 and 4 work (`AMD_COV=4` builds version 4).

Proper fix: use one LLVM version end to end, or link where the versions match.

## What to borrow from IREE

IREE is a complete ML compiler with strong AMD support. Its scope differs from MARSIR's, but
several of its design choices fit our problem. What follows is from general knowledge of
IREE's design, not checked against its current source.

- **Matrix-instruction layouts as data.** IREE describes each MFMA/MMA instruction (shape and
  which lane holds which entry) as an attribute, and generic passes derive the data movement
  from it. Our tensor-core pass writes the lane formulas by hand, separately for NVIDIA and
  AMD. With layouts as data, a new GPU would be a table entry instead of a hand derivation.
- **A target description.** IREE attaches to each target its wave size, shared-memory limit
  and available matrix instructions, and the passes read them. Ours hard-code 32 lanes
  (`WarpPerElement`, `WarpDistribute`), which is what blocks the 64-lane tensor-core mapping
  on AMD.
- **Pipelines in C++, chosen by target.** Our pipelines are Python scripts that call
  `mir-opt` several times, and `amd_target.py` patches them at run time. Registered C++
  pipelines selected by the `gpu.module` target attribute (`#nvvm.target`, `#rocdl.target`)
  would replace that.
- **Per-target tuning.** IREE keeps tuning choices per target. Our measured schedule choices
  (lanes per element, elements per block, register caps) were all made on GH200 and need
  their own AMD values.

## Tensor cores on AMD (in progress)

`--mir-chain-contracts=target=amd` emits `v_mfma_f64_4x4x4_4b_f64`: four 4x4 blocks as the
quadrants of our 8x8 tiles, two instructions per contraction, as on NVIDIA. The lane layout
comes from AMD's matrix instruction calculator and is checked by
`../paper1-work/amd/mfma/layout_model.py`. A result feeding the next contraction's B operand
needs only a swap with lane `L xor 8`; feeding the A operand needs a general shuffle.

The production p=7 CVFEM kernel compiles with 264 MFMA instructions, one per NVIDIA m8n8k4,
but cannot run yet: the rest of the schedule still assigns 32 lanes per element, and the AMD
instruction needs the element on all 64 lanes of a wavefront.

## Target layer plan (best practice)

Goal: one MARSIR pipeline per kernel family, written in C++, that reads the GPU's properties
from a target description; NVIDIA and AMD differ by one option. Every step must leave the
committed NVIDIA PTX byte-identical, which step 0 checks.

0. **Golden check.** `test/check_generated.py`: a list of every committed kernel file with the
   command that produces it; regenerate all into a temporary directory and compare bytes. Run
   before and after each step.
1. **Target description.** An attribute `#mir.target<wave = 32, shared = 49152, mma = ...>` on
   the module, set by `--mir-set-target=nvidia|amd`, and one C++ helper the passes call. When
   the attribute is absent the helper returns today's NVIDIA values, so nothing changes.
2. **Wave size from the target.** Replace the hard-coded 32 in WarpPerElement, WarpWrap,
   WarpDistribute, RegisterResidentChain, DistributeFills, LowerCellwise, CellwiseSchedule
   and EmulateWarp. This unblocks 64 lanes per element for the AMD tensor-core kernels.
3. **Text edits become passes.** The Python scripts edit the IR as text between passes (for
   example the `mir.constant` / `mir.reflect` argument attributes); each becomes a small pass.
4. **Registered C++ pipelines** for bp, hl, dss, block_dss and tet, with typed options. The
   "halve the warps until shared memory fits" retry becomes arithmetic in the schedule pass
   against the target's limit; the measured tuning tables become per-target defaults. The
   Python scripts shrink to one `mir-opt` call plus the gates.
5. **AMD lowering in MLIR.** Private-space allocas, `gpu.grid_dim` read from the dispatch
   packet, and the ROCDL tail as the `target=amd` branch of each pipeline; `amd_target.py`
   and its patches go away. Checked on Beverin: same correctness and timing as today.
6. **Matrix-instruction layouts as data.** A table per instruction (shape, which lane holds
   which A/B/C entry); the tensor-core pass derives its shuffles from it for both vendors.

The AMD performance work (the operator is 3.2x to 5x above MI300A's bandwidth floor) builds
on steps 1 and 2.

## Next steps

1. Find where the time goes inside the kernel on MI300A. A streaming kernel moving the same
   bytes reaches 3.5-3.6 TB/s there (Beverin 757094; 8 B and 16 B loads within 3%), so the
   operator (1.1 / 0.7 TB/s) loses 3.2x / 5x inside its own schedule, not in the platform.
   Hardware counters next (rocprofv3 from the image crashed on counter collection).
2. Tune for AMD with the same knobs as on GH200, measured as an AMD ablation.
3. Target layer: see "Target layer plan" above.
4. Port MARS's hand-written cell-wise kernels to HIP (32-bit warp masks), for full solves.
5. 64 lanes per element for the AMD tensor-core kernels.
