# MARSIR on AMD (MI300A, gfx942)

Status, 20 November 2026. The generated BP3 operator runs correctly on MI300A, built from
the same MARSIR passes as the NVIDIA kernels. It is not yet tuned for AMD.

## Where things stand

| Kernel | GH200 (NVIDIA) | MI300A (AMD) |
|---|---|---|
| BP3 p=3, 262,144 elements | 0.661 ms, 25.4 GDoF/s | 1.184 ms, 14.2 GDoF/s |
| BP3 p=7, 32,768 elements | 0.578 ms, 29.0 GDoF/s | 1.556 ms, 10.8 GDoF/s |

Same gate (`mars_marsir_galerkin_gate`), same source, median of 20 warm applies; both pass the
check against host quadrature. GDoF/s counts element-local values. Jobs: Daint 5410396,
Beverin 757106 (with opt, rough edge 5; before it 1.638 / 1.986 ms, Beverin 757090).

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

## Five rough edges, and what we did about them

Upstream MLIR (LLVM 19) can target AMD, but in five places the AMD path is less smooth than
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

Fix (done in MLIR): `--mir-amdgpu-grid-dim` replaces `gpu.grid_dim` on AMD targets, before
the ROCDL conversion, by a read of the *dispatch packet*: a small record the GPU receives with
every kernel launch that holds the grid size and the block size; the block count is their
ceiling quotient. It first worked as a patch on the LLVM IR text; the pass gives the same
machine code, byte for byte.

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
and convert the pointer to a generic one for its uses. Since `opt -O3` runs before `llc`
(rough edge 5), LLVM promotes them to registers through that cast: the p=7 tensor-core kernel
went from 2,568 B of scratch per thread to none (116 VGPRs, 264 MFMA instructions).

Why it stays a patch: these arrays are not ours. Upstream `--convert-vector-to-scf` creates
them (temporaries of vector transfers) after `--mir-gpu-wrap`, in the default space; a
MARSIR pass that retypes them would have to fix every view of them. The clean fix belongs
upstream: allocas created during lowering should take the target's alloca address space
(the DLTI `alloca_memory_space` entry exists for this).

### 4. Our newer LLVM and Beverin's ROCm disagree on a directive

The kernel's assembly starts with `.amdhsa_code_object_version 5`, which names the format of
the binary. LLVM 19 writes it; the older assembler in Beverin's ROCm does not recognize it, so
assembling on Beverin fails. That was the first Beverin job.

Workaround: assemble on the laptop with LLVM 19, and only link on Beverin. Both code object
versions 5 and 4 work (`AMD_COV=4` builds version 4).

Proper fix: use one LLVM version end to end, or link where the versions match.

### 5. llc alone reads global memory with flat loads

LLVM has two halves: `opt`, the optimizer (LLVM IR to LLVM IR), and `llc`, the code
generator (LLVM IR to assembly). clang runs both; our first AMD path ran only `llc`. The
kernel's arrays reach LLVM as generic pointers inside MLIR's memref descriptor structs, and
`llc` alone cannot prove they point to global memory, so every access became `flat_load`: the
hardware resolves the address space per access, and each one also counts as an LDS access, so
waits on shared memory also wait for HBM. NVIDIA's backend treats kernel pointer arguments as
global by itself, so the PTX was never affected.

Fix: run `opt -O3` for gfx942 before `llc`, as clang does for HIP. AMD's part of that pipeline
infers the address spaces: all flat accesses become `global_load`/`global_store`, and the
cleanup saves about 15 registers. Measured on MI300A (Beverin 757106, two rounds, same
kernels): BP3 p=3 1.63 -> 1.18 ms, p=7 2.07 -> 1.56 ms.

Proper fix in MLIR as well: `--mir-gpu-wrap` should give the kernel's array arguments
`#gpu.address_space<global>`, as it already does `workgroup` for shared memory, so the IR
states the address space instead of LLVM inferring it.

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
   bytes reaches 3.5-3.6 TB/s there (Beverin 757094; 8 B and 16 B loads within 3%). With
   rough edge 5 fixed the operator reaches 1.6 / 0.9 TB/s, still 2.3x / 3.9x above that
   floor. Hardware counters (Beverin 757110, /opt/rocm-6.3.0/bin/rocprof; the uenv rocprof
   looks for rocm_agent_enumerator only in its own bin/ and silently runs nothing):
   - global traffic is ideal: FETCH_SIZE = G + u bytes, WRITE_SIZE = y bytes, no 32 B requests;
   - LDS bank-conflict cycles are 36% (p=3) / 44% (p=7) of LDS active cycles: the scratch
     padding follows NVIDIA's bank model (WarpPerElement kBanks), not CDNA3's;
   - the measured shapes keep one wavefront per SIMD (64 KB / 46 KB of LDS per block).
   The occupancy ablation (757093) predates the flat_load fix, so it is repeated with opt.
2. Tune for AMD with the same knobs as on GH200, measured as an AMD ablation.
3. Target layer: see "Target layer plan" above.
4. Port MARS's hand-written cell-wise kernels to HIP (32-bit warp masks), for full solves.
5. 64 lanes per element for the AMD tensor-core kernels.
