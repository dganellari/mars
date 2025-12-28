# CVFEM/EBVC Assembly Bench

Test bench for CVFEM/EBVC assembly operators using STK NGP utilities.

## Build

Before the first build run `git submodule update --init --recursive`. Build
utilities for several platforms are located in `scripts/`.  Some hard-coded
paths in these scripts may need to be adjusted.  Building is based on `cmake`
and executing the provided helper script should build the source (some
environment may need to be loaded first).  Alternatively, standard `cmake`
workflow may be used to build the source.

## Running

Running the test executable requires an input mesh argument.  A small test mesh
is provided in `test/mesh/pitz_daly.tar.gz` and must be extracted first using
`tar -C test/mesh -xzf test/mesh/pitz_daly.tar.gz`.

Building the code will create the `AssemblyTest.exe` executable in the build
directory.  If one of the helper scripts is used, a copy of the built executable
is placed in the root of the project according to the value of `target_name` in
the script.  For example, an OpenMP build using the build script on Alps can be
executed by:

```sh
./01_test.rel.O3.omp test/mesh/pitz_daly
```

## Assembler versions
### CvFem
Basic control-volume-finite-element assembler corresponding to a naive implementation very close to the 
original Nalu method. Uses master-element implementation from Nalu-Wind. Constructor
prepares shape function objects for all element-types. All scratch spaces are 
allocated in Kokkos' scratch level 1, i.e. global memory. 
### CvFemTeamScratch
Same as `CvFem`, except now all scratch spaces are allocated per team instead of 
per thread. This reduces non-coalesced accesses to the global scratch spaces. 
In addition, scratch spaces for index permutations and connected nodes are replaced by 
stack arrays. Minor performance improvement over original version. 
### CvFemInlineSF
Same as `CvFemTeamScratch` with the following adaptations: 
- Assembler is templated on an element type and precomputes some shape function data. 
- Shape function gradient data is calculated inline per integration point per node, 
thus reducing the number of required registers. 
- Uses pre-computed integration point area vectors instead of recalculating them every
time, further reducing registers but increasing overall memory footprint and bandwidth 
pressure.
- Node coordinates are stored in a team-scratch space.
- Currently only hex elements are supported and a way to loop over a mesh by element 
type is still needed. 
- This version is the fastest CvFem approach so far (speedup ~8x vs. the original `CvFem`).
### Edge
Edge-based assembly. Similar overall assembly approach, but ignoring shape functions
entirely.  Naturally significantly faster than the original `CvFem`, but lacking the 
added benefit of shape function approach, e.g. will require non-orthogonality 
correction for actual use, which is currently off. 
Currently does not support block sizes > 1.