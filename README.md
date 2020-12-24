# M.A.R.S #
## Mesh Adaptive Refinement for Supercomputing ##

MARS is an open-source mesh management library designed to handle N-dimensional elements (N <= 4). 
MARS is developed in C++ and makes use of template meta-programming to have compile time dimensions of elements and vectors, thus allowing for both compile time performance optimizations and concise and reusable code.

The main features of MARS consist of:

1. Parallel mesh generation

2. Adaptive mesh refinement using bisection algorithms

3. Conforming mesh data-structure

4. Mesh quality estimators to study the output of different mesh-refinement strategies

5. Performance portable algorithms and data-structures targetting different accelerators.

MARS targets multi-core CPUs and GPUs using the C++ Kokkos programming model. The mesh is entirely constructed and stored on the device (GPUs). This enables libraries using MARS to perform further operations directly on the device, avoiding going through the host. Currently, MARS supports as its performance portable, parallel, adaptive refinement based algorithm the LEPP (Longest edge propagation path) from Rivara. Mesh generation is fully supported in parallel. As an example, 143 million Hex8 elements can be generated on a single node in just 0.86 sec.

A distributed memory, parallel implementation based on MPI is ongoing work, and forest of octrees and space filling curves algorithms for efficient mesh partitioning are being planned.

## Downloading MARS and its dependencies ##

Clone the repository and its submodules. MARS relies on googletest and google/benchmark.

`git clone --recurse-submodules https://bitbucket.org/zulianp/mars.git`
or for older git versions

`git clone https://bitbucket.org/zulianp/mars.git && cd mars && git submodule update --init --recursive`

Compiling M.A.R.S for serial usage:

	- cd mars/
	- mkdir build
	- cd build
	- cmake ..
	- make


MARS will automatically find Kokkos if installed into your system. It looks for KOKKOS_DIR or TRILINOS_DIR into the environment variables. It can work with kokkos standalone or with kokkos from the Trilinos library. Use -DTRY_WITH_KOKKOS=ON to use the feature. For more details check CMakeLists.txt.

The default when compiling MARS with Kokkos without specifing any other CMAKE flag is the Kokkos/OpenMP execution space. Kokkos should also be compiled with OpenMP support. Otherwise the default is the serial execution space.

To compile for CUDA the Cmake flag needs to be set: MARS_USE_CUDA=ON. An example would be: 
cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_BUILD_TYPE=Release -DTRY_WITH_KOKKOS=ON -DMARS_USE_CUDA=ON ..
In this case Kokkos should also be compiled with CUDA support.
To disable RDMA set MARS_NO_RDMA=ON.


## Cite MARS ##

If you use MARS CPU please use the following bibliographic entry


```
#!bibtex

@misc{marscpu,
    author = {Zulian, Patrick and Ganellari, Daniel and Rovi, Gabriele},
    title = {{MARS} - {M}esh {A}daptive {R}efinement for {S}upercomputing. {G}it repository},
    url = {https://bitbucket.org/zulianp/mars},
    year = {2018}
}
```

If you use MARS GPU please use the following bibliographic entry


```
#!bibtex

@misc{marsgpu,
    author = {Ganellari, Daniel and Zulian, Patrick and Rovi, Gabriele},
    title = {{MARS} - {M}esh {A}daptive {R}efinement for {S}upercomputing. {G}it repository},
    url = {https://bitbucket.org/zulianp/mars},
    year = {2018}
}
```




