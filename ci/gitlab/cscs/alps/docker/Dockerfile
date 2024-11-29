FROM jfrog.svc.cscs.ch/docker-ci-ext/base-containers/public/spack-build:spack0.21.0-ubuntu22.04-cuda12.4.1 as builder

ARG NUM_PROCS

ARG KOKKOS_SPEC

ARG ALPS_SPEC

ARG MARS_CXX_SPEC

ARG MARS_CUDA_SPEC

RUN spack-install-helper --target $ALPS_SPEC \
    "git" "cmake" "python@3.11" "vim +python +perl +lua" \
    "kokkos +aggressive_vectorization $KOKKOS_SPEC cxxstd=17 +pic +serial +shared +tuning" \
    "kokkos-kernels"

# copy only relevant parts to the final container
FROM jfrog.svc.cscs.ch/docker-ci-ext/base-containers/public/spack-runtime:ubuntu22.04-cuda12.4.1

# it is important to keep the paths, otherwise your installation is broken
# all these paths are created with the above `spack-install-helper` invocation
COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software
COPY --from=builder /opt/._view /opt/._view
COPY --from=builder /etc/profile.d/z10_spack_environment.sh /etc/profile.d/z10_spack_environment.sh

# Some boilerplate to get all paths correctly - fix_spack_install is part of the base image
# and makes sure that all important things are being correctly setup
RUN fix_spack_install

# Finally install software that is needed, e.g. compilers
# It is also possible to build compilers via spack and let all dependencies be handled by spack
RUN apt-get -yqq update && apt-get -yqq upgrade \
 && apt-get -yqq install build-essential gfortran \
 && rm -rf /var/lib/apt/lists/*


ENV MARS_DIR=/mars.install

COPY . /mars.src

RUN cd /mars.src && git submodule update --init --recursive
RUN mkdir -p /mars.src/build
RUN cd /mars.src/build \
  && cmake \
    -DCMAKE_BUILD_TYPE=Release \
    $MARS_CXX_SPEC \
    -DCMAKE_INSTALL_PREFIX=$MARS_DIR \
    -DCMAKE_VERBOSE_MAKEFILE=ON \
    -DMARS_ENABLE_KOKKOS=ON \
    -DMARS_ENABLE_KOKKOS_KERNELS=ON \
    $MARS_CUDA_SPEC \
    -DMARS_ENABLE_TESTS=ON \
    /mars.src |& tee cmake_output.log \
  && ! grep 'Manually-specified variables were' cmake_output.log \
  && make -j$NUM_PROCS \
  && make install
RUN rm -rf /mars.src/build

