FROM jfrog.svc.cscs.ch/docker-ci-ext/base-containers/public/spack-build:spack0.21.0-ubuntu22.04-cuda12.4.1 as builder

ARG NUM_PROCS

ARG KOKKOS_SPEC

ARG SPACK_ENV

RUN spack-install-helper \
    alps-zen2 \
    "kokkos +aggressive_vectorization ~alloc_async +pic +serial +shared +tuning $KOKKOS_SPEC" \
    # "kokkos-kernels" \
    "cmake" \
    "git"

ENV MARS_DIR=/mars.install

RUN git clone https://github.com/dganellari/mars.git --branch development /mars.src
RUN cd /mars.src && git submodule update --init --recursive
RUN mkdir -p /mars.src/build
RUN spack env activate /opt/spack-environment \
  && cd /mars.src/build \
  && cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=/opt/view/bin/nvcc_wrapper \
    -DCMAKE_INSTALL_PREFIX=$MARS_DIR \
    -DCMAKE_VERBOSE_MAKEFILE=ON \
    -DMARS_ENABLE_KOKKOS=ON \
    -DMARS_ENABLE_KOKKOS_KERNELS=OFF \
    -DMARS_ENABLE_CUDA=ON \
    -DMARS_ENABLE_TESTS=ON \
    -DMARS_ENABLE_EXAMPLES=ON \
    /mars.src |& tee cmake_output.log \
  && ! grep 'Manually-specified variables were' cmake_output.log \
  && make -j$NUM_PROCS \
  && make install \
  && spack env deactivate
RUN rm -rf /mars.src/build

# end of builder container, now we are ready to copy necessary files

# copy only relevant parts to the final container
FROM jfrog.svc.cscs.ch/docker-ci-ext/base-containers/public/spack-base:spack0.21.0-ubuntu22.04-cuda12.4.1

# copy mars installation to the final container
COPY --from=builder $MARS_DIR $MARS_DIR

# it is important to keep the paths, otherwise your installation is broken
# all these paths are created with the above `spack-install-helper` invocation
COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software
COPY --from=builder /opt/._view /opt/._view
COPY --from=builder /etc/profile.d/z10_spack_environment.sh /etc/profile.d/z10_spack_environment.sh

RUN fix_spack_install

RUN apt-get -yqq update && apt-get -yqq upgrade \
 && apt-get -yqq install build-essential gdb gfortran \
 && rm -rf /var/lib/apt/lists/*
