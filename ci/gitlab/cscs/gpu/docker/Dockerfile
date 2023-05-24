ARG BASE_IMAGE
FROM $BASE_IMAGE

ARG NUM_PROCS

COPY . /mars.src

ENV MARS_DIR=/mars.install

# scripting broken, hence we disable it.
# We are grepping in the output of cmake for the string 'Manually-specified variables were' to catch early error were variables are being renamed.
# It is an error to specify variables on the commandline that are unused.
RUN mkdir -p /mars.src/build \
  && cd /mars.src/build \
  && cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_COMPILER=/opt/view/bin/nvcc_wrapper \
    -DCMAKE_INSTALL_PREFIX=$MARS_DIR \
    -DCMAKE_VERBOSE_MAKEFILE=ON \
    -DMARS_ENABLE_KOKKOS=ON \
    -DMARS_ENABLE_CUDA=ON \
    -DMARS_ENABLE_BENCHMARK=ON \
    -DMARS_ENABLE_TESTING=ON \
    .. |& tee cmake_output.log \
  && ! grep 'Manually-specified variables were' cmake_output.log \
  && make -j$NUM_PROCS \
  && make install

RUN mkdir -p /utopia.src/utopia_fe/build \
  && cd /utopia.src/utopia_fe/build \
  && cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$UTOPIA_DIR \
    -DCMAKE_VERBOSE_MAKEFILE=OFF \
    -DCMAKE_CXX_COMPILER=/opt/view/bin/nvcc_wrapper \
    -DCMAKE_PREFIX_PATH=$MARS_DIR/lib/cmake \
    -DBUILD_SHARED_LIBS=OFF \
    -DUTOPIA_ENABLE_CXX14_FEATURES=ON \
    -DUTOPIA_ENABLE_LIBMESH=OFF \
    -DUTOPIA_ENABLE_MARS=ON \
    -DUTOPIA_ENABLE_MOONOLITH=OFF \
    -DUTOPIA_ENABLE_TRILINOS_ALGEBRA=ON \
    .. |& tee cmake_output.log \
  && ! grep 'Manually-specified variables were' cmake_output.log \
  && make -j$NUM_PROCS complete \
  && make install
