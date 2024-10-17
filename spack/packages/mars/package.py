# Copyright 2013-2022 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.package import *


class Mars(CMakePackage, CudaPackage,  ROCmPackage):
    """Mesh Adaptive Refinement for Supercomputing."""

    homepage = "https://bitbucket.org/zulianp/mars"
    git = "https://dganellari@bitbucket.org/zulianp/mars.git"
    #  url = "https://dganellari@bitbucket.org/zulianp/mars.git"

    maintainers = ["danielganellari"]

    version("develop", branch="development")

    variant("kokkos", default=True)
    variant("openmp", default=True)
    variant("cxxopts", default=True)
    variant("tests", default=False)
    variant("build_type",
            default="Release",
            description="CMake build type",
            values=("Debug", "Release", "RelWithDebInfo"),
            )

#      with when('+rocm'):
        #  variant("magma", default=True, description="Use magma eigenvalue solver (AMDGPU)")
#          depends_on("magma+rocm", when="+magma+rocm")

    depends_on("cmake@3.18:", type="build")
    depends_on("mpi")
    #  depends_on("lapack")

    depends_on("kokkos")
    depends_on("kokkos+openmp", when="+openmp")

    depends_on("kokkos+cuda+cuda_lambda+wrapper", when="+cuda%gcc")
    depends_on("kokkos+cuda", when="+cuda")

    # rocm dependencies
    depends_on("kokkos+rocm", when="+rocm")

    depends_on("kokkos-kernels")
    #  depends_on("rocblas", when="+rocm")
    #  depends_on("rocsolver", when="+rocm")

    depends_on("googletest", type="build", when="+tests")
    depends_on("cxxopts")

    def cmake_args(self):
        options = [
            self.define_from_variant("MARS_ENABLE_TESTS", "tests"),
            self.define_from_variant("MARS_ENABLE_KOKKOS", "kokkos"),
            #  self.define_from_variant("MARS_ENABLE_OPENMP", "openmp"),
            self.define_from_variant("MARS_ENABLE_HIP", "rocm"),
            self.define_from_variant("MARS_ENABLE_CUDA", "cuda"),
            self.define_from_variant("MARS_ENABLE_CXXOPTS", "cxxopts"),
        ]

        if "+cuda%gcc" in self.spec:
            options.append(
                "-DCMAKE_CXX_COMPILER=%s" % self.spec["kokkos-nvcc-wrapper"].kokkos_cxx
            )
            options.append("-DMARS_ENABLE_CUDA=ON")

        if "+cuda" in self.spec:
            cuda_arch = self.spec.variants["cuda_arch"].value
            if cuda_arch[0] != "none":
                options += ["-DCMAKE_CUDA_FLAGS=-arch=sm_{0}".format(cuda_arch[0])]

        if "+rocm" in self.spec:
            options.append(self.define(
                "CMAKE_CXX_COMPILER", self.spec["hip"].hipcc))
            archs = ",".join(self.spec.variants['amdgpu_target'].value)
            options.append("-DHIP_HCC_FLAGS=--amdgpu-target={0}".format(archs))
            options.append("-DCMAKE_CXX_FLAGS=--amdgpu-target={0} --offload-arch={0}".format(archs))
            options.append("-DMARS_ENABLE_HIP=ON")

        return options
