# Copyright 2024 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import os
import reframe as rfm
import reframe.utility.sanity as sn
import reframe.utility.udeps as udeps
from collections import defaultdict
import uenv

mars_references = {
    'gh200': {
        'serial': {
            'small_mesh': {
                'time_run':  (9.25, -0.1, 0.1, 's'),
            },
            'medium_mesh': {
                'time_run':  (35.0, -0.05, 0.05, 's'),
            }
        },
        'distributed_mesh': {
            'medium': {
                'time_run':  (9.2, -0.05, 0.05, 's'),
            }
        },
    }
}

#  class cscs_reframe_tests_download(rfm.RunOnlyRegressionTest):
    #  descr = 'Fetch MARS sources code'
    #  sourcesdir = None
    #  executable = 'git'
    #  executable_opts = [
    #      'clone', 'git@github.com:dganellari/mars.git'
    #      '&&', 'cd', 'mars'
    #      '&&', 'git', 'checkout', 'master'
    #  ]
    #  local = True
    #
    #  @sanity_function
    #  def validate_download(self):
#          return sn.assert_eq(self.job.exitcode, 0)


class mars_download(rfm.RunOnlyRegressionTest):
    descr = 'Fetch MARS sources code'
    sourcesdir = None
    executable = 'git'
    executable_opts = [
        'clone', 'git@github.com:dganellari/mars.git'
        '&&', 'cd', 'mars'
        '&&', 'git', 'checkout', 'master'
    ]
    local = True

    @sanity_function
    def validate_download(self):
        return sn.assert_eq(self.job.exitcode, 0)

class mars_build(rfm.CompileOnlyRegressionTest):
    descr = 'Build MARS'
    valid_systems = ['+uenv']
    #  valid_prog_environs = ['+mpi +kokkos']
    valid_prog_environs = ['+mpi']
    build_system = 'CMake'
    sourcedir = None
    maintainers = ['dganellari']
    mars_sources = fixture(mars_download, scope='session')
    # NOTE: required so that the build stage is performed on
    # a compute node using an sbatch job.
    # This will force the uenv and view to be loaded using
    # "#SBATCH --uenv=" etc
    build_locally = False

    @run_before('compile')
    def prepare_build(self):
        self.uarch = uenv.uarch(self.current_partition)
        self.build_system.builddir = os.path.join(self.stagedir, 'build')

        self.prebuild_cmds = [
            f'rsync -a {self.mars_sources.stagedir}/mars/ .'
        ]

        self.build_system.config_opts = [
            ' -DCMAKE_BUILD_TYPE=Release',
            ' -DCMAKE_CXX_COMPILER=nvcc_wrapper',
            ' -DMARS_ENABLE_KOKKOS=ON',
            ' -DMARS_ENABLE_KERNELS=ON',
            ' -DMARS_ENABLE_CUDA=ON',
            ' -DMARS_ENABLE_TESTS=ON',
        ]
        # set architecture-specific flags
        if self.uarch == 'gh200':
            self.build_system.config_opts += [
                ' -DCMAKE_CUDA_ARCHITECTURES=90',
                ' -DCMAKE_CUDA_FLAGS_RELEASE="-expt-extended-lambda"',
            ]
        elif self.uarch == 'a100':
            self.build_system.config_opts += [
                ' -DCMAKE_CUDA_ARCHITECTURES=80',
                ' -DCMAKE_CUDA_FLAGS_RELEASE="-expt-extended-lambda"',
            ]
        elif self.uarch == 'zen2':
            self.build_system.config_opts += []

        self.build_system.max_concurrency = 64

        self.build_system.make_opts = []


@rfm.simple_test
class mars_unit(rfm.RunOnlyRegressionTest):
    descr = 'Run the mars unit tests'
    valid_systems = ['+uenv ']
    valid_prog_environs = ['+mpi']
    target_executable = variable(str, value='ctest')
    time_limit = '5m'
    maintainers = ['dganellari']

    mars_build = fixture(mars_build, scope='environment')

    @run_before('run')
    def prepare_run(self):
        self.executable = self.target_executable
        self.executable_opts = [
            '--test-dir', f'{self.mars_build.stagedir}/build', '-V', '-R', self.test
        ]

    @sanity_function
    def validate_test(self):
        return sn.assert_found(r'PASSED', self.stdout)

@rfm.simple_test
class mars_discretization(rfm.RunOnlyRegressionTest):
    descr = 'Run the mars FEM discretization tests small and medium mesh'
    valid_systems = ['uenv']
    valid_prog_environs = ['+mpi']
    target_executable = variable(str, value='ctest')
    maintainers = ['dganellari']
    model_size = parameter(['small_mesh', 'medium_mesh'])

    mars_build = fixture(mars_build, scope='environment')


    @run_before('run')
    def prepare_run(self):
        self.executable = self.target_executable
        self.executable_opts = [
            '--test-dir', f'{self.mars_build.stagedir}/build', '-V', '-R', self.test
        ]

        # Instead of explicitly listing performance targets for all possible
        # system:partition combinations, set the reference targets to those
        # for the uarch of the current partition.
        # * self.uarch is one of the alps arch: gh200, zen2, a100, ... or None
        # * self.current_partition.fullname is the vcluster:partition string,
        #   for example "daint:normal" or "todi:debug".
        self.uarch = uenv.uarch(self.current_partition)
        if (self.uarch is not None) and (self.uarch in mars_references):
            self.reference = {
                self.current_partition.fullname:
                    arbor_references[self.uarch]['serial'][self.model_size]
            }

    @sanity_function
    def assert_sanity(self):
        return sn.assert_found(r'All tests passed', self.stdout)

    @performance_function('s')
    def time_run(self):
        return sn.extractsingle(r'model-run\s+(\S+)', self.stdout, 1, float)

slurm_config = {
        'gh200': {"ranks": 4, "cores": 64, "gpu": True},
        'zen2': {"ranks": 2, "cores": 64, "gpu": False},
}


@rfm.simple_test
class mars_discretization_mpi(mars_discretization):
    """
    adapt the mars discretization test to run with MPI
    """

    descr = 'arbor busyring model MPI on a single node'
    model_size = parameter(['medium'])

    @run_before('run')
    def prepare_run(self):
        self.uarch = uenv.uarch(self.current_partition)
        self.executable = self.target_executable
        self.executable_opts = [
            '--test-dir', f'{self.mars_build.stagedir}/build', '-V', '-R', self.test
        ]

        self.num_tasks = slurm_config[self.uarch]["ranks"]
        self.num_cpus_per_task = slurm_config[self.uarch]["cores"]
        if slurm_config[self.uarch]["gpu"]:
            self.job.options = ['--gpus-per-task=1']

        self.uarch = uenv.uarch(self.current_partition)
        if (self.uarch is not None) and (self.uarch in mars_references):
            self.reference = {
                self.current_partition.fullname:
                    arbor_references[self.uarch]['distributed'][self.model_size]
            }
