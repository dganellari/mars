# Copyright 2024 Swiss National Supercomputing Centre (CSCS/ETH Zurich)
# ReFrame Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: BSD-3-Clause

import os
import reframe as rfm
import reframe.utility.sanity as sn
import reframe.utility.udeps as udeps
import mars.reframe.utils as mars_utils
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

class mars_download(rfm.RunOnlyRegressionTest):
    descr = 'Fetch MARS sources code'
    sourcesdir = None
    executable = 'git'
    executable_opts = [
        'clone', 'git@github.com:dganellari/mars.git'
        '&&', 'cd', 'mars'
        '&&', 'git', 'checkout', 'hilbert'
    ]
    local = True

    @sanity_function
    def validate_download(self):
        return sn.assert_eq(self.job.exitcode, 0)

class mars_build(rfm.CompileOnlyRegressionTest):
    descr = 'Build MARS'
    valid_systems = ['+uenv']
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
            ' -DMARS_ENABLE_CUDA=ON',
            ' -DMARS_ENABLE_TESTING=ON',
            ' -DMARS_ENABLE_BENCHMARKS=OFF',
            '-DMARS_ENABLE_CXXOPTS=ON',
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
    valid_systems = ['*']
    valid_prog_environs = ['+mars-dev']
    time_limit = '5m'
    maintainers = ['dganellari']
    mars_build = fixture(mars_build, scope='environment')

    @run_before('run')
    def prepare_run(self):
        self.executable = os.path.join(self.arbor_build.stagedir,
                                       'build', 'bin', 'unit')
        self.executable_opts = []

    @sanity_function
    def validate_test(self):
        return sn.assert_found(r'PASSED', self.stdout)

class testBase(rfm.RunOnlyRegressionTest):
    valid_systems = ['+uenv']
    valid_prog_environs = ['+mpi']
    target_executable = variable(str, value='ctest')
    maintainers = ['dganellari']

    test = ''
    region = []
    is_serial_test = True
    use_multithreading = False
    refs = {}

    quicc_build = fixture(quicc_build, scope='environment')

    @run_before('run')
    def set_num_tasks(self):
        """Set num tasks based on machine"""
        proc = self.current_partition.processor
        self.num_tasks_per_node = proc.num_cores
        self.num_tasks = 1
        self.time_limit = '30m'

    @run_before('run')
    def set_exec(self):
        self.executable = self.target_executable
        self.executable_opts = [
            '--test-dir', f'{self.mars_build.stagedir}/build', '-V', '-R', self.test
        ]

    @sanity_function
    def assert_sanity(self):
        return sn.assert_found(r'All tests passed', self.stdout)


