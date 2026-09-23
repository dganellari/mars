"""Exercise the real Cornerstone CMake setup twice without requiring a GPU toolkit."""
import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
CMAKE_FILE = ROOT / 'CMakeLists.txt'


def cornerstone_setup():
    text = CMAKE_FILE.read_text()
    start = text.index('# Find Cornerstone - centralize')
    end = text.index('set(MARS_MODULES', start)
    return text[start:end]


class CornerstoneReconfigure(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='mars-cmake-link-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / 'mars'
        self.build = self.root / 'build'
        self.cornerstone = self.root / 'cornerstone'
        for directory in (self.source / 'cmake', self.cornerstone / 'include/cstone/domain',
                          self.cornerstone / 'include/cstone', self.source / 'cuda-headers'):
            directory.mkdir(parents=True, exist_ok=True)
        module = ROOT / 'cmake/Findcornerstone.cmake'
        if not module.exists():
            module = ROOT / 'base/cmake/Findcornerstone.cmake'
        shutil.copyfile(str(module), str(self.source / 'cmake/Findcornerstone.cmake'))
        (self.source / 'cmake/FindCUDAToolkit.cmake').write_text('set(CUDAToolkit_FOUND TRUE)\n')
        (self.cornerstone / 'include/cstone/domain/domain.hpp').write_text('int cornerstone_value();\n')
        (self.cornerstone / 'CMakeLists.txt').write_text('''cmake_minimum_required(VERSION 3.22)
project(cornerstone-octree LANGUAGES CXX)
add_subdirectory(include/cstone)
''')
        (self.cornerstone / 'include/cstone/CMakeLists.txt').write_text('''add_library(cstone_gpu STATIC vector.cpp)
''')
        self.gpu_source = self.cornerstone / 'include/cstone/vector.cpp'
        self.gpu_source.write_text('int cornerstone_value() { return 17; }\n')
        (self.source / 'domain.cpp').write_text('int cornerstone_value(); int domain_value() { return cornerstone_value(); }\n')
        (self.source / 'mars.cpp').write_text('int mars_version() { return 1; }\n')
        (self.source / 'main.cpp').write_text('''#include <cstdlib>
int domain_value();
int main(int argc,char** argv) { return argc==2 && domain_value()==std::atoi(argv[1]) ? 0 : 1; }
''')
        self.write_project()

    def write_project(self, setup=None):
        # Use the actual common links, not a per-executable workaround in the fixture.
        backend = (ROOT / 'backend/distributed/unstructured/CMakeLists.txt').read_text()
        backend_link = backend[backend.index('# BUILD_INTERFACE only'):backend.index('# Export targets')]
        top = CMAKE_FILE.read_text()
        top_link = top[top.index('# Ensure DeviceVector instantiations'):top.index('add_custom_target(mars_complete)')]
        prefix = '''cmake_minimum_required(VERSION 3.22)
project(MarsLinkGate LANGUAGES CXX)
set(MARS_ENABLE_UNSTRUCTURED ON)
option(MARS_ENABLE_CUDA "Test CUDA dependency selection" ON)
option(MARS_ENABLE_HIP "Test HIP dependency selection" OFF)
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CUDAToolkit_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/cuda-headers")
set(HIP_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/cuda-headers")
foreach(library cudart cublas cusparse cusolver)
    add_library(CUDA::${library} INTERFACE IMPORTED GLOBAL)
endforeach()
add_library(mars::rocmlibs INTERFACE IMPORTED GLOBAL)
'''
        suffix = '''
if(TARGET cstone_gpu)
    get_target_property(imported cstone_gpu IMPORTED)
    file(WRITE "${CMAKE_BINARY_DIR}/gpu-target.txt" "${imported}")
else()
    file(WRITE "${CMAKE_BINARY_DIR}/gpu-target.txt" "missing")
endif()
if(MARS_ENABLE_CUDA OR MARS_ENABLE_HIP)
    add_library(mars_unstructured STATIC domain.cpp)
    add_library(mars STATIC mars.cpp)
    target_link_libraries(mars PUBLIC mars_unstructured)
'''+backend_link+top_link+'''
    add_executable(full_mars main.cpp)
    target_link_libraries(full_mars PRIVATE mars)
    add_executable(direct_backend main.cpp)
    target_link_libraries(direct_backend PRIVATE mars_unstructured)
endif()
'''
        (self.source / 'CMakeLists.txt').write_text(prefix+(setup if setup is not None else cornerstone_setup())+suffix)

    def run_command(self, command, success=True):
        result = subprocess.run([str(arg) for arg in command], stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, universal_newlines=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
        return result.stdout

    def configure(self, *options, **kwargs):
        return self.run_command(['cmake', '-S', self.source, '-B', self.build,
                                '-DFETCHCONTENT_SOURCE_DIR_CORNERSTONE_FETCH='+str(self.cornerstone)]
                               +list(options), **kwargs)

    def build_and_run_all(self, expected):
        self.run_command(['cmake', '--build', self.build, '--parallel', '2'])
        for executable in ('full_mars', 'direct_backend'):
            self.run_command([self.build / executable, str(expected)])

    def repeated_source_build(self, *options):
        self.configure(*options)
        self.build_and_run_all(17)
        for value in (29, 43):
            # Older Make versions compare timestamps at whole-second resolution.
            time.sleep(1.05)
            self.gpu_source.write_text('int cornerstone_value() { return %d; }\n' % value)
            self.configure()
            self.assertEqual((self.build / 'gpu-target.txt').read_text(), 'FALSE')
            self.build_and_run_all(value)

    def test_reconfigure_rebuilds_all_consumers(self):
        self.repeated_source_build()

    def test_reconfigure_hip_selection(self):
        self.repeated_source_build('-DMARS_ENABLE_CUDA=OFF', '-DMARS_ENABLE_HIP=ON')

    def test_default_fetch_layout_recovers_existing_cache(self):
        # No source override on the repair configure: mirrors the affected Daint cache.
        source = self.build / '_deps/cornerstone_fetch-src'
        source.parent.mkdir(parents=True)
        shutil.copytree(str(self.cornerstone), str(source))
        self.cornerstone = source
        self.gpu_source = source / 'include/cstone/vector.cpp'
        self.configure('-DFETCHCONTENT_SOURCE_DIR_CORNERSTONE_FETCH=',
                       '-Dcornerstone-octree_SOURCE_DIR='+str(source),
                       '-Dcornerstone-octree_BINARY_DIR='+str(self.build / '_deps/cornerstone_fetch-build'),
                       '-DCORNERSTONE_INCLUDE_DIRS='+str(source / 'include'))
        self.build_and_run_all(17)
        time.sleep(1.05)
        self.gpu_source.write_text('int cornerstone_value() { return 29; }\n')
        self.configure('-DFETCHCONTENT_SOURCE_DIR_CORNERSTONE_FETCH=')
        self.assertEqual((self.build / 'gpu-target.txt').read_text(), 'FALSE')
        self.build_and_run_all(29)

    def test_installed_library_takes_precedence(self):
        self.configure()
        self.build_and_run_all(17)
        prefix = self.root / 'installed'
        shutil.copytree(str(self.cornerstone / 'include'), str(prefix / 'include'))
        (prefix / 'lib').mkdir()
        libraries = list(self.build.glob('**/libcstone_gpu.a'))
        self.assertEqual(len(libraries), 1)
        library = prefix / 'lib/libcstone_gpu.a'
        shutil.copyfile(str(libraries[0]), str(library))
        self.gpu_source.write_text('#error must use the explicitly installed library\n')
        self.configure('-DCORNERSTONE_INSTALL_DIR='+str(prefix), '-DCORNERSTONE_GPU_LIBRARY='+str(library))
        self.assertEqual((self.build / 'gpu-target.txt').read_text(), 'TRUE')
        self.build_and_run_all(17)
        self.configure()
        self.build_and_run_all(17)

    def test_headers_without_gpu_library_fail_during_configure(self):
        output = self.configure('-DCORNERSTONE_INSTALL_DIR='+str(self.cornerstone), success=False)
        self.assertIn('requires the Cornerstone GPU library', output)

    def test_cpu_headers_do_not_require_gpu_library(self):
        self.configure('-DMARS_ENABLE_CUDA=OFF', '-DCORNERSTONE_INSTALL_DIR='+str(self.cornerstone))
        self.assertEqual((self.build / 'gpu-target.txt').read_text(), 'missing')
        self.run_command(['cmake', '--build', self.build])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cmake-file', type=Path, default=CMAKE_FILE)
    args, remaining = parser.parse_known_args()
    CMAKE_FILE = args.cmake_file
    unittest.main(argv=[__file__]+remaining)
