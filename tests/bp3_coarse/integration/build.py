"""Build explicit OFF and automatic targets through MARS's top-level CMake."""
import argparse
import json
from pathlib import Path
import re
import shutil
from common import ACCEPTED, ACCEPTED_SHA, CALLBACKS, CONFIG_SHA, HYPRE, HYPRE_SHA, command, digest, save, step, verify_bundle


def mpi_args():
    flags = [arg for arg in step(1, 3)[1:] if arg != '-n1']
    return ['-DMPIEXEC_EXECUTABLE=/usr/bin/srun', '-DMPIEXEC_NUMPROC_FLAG=-n',
            '-DMPIEXEC_PREFLAGS=' + ';'.join(flags)]


def configure_args(bundle, build, enabled):
    prefix = HYPRE.parent.parent
    return ['cmake', '-S', str(bundle / 'source'), '-B', str(build),
            '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON',
            '-DCMAKE_CXX_COMPILER=' + str(shutil.which('mpicxx')), '-DCMAKE_CUDA_HOST_COMPILER=' + str(shutil.which('mpicxx')),
            '-DCMAKE_CUDA_COMPILER=' + str(shutil.which('nvcc')), '-DCMAKE_CUDA_ARCHITECTURES=90',
            '-DCMAKE_CUDA_RUNTIME_LIBRARY=Shared',
            '-DMARS_ENABLE_KOKKOS=OFF', '-DMARS_ENABLE_CUDA=ON', '-DMARS_ENABLE_HIP=OFF',
            '-DMARS_ENABLE_UNSTRUCTURED=ON', '-DMARS_ENABLE_HYPRE=ON', '-DMARS_ENABLE_MARSIR=ON',
            '-DMARS_ENABLE_FEM_EXAMPLES=ON', '-DMARS_ENABLE_MPI=ON', '-DMARS_ENABLE_TESTS=ON',
            '-DMARS_ENABLE_DEV_EXAMPLES=OFF',
            '-DHYPRE_INCLUDE_DIR=' + str(prefix / 'include'), '-DHYPRE_LIBRARY=' + str(HYPRE),
            '-DFETCHCONTENT_SOURCE_DIR_CORNERSTONE_FETCH=' + str(bundle / 'dependencies/cornerstone'),
            '-DFETCHCONTENT_SOURCE_DIR_GOOGLETEST=' + str(bundle / 'dependencies/googletest'),
            '-DCMAKE_DISABLE_FIND_PACKAGE_GTest=ON', '-DFETCHCONTENT_FULLY_DISCONNECTED=ON'] + mpi_args() + (
                [] if enabled else ['-DMARS_ENABLE_BP3_COARSE_CONSOLIDATION=OFF'])


def build(bundle, out, extra_on=()):
    verify_bundle(bundle)
    if not shutil.which('mpicxx') or not shutil.which('nvcc'):
        raise ValueError('The pinned CUDA uenv must expose mpicxx and nvcc')
    if digest(HYPRE) != HYPRE_SHA or digest(HYPRE.parent.parent / 'include/HYPRE_config.h') != CONFIG_SHA:
        raise ValueError('Hypre library or headers changed')
    out.mkdir(parents=True, exist_ok=False)
    control = ACCEPTED / 'build/bp3_solver'
    if digest(control) != ACCEPTED_SHA:
        raise ValueError('Accepted executable changed')
    link = command(['ldd', '-r', control], out / 'accepted-linkage.txt')
    libraries = re.findall(r'libHYPRE[^\s]*\s+=>\s+(\S+)', link)
    if 'not found' in link or 'undefined symbol:' in link or len(libraries) != 1 or Path(libraries[0]).resolve() != HYPRE.resolve():
        raise ValueError('Accepted control resolves a different Hypre library')
    binaries = {}
    for mode in ('off', 'on'):
        directory = out / mode
        options = list(extra_on) if mode == 'on' else []
        command(configure_args(bundle, directory, mode == 'on') + options, out / (mode + '-configure.log'), timeout=600)
        targets = ['mars_marsir_bp', 'mars_marsir_cellwise_solve']
        if mode == 'on':
            targets.append('mars_bp3_coarse_test')
        command(['cmake', '--build', directory, '--parallel', '4', '--target'] + targets,
                out / (mode + '-compile.log'), timeout=6000)
        database = json.loads((directory / 'compile_commands.json').read_text())
        for entry in database:
            flags = entry.get('command', ' '.join(entry.get('arguments', [])))
            opted_in = 'MARS_BP3_COARSE_CONSOLIDATION' in flags
            expected = mode == 'on' and Path(entry['file']).name == 'mars_marsir_bp.cu'
            if opted_in != expected:
                raise ValueError('Coarse macro escaped its target: ' + entry['file'])
        for target in targets:
            binary = directory / 'examples/distributed/unstructured' / target
            link = command(['ldd', '-r', binary], out / (mode + '-' + target + '-linkage.txt'))
            libraries = re.findall(r'libHYPRE[^\s]*\s+=>\s+(\S+)', link)
            if 'not found' in link or 'undefined symbol:' in link or len(libraries) != 1 or Path(libraries[0]).resolve() != HYPRE.resolve():
                raise ValueError('Wrong runtime linkage: ' + target)
            symbols = command(['nm', '-D', '--defined-only', binary], out / (mode + '-' + target + '-symbols.txt'))
            expected = mode == 'on' and target in ('mars_marsir_bp', 'mars_bp3_coarse_test')
            for symbol in CALLBACKS:
                if bool(re.search(r'\bT ' + symbol + r'$', symbols, re.M)) != expected:
                    raise ValueError('Wrong callback scope: ' + target + ':' + symbol)
            binaries[str(binary.relative_to(out))] = digest(binary)
    verify_bundle(bundle)
    save(out / 'build.json', dict(bundle_sha256=digest(bundle / 'manifest.json'), binaries=binaries,
                                  hypre_sha256=HYPRE_SHA, config_sha256=CONFIG_SHA))
    print('PASS standard MARS builds, target isolation and runtime linkage', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bundle', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    build(args.bundle.resolve(), args.output.resolve())
