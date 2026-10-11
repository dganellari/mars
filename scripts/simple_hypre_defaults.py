#!/usr/bin/env python3
"""Query fresh Hypre objects using the library recorded in an OpenAccel capture."""

import argparse
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

from simple_pressure_settings import saved_launch
from simple_public_diagnostics import SafeParser
from simple_snapshot_compare import digest, require
import simple_startup_probe as startup


SCHEMA = 'mars-simple-hypre-defaults-v1'
RAW_SCHEMA = 'mars-hypre-defaults-raw-v1'
REAL_KEYS = frozenset(('amg_strong_threshold', 'amg_trunc_factor',
                       'amg_jacobi_trunc_threshold', 'amg_max_row_sum', 'amg_agg_trunc_factor'))
INT_KEYS = frozenset(('amg_coarsen_type', 'amg_interp_type', 'amg_relax_order',
    'amg_p_max_elmts', 'amg_max_levels', 'amg_min_coarse_size', 'amg_max_coarse_size',
    'amg_num_functions', 'amg_coarsen_cut_factor', 'amg_cycle_type', 'amg_fcycle',
    'amg_relax_down', 'amg_relax_up', 'amg_relax_coarse',
    'amg_sweeps_down', 'amg_sweeps_up', 'amg_sweeps_coarse',
    'amg_agg_num_levels', 'amg_agg_interp_type', 'amg_num_paths', 'amg_keep_transpose',
    'amg_nodal', 'amg_nodal_diag', 'gmres_restart_dimension', 'gmres_minimum_iterations',
    'flexgmres_restart_dimension', 'flexgmres_minimum_iterations'))


def library_kind(path):
    name = Path(path).name.lower()
    if re.fullmatch(r'libhypre(?:-[0-9.]+)?\.so(?:\.[0-9]+)*', name):
        return 'hypre'
    if re.match(r'lib(?:mpi|mpich)(?:[_.-]|$)', name):
        return 'mpi'
    return None


def matching_install(libraries):
    matches = [Path(path) for path in libraries if library_kind(path) == 'hypre']
    require(len(matches) == 1)
    library = matches[0].resolve(strict=True)
    require(library.parent.name in ('lib', 'lib64'))
    include = library.parent.parent / 'include'
    for name in ('HYPRE.h', 'HYPRE_config.h', '_hypre_parcsr_ls.h'):
        require((include / name).is_file())
    return library, include


def unchanged(files):
    require(all(Path(path).is_file() and digest(path) == expected
                for path, expected in files.items()))


def find_compiler(requested):
    # Cray environments can provide CC without the mpicxx alias.
    for name in ([requested] if requested is not None else ('mpicxx', 'mpic++', 'mpiCC', 'CC')):
        found = shutil.which(name)
        if found is not None:
            return found
    return None


def cached_toolchain(path, parallel):
    cache = {}
    for line in path.read_text().splitlines():
        if line.startswith(('#', '//')) or '=' not in line or ':' not in line.split('=', 1)[0]:
            continue
        key, value = line.split('=', 1)
        cache[key.split(':', 1)[0]] = value
    compiler = cache['CMAKE_CXX_COMPILER']
    require(Path(compiler).is_absolute())
    inputs = {str(path): digest(path)}
    compile_flags, link_flags = [], []
    if parallel:
        includes = set()
        for key in ('MPI_CXX_COMPILER_INCLUDE_DIRS', 'MPI_CXX_ADDITIONAL_INCLUDE_DIRS',
                    'MPI_CXX_INCLUDE_DIRS', 'MPI_CXX_HEADER_DIR'):
            includes.update(value for value in cache.get(key, '').split(';') if value)
        require(includes and any((Path(value) / 'mpi.h').is_file() for value in includes))
        for value in sorted(includes):
            require(Path(value).is_absolute() and Path(value).is_dir())
            compile_flags.append('-I' + value)
            for name in ('mpi.h', 'mpio.h'):
                header = Path(value) / name
                if header.is_file():
                    inputs[str(header)] = digest(header)
        compile_flags += ['-D' + value for value in cache.get('MPI_CXX_COMPILE_DEFINITIONS', '').split(';') if value]
        compile_flags += [value for value in cache.get('MPI_CXX_COMPILE_OPTIONS', '').split(';') if value]
        link_flags += shlex.split(cache.get('MPI_CXX_LINK_FLAGS', ''))
        names = cache['MPI_CXX_LIB_NAMES'].split(';')
        require(names and all(re.fullmatch('[A-Za-z0-9_]+', name) for name in names))
        for name in names:
            library = Path(cache['MPI_' + name + '_LIBRARY'])
            require(library.is_absolute() and library.is_file())
            link_flags.append(str(library))
            inputs[str(library)] = digest(library)
        require(not any('$<' in flag or flag.startswith('SHELL:') for flag in compile_flags + link_flags))
    return compiler, compile_flags, link_flags, inputs


def probe_dependencies(libraries, captured, parallel):
    # Compiler/loader libraries may differ; Hypre and MPI must match the capture.
    expected = {kind: {value for path, value in captured.items() if library_kind(path) == kind}
                for kind in ('hypre', 'mpi')}
    found = {kind: {value for path, value in libraries.items() if library_kind(path) == kind}
             for kind in ('hypre', 'mpi')}
    require(len(found['hypre']) == 1 and found['hypre'] == expected['hypre'])
    require(found['mpi'] <= expected['mpi'] and (found['mpi'] or not parallel))


def raw_output(path):
    # Launchers and binding scripts may also write to stdout. Keep their text private.
    lines = [line for line in path.read_text().splitlines()
             if line.startswith('{"schema":"' + RAW_SCHEMA + '"')]
    require(len(lines) == 1)
    return json.loads(lines[0])


def public_defaults(raw):
    require(raw['schema'] == RAW_SCHEMA)
    require(raw['header_version_matches'] is True and raw['getter_layout_checks_passed'] is True)
    require(isinstance(raw['version'], str) and re.fullmatch(r'[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}', raw['version']))
    require(type(raw['build_cuda']) is bool)
    require(type(raw['memory_location']) is int and raw['memory_location'] in (0, 1, 2))
    require(type(raw['execution_policy']) is int and raw['execution_policy'] in (-1, 0, 1))
    values = raw['defaults']
    require(isinstance(values, dict) and set(values) == INT_KEYS | REAL_KEYS)
    result = {}
    for key in sorted(values):
        value = values[key]
        require(type(value) in (int, float) and math.isfinite(value))
        if key in INT_KEYS:
            require(value == int(value) and -2**31 <= value < 2**31)
            value = int(value)
        result[key] = value
    # Nothing else from the raw log is eligible for publication.
    return dict(version=raw['version'], build_cuda=raw['build_cuda'],
                memory_location=raw['memory_location'], execution_policy=raw['execution_policy'],
                header_version_matches=True, getter_layout_checks_passed=True, defaults=result)


def inspect(pair, output, compiler, launcher, include_dirs=(), build_cache=None):
    public = dict(schema=SCHEMA, comparison_status='invalid_evidence', failed_check='saved_capture',
                  scope='fresh_objects_in_captured_hypre_library_no_matrix_or_solve',
                  application_effective_settings_verified=False,
                  matrix_dependent_amg_setup_verified=False, identical_linear_solvers_verified=False)
    try:
        pair = pair.resolve()
        record = saved_launch(pair, 'openaccel')
        require(record['status'] == 'finished')
        public['failed_check'] = 'capture_libraries'
        captured = record['libraries']
        identity = dict(captured, **{record['executable']: record['executable_sha256']})
        unchanged(identity)
        public['failed_check'] = 'installed_headers'
        library, include = matching_install(captured)
        headers = {str(path): digest(path) for path in include.glob('*.h')}
        config = (include / 'HYPRE_config.h').read_text()
        parallel = not re.search(r'^\s*#\s*define\s+HYPRE_SEQUENTIAL\b', config, re.M)
        compile_flags, link_flags, build_inputs = [], [], {}
        public['compiler_configuration'] = 'wrapper_lookup'
        if build_cache is not None:
            public['failed_check'] = 'recorded_build_configuration'
            require(compiler is None)
            compiler, compile_flags, link_flags, build_inputs = cached_toolchain(build_cache.resolve(), parallel)
            public['compiler_configuration'] = 'reference_build_cache'
        public['failed_check'] = 'compiler_unavailable'
        compiler = find_compiler(compiler)
        public['compiler_available'] = compiler is not None
        require(compiler is not None)
        public['failed_check'] = 'launcher_unavailable'
        public['launcher_available'] = bool(launcher and shutil.which(launcher[0]) is not None)
        require(public['launcher_available'])
        public['failed_check'] = 'probe_compile'
        source = Path(__file__).resolve().parent.parent / 'tests/reference/openaccel/simple_performance/hypre_defaults.cpp'
        binary = output / 'probe'
        command = [compiler, '-std=c++11', '-fPIC', '-I' + str(include)]
        command += compile_flags
        command += ['-I' + str(path.resolve(strict=True)) for path in include_dirs]
        command += [str(source), str(library)] + link_flags
        command += ['-Wl,-rpath,' + str(library.parent), '-ldl', '-o', str(binary)]
        launch = list(launcher) + [str(binary)]
        metadata = dict(capture=record, compiler_command=command, launch_command=launch,
                        source_sha256=digest(source), headers=headers, build_inputs=build_inputs)
        startup.write_json(output / 'private-provenance.json', metadata)
        with (output / 'compile.log').open('xb') as log:
            built = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        (output / 'compile.exit').write_text(str(built.returncode) + '\n')
        require(built.returncode == 0 and binary.is_file())
        public['failed_check'] = 'probe_libraries'
        libraries = startup.runtime_libraries(binary)
        probe_dependencies(libraries, captured, parallel)
        probe_hash = digest(binary)
        startup.write_json(output / 'private-probe-libraries.json', libraries)
        public['failed_check'] = 'probe_launch'
        with (output / 'probe.log').open('xb') as log:
            ran = subprocess.run(launch, cwd=str(output), stdout=log, stderr=subprocess.STDOUT)
        (output / 'probe.exit').write_text(str(ran.returncode) + '\n')
        require(ran.returncode == 0)
        public['failed_check'] = 'probe_output'
        raw = raw_output(output / 'probe.log')
        projection = public_defaults(raw)
        public['failed_check'] = 'loaded_library_identity'
        require(Path(raw['library_path']).is_absolute())
        require(digest(raw['library_path']) == digest(library))
        require(isinstance(raw['mpi_library_path'], str))
        if parallel:
            require(Path(raw['mpi_library_path']).is_absolute())
            require(digest(raw['mpi_library_path']) in
                    {value for path, value in captured.items() if library_kind(path) == 'mpi'})
        else:
            require(raw['mpi_library_path'] == '')
        public['failed_check'] = 'inputs_changed'
        require(record == saved_launch(pair, 'openaccel'))
        unchanged(identity)
        unchanged(headers)
        unchanged(build_inputs)
        unchanged(libraries)
        require(probe_hash == digest(binary) and metadata['source_sha256'] == digest(source))
        public.update(projection)
        public.update(comparison_status='completed', failed_check='none',
                      saved_capture_identity_verified=True, captured_libraries_unchanged=True,
                      loaded_hypre_library_matches_capture=True, mpi_library_identity_verified=bool(parallel),
                      mpi_enabled=bool(parallel))
    except Exception:
        # Paths, captured arguments and build errors stay in the private directory.
        pass
    return public


def main(argv=None):
    parser = SafeParser(description=__doc__)
    parser.add_argument('--pair', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    toolchain = parser.add_mutually_exclusive_group()
    toolchain.add_argument('--cxx')
    toolchain.add_argument('--build-cache', type=Path)
    parser.add_argument('--include-dir', type=Path, action='append', default=[])
    parser.add_argument('--launcher', nargs=argparse.REMAINDER, required=True)
    args = parser.parse_args(argv)
    os.umask(0o077)
    try:
        args.output_dir.mkdir(mode=0o700)
    except OSError:
        print('ERROR: a new output directory is required; existing files were not replaced.', file=sys.stderr)
        return 1
    output = args.output_dir.resolve()
    result = inspect(args.pair, output, args.cxx, args.launcher, args.include_dir, args.build_cache)
    startup.write_json(output / 'public.json', result)
    print('Hypre defaults inspection written. Share only public.json; no flow solver was launched.')
    return 0 if result['comparison_status'] == 'completed' else 1


if __name__ == '__main__':
    sys.exit(main())
