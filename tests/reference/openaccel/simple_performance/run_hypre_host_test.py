#!/usr/bin/env python3
"""Execute the GMRES wrapper with host-emulated CUDA and a real CPU Hypre library.

Only kernel launch syntax/includes are adapted in a generated, ignored header.
The wrapper's packing, cache, update, setup, solve and destruction code is used.
This is a lifecycle/numeric regression, not CUDA, distributed MPI or speed evidence.
"""
import argparse
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hypre-prefix', type=Path, required=True)
    parser.add_argument('--cxx', default='c++')
    args = parser.parse_args()
    build = ROOT / '.local-worktrees/simple-performance/hypre-host'
    build.mkdir(parents=True, exist_ok=True)
    directory = ROOT / 'backend/distributed/unstructured/solvers'
    pcg = (directory / 'mars_hypre_pcg_solver.hpp').read_text()
    kernels = pcg[pcg.index('// Pass 1:'):pcg.index('// GPU-resident Hypre PCG')]
    wrapper = (directory / 'mars_hypre_gmres_solver.hpp').read_text()
    # Observe the selection while still calling the installed Hypre API.
    wrapper = wrapper.replace('HYPRE_SetSpMVUseVendor(', 'host_set_spmv_use_vendor(')
    wrapper = wrapper.replace('hypre_ForceSyncComputeStream()', 'host_hypre_stream_sync()')
    # Keep Hypre includes so version-specific declarations are compiled too.
    wrapper = re.sub(r'^#include(?!\s+[<"](?:HYPRE|_hypre|mars_hypre_pressure_recovery\.hpp))[^\n]*\n', '', wrapper, flags=re.M)
    # Sequential internal headers alias MPI names; retain our instrumented stubs.
    mpi_names = ('Comm', 'COMM_WORLD', 'INT', 'DOUBLE', 'MAX', 'SUM', 'Comm_rank',
                 'Allreduce', 'Barrier', 'Abort', 'Wtime')
    mpi_restore = '\n'.join('#undef MPI_' + name for name in mpi_names)
    wrapper = wrapper.replace('\nnamespace mars {', '\n' + mpi_restore + '\nnamespace mars {', 1)
    source = ('#include "hypre_host_shim.hpp"\nnamespace mars::fem {\n' + kernels
              + '\n}\n' + wrapper)
    # Each launch is executed with the same grid, block and thread indices on CPU.
    pattern = r'(\w+(?:<[^;{}\n]*>)?)<<<(.*?)>>>\((.*?)\);'
    source, launches = re.subn(pattern, r'launch_host(\2, [&] { \1(\3); });', source, flags=re.S)
    if launches != 7 or '<<<' in source:
        raise RuntimeError(f'Kernel launch adaptation needs review: {launches} launches')
    (build / 'host_gmres.hpp').write_text(source)
    exe = build / 'hypre_numeric_refresh'
    command = [args.cxx, '-std=c++20', '-O1', '-g', '-fsanitize=address,undefined',
               '-fno-omit-frame-pointer', '-I' + str(HERE), '-I' + str(build), '-I' + str(directory),
               '-I' + str(args.hypre_prefix / 'include'),
               str(HERE / 'hypre_numeric_refresh.cpp'),
               str(args.hypre_prefix / 'lib/libHYPRE.a'), '-o', str(exe)]
    subprocess.run(command, check=True)
    subprocess.run([str(exe)], check=True)


if __name__ == '__main__':
    main()
