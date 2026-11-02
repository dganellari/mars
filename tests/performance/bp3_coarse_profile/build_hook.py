#!/usr/bin/env python3
"""Build annotations against the exact shared Hypre loaded by the solver."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess


def output(*args):
    return subprocess.check_output(args, universal_newlines=True, stderr=subprocess.STDOUT)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('executable', type=Path)
    parser.add_argument('output_dir', type=Path)
    args = parser.parse_args()
    dest = args.output_dir.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    binary = args.executable.resolve(strict=True)
    linkage = output('ldd', str(binary))
    (dest / 'linkage.txt').write_text(linkage)
    found = re.findall(r'libHYPRE[^\s]*\s+=>\s+(\S+)', linkage)
    if len(found) != 1 or not Path(found[0]).is_file():
        raise SystemExit('Need exactly one dynamically linked Hypre; static linkage cannot use this probe.')
    lib = Path(found[0]).resolve(strict=True)
    include = lib.parent.parent / 'include'
    config = (include / 'HYPRE_config.h').read_text()
    if not re.search(r'^#define HYPRE_RELEASE_VERSION "2\.33\.0"', config, re.M):
        raise SystemExit('This experiment is checked for Hypre 2.33.0 only; inspect a different version first.')
    if re.search(r'^#define HYPRE_USING_(NVTX|ROCTX)\b', config, re.M):
        raise SystemExit('Native Hypre profiling adds synchronization; use the original profiling-disabled build.')
    relocations = output('readelf', '-rW', str(lib))
    for name in ('hypre_GpuProfilingPushRange', 'hypre_GpuProfilingPopRange'):
        if not any(name in line and ('JUMP_SLOT' in line or 'GLOB_DAT' in line) for line in relocations.splitlines()):
            raise SystemExit(f'{name} is not dynamically interposable in this library; no trace will be claimed.')
    roots = [Path(os.environ[k]) for k in ('NVTX_INCLUDE',) if k in os.environ]
    for key in ('CUDA_HOME', 'CUDA_PATH'):
        if key in os.environ: roots.append(Path(os.environ[key]) / 'include')
    for tool in ('nvcc', 'nsys'):
        location = shutil.which(tool)
        if location:
            root = Path(location).resolve().parent.parent
            roots += [root / 'include', root / 'targets/sbsa-linux/include',
                      root / 'targets/aarch64-linux/include', root / 'target-linux-sbsa/nvtx/include']
    nvtx = next((p for p in roots if (p / 'nvtx3/nvToolsExt.h').is_file()), None)
    if nvtx is None:
        raise SystemExit('Set NVTX_INCLUDE to the directory containing nvtx3/nvToolsExt.h from CUDA or Nsight.')
    source = Path(__file__).with_name('hypre_ranges.cpp')
    command = [os.environ.get('MPICXX', 'mpicxx'), '-std=c++17', '-O2', '-fPIC', '-shared',
               '-Wall', '-Wextra', '-Werror', '-I'+str(include), '-I'+str(nvtx),
               str(source), '-ldl', '-o', str(dest / 'libbp3_ranges.so')]
    subprocess.run(command, check=True)
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    manifest = {'binary': str(binary), 'binary_sha256': digest(binary), 'hypre': str(lib),
                'hypre_sha256': digest(lib), 'hypre_config_sha256': digest(include / 'HYPRE_config.h'),
                'hook_source_sha256': digest(source), 'compile_command': command,
                'nsys_version': output('nsys', '--version')}
    (dest / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    (dest / 'HYPRE_config.h').write_text(config)
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__': main()
