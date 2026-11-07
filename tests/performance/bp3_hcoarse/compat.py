"""Check the accepted h hierarchy with the newer constant-memory operator ABI."""
import json
import os
from pathlib import Path
import re
import sys

import integrate as i
from common import CONFIG_SHA, HYPRE, HYPRE_SHA, clean_environment, command, digest, save, step, verify_bundle, verify_files

HERE = Path(__file__).resolve().parent


def build(bundle, root):
    if digest(HYPRE) != HYPRE_SHA or digest(HYPRE.parent.parent / 'include/HYPRE_config.h') != CONFIG_SHA:
        raise ValueError('Runtime Hypre changed')
    out = root / 'build'
    out.mkdir()
    directory = out / 'on'
    command(i.standard_build.configure_args(bundle, directory, True) +
            ['-DCMAKE_PROJECT_TOP_LEVEL_INCLUDES=' + str(bundle / 'hcoarse/extra.cmake')],
            out / 'configure.log', timeout=600)
    command(['cmake', '--build', directory, '--parallel', '4', '--target', 'mars_marsir_bp', 'bp3_h_solver'],
            out / 'compile.log', timeout=5400)
    files = {}
    for mode in ('on', 'check'):
        executable = i.binary(root, mode)
        link = command(['ldd', '-r', executable], out / (mode + '-linkage.txt'))
        libs = re.findall(r'libHYPRE[^\s]*\s+=>\s+(\S+)', link)
        if 'not found' in link or 'undefined symbol:' in link or len(libs) != 1 or Path(libs[0]).resolve() != HYPRE.resolve():
            raise ValueError('Unexpected runtime linkage')
        files[str(executable.relative_to(out))] = digest(executable)
    db = json.loads((directory / 'compile_commands.json').read_text())
    for entry in db:
        flags = entry.get('command', ' '.join(entry.get('arguments', [])))
        private = Path(entry['file']) == bundle / 'hcoarse/private/mars_marsir_bp.cu'
        if ('MARS_BP3_HCOARSE_CHECK' in flags) != private:
            raise ValueError('Check instrumentation escaped its target')
    save(out / 'build.json', dict(binaries=files, manifest_sha256=digest(bundle / 'manifest.json'),
                                  hypre_sha256=HYPRE_SHA, config_sha256=CONFIG_SHA))


def check(bundle, root):
    command(step(1, 100) + ['/usr/bin/python3', HERE / 'compat.py', 'build', bundle, root], root / 'build.log')
    rows = []
    for ranks, p, ne in ((1, 1, 8), (4, 1, 8), (1, 3, 16), (1, 7, 16), (4, 3, 32), (4, 7, 32)):
        rows.append(i.solve(bundle, root, root, 'fallback_p%d_r%d' % (p, ranks), 'on', ranks, p, ne, automatic=True))
    for p in (3, 7):
        rows.append(i.solve(bundle, root, root, 'single_check_p%d' % p, 'check', 1, p, 32, 2, checked=True))
        for depth in (1, 2, 3):
            normal = i.solve(bundle, root, root, 'normal_p%d_h%d' % (p, depth), 'on', 4, p, 32, depth)
            checked = i.solve(bundle, root, root, 'check_p%d_h%d' % (p, depth), 'check', 4, p, 32, depth, checked=True)
            if not i.checked_equivalent(normal, checked):
                raise ValueError('Diagnostic and normal hierarchy disagree')
            rows += [normal, checked]
    for p, ne, depth, solves in ((3, 128, 2, 3), (7, 64, 1, 5)):
        normal = i.solve(bundle, root, root, 'auto_p%d_r8' % p, 'on', 8, p, ne, depth, automatic=True, solves=solves)
        explicit = i.solve(bundle, root, root, 'explicit_p%d_r8' % p, 'on', 8, p, ne, depth, solves=solves)
        if not i.equivalent([normal, explicit]):
            raise ValueError('Automatic and explicit choices disagree')
        rows += [normal, explicit]
    verify_bundle(bundle)
    verify_files(root / 'build', json.loads((root / 'build/build.json').read_text())['binaries'])
    files = {str(p.relative_to(root)): digest(p) for p in root.rglob('*') if p.is_file() and
             p.suffix in ('.json', '.log', '.txt') and 'CMakeFiles' not in p.parts and '_deps' not in p.parts}
    save(root / 'passed.json', dict(passed=True, solves=len(rows), direct_checks=8,
        max_true=max(r['true_residual'] for r in rows), max_relative=max(r['relative'] for r in rows),
        max_direct=max(r.get('h_check_worst') or 0 for r in rows), evidence=files,
        manifest_sha256=digest(bundle / 'manifest.json'),
        scope='Combined-source correctness and ABI compatibility; not a new scaling measurement'))
    print('PASS combined operator/h hierarchy: 24 solves and eight direct GPU checks', flush=True)


if __name__ == '__main__':
    if len(sys.argv) != 4 or sys.argv[1] not in ('build', 'check'):
        raise SystemExit('usage: compat.py build|check bundle results')
    clean_environment()
    bundle, root = map(lambda x: Path(x).resolve(), sys.argv[2:])
    verify_bundle(bundle)
    if sys.argv[1] == 'build':
        build(bundle, root)
    else:
        if int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) < 2:
            raise ValueError('Compatibility check needs two nodes')
        root.mkdir(parents=True, exist_ok=False)
        save(root / 'inputs.json', dict(manifest_sha256=digest(bundle / 'manifest.json'), phase='compatibility'))
        check(bundle, root)
