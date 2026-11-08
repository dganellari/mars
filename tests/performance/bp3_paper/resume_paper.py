"""Reuse the verified MARS checks and MFEM installation; build only the missing driver."""
import json
from pathlib import Path
import shutil
from common import CONFIG_SHA, HYPRE, HYPRE_SHA, digest, save, verify_bundle, verify_files


def protected(name):
    if name.startswith('source/tests/performance/bp3_paper/') or name == 'source/tests/performance/mfem_bp3/CMakeLists.txt':
        return False
    return name.startswith(('source/', 'generated/', 'dependencies/', 'hcoarse/private/')) or name in (
        'hcoarse/extra.cmake', 'hcoarse/check.hpp', 'kit/build.py', 'kit/rank.py')


def restore(bundle, previous, root):
    receipt = json.loads((bundle / 'paper/resume-receipt.json').read_text())
    if str(previous) != receipt['previous']:
        raise ValueError('Wrong previous preparation')
    original = Path(receipt['bundle'])
    if digest(original / 'manifest.json') != receipt['manifest_sha256']:
        raise ValueError('Original source manifest changed')
    old, new = verify_bundle(original), verify_bundle(bundle)
    if {n:s for n,s in old['files'].items() if protected(n)} != {n:s for n,s in new['files'].items() if protected(n)}:
        raise ValueError('Reusing binaries requires unchanged sources, dependencies and operators')
    if digest(HYPRE) != HYPRE_SHA or digest(HYPRE.parent.parent / 'include/HYPRE_config.h') != CONFIG_SHA:
        raise ValueError('Runtime Hypre changed')
    verify_files(previous, receipt['files'])
    verify_files(previous / 'build', receipt['binaries'])
    passed = json.loads((previous / 'passed.json').read_text())
    if not passed['passed'] or passed['solves'] != 24 or passed['direct_checks'] != 8:
        raise ValueError('Missing completed operator checks')
    verify_files(previous, passed['evidence'])
    for name in receipt['files']:
        target = root / ('previous-' + name if name in ('inputs.json', 'passed.json') else name)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(previous / name), str(target))
    for name in receipt['binaries']:
        target = root / 'build' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(previous / 'build' / name), str(target))
    copied = {('previous-' + n if n in ('inputs.json', 'passed.json') else n): s for n,s in receipt['files'].items()}
    verify_files(root, copied)
    verify_files(root / 'build', receipt['binaries'])
    save(root / 'reused.json', dict(receipt, resumed_manifest_sha256=digest(bundle / 'manifest.json'),
                                  rerun_completed_solves=False, rebuilt_mars=False, rebuilt_mfem_library=False))
    print('Verified reuse: two MARS binaries, 24 solves and installed MFEM library', flush=True)
