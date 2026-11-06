"""Reuse the audited preparation's binaries and successful solves without rerunning them."""
import json
from pathlib import Path
import shutil
from common import CONFIG_SHA, HYPRE, HYPRE_SHA, digest, save, verify_bundle, verify_files

HERE = Path(__file__).resolve().parent


def protected(name):
    return name.startswith(('source/', 'dependencies/', 'generated/', 'hcoarse/private/')) or name in (
        'hcoarse/check.hpp', 'hcoarse/extra.cmake', 'kit/build.py', 'kit/rank.py')


def check_sources(old, bundle, expected):
    before, after = verify_bundle(old), verify_bundle(bundle)
    if digest(old / 'manifest.json') != expected:
        raise ValueError('Original source snapshot changed')
    first = {k:v for k,v in before['files'].items() if protected(k)}
    second = {k:v for k,v in after['files'].items() if protected(k)}
    if first != second:
        raise ValueError('Reusing binaries requires identical source, dependencies and PTX')


def verify_reused(root, name):
    record = json.loads((root / 'reused.json').read_text())
    if name not in record['solves']:
        raise ValueError('An existing result was not approved for reuse')
    files = {n:v for n,v in record['files'].items() if n.startswith(name + '/')}
    if not files:
        raise ValueError('Missing reused evidence')
    verify_files(root, files)


def restore(bundle, previous, root):
    record = json.loads((HERE / 'resume_receipt.json').read_text())
    if str(previous) != record['previous']:
        raise ValueError('Resume requires the audited preparation5279173')
    check_sources(Path(record['bundle']), bundle, record['manifest_sha256'])
    if digest(previous / 'build/build.json') != record['build_sha256']:
        raise ValueError('Original build receipt changed')
    build = json.loads((previous / 'build/build.json').read_text())
    if build['bundle_sha256'] != record['manifest_sha256'] or build['binaries'] != record['binaries']:
        raise ValueError('Build does not match the audited input snapshot')
    if digest(HYPRE) != HYPRE_SHA or digest(HYPRE.parent.parent / 'include/HYPRE_config.h') != CONFIG_SHA:
        raise ValueError('Runtime Hypre library or headers changed')
    verify_files(previous, record['files'])
    verify_files(previous / 'build', record['binaries'])
    for name in record['files']:
        if name == 'inputs.json':
            continue
        target = root / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(previous / name), str(target))
    for name in record['binaries']:
        target = root / 'build' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(previous / 'build' / name), str(target))
    verify_files(root, {n:v for n,v in record['files'].items() if n != 'inputs.json'})
    verify_files(root / 'build', record['binaries'])
    save(root / 'previous-inputs.json', json.loads((previous / 'inputs.json').read_text()))
    save(root / 'reused.json', dict(record, recompiled=False, rerun_completed_solves=False,
                                  resumed_manifest_sha256=digest(bundle / 'manifest.json')))
    print('Verified and copied six binaries and 30 successful solves; no rebuild or repeat', flush=True)
