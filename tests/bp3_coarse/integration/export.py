"""Freeze code and pinned public dependencies without reading meshes or private notes."""
import argparse
import json
from pathlib import Path
import shutil
import subprocess
from common import digest, save, verify_files

ROOTS = ('backend/', 'base/', 'cmake/', 'core/', 'examples/', 'io/', 'kokkos/', 'scripts/', 'tests/')
SUFFIXES = {'.cmake', '.in', '.hpp', '.h', '.hxx', '.cpp', '.cc', '.c', '.cu', '.cuh', '.py', '.sh'}
NEW = ('backend/distributed/unstructured/solvers/mars_cellwise_amg_consolidate.cu',
       'backend/distributed/unstructured/solvers/mars_cellwise_amg_consolidate.hpp',
       'cmake/MarsBP3Coarse.cmake', 'cmake/MarsBP3CoarseSupport.cmake', 'tests/bp3_coarse/fixture.cu', 'tests/bp3_coarse/register.cmake')
REVISIONS = {'cornerstone': '2a95e63b0ca5454712ab83f27505055c078887dc',
             'googletest': '15460959cbbfa20e66ef0b5ab497367e47fc0a04'}


def eligible(name):
    path = Path(name)
    if path.is_absolute() or '..' in path.parts or any(p.startswith('.') for p in path.parts):
        return False
    if 'pump' in name.lower() or 'data' in path.parts:
        return False
    if name in ('CMakeLists.txt', 'LICENSE.md'):
        return True
    return name.startswith(ROOTS) and (path.name == 'CMakeLists.txt' or path.suffix in SUFFIXES)


def copy(source, target):
    if source.is_symlink() or not source.is_file():
        raise ValueError('Expected a regular file: ' + str(source))
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(source), str(target))


def export(repo, deps, accepted, target):
    revision = subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], universal_newlines=True).strip()
    base = 'df67598afa97c8f80ba288e7e59095d94d8cad19'
    lineage = subprocess.run(['git', '-C', str(repo), 'merge-base', '--is-ancestor', base, revision],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if lineage.returncode:
        raise ValueError('Rebase needs a new integration review')
    target.mkdir(parents=True, exist_ok=False)
    names = subprocess.check_output(['git', '-C', str(repo), 'ls-files', '-z']).decode().split('\0')
    names = sorted({n for n in names if eligible(n)} | set(NEW))
    for name in names:
        for parent in (repo / name).parents:
            if parent == repo:
                break
            if parent.is_symlink():
                raise ValueError('Symlinked source parent')
        copy(repo / name, target / 'source' / name)
    for name, expected in REVISIONS.items():
        source = deps / name
        record = json.loads((source / 'dependency-manifest.json').read_text())
        if record['revision'] != expected:
            raise ValueError('Wrong dependency revision: ' + name)
        verify_files(source, record['files'])
        for path in record['files']:
            copy(source / path, target / 'dependencies' / name / path)
        copy(source / 'dependency-manifest.json', target / 'dependencies' / name / 'dependency-manifest.json')
    manifest = json.loads((accepted / 'source/source.json').read_text())
    for name, sha in manifest['files'].items():
        if name.startswith('generated/'):
            path = accepted / 'source' / name
            if digest(path) != sha:
                raise ValueError('Accepted PTX changed')
            copy(path, target / name)
    here = Path(__file__).resolve().parent
    for source in sorted(here.iterdir()):
        if source.suffix in ('.py', '.sh') and source.name not in ('test_integration.py', 'export.py'):
            copy(source, target / 'kit' / source.name)
    files = {str(p.relative_to(target)): digest(p) for p in target.rglob('*') if p.is_file()}
    save(target / 'manifest.json', dict(version=1, base_revision=revision, dependencies=REVISIONS, files=files,
         scope='Code-only standard MARS build, coarse integration changes, accepted PTX; no meshes or private notes'))
    print('Exported {} files'.format(len(files)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('repo', 'deps', 'accepted', 'target'):
        p.add_argument('--' + name, type=Path, required=True)
    a = p.parse_args()
    export(a.repo.resolve(), a.deps.resolve(), a.accepted.resolve(), a.target.resolve())
