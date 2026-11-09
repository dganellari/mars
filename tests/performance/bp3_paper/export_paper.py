"""Freeze the final paper sources, current PTX and sealed public dependencies."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'bp3_hcoarse'))
from export_compat import export_compat
from common import digest, save, verify_bundle, verify_files
import export as base_export

MARS_REVISION = 'ea1c0ef923246f7489ece4b847a11e43f78ebe4c'
MFEM_REVISION = 'a01719101027383954b69af1777dc828bf795d62'
REFERENCE_SHA = '97c60edb7e9a6cbd82323cbf1a54a0a496717bd7e627a2b190bc5740b9c0a2e8'


def require_source(repo):
    def git(*args):
        return subprocess.check_output(['git', '-C', str(repo)] + list(args), text=True).strip()
    if git('rev-parse', 'HEAD') != MARS_REVISION:
        raise ValueError('Paper source must be at ' + MARS_REVISION)
    if git('status', '--porcelain', '--untracked-files=all'):
        raise ValueError('Paper source checkout must be clean')


def freeze(repo, reference, target):
    require_source(repo)
    if digest(reference / 'manifest.json') != REFERENCE_SHA:
        raise ValueError('Wrong dependency reference manifest')
    previous = verify_bundle(reference)
    if previous['dependencies'] != dict(base_export.REVISIONS, mfem=MFEM_REVISION):
        raise ValueError('Wrong public dependency revisions')

    names = subprocess.check_output(['git', '-C', str(repo), 'ls-files', '-z',
                                     'marsir-mlir/generated'], text=True).split('\0')
    generated = {Path(n).relative_to('marsir-mlir').as_posix(): digest(repo / n) for n in names if n}
    if not generated:
        raise ValueError('Missing tracked generated operators')
    target.parent.mkdir(parents=True, exist_ok=True)
    # The shared exporter takes an operator manifest; supply the pinned checkout,
    # not historical PTX from the dependency reference.
    with tempfile.TemporaryDirectory(prefix='paper-operators-', dir=target.parent) as temporary:
        accepted = Path(temporary)
        for name in generated:
            base_export.copy(repo / 'marsir-mlir' / name, accepted / 'source' / name)
        save(accepted / 'source/source.json', dict(files=generated))
        export_compat(repo, reference / 'dependencies', accepted, target)

    base_export.copy(repo / 'tests/performance/mfem_bp3/Makefile',
                     target / 'source/tests/performance/mfem_bp3/Makefile')
    for name in previous['files']:
        if name.startswith('dependencies/mfem/'):
            base_export.copy(reference / name, target / name)
    for source in HERE.iterdir():
        if source.suffix in ('.py', '.sh', '.md'):
            base_export.copy(source, target / 'paper' / source.name)
    # This campaign rebuilds everything; no historical binary receipt applies.
    (target / 'hcoarse/resume_receipt.json').unlink()

    require_source(repo)
    verify_files(target, generated)
    record = json.loads((target / 'manifest.json').read_text())
    record.update(scope='Final BP3/MFEM comparison and shared CVFEM build; native validation pending',
                  dependency_reference_sha256=REFERENCE_SHA, generated_revision=MARS_REVISION,
                  generated_files=generated, preparation='fresh',
                  harness_base_revision=subprocess.check_output(
                      ['git', '-C', str(HERE), 'rev-parse', 'HEAD'], text=True).strip())
    record['dependencies']['mfem'] = MFEM_REVISION
    record['files'] = {str(p.relative_to(target)): digest(p) for p in target.rglob('*') if p.is_file() and p != target / 'manifest.json'}
    # Runner changes are reviewed separately from the immutable production source.
    record['harness_files'] = {n: sha for n, sha in record['files'].items()
                               if n.startswith(('paper/', 'hcoarse/', 'kit/'))}
    save(target / 'manifest.json', record)
    verify_bundle(target)
    print('Frozen paper kit: %d files, %s' % (len(record['files']), digest(target / 'manifest.json')))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('repo', 'reference', 'target'):
        parser.add_argument('--' + name, type=Path, required=True)
    a = parser.parse_args()
    freeze(*(getattr(a, n).resolve() for n in ('repo', 'reference', 'target')))
