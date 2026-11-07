"""Freeze current MARS operators and the immutable MFEM 4.8 source dependency."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'bp3_hcoarse'))
from export_compat import export_compat
from common import digest, save, verify_bundle
import export as base_export

MFEM_REVISION = 'a01719101027383954b69af1777dc828bf795d62'


def freeze(repo, deps, accepted, mfem, target, resume_receipt=None):
    revision = subprocess.check_output(['git', '-C', str(mfem), 'rev-parse', 'v4.8^{commit}'], universal_newlines=True).strip()
    if revision != MFEM_REVISION:
        raise ValueError('Wrong MFEM release')
    export_compat(repo, deps, accepted, target)
    # The ordinary source filter omits Makefiles; MFEM's installed config.mk drives this one.
    base_export.copy(repo / 'tests/performance/mfem_bp3/Makefile', target / 'source/tests/performance/mfem_bp3/Makefile')
    library = target / 'dependencies/mfem'
    library.mkdir()
    archive = target / 'mfem.tar'
    with archive.open('wb') as stream:
        subprocess.run(['git', '-C', str(mfem), 'archive', MFEM_REVISION], stdout=stream, check=True)
    import tarfile
    with tarfile.open(archive) as tar:
        for item in tar.getmembers():
            path = Path(item.name)
            if path.is_absolute() or '..' in path.parts or not (item.isfile() or item.isdir()):
                raise ValueError('Unsafe dependency archive member')
        tar.extractall(library, filter='data')
    archive.unlink()
    for source in HERE.iterdir():
        if source.suffix in ('.py', '.sh', '.md'):
            base_export.copy(source, target / 'paper' / source.name)
    if resume_receipt is not None:
        base_export.copy(resume_receipt, target / 'paper/resume-receipt.json')
    record = json.loads((target / 'manifest.json').read_text())
    record['scope'] = 'Current-operator BP3/MFEM paper comparison; 256 ranks are a pilot, not validated scaling'
    record['dependencies']['mfem'] = MFEM_REVISION
    record['files'] = {str(p.relative_to(target)): digest(p) for p in target.rglob('*') if p.is_file() and p != target / 'manifest.json'}
    save(target / 'manifest.json', record)
    verify_bundle(target)
    print('Frozen paper kit: %d files, %s' % (len(record['files']), digest(target / 'manifest.json')))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('repo', 'deps', 'accepted', 'mfem', 'target'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--resume-receipt', type=Path)
    a = parser.parse_args()
    freeze(*(getattr(a, n).resolve() for n in ('repo', 'deps', 'accepted', 'mfem', 'target')),
           resume_receipt=a.resume_receipt.resolve() if a.resume_receipt else None)
