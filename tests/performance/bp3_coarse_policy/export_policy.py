"""Freeze the automatic-default build and its isolated larger-root comparison."""
import argparse
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'bp3_coarse/integration'))
import export
from common import digest, save, verify_bundle
from candidate import generate


def freeze(repo, deps, accepted, target):
    export.export(repo, deps, accepted, target)
    # Old integration submitters are not valid entry points for this experiment.
    for name in ('job.sh', 'submit.sh'):
        (target / 'kit' / name).unlink()
    for source in HERE.iterdir():
        if source.suffix in ('.py', '.sh', '.cmake', '.md'):
            export.copy(source, target / 'policy' / source.name)
    generate(target / 'source', target / 'policy/private')
    record = verify_source_manifest(target)
    save(target / 'manifest.json', record)
    verify_bundle(target)
    print('Frozen policy bundle: %d files, %s' % (len(record['files']), digest(target / 'manifest.json')))


def verify_source_manifest(target):
    import json
    record = json.loads((target / 'manifest.json').read_text())
    record['scope'] = 'Automatic BP3 default and isolated fixed 500000-row root; adaptive policy unchanged'
    record['files'] = {str(p.relative_to(target)): digest(p) for p in target.rglob('*')
                       if p.is_file() and p != target / 'manifest.json'}
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('repo', 'deps', 'accepted', 'target'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    freeze(args.repo.resolve(), args.deps.resolve(), args.accepted.resolve(), args.target.resolve())
