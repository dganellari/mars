"""Freeze the current matching operator loader and PTX for the combined-source check."""
import argparse
import json
from pathlib import Path
from export_hcoarse import freeze
from common import digest, save, verify_bundle


def export_compat(repo, deps, accepted, target):
    freeze(repo, deps, accepted, target, current_operators=True)
    replacements = {}
    for spec in ('bp1', 'bp3', 'bp1_tri', 'bp3_tri'):
        for p in (1, 3, 7):
            name = '%s_p%d_sm90.ptx' % (spec, p)
            source = repo / 'marsir-mlir/generated' / name
            text = source.read_text()
            if '%s_apply_p%d__m0' % (spec, p) not in text:
                raise ValueError('Expected the constant-memory operator ABI: ' + name)
            (target / 'generated' / name).write_bytes(source.read_bytes())
            replacements[name] = digest(source)
    record = json.loads((target / 'manifest.json').read_text())
    record['scope'] = 'Combined-source compatibility after bf087b34; current operator loader and its matching PTX'
    record['operator_ptx'] = replacements
    record['files'] = {str(p.relative_to(target)): digest(p) for p in target.rglob('*') if p.is_file() and p != target / 'manifest.json'}
    save(target / 'manifest.json', record)
    verify_bundle(target)
    print('Frozen compatibility kit: %d files, %s' % (len(record['files']), digest(target / 'manifest.json')))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for n in ('repo', 'deps', 'accepted', 'target'):
        parser.add_argument('--' + n, type=Path, required=True)
    a = parser.parse_args()
    export_compat(a.repo.resolve(), a.deps.resolve(), a.accepted.resolve(), a.target.resolve())
