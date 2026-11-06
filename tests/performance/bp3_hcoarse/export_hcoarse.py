"""Freeze normal-build integration with the already validated image, PTX and reference."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'bp3_coarse/integration'))
import export
from common import digest, save, verify_bundle


def freeze(repo, deps, accepted, target):
    export.NEW += tuple('backend/distributed/unstructured/solvers/' + n for n in
                        ('mars_cellwise_h_transfer.hpp', 'mars_cellwise_h_coarse.hpp'))
    export.NEW += ('tests/bp3_coarse/integration/cli_stubs.hpp',)
    export.export(repo, deps, accepted, target)
    for name in ('job.sh', 'submit.sh'):
        (target / 'kit' / name).unlink()
    for source in HERE.iterdir():
        if source.suffix in ('.hpp', '.py', '.sh', '.cmake', '.md', '.json'):
            export.copy(source, target / 'hcoarse' / source.name)
    export.copy(repo / 'examples/distributed/unstructured/mars_marsir_bp.cu', target / 'hcoarse/private/mars_marsir_bp.cu')
    export.copy(HERE / 'check.hpp', target / 'hcoarse/private/check.hpp')
    record = json.loads((target / 'manifest.json').read_text())
    record['scope'] = 'Normal BP3 h hierarchy integration with scoped automatic defaults; unchanged numerical candidate'
    record['files'] = {str(p.relative_to(target)): digest(p) for p in target.rglob('*')
                       if p.is_file() and p != target / 'manifest.json'}
    save(target / 'manifest.json', record)
    verify_bundle(target)
    print('Frozen integration kit: %d files, %s' % (len(record['files']), digest(target / 'manifest.json')))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('repo', 'deps', 'accepted', 'target'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    freeze(args.repo.resolve(), args.deps.resolve(), args.accepted.resolve(), args.target.resolve())
