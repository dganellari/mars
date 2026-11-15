"""Summarize only complete three-round comparisons; preserve failed and missing cases."""
import argparse
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
KIT = HERE.parent / 'kit' if (HERE.parent / 'kit').is_dir() else HERE.parents[1] / 'bp3_coarse/integration'
sys.path.insert(0, str(KIT))
from common import verify_files


def summarize(paths):
    records = []
    manifest = None
    failures = []
    for path in paths:
        record = json.loads((path / 'results.json').read_text())
        evidence = json.loads((path / 'evidence.json').read_text())
        verify_files(path, evidence['evidence'])
        sha = record['manifest_sha256']
        if evidence['manifest_sha256'] != sha or manifest is not None and manifest != sha:
            raise ValueError('Do not combine different frozen revisions')
        manifest = sha
        records += record['cases']
        failures += record['failures']
    keys = [(c['workload'], c['ranks'], c['p']) for c in records]
    if len(keys) != len(set(keys)):
        raise ValueError('Duplicate cases need an explicit new-round analysis')
    baseline = {(c['workload'], c['p'], v): row['median_ms'] for c in records if c['ranks'] == 1
                for v, row in c['summary'].items()}
    lines = ['# BP3 comparison', '', 'Manifest: `' + str(manifest) + '`.', '',
             'Times are median fifth-solve host wall times across three rounds. Efficiency uses',
             'the same configuration on one GPU; it includes iteration changes. It is not',
             'a comparison against the fastest one-GPU configuration.', '',
             '| Workload | p | GPUs | Configuration | Solve ms [min, max] | Iterations | Efficiency |',
             '|---|---:|---:|---|---|---|---|']
    ratios = []
    for c in records:
        for v, row in c['summary'].items():
            base = baseline.get((c['workload'], c['p'], v))
            efficiency = base / row['median_ms'] / (c['ranks'] if c['workload'] == 'strong' else 1) if base else None
            lines.append('| %s | %d | %d | %s | %.3f [%.3f, %.3f] | %s | %s |' % (
                c['workload'], c['p'], c['ranks'], v, row['median_ms'], row['min_ms'], row['max_ms'],
                ', '.join(str(x) for x in row['iterations']), '%.1f%%' % (100 * efficiency) if efficiency is not None else 'unavailable'))
        if 'mfem' in c['summary']:
            mfem = c['summary']['mfem']['median_ms']
            for v in ('matched', 'best_h0', 'best_h'):
                if v in c['summary']:
                    ratios.append('%s p%d, %d GPUs: MFEM / %s = %.3fx.' % (
                        c['workload'], c['p'], c['ranks'], v, mfem / c['summary'][v]['median_ms']))
    lines += [''] + ratios + ['', 'Failures retained: %d.' % len(failures)]
    for f in failures:
        lines.append('- %s: %s' % (f['name'], f['error']))
    lines += ['', 'Inspect setup/first-solve costs, per-iteration costs and iteration counts in results.json.',
              'Missing configurations have no accepted timing or inferred speedup.']
    return '\n'.join(lines) + '\n'


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('results', type=Path, nargs='+')
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    a.output.write_text(summarize(a.results))
