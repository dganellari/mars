"""Validate the standard build, then compare it with the accepted executable."""
import argparse
import json
import os
from pathlib import Path
import statistics
from common import (ACCEPTED, ACCEPTED_SHA, clean_environment, command, digest, parse,
                    save, solver_args, step, verify_bundle, verify_files)

HERE = Path(__file__).resolve().parent


def timing_cases():
    cases = [('weak', p, ranks, {1: 64, 8: 128, 64: 256}[ranks] // (2 if p == 7 else 1))
             for ranks in (1, 8, 64) for p in (3, 7)]
    return cases + [('strong', 3, 64, 128), ('strong', 7, 64, 64)]


def binary(prepared, mode):
    if mode == 'accepted':
        return ACCEPTED / 'build/bp3_solver'
    return prepared / 'build' / mode / 'examples/distributed/unstructured/mars_marsir_bp'


def verify_accepted(bundle):
    if digest(ACCEPTED / 'build/bp3_solver') != ACCEPTED_SHA:
        raise ValueError('Accepted control executable changed')
    manifest = json.loads((ACCEPTED / 'source/source.json').read_text())
    verify_files(ACCEPTED / 'source', manifest['files'])
    exported = json.loads((bundle / 'manifest.json').read_text())
    ptx = {k[10:]: v for k, v in manifest['files'].items() if k.startswith('generated/')}
    if ptx != {k[10:]: v for k, v in exported['files'].items() if k.startswith('generated/')}:
        raise ValueError('Generated kernels differ from the accepted control')


def solve(bundle, prepared, root, name, mode, ranks, p, ne, cutoff):
    out = root / name
    out.mkdir()
    args = step(ranks) + ['/usr/bin/python3', str(HERE / 'rank.py'), str(out),
                         'no' if mode == 'off' else 'yes', str(binary(prepared, mode))]
    args += solver_args(bundle / 'generated', p, ne, cutoff)
    command(args, out / 'solver.log')
    result = parse(out, ranks, cutoff, mode != 'off')
    expected_sha = digest(binary(prepared, mode))
    for rank in range(ranks):
        receipt = json.loads((out / ('process%d.json' % rank)).read_text())
        if receipt['rank'] != rank or receipt['returncode'] != 0 or receipt['binary_sha256'] != expected_sha:
            raise ValueError('Wrong per-rank execution receipt')
        if mode != 'off' and receipt['audit_sha256'] != digest(out / ('rank%d.json' % rank)):
            raise ValueError('Audit export changed')
    print('PASS ' + name, flush=True)
    return result


def verify_prepared(prepared, bundle):
    record = json.loads((prepared / 'prepared.json').read_text())
    if record['passed'] is not True or record['manifest_sha256'] != digest(bundle / 'manifest.json'):
        raise ValueError('Preparation or source manifest changed')
    verify_files(prepared, record['evidence'])
    build = json.loads((prepared / 'build/build.json').read_text())
    verify_files(prepared / 'build', build['binaries'])
    return record


def prepare(bundle, root):
    command(step(1, 115) + ['/usr/bin/python3', HERE / 'build.py', bundle, root / 'build'], root / 'build.log')
    validate(bundle, root, root / 'build/on')


def validate(bundle, root, test_build):
    fixture = root / 'build/on/examples/distributed/unstructured/mars_bp3_coarse_test'
    # CTest coordinates in the batch process; its test launches a fresh MPI step.
    inventory = command(['ctest', '--test-dir', test_build, '-R', '^marsBP3CoarseSingle$',
                         '--show-only=json-v1'], root / 'ctest-inventory.json', timeout=30)
    if [test['name'] for test in json.loads(inventory)['tests']] != ['marsBP3CoarseSingle']:
        raise ValueError('Expected integration CTest is missing')
    ctest = command(['ctest', '--test-dir', test_build, '-R', '^marsBP3CoarseSingle$',
                     '--no-tests=error', '--verbose', '--output-on-failure', '--timeout', '180'],
                    root / 'ctest.log', timeout=240)
    if '100% tests passed, 0 tests failed out of 1' not in ctest:
        raise ValueError('The integration CTest did not pass')
    for ranks, partition, n in ((1, 'balanced', 511), (4, 'balanced', 511),
                                (4, 'empty-root', 511), (4, 'empty-tail', 511), (4, 'balanced', 32767)):
        for cutoff in (0, 100 if n == 511 else 5000, 50000):
            log = root / 'fixture_r{}_{}_n{}_c{}.log'.format(ranks, partition, n, cutoff)
            env = dict(os.environ, BP3_CONSOLIDATE_ROWS=str(cutoff), BP3_FIXTURE_PROFILE='mars')
            text = command(step(ranks) + [fixture, partition, n], log, env)
            if 'PASS {} ranks={}'.format(partition, ranks) not in text:
                raise ValueError('Fixture failed')
    # The disabled build and the enabled-but-unused adapter must match stock.
    for ranks in (1, 4):
        for p in (1, 3, 7):
            records = {}
            for label, mode, cutoff in (('off', 'off', 0), ('on0', 'on', 0), ('stock', 'accepted', 0),
                                         ('on', 'on', -1), ('accepted', 'accepted', -1)):
                name = 'smoke_p{}_r{}_{}'.format(p, ranks, label)
                records[label] = solve(bundle, root, root, name, mode, ranks, p, 8, cutoff)
            if (len({records[k]['iterations'] for k in ('off', 'on0', 'stock')}) != 1 or
                    records['on']['iterations'] != records['accepted']['iterations']):
                raise ValueError('Smoke iteration counts changed')
            if records['on']['hierarchy'] != records['accepted']['hierarchy'] or records['on0']['hierarchy'] != records['stock']['hierarchy']:
                raise ValueError('Smoke coarse hierarchy changed')
    for mode, cutoff, needle in (('off', -1, 'MARS_ENABLE_BP3_COARSE_CONSOLIDATION=ON'),
                                  ('on', -2, 'invalid --amg-consolidate'), ('on', 9, 'expects -1')):
        log = root / ('reject_{}_{}.log'.format(mode, cutoff))
        text = command(step(1) + [binary(root, mode)] + solver_args(bundle / 'generated', 3, 8, cutoff), log, expected=1)
        if needle not in text:
            raise ValueError('CLI did not reject the unsupported mode')
    build = json.loads((root / 'build/build.json').read_text())
    verify_files(root / 'build', build['binaries'])
    verify_bundle(bundle)
    files = {str(p.relative_to(root)): digest(p) for p in root.rglob('*')
             if p.is_file() and (p.suffix in ('.log', '.json') or p.name.endswith('.txt'))
             and 'CMakeFiles' not in p.parts and '_deps' not in p.parts}
    # The outer build.log is complete before this receipt is written.
    save(root / 'prepared.json', dict(passed=True, manifest_sha256=digest(bundle / 'manifest.json'), evidence=files))
    print('PASS preparation; dependent integration timing may run', flush=True)


def trial(bundle, prepared, root):
    verify_prepared(prepared, bundle)
    rows = []
    for family, p, ranks, ne in timing_cases():
        samples = {mode: [] for mode in ('accepted', 'on')}
        for number, order in enumerate((('accepted', 'on'), ('on', 'accepted'), ('accepted', 'on')), 1):
            for mode in order:
                name = '{}_p{}_r{}_{}{}'.format(family, p, ranks, mode, number)
                samples[mode].append(solve(bundle, prepared, root, name, mode, ranks, p, ne, -1))
        old = [r['wall_ms'] for r in samples['accepted']]
        new = [r['wall_ms'] for r in samples['on']]
        ratio = statistics.median(new) / statistics.median(old)
        equal = len({(r['iterations'], tuple(r['hierarchy'])) for records in samples.values() for r in records}) == 1
        rows.append(dict(family=family, p=p, ranks=ranks, ne=ne, accepted_ms=old, integrated_ms=new,
                         median_ratio=ratio, matching_iterations_and_hierarchy=equal,
                         passed=equal and ratio <= 1.03))
    verify_prepared(prepared, bundle)
    verify_bundle(bundle)
    verify_accepted(bundle)
    report = dict(passed=all(r['passed'] for r in rows), cases=rows,
                  metric='Last of three solves, host wall maximum over ranks; three alternating rounds',
                  native_integration_only=True, default_changed=False)
    save(root / 'verdict.json', report)
    print(json.dumps(report, indent=2), flush=True)
    if not report['passed']:
        raise RuntimeError('Integrated-build equivalence gate failed; retain all measurements')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('prepare', 'resume', 'trial'))
    parser.add_argument('bundle', type=Path)
    parser.add_argument('results', type=Path)
    parser.add_argument('--prepared', type=Path)
    args = parser.parse_args()
    clean_environment()
    needed = 16 if args.phase == 'trial' else 1
    if int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) < needed:
        raise ValueError('Needs {} allocated nodes'.format(needed))
    bundle, root = args.bundle.resolve(), args.results.resolve()
    verify_bundle(bundle)
    verify_accepted(bundle)
    root.mkdir(parents=True, exist_ok=False)
    save(root / 'inputs.json', dict(manifest_sha256=digest(bundle / 'manifest.json'),
                                  accepted_sha256=ACCEPTED_SHA, phase=args.phase))
    if args.phase == 'prepare':
        prepare(bundle, root)
    elif args.phase == 'resume':
        from reuse import resume
        if not args.prepared:
            raise ValueError('Resume requires a verified native build')
        test_build = resume(bundle, args.prepared.resolve(), root)
        validate(bundle, root, test_build)
    else:
        if not args.prepared:
            raise ValueError('Trial requires passing preparation')
        trial(bundle, args.prepared.resolve(), root)


if __name__ == '__main__':
    main()
