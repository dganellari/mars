"""Check normal-build selection and equivalence to the validated h hierarchy."""
import argparse
import json
import os
from pathlib import Path
import re
import statistics
import sys

HERE = Path(__file__).resolve().parent
KIT = HERE.parent / 'kit' if (HERE.parent / 'kit').is_dir() else HERE.parents[1] / 'bp3_coarse/integration'
sys.path.insert(0, str(KIT))
import build as standard_build
import run as standard_run
from common import (CALLBACKS, clean_environment, command, digest, parse, save, solver_args,
                    step, verify_bundle, verify_files)
from receipts import parse_h

REFERENCE = Path('/capstor/scratch/cscs/gandanie/bp3-hcoarse/5278195/build/native/hcoarse/bp3_h_solver')
REFERENCE_SHA = '557f2b0b6c3667912aebec4ac906ab4956f7712adbc9750de4605fea0b75fc6a'


def verify_reference(bundle):
    standard_run.verify_accepted(bundle)
    if digest(REFERENCE) != REFERENCE_SHA:
        raise ValueError('Validated h executable changed')


def binary(prepared, mode):
    if mode == 'reference':
        return REFERENCE
    if mode == 'check':
        return prepared / 'build/on/hcoarse/bp3_h_solver'
    return prepared / 'build' / mode / 'examples/distributed/unstructured/mars_marsir_bp'


def build(bundle, root):
    standard_build.build(bundle, root, ['-DCMAKE_PROJECT_TOP_LEVEL_INCLUDES=' + str(bundle / 'hcoarse/extra.cmake')])
    command(['cmake', '--build', root / 'on', '--parallel', '4', '--target', 'bp3_h_solver'],
            root / 'check-compile.log', timeout=1800)
    executable = binary(root.parent, 'check')
    link = command(['ldd', '-r', executable], root / 'check-linkage.txt')
    if 'not found' in link or 'undefined symbol:' in link:
        raise ValueError('Check executable has unresolved dependencies')
    symbols = command(['nm', '-D', '--defined-only', executable], root / 'check-symbols.txt')
    if not all(re.search(r'\bT ' + s + r'$', symbols, re.M) for s in CALLBACKS):
        raise ValueError('Check target lost the adapter')
    db = json.loads((root / 'on/compile_commands.json').read_text())
    for entry in db:
        flags = entry.get('command', ' '.join(entry.get('arguments', [])))
        expected = Path(entry['file']) == bundle / 'hcoarse/private/mars_marsir_bp.cu'
        if ('MARS_BP3_HCOARSE_CHECK' in flags) != expected:
            raise ValueError('Validation instrumentation escaped its target')
    record = json.loads((root / 'build.json').read_text())
    record['binaries'][str(executable.relative_to(root))] = digest(executable)
    save(root / 'build.json', record)


def solve(bundle, prepared, root, name, mode, ranks, p, ne, depth=0, automatic=False, checked=False, solves=3):
    out = root / name
    reused = out.exists()
    if reused:
        from resume import verify_reused
        verify_reused(root, name)
    else:
        out.mkdir()
    executable = binary(prepared, mode)
    cutoff = 0 if mode == 'off' else -1
    flags = solver_args(bundle / 'generated', p, ne, cutoff)
    flags[flags.index('--solves') + 1] = str(solves)
    if automatic:
        i = flags.index('--amg-consolidate')
        del flags[i:i + 2]
    else:
        flags += ['--h-levels', str(depth)]
    if checked:
        flags.append('--h-check')
    if not reused:
        command(step(ranks) + ['/usr/bin/python3', KIT / 'rank.py', out,
                'no' if mode == 'off' else 'yes', executable] + flags, out / 'solver.log')
    result = parse_h(out, ranks, ne, depth, checked, solves) if depth else parse(out, ranks, cutoff, mode != 'off', solves)
    text = (out / 'solver.log').read_text()
    if not depth and 'BP3 geometric levels:' in text:
        raise ValueError('Fallback unexpectedly constructed a hierarchy')
    if mode != 'reference':
        expected = [('auto' if automatic else 'explicit', str(depth))]
        if re.findall(r'BP3 h selection: requested=(\w+) levels=(\d+)', text) != expected:
            raise ValueError('Wrong automatic or explicit selection')
    if mode in ('on', 'off') and ('BP3 h memory after setup:' in text or 'PASS h validation:' in text):
        raise ValueError('Production binary ran test instrumentation')
    expected_sha = digest(executable)
    for rank in range(ranks):
        receipt = json.loads((out / ('process%d.json' % rank)).read_text())
        if receipt['rank'] != rank or receipt['returncode'] or receipt['binary_sha256'] != expected_sha:
            raise ValueError('Wrong process receipt')
        if mode != 'off' and receipt['audit_sha256'] != digest(out / ('rank%d.json' % rank)):
            raise ValueError('Wrong rank audit')
    print(('PASS reused ' if reused else 'PASS ') + name, flush=True)
    return result


def equivalent(rows):
    return len({(r['iterations'], tuple(r['hierarchy'])) for r in rows}) == 1


def checked_equivalent(normal, checked):
    # validate() applies the hierarchy four times, each with three bottom AMG cycles.
    # The adapter counts none of those calls when consolidation is inactive.
    expected_extra = 12 if normal['hierarchy'][0] >= 0 else 0
    return (checked.get('h_check_worst') is not None and
            checked['iterations'] == normal['iterations'] and
            checked['hierarchy'][:3] == normal['hierarchy'][:3] and
            checked['hierarchy'][3] == normal['hierarchy'][3] + expected_extra)


def prepare(bundle, root, previous=None):
    if previous is not None:
        from resume import restore
        restore(bundle, previous, root)
    else:
        command(step(1, 170) + ['/usr/bin/python3', HERE / 'integrate.py', 'build', bundle, root / 'build'],
                root / 'build.log')
        inventory = command(['ctest', '--test-dir', root / 'build/on', '-R', '^marsBP3CoarseSingle$', '--show-only=json-v1'],
                            root / 'ctest-inventory.json', timeout=30)
        if [t['name'] for t in json.loads(inventory)['tests']] != ['marsBP3CoarseSingle']:
            raise ValueError('Missing adapter CTest')
        text = command(['ctest', '--test-dir', root / 'build/on', '-R', '^marsBP3CoarseSingle$', '--no-tests=error',
                        '--verbose', '--output-on-failure', '--timeout', '180'], root / 'ctest.log', timeout=240)
        if '100% tests passed, 0 tests failed out of 1' not in text:
            raise ValueError('Adapter CTest failed')
    for ranks in (1, 4):
        for p in (1, 3, 7):
            solve(bundle, root, root, 'fallback_off_p%d_r%d' % (p, ranks), 'off', ranks, p, 8, automatic=True)
            normal = solve(bundle, root, root, 'fallback_on_p%d_r%d' % (p, ranks), 'on', ranks, p, 8, automatic=True)
            previous = solve(bundle, root, root, 'fallback_ref_p%d_r%d' % (p, ranks), 'reference', ranks, p, 8)
            if not equivalent([normal, previous]):
                raise ValueError('Fallback changed the solver')
    # Exercise every public explicit depth on one and several GPUs, with direct algebra checks.
    for ranks in (1, 4):
        for depth in (1, 2, 3):
            rows = [solve(bundle, root, root, 'explicit_r%d_h%d_%s' % (ranks, depth, mode), mode,
                          ranks, 7, 32, depth, checked=mode == 'check') for mode in ('on', 'reference', 'check')]
            if not equivalent(rows[:2]) or not checked_equivalent(rows[0], rows[2]):
                raise ValueError('Integrated hierarchy differs from the validated one')
    for suffix, flags, needle in (
        ('range', ['--h-levels', '4'], 'expects -1'),
        ('negative', ['--h-levels', '-2'], 'expects -1'),
        ('p1', ['--h-levels', '1', '--p', '1'], 'require unblocked'),
        ('rotated', ['--h-levels', '1', '--rotate'], 'require unblocked'),
        ('alignment', ['--h-levels', '3', '--ne', '12'], 'need aligned'),
        ('check', ['--h-check'], 'an option needs a value')):
        text = command(step(1) + [binary(root, 'on')] + solver_args(bundle / 'generated', 3, 32, -1) + flags,
                       root / ('reject_' + suffix + '.log'), expected=1)
        if needle not in text:
            raise ValueError('Invalid command not rejected as expected')
    verify_bundle(bundle)
    verify_reference(bundle)
    verify_files(root / 'build', json.loads((root / 'build/build.json').read_text())['binaries'])
    files = {str(p.relative_to(root)): digest(p) for p in root.rglob('*') if p.is_file() and
             (p.suffix in ('.log', '.json') or p.name.endswith('.txt')) and
             'CMakeFiles' not in p.parts and '_deps' not in p.parts}
    save(root / 'prepared.json', dict(passed=True, manifest_sha256=digest(bundle / 'manifest.json'), evidence=files))
    print('PASS preparation: six binaries, 36 solver invocations and six rejected commands', flush=True)


def cases():
    # One fallback per class and every automatic winner; no repeat search over depths.
    return [('single_p3', 1, 3, 64, 0, 3), ('single_p7', 1, 7, 32, 0, 3),
            ('weak_p3_8', 8, 3, 128, 2, 3), ('short_p7_8', 8, 7, 64, 0, 3),
            ('repeated_p7_8', 8, 7, 64, 1, 5), ('weak_p3_64', 64, 3, 256, 2, 3),
            ('weak_p7_64', 64, 7, 128, 1, 3), ('strong_p3_64', 64, 3, 128, 2, 3),
            ('small_p7_64', 64, 7, 64, 0, 3)]


def verdict(samples):
    if any(len(samples[m]) != 3 for m in ('on', 'reference')):
        return dict(passed=False, incomplete=True)
    ratio = statistics.median(r['wall_ms'] for r in samples['on']) / statistics.median(r['wall_ms'] for r in samples['reference'])
    stable = equivalent(samples['on'] + samples['reference'])
    return dict(passed=stable and abs(ratio - 1) <= .03, equivalent=stable, warm_median_ratio=ratio,
                cold_ms={m: [r['setup_ms'] + r['wall_solves_ms'][0] for r in rs] for m, rs in samples.items()})


def trial(bundle, prepared, root):
    standard_run.verify_prepared(prepared, bundle)
    completed, failures = [], []
    for name, ranks, p, ne, depth, solves in cases():
        samples = dict(on=[], reference=[])
        for round_ in range(3):
            for mode in (('reference', 'on') if round_ % 2 == 0 else ('on', 'reference')):
                label = '%s_%s_%d' % (name, mode, round_ + 1)
                try:
                    samples[mode].append(solve(bundle, prepared, root, label, mode, ranks, p, ne,
                                               depth, automatic=mode == 'on', solves=solves))
                except (ValueError, RuntimeError) as error:
                    failures.append(dict(name=label, error=str(error)))
                    print('FAIL %s: %s; retaining output and continuing independent checks' % (label, error), flush=True)
                save(root / 'samples.json', dict(completed=completed, current=name, samples=samples, failures=failures))
        completed.append(dict(name=name, ranks=ranks, p=p, ne=ne, depth=depth, solves=solves,
                              samples=samples, verdict=verdict(samples)))
        save(root / 'cases.json', completed)
    standard_run.verify_prepared(prepared, bundle)
    verify_reference(bundle)
    verify_bundle(bundle)
    report = dict(cases=completed, failures=failures, passed=not failures and all(c['verdict']['passed'] for c in completed))
    save(root / 'verdict.json', report)
    if not report['passed']:
        raise RuntimeError('Integration gate rejected; evidence retained')
    print('PASS normal-build integration: all automatic choices and fallbacks match the validated reference', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('build', 'prepare', 'resume', 'trial'))
    parser.add_argument('bundle', type=Path)
    parser.add_argument('results', type=Path)
    parser.add_argument('--prepared', type=Path)
    args = parser.parse_args()
    clean_environment()
    needed = 16 if args.phase == 'trial' else 1
    if int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) < needed:
        raise ValueError('Needs %d allocated nodes' % needed)
    bundle, root = args.bundle.resolve(), args.results.resolve()
    verify_bundle(bundle)
    verify_reference(bundle)
    if args.phase == 'build':
        build(bundle, root)
        return
    root.mkdir(parents=True, exist_ok=False)
    save(root / 'inputs.json', dict(manifest_sha256=digest(bundle / 'manifest.json'), phase=args.phase, reference_sha256=REFERENCE_SHA))
    if args.phase == 'prepare':
        prepare(bundle, root)
    elif args.phase == 'resume':
        if not args.prepared:
            raise ValueError('Resume requires the reviewed previous run')
        prepare(bundle, root, args.prepared.resolve())
    elif args.prepared:
        trial(bundle, args.prepared.resolve(), root)
    else:
        raise ValueError('Trial requires sealed preparation')


if __name__ == '__main__':
    main()
