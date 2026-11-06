"""Compare automatic selection and a larger root against the accepted policy."""
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
import build as integration_build
import run as integration
from common import (CALLBACKS, HYPRE, clean_environment, command, digest, save, step,
                    verify_bundle, verify_files)


def build(bundle, root):
    options = ['-DCMAKE_PROJECT_TOP_LEVEL_INCLUDES=' + str(bundle / 'policy/extra.cmake')]
    integration_build.build(bundle, root, options)
    cache = (root / 'on/CMakeCache.txt').read_text()
    if 'MARS_ENABLE_BP3_COARSE_CONSOLIDATION:STRING=AUTO\n' not in cache:
        raise ValueError('Normal build did not use the fresh AUTO default')
    command(['cmake', '--build', root / 'on', '--parallel', '4', '--target',
             'bp3_policy_solver', 'bp3_policy_fixture'], root / 'policy-compile.log', timeout=2400)
    record = json.loads((root / 'build.json').read_text())
    for target in ('bp3_policy_solver', 'bp3_policy_fixture'):
        binary = root / 'on/policy' / target
        links = command(['ldd', '-r', binary], root / (target + '-linkage.txt'))
        libraries = re.findall(r'libHYPRE[^\s]*\s+=>\s+(\S+)', links)
        if ('not found' in links or 'undefined symbol:' in links or len(libraries) != 1 or
                Path(libraries[0]).resolve() != HYPRE.resolve()):
            raise ValueError('Wrong policy Hypre linkage')
        symbols = command(['nm', '-D', '--defined-only', binary], root / (target + '-symbols.txt'))
        if not all(re.search(r'\bT ' + symbol + r'$', symbols, re.M) for symbol in CALLBACKS):
            raise ValueError('Policy callbacks are missing')
        record['binaries'][str(binary.relative_to(root))] = digest(binary)
    save(root / 'build.json', record)
    verify_bundle(bundle)


def private_binary(prepared):
    return prepared / 'build/on/policy/bp3_policy_solver'


def prepare(bundle, root):
    command(step(1, 160) + ['/usr/bin/python3', HERE / 'policy.py', 'build', bundle, root / 'build'],
            root / 'build.log')
    integration.validate(bundle, root, root / 'build/on')
    # Replace the standard seal only after all extra checks complete.
    (root / 'prepared.json').unlink()
    for ranks in (1, 4):
        for p in (1, 3, 7):
            auto = integration.solve(bundle, root, root, 'auto_p%d_r%d' % (p, ranks), 'on', ranks, p, 8, -1,
                                     automatic=True)
            expected = json.loads((root / ('smoke_p%d_r%d_on' % (p, ranks)) / 'summary.json').read_text())
            if (auto['iterations'], auto['hierarchy']) != (expected['iterations'], expected['hierarchy']):
                raise ValueError('Omitting the flag changed the accepted policy')
    fixture = root / 'build/on/policy/bp3_policy_fixture'
    for ranks, partition, cutoff in ((1, 'balanced', 500000), (4, 'balanced', 50000),
                                      (4, 'balanced', 500000), (4, 'empty-root', 500000),
                                      (4, 'empty-tail', 500000)):
        env = dict(os.environ, BP3_CONSOLIDATE_ROWS=str(cutoff), BP3_FIXTURE_PROFILE='mars')
        text = command(step(ranks, 5) + [fixture, partition, '262143'],
                       root / ('large_r%d_%s_c%d.log' % (ranks, partition, cutoff)), env)
        selected = [int(x) for x in re.findall(r'BP3 consolidation: level \d+, (\d+) rows', text)]
        if ('PASS %s ranks=%d' % (partition, ranks)) not in text:
            raise ValueError('Large fixture failed')
        if ranks > 1 and cutoff > 50000 and (not selected or min(selected) <= 50000):
            raise ValueError('Fixture did not exercise the larger root')
    for ranks in (1, 4):
        for p in (1, 3, 7):
            for cutoff in (-1, 500000):
                result = integration.solve(bundle, root, root, 'private_p%d_r%d_c%d' % (p, ranks, cutoff),
                                           'on', ranks, p, 8, cutoff, executable=private_binary(root))
                if cutoff == -1:
                    expected = json.loads((root / ('auto_p%d_r%d' % (p, ranks)) / 'summary.json').read_text())
                    if (result['iterations'], result['hierarchy']) != (expected['iterations'], expected['hierarchy']):
                        raise ValueError('Private control changed the accepted policy')
    verify_bundle(bundle)
    verify_files(root / 'build', json.loads((root / 'build/build.json').read_text())['binaries'])
    files = {str(p.relative_to(root)): digest(p) for p in root.rglob('*')
             if p.is_file() and (p.suffix in ('.log', '.json') or p.name.endswith('.txt'))
             and 'CMakeFiles' not in p.parts and '_deps' not in p.parts}
    save(root / 'prepared.json', dict(passed=True, manifest_sha256=digest(bundle / 'manifest.json'), evidence=files))
    print('PASS automatic selection and larger-root preparation', flush=True)


def auto_verdict(rows):
    old, new = rows['accepted'], rows['auto']
    ratio = statistics.median([r['wall_ms'] for r in new]) / statistics.median([r['wall_ms'] for r in old])
    equal = len({(r['iterations'], tuple(r['hierarchy'])) for r in old + new}) == 1
    return dict(passed=equal and ratio <= 1.03, median_ratio=ratio, matching_iterations_and_hierarchy=equal)


def root_verdict(rows):
    control, large, auto = rows['private'], rows['large'], rows['auto']
    control_ms, large_ms, auto_ms = ([r['wall_ms'] for r in items] for items in (control, large, auto))
    control_equal = len({(r['iterations'], tuple(r['hierarchy'])) for r in control + auto}) == 1
    control_ratio = statistics.median(control_ms) / statistics.median(auto_ms)
    ratio = statistics.median(large_ms) / statistics.median(control_ms)
    earlier = all(r['hierarchy'][0] < control[0]['hierarchy'][0] and r['hierarchy'][1] > 50000 for r in large)
    unchanged = all(r['hierarchy'][:3] == control[0]['hierarchy'][:3] for r in large)
    timing = (ratio <= 0.95 and max(large_ms) < min(control_ms) and max(large_ms) < min(auto_ms)) if earlier else (unchanged and ratio <= 1.03)
    cold = [r['setup_ms'] + r['wall_solves_ms'][0] for r in large]
    cold_control = [r['setup_ms'] + r['wall_solves_ms'][0] for r in control]
    return dict(control_passed=control_equal and control_ratio <= 1.03, control_ratio=control_ratio,
                control_matching_iterations_and_hierarchy=control_equal, selected_earlier_level=earlier,
                unchanged_level=unchanged, candidate_median_ratio=ratio,
                warm_candidate_passed=timing and control_equal and control_ratio <= 1.03,
                cold_setup_plus_first_ms=cold, cold_control_ms=cold_control,
                promotion_pending_cold_cost_and_memory_review=True)


def root_memory(log):
    matches = re.findall(r'BP3 root memory: (\d+) used of (\d+) bytes after setup', log)
    if len(matches) != 1:
        raise ValueError('Missing root memory snapshot')
    used, total = map(int, matches[0])
    if not 0 < used < total:
        raise ValueError('Invalid root memory snapshot')
    return dict(used_bytes=used, total_bytes=total, peak_measured=False)


def trial(bundle, prepared, root):
    integration.verify_prepared(prepared, bundle)
    cases, failures = [], []
    for family, p, ranks, ne in integration.timing_cases():
        modes = ['accepted', 'auto'] + (['private', 'large'] if ranks == 64 else [])
        rows = {mode: [] for mode in modes}
        for number in range(1, 4):
            order = modes if number % 2 else list(reversed(modes))
            for mode in order:
                name = '%s_p%d_r%d_%s%d' % (family, p, ranks, mode, number)
                cutoff = 500000 if mode == 'large' else -1
                executable = private_binary(prepared) if mode in ('private', 'large') else None
                try:
                    result = integration.solve(bundle, prepared, root, name, 'accepted' if mode == 'accepted' else 'on',
                                               ranks, p, ne, cutoff, executable=executable, automatic=mode == 'auto')
                    if executable:
                        result['root_memory'] = root_memory((root / name / 'solver.log').read_text())
                    rows[mode].append(result)
                except (ValueError, RuntimeError) as error:
                    failures.append(dict(name=name, error=str(error)))
                    print('FAIL %s: %s; continuing independent runs' % (name, error), flush=True)
                save(root / 'samples.json', dict(completed_cases=cases, current_case=[family, p, ranks, ne],
                                                 samples=rows, failures=failures))
        auto = auto_verdict(rows) if len(rows['accepted']) == len(rows['auto']) == 3 else dict(passed=False, incomplete=True)
        case = dict(family=family, p=p, ranks=ranks, ne=ne, samples=rows, auto=auto)
        if ranks == 64:
            case['larger_root'] = root_verdict(rows) if all(len(rows[k]) == 3 for k in ('auto', 'private', 'large')) else dict(
                warm_candidate_passed=False, incomplete=True)
        cases.append(case)
        save(root / 'cases.json', cases)
    integration.verify_prepared(prepared, bundle)
    integration.verify_accepted(bundle)
    verify_bundle(bundle)
    report = dict(auto_passed=all(c['auto']['passed'] for c in cases), cases=cases, failures=failures,
                  larger_root_promoted=False,
                  metric='Third solve, host wall maximum over ranks; three alternating rounds. Cold times retained separately.')
    save(root / 'verdict.json', report)
    print(json.dumps(report, indent=2), flush=True)
    if not report['auto_passed'] or failures:
        raise RuntimeError('One or more gates rejected; all completed measurements retained')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('build', 'prepare', 'trial'))
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
    integration.verify_accepted(bundle)
    if args.phase == 'build':
        build(bundle, root)
        return
    root.mkdir(parents=True, exist_ok=False)
    save(root / 'inputs.json', dict(manifest_sha256=digest(bundle / 'manifest.json'), phase=args.phase))
    if args.phase == 'prepare':
        prepare(bundle, root)
    else:
        if not args.prepared:
            raise ValueError('Trial requires passing preparation')
        trial(bundle, args.prepared.resolve(), root)


if __name__ == '__main__':
    main()
