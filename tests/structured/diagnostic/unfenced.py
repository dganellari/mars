"""Observe setup stalls without changing the original binaries or adding GPU fences."""
import argparse
import json
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import diagnose
import performance
import run
import setup_trace

CONTEXTS = (('direct', 'single', 1), ('direct', 'multi', 2), ('direct', 'multi', 4),
            ('staged', 'single', 1), ('staged', 'single', 2), ('staged', 'single', 4),
            ('staged', 'multi', 4))


def schedule():
    tasks = []
    for observer, rounds in (('control', 5), ('mpi-quiet', 2), ('host', 2)):
        for transport, layout, ranks in CONTEXTS:
            for repeat in range(rounds):
                variants = ('reference', 'candidate') if repeat % 2 == 0 else ('candidate', 'reference')
                for variant in variants:
                    tasks.append(dict(transport=transport, layout=layout, ranks=ranks,
                                      observer=observer, repeat=repeat, variant=variant))
    return tasks


def comparisons(records):
    result = []
    for transport, layout, ranks in CONTEXTS:
        group = [r for r in records if r['observer'] == 'control' and
                 (r['transport'], r['layout'], r['ranks']) == (transport, layout, ranks)]
        rounds = {(r['variant'], r['repeat']) for r in group}
        context = dict(transport=transport, layout=layout, ranks=ranks)
        if len(group) != 10 or rounds != {(v, i) for v in ('reference', 'candidate') for i in range(5)} or any(
                r['status'] != 'PASS' for r in group):
            result.append(dict(context, status='INCOMPLETE'))
            continue
        for case in performance.CASES:
            for metric, allowance in performance.METRICS.items():
                values = {v: [r['audit']['cases'][case][metric] for r in group if r['variant'] == v]
                          for v in ('reference', 'candidate')}
                reference, candidate = (statistics.median(values[v]) for v in ('reference', 'candidate'))
                result.append(dict(context, case=case, metric=metric, samples=values,
                                   reference_seconds=reference, candidate_seconds=candidate,
                                   absolute_allowance_seconds=allowance,
                                   passed=candidate <= 1.10 * reference + allowance))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', required=True, type=Path)
    args = parser.parse_args()
    kit = Path(__file__).resolve().parent
    if args.work.exists():
        parser.error('Use a fresh diagnostic directory')
    args.work.mkdir(parents=True)
    result = {'status': 'INCOMPLETE', 'records': [], 'preflights': [],
              'scope': 'Focused diagnostics only. Host callbacks measure host time; overlapping calls are not additive. '
                       'Profiled samples are never release acceptance samples.'}

    def save():
        (args.work / 'result.json').write_text(json.dumps(result, indent=2) + '\n')

    try:
        pins = json.loads((kit / 'pins.json').read_text())
        manifest = json.loads((kit / 'manifest.json').read_text())
        kit_files = {str(kit / p): h for p, h in manifest['files'].items()}
        result['kit_manifest_sha256'] = run.digest(kit / 'manifest.json')
        diagnose.verify_files(kit_files)
        diagnose.verify_files(pins['cuda']['files'])
        result['input_files'] = pins['cuda']['files']
        configs = {}
        for transport in ('direct', 'staged'):
            original = Path(pins['prepared']) / ('cuda-' + transport)
            checked = json.loads((original / 'result.json').read_text())
            if checked['status'] != 'PASS':
                raise ValueError('Original capability tests did not pass')
            templates = {int(t['name'].split('mpi')[-1]) if 'mpi' in t['name'] else 1: t['command']
                         for t in json.loads((original / 'consumer-list.log').read_text())['tests']}
            built = json.loads((original.parent / ('cuda-' + transport + '-performance-build/result.json')).read_text())
            if built['status'] != 'BUILT_NOT_RUN':
                raise ValueError('Original benchmark build did not pass')
            configs[transport] = dict(environment=checked['environment'], templates=templates,
                                      binaries=built['binary_paths'])
        build = args.work / 'observer-build'
        commands = [('configure', ['cmake', '-S', str(kit), '-B', str(build), '-DCMAKE_BUILD_TYPE=Release',
                    '-DCMAKE_CXX_COMPILER=/user-environment/env/default/bin/mpicxx',
                    '-DKOKKOS_TOOLS_INCLUDE=' + str(kit / 'include'),
                    '-DMARS_DIAG_CUDA_INCLUDE=' + pins['cuda_include']]),
                    ('build', ['cmake', '--build', str(build), '--parallel', '8'])]
        result['build_commands'] = commands
        for name, command in commands:
            if diagnose.execute(command, args.work / (name + '.log'), run.execution_environment({}), 600):
                raise ValueError('Observer ' + name + ' failed')
        library = build / 'libmars_structured_trace.so'
        probe = build / 'mars_mpi_trace_probe'
        result['observer_files'] = {str(p): run.digest(p) for p in (library, probe)}

        def execute(task, name, probe=False):
            config = configs[task['transport']]
            prefix = args.work / name
            env = run.execution_environment(config['environment'])
            env.update(MARS_GATE_LAYOUT=task['layout'], OMP_NUM_THREADS='1')
            binary = build / 'mars_mpi_trace_probe' if probe else config['binaries'][task['variant']]
            command = diagnose.command(config['templates'][task['ranks']], kit / 'rank.sh', library,
                                       prefix, binary, task['observer'])
            record = dict(task, name=name, command=command, environment=run.environment_receipt(env), status='FAILED')
            print('Observing ' + name, flush=True)
            record['returncode'] = diagnose.execute(command, prefix.with_suffix('.log'), env)
            try:
                if record['returncode']:
                    raise ValueError('Invocation exited ' + str(record['returncode']))
                if not probe:
                    record['audit'] = performance.audit(prefix.with_suffix('.log').read_text(), task['ranks'],
                        'cuda', task['transport'], 2 if task['layout'] == 'multi' else 1, 1)
                if task['observer'] != 'control':
                    record['observer_audit'] = diagnose.observer(prefix, task['ranks'], task['observer'], probe)
                    if not probe:
                        record['phases'] = setup_trace.audit(prefix, task['ranks'], record['audit'], task['observer'] == 'host')
                record['status'] = 'PASS'
            except Exception as error:
                record['error'] = str(error)
            return record

        for observer, ranks, probe in (('mpi-quiet', 2, True), ('mpi-quiet', 2, False), ('host', 1, False)):
            task = dict(transport='direct', layout='single', ranks=ranks, observer=observer,
                        repeat=0, variant='candidate')
            result['preflights'].append(execute(task, 'preflight-' + observer + ('-error' if probe else ''), probe))
            save()
        # Failed diagnostics remain visible and do not suppress the other observation modes.
        for task in schedule():
            name = '{transport}-{layout}-np{ranks}-{observer}-{repeat}-{variant}'.format(**task)
            result['records'].append(execute(task, name))
            save()
        result['control_comparisons'] = comparisons(result['records'])
        diagnose.verify_files(kit_files)
        diagnose.verify_files(pins['cuda']['files'])
        diagnose.verify_files(result['observer_files'])
        result['status'] = ('DIAGNOSTIC_COLLECTED' if all(r['status'] == 'PASS' for r in
                            result['records'] + result['preflights']) else 'COLLECTED_WITH_FAILURES')
    except Exception as error:
        result['error'] = str(error)
    result['artifacts'] = {str(p): run.digest(p) for p in args.work.glob('*') if p.is_file() and p.name != 'result.json'}
    save()
    print(result['status'] + ': ' + str(args.work / 'result.json'))
    return 0 if result['status'] == 'DIAGNOSTIC_COLLECTED' else 1


if __name__ == '__main__':
    sys.exit(main())
