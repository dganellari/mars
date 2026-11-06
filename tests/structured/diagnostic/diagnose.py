#!/usr/bin/env python3
"""Inspect the remaining CPU/CUDA failures using the original benchmark binaries."""
import argparse
import collections
import json
import os
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
import run
import performance


def schedule(backend):
    tasks = []
    configurations = ('openmp-direct',) if backend == 'cpu' else ('cuda-direct', 'cuda-staged')
    for config in configurations:
        if config == 'cuda-direct':
            for ranks in (2, 4):
                for variant in ('reference', 'candidate'):
                    tasks.append(dict(config=config, mode='mpi', layout='multi', ranks=ranks, repeat=0, variant=variant))
        layouts = ('single', 'multi') if config == 'cuda-staged' else ('single',)
        for layout in layouts:
            for repeat in range(5):
                order = ('reference', 'candidate') if repeat % 2 == 0 else ('candidate', 'reference')
                for variant in order:
                    tasks.append(dict(config=config, mode='control', layout=layout, ranks=4, repeat=repeat, variant=variant))
        for layout in (('single',) if backend == 'cpu' else ('single', 'multi')):
            for variant in ('reference', 'candidate'):
                tasks.append(dict(config=config, mode='profile', layout=layout, ranks=4, repeat=0, variant=variant))
    return tasks


def verify_files(files):
    for filename, expected in files.items():
        path = Path(filename)
        if run.digest(path) != expected:
            raise ValueError('Input changed: ' + filename)


def observer(prefix, ranks, mode, probe=False):
    records = []
    for rank in range(ranks):
        data = json.loads(Path(str(prefix) + '.rank' + str(rank) + '.json').read_text())
        if data['rank'] != rank or data['status'] != 'complete' or data['dropped']:
            raise ValueError('Incomplete observer receipt')
        if not any(e['kind'] == 'mpi' and e['name'] == 'MPI_Allreduce' for e in data['events']):
            raise ValueError('MPI observer was not active')
        if probe and not any(e['kind'] == 'mpi_error' and e['name'].startswith('MPI_Isend:') for e in data['events']):
            raise ValueError('MPI returned-error probe was not observed')
        if mode in ('profile', 'host') and not (data['tools_loaded'] and data['tools_finalized']):
            raise ValueError('Kokkos profiler was not active')
        if mode in ('profile', 'host') and data['global_fencing'] != (mode == 'profile'):
            raise ValueError('Unexpected profiler synchronization')
        if mode in ('profile', 'host') and not any(e['kind'] == 'parallel_for' for e in data['events']):
            raise ValueError('No Kokkos kernels were observed')
        if mode == 'mpi-quiet' and (data['tools_loaded'] or data['global_fencing']):
            raise ValueError('MPI-only trace enabled Kokkos profiling')
        totals = collections.defaultdict(lambda: {'calls': 0, 'seconds': 0.0})
        for event in data['events']:
            key = event['kind'] + ':' + event['name']
            totals[key]['calls'] += 1
            totals[key]['seconds'] += event['seconds']
        records.append({'rank': rank, 'hostname': data['hostname'], 'event_count': len(data['events']),
                        'totals': dict(totals)})
    return records


def command(template, wrapper, library, prefix, binary, mode):
    indices = [i for i, word in enumerate(template) if Path(word).name == 'mars_structured_consumer']
    if len(indices) != 1 or indices[0] != len(template) - 1:
        raise ValueError('Unexpected prepared launch command')
    index = indices[0]
    return template[:index] + ['--launcher-arg=--time=00:05:00', '/bin/bash', str(wrapper), mode,
                               str(library), str(prefix), str(binary)]


def execute(argv, log, env, timeout=330):
    with log.open('w') as output:
        try:
            process = run.run_child(argv, stdout=output, stderr=subprocess.STDOUT, env=env, timeout=timeout)
            return process.returncode
        except subprocess.TimeoutExpired:
            return 124


def compare_controls(records):
    comparisons = []
    groups = {(r['config'], r['layout']) for r in records if r['mode'] == 'control'}
    for config, layout in sorted(groups):
        group = [r for r in records if r['mode'] == 'control' and r['config'] == config and r['layout'] == layout]
        rounds = {(r['variant'], r['repeat']) for r in group}
        expected = {(v, i) for v in ('reference', 'candidate') for i in range(5)}
        if len(group) != 10 or rounds != expected or any(r.get('status') != 'PASS' for r in group):
            comparisons.append({'config': config, 'layout': layout, 'status': 'INCOMPLETE'})
            continue
        import statistics
        for case in performance.CASES:
            for metric, allowance in performance.METRICS.items():
                values = {v: [r['audit']['cases'][case][metric] for r in group if r['variant'] == v]
                          for v in ('reference', 'candidate')}
                reference, candidate = (statistics.median(values[v]) for v in ('reference', 'candidate'))
                comparisons.append({'config': config, 'layout': layout, 'ranks': 4, 'case': case, 'metric': metric,
                    'samples': values, 'reference_seconds': reference, 'candidate_seconds': candidate,
                    'pass': candidate <= 1.10 * reference + allowance})
    return comparisons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--backend', choices=('cpu', 'cuda'), required=True)
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()
    args.work = args.work.resolve()
    if args.work.exists():
        parser.error('Use a new result directory')
    args.work.mkdir(parents=True)
    kit = Path(__file__).resolve().parent
    manifest = json.loads((kit / 'manifest.json').read_text())
    kit_files = {str(kit / name): checksum for name, checksum in manifest['files'].items()}
    pins = json.loads((kit / 'pins.json').read_text())
    pinned = pins[args.backend]
    result = {'status': 'INCOMPLETE', 'scope': 'Diagnostic only; profiled timings are not release acceptance evidence',
              'backend': args.backend, 'commands': [], 'records': [], 'kit_manifest_sha256': run.digest(kit / 'manifest.json')}
    result_path = args.work / 'result.json'
    def save():
        result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    try:
        verify_files(kit_files)
        verify_files(pinned['files'])
        for name in ('source', 'reference'):
            if run.verify_snapshot(Path(pins[name])) != pins[name + '_manifest_sha256']:
                raise ValueError('Source manifest differs')
        result['input_files'] = pinned['files']
        build = args.work / 'observer-build'
        env = run.execution_environment({})
        argv = ['cmake', '-S', str(kit), '-B', str(build), '-DCMAKE_BUILD_TYPE=Release',
                '-DCMAKE_CXX_COMPILER=' + pinned['compiler'], '-DKOKKOS_TOOLS_INCLUDE=' + pins['kokkos_tools_include']]
        for name, cmd in (('configure', argv), ('build', ['cmake', '--build', str(build), '--parallel', '8'])):
            result['commands'].append(cmd)
            if execute(cmd, args.work / (name + '.log'), env, 600):
                raise ValueError('Observer ' + name + ' failed')
        library = build / 'libmars_structured_trace.so'
        probe = build / 'mars_mpi_trace_probe'
        if not library.is_file() or not probe.is_file():
            raise ValueError('Missing native observer library or probe')
        result['observer_files'] = {str(path): run.digest(path) for path in (library, probe)}
        first = next(iter(pinned['configurations'].values()))
        probe_env = run.execution_environment(first['environment'])
        probe_env.update(MARS_GATE_LAYOUT='single', OMP_PROC_BIND='false')
        prefix = args.work / 'probe'
        argv = command(first['templates']['2'], kit / 'rank.sh', library, prefix, probe, 'mpi')
        result['commands'].append(argv)
        if execute(argv, args.work / 'probe.log', probe_env):
            raise ValueError('Native observer probe failed')
        result['probe'] = observer(prefix, 2, 'mpi', probe=True)
        first_name = next(iter(pinned['configurations']))
        first_backend, first_transport = first_name.split('-')
        probe_env.update(OMP_NUM_THREADS=str(first['threads']), KOKKOS_NUM_THREADS=str(first['threads']))
        prefix = args.work / 'profile-probe'
        argv = command(first['templates']['1'], kit / 'rank.sh', library, prefix, first['binaries']['candidate'], 'profile')
        result['commands'].append(argv)
        if execute(argv, args.work / 'profile-probe.log', probe_env):
            raise ValueError('Native Kokkos profiler activation failed')
        result['profile_probe'] = {
            'audit': performance.audit((args.work / 'profile-probe.log').read_text(), 1,
                                       first_backend, first_transport, 1, first['threads']),
            'observer': observer(prefix, 1, 'profile')}
        save()
        for task in schedule(args.backend):
            config = pinned['configurations'][task['config']]
            backend, transport = task['config'].split('-')
            name = '{config}-{mode}-{layout}-np{ranks}-{repeat}-{variant}'.format(**task)
            prefix = args.work / name
            env = run.execution_environment(config['environment'])
            env.update(MARS_GATE_LAYOUT=task['layout'], OMP_NUM_THREADS=str(config['threads']),
                       KOKKOS_NUM_THREADS=str(config['threads']), OMP_PROC_BIND='false')
            argv = command(config['templates'][str(task['ranks'])], kit / 'rank.sh', library, prefix,
                           config['binaries'][task['variant']], task['mode'])
            record = dict(task, name=name, command=argv, environment=run.environment_receipt(env), status='FAILED')
            print('Running ' + name, flush=True)
            record['returncode'] = execute(argv, prefix.with_suffix('.log'), env)
            text = prefix.with_suffix('.log').read_text()
            try:
                if record['returncode']:
                    raise ValueError('Benchmark exited ' + str(record['returncode']))
                record['audit'] = performance.audit(text, task['ranks'], backend, transport,
                                                   2 if task['layout'] == 'multi' else 1, config['threads'])
                if task['mode'] != 'control':
                    record['observer'] = observer(prefix, task['ranks'], task['mode'])
                record['status'] = 'PASS'
            except (ValueError, OSError, KeyError) as error:
                record['error'] = str(error)
            record['mpi_error_lines'] = [line for line in text.splitlines() if 'MARS_MPI_ERROR ' in line]
            result['records'].append(record)
            save()
        result['control_comparisons'] = compare_controls(result['records'])
        verify_files(pinned['files'])
        verify_files(kit_files)
        verify_files(result['observer_files'])
        for name in ('source', 'reference'):
            if run.verify_snapshot(Path(pins[name])) != pins[name + '_manifest_sha256']:
                raise ValueError('Source changed during diagnostics')
        result['status'] = 'DIAGNOSTIC_COLLECTED'
    except Exception as error:
        result['error'] = str(error)
        print('INCOMPLETE: ' + str(error), file=sys.stderr)
    finally:
        result['artifacts'] = {str(p): run.digest(p) for p in args.work.glob('*')
                               if p.is_file() and p != result_path and p.suffix in ('.json', '.log', '.affinity')}
        save()
    print(result['status'] + ': ' + str(result_path))
    return 0 if result['status'] == 'DIAGNOSTIC_COLLECTED' else 1


if __name__ == '__main__':
    sys.exit(main())
