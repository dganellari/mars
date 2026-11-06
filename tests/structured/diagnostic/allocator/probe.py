#!/usr/bin/env python3
"""Compare CUDA allocation mechanisms without MARS mesh or communication code."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'cluster'))
import run
import mpi_launch


def schedule():
    for layout in ('single', 'multi'):
        for ranks in (2, 4):
            for words, count, offset in ((4096, 1024, 0), (6144, 1024, 128), (8192, 2048, 0)):
                for build, mode in (('original', 'host'), ('original', 'cuda'), ('original', 'async'),
                                    ('original', 'view'), ('compatible', 'view')):
                    yield dict(layout=layout, ranks=ranks, words=words, count=count,
                               offset=offset, build=build, mode=mode)


def audit(text, task):
    def rows(prefix):
        return [json.loads(line.split(prefix, 1)[1]) for line in text.splitlines() if prefix in line]
    passed = rows('MARS_ALLOCATOR_PASS ')
    buffers = rows('MARS_ALLOCATOR_BUFFER ')
    ranks = task['ranks']
    if len(passed) != ranks or {r['rank'] for r in passed} != set(range(ranks)):
        raise ValueError('Missing or duplicate completed ranks')
    nodes = 2 if task['layout'] == 'multi' else 1
    if len({r['hostname'] for r in passed}) != nodes:
        raise ValueError('Wrong node count')
    if len({(r['hostname'], r['device_id']) for r in passed}) != ranks:
        raise ValueError('Ranks share a GPU')
    for receipt in passed + buffers:
        if any(receipt[k] != task[k] for k in ('mode', 'ranks', 'count', 'offset', 'words')):
            raise ValueError('Wrong allocation or message shape')
        if receipt['kokkos_async'] != (task['build'] == 'original'):
            raise ValueError('Wrong Kokkos allocator configuration')
    if any(r['errors'] != 0 or r['exchanges'] != 9 for r in passed):
        raise ValueError('Incorrect exchange')
    expected = {(r, c, label) for r in range(ranks) for c in range(3) for label in ('send', 'recv')}
    if len(buffers) != len(expected) or {(b['rank'], b['cycle'], b['label']) for b in buffers} != expected:
        raise ValueError('Missing or duplicate allocation receipts')
    pooled = task['mode'] == 'async' or (task['mode'] == 'view' and task['build'] == 'original' and task['words'] * 8 + 128 >= 40000)
    for b in buffers:
        if b['pooled'] != pooled or (task['mode'] != 'host' and (b['pool_query'] != 0 or b['memory_type'] != 2)):
            raise ValueError('Allocation provenance differs from requested mode')
    return {'passed': passed, 'buffers': buffers}


def execute(argv, log, env, timeout):
    with log.open('w') as output:
        try:
            return run.run_child(argv, stdout=output, stderr=subprocess.STDOUT, env=env, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            return 124


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('original', 'compatible', 'work', 'cxx', 'gpu-compiler', 'launcher'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--launcher-arg', action='append', default=[])
    args = parser.parse_args()
    if args.work.exists():
        parser.error('Use a fresh probe directory')
    args.work.mkdir(parents=True)
    source = Path(__file__).resolve().parent
    result = {'status': 'INCOMPLETE', 'records': [], 'commands': [], 'inputs': {}}
    for prefix in (args.original, args.compatible):
        for path in prefix.rglob('*'):
            if path.is_file(): result['inputs'][str(path)] = run.digest(path)
    env = run.execution_environment({})
    for name in ('LD_PRELOAD', 'KOKKOS_TOOLS_LIBS', 'KOKKOS_PROFILE_LIBRARY'):
        env.pop(name, None)
    try:
        for variant in ('original', 'compatible'):
            build = args.work / variant
            configure = ['cmake', '-S', str(source), '-B', str(build), '-G', 'Ninja',
                '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_CXX_COMPILER=' + str(args.cxx),
                '-DCMAKE_CUDA_COMPILER=' + str(args.gpu_compiler), '-DCMAKE_CUDA_HOST_COMPILER=' + str(args.cxx),
                '-DCMAKE_CUDA_ARCHITECTURES=90', '-DKokkos_DIR=' + str(getattr(args, variant) / 'lib/cmake/Kokkos')]
            for phase, argv in (('configure', configure), ('build', ['cmake', '--build', str(build), '--parallel', '8'])):
                result['commands'].append(argv)
                if execute(argv, args.work / (variant + '-' + phase + '.log'), env, 600):
                    raise ValueError(variant + ' probe ' + phase + ' failed')
        result['binaries'] = {str(args.work / v / 'mars_mpi_allocator_probe'):
            run.digest(args.work / v / 'mars_mpi_allocator_probe') for v in ('original', 'compatible')}
        for index, task in enumerate(schedule()):
            binary = args.work / task['build'] / 'mars_mpi_allocator_probe'
            step = argparse.Namespace(launcher=str(args.launcher), ranks=task['ranks'],
                launcher_arg=[a for a in args.launcher_arg if not a.startswith('--time=')] + ['--time=00:02:00'],
                binary=[str(binary), task['mode'], str(task['words']), str(task['count']), str(task['offset'])])
            argv = mpi_launch.command(step, task['layout'])
            name = '{:02d}-{build}-{mode}-{layout}-np{ranks}-w{words}'.format(index, **task)
            print('Running ' + name, flush=True)
            record = dict(task, name=name, command=argv, status='FAILED')
            log = args.work / (name + '.log')
            record['returncode'] = execute(argv, log, env, 150)
            if not record['returncode']:
                try:
                    record['audit'] = audit(log.read_text(), task)
                    record['status'] = 'PASS'
                except (ValueError, KeyError) as error:
                    record['error'] = str(error)
            result['records'].append(record)
            (args.work / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
        for path, checksum in dict(result['inputs'], **result['binaries']).items():
            if run.digest(path) != checksum: raise ValueError('Probe input changed: ' + path)
        result['status'] = 'COLLECTED'
        result['compatible_controls_pass'] = all(r['status'] == 'PASS' for r in result['records']
            if r['mode'] in ('host', 'cuda') or r['build'] == 'compatible')
    except Exception as error:
        result['error'] = str(error)
    result['artifacts'] = {str(p): run.digest(p) for p in args.work.glob('*.log')}
    (args.work / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'] + ': ' + str(args.work / 'result.json'))
    return 0 if result['status'] == 'COLLECTED' else 1


if __name__ == '__main__':
    sys.exit(main())
