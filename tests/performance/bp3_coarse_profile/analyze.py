#!/usr/bin/env python3
"""Validate solves and report level-local host, GPU and MPI measurements separately."""
import argparse
from bisect import bisect_right
import csv
import io
import json
import math
from pathlib import Path
import re

LEVEL = re.compile(r'bp3\.amg/L(\d+)$')


def result(path):
    text = Path(path).read_text()
    found = re.findall(r'FCG \+ V-cycle: (\d+) iterations, sqrt\(<r,z>/<r0,z0>\) = (\S+), true assembled residual (\S+)', text)
    if len(found) != 1 or not re.search(r'\bPASS\b', text) or re.search(r'\bFAIL\b', text):
        raise ValueError(f'{path}: missing or failed solver gate')
    iterations, reported, true = found[0]
    iterations, reported, true = int(iterations), float(reported), float(true)
    if not (0 < iterations < 2000 and math.isfinite(reported) and 0 <= reported <= 1e-13 and
            math.isfinite(true) and 0 <= true <= 1e-10):
        raise ValueError(f'{path}: convergence or true residual gate failed')
    timing = re.search(r'FCG iteration\s+([\d.]+) ms\s+\(solve ([\d.]+) ms', text)
    if not timing: raise ValueError(f'{path}: no solve timing')
    return dict(iterations=iterations, relative=reported, true_residual=true,
                iteration_ms=float(timing[1]), solve_ms=float(timing[2]),
                history=re.findall(r'^\s+it\s+(\d+)\s+(\S+)', text, re.M))


def check_audits(directory, ranks):
    paths = sorted(directory.glob('rank*.audit.json'))
    if len(paths) != ranks: raise ValueError(f'Expected {ranks} hook audits, got {len(paths)}')
    for rank in range(ranks):
        data = json.loads((directory / f'rank{rank}.audit.json').read_text())
        if not (data['pushes'] == data['pops'] and data['cycles'] > 0 and data['errors'] == 0
                and data['levels'][0] > 0 and data['levels'][1] > 0):
            raise ValueError(f'Hook missing or unbalanced on rank {rank}: {data}')


def key(name):
    return re.sub(r'[^a-z0-9]', '', re.sub(r'\([^)]*\)', '', name).lower())


def table(path, required):
    # nsys prints a short preamble before CSV; find the actual header explicitly.
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        columns = [key(x) for x in next(csv.reader([line]))]
        if set(required) <= set(columns):
            records = list(csv.DictReader(io.StringIO('\n'.join(lines[i:]))))
            return [{key(k): v for k, v in row.items() if k is not None} for row in records]
    raise ValueError(f'{path}: report absent or unrecognized CSV columns')


def trace_summary(directory, rank, require_mpi=True):
    prefix = directory / f'rank{rank}'
    read = lambda kind, columns: table(Path(str(prefix)+'.'+kind+'.csv'), columns)
    ranges = read('nvtx_pushpop_trace', ('name', 'start', 'end', 'pid', 'tid'))
    levels, intervals, cycles = {}, {}, 0
    for row in ranges:
        name = row['name']
        if name.endswith('bp3.amg/cycle'): cycles += 1
        match = LEVEL.search(name)
        if not match: continue
        level = int(match[1])
        entry = levels.setdefault(level, dict(host_ns=0, kernel_ns=0, projected_span_ns=0,
                                              mpi_ns=0, mpi_wait_ns=0, visits=0, kernels=0))
        start, end = float(row['start']), float(row['end'])
        entry['host_ns'] += end-start
        entry['visits'] += 1
        intervals.setdefault((row['pid'], row['tid']), []).append((start, end, level))
    if cycles == 0 or len(levels) < 2: raise ValueError(f'{prefix}: no traced multilevel AMG cycles')
    for ranges_on_thread in intervals.values():
        ranges_on_thread.sort()
        if any(a[1] > b[0] for a, b in zip(ranges_on_thread, ranges_on_thread[1:])):
            raise ValueError('Unexpected nested AMG level roots; cannot sum them')
    for row in read('nvtx_kern_sum', ('nvtxrange', 'totaltime', 'kerninst')):
        match = LEVEL.search(row['nvtxrange'])
        if match:
            entry = levels[int(match[1])]
            entry['kernel_ns'] += float(row['totaltime'])
            entry['kernels'] += int(row['kerninst'])
    if not sum(v['kernels'] for v in levels.values()):
        raise ValueError(f'{prefix}: no GPU kernels correlated to level ranges')
    for row in read('nvtx_gpu_proj_sum', ('range', 'totalprojtime')):
        match = LEVEL.search(row['range'])
        if match: levels[int(match[1])]['projected_span_ns'] += float(row['totalprojtime'])
    # Dedupe split send/receive report rows. MPI call time is not network-only time.
    try:
        mpi = read('mpi_event_trace', ('start', 'end', 'event', 'pid', 'tid'))
    except ValueError:
        stderr = (Path(str(prefix)+'.mpi_event_trace.stderr').read_text() +
                  Path(str(prefix)+'.mpi_event_trace.csv').read_text())
        if require_mpi or not re.search(r'(does not contain.*MPI|no MPI.*events)', stderr, re.I): raise
        mpi = []
    if require_mpi and not mpi: raise ValueError(f'{prefix}: no MPI events on a multi-rank run')
    starts = {thread: [x[0] for x in spans] for thread, spans in intervals.items()}
    seen, outside = set(), 0
    for row in mpi:
        identity = tuple(row[k] for k in ('pid', 'tid', 'start', 'end', 'event'))
        if identity in seen: continue
        seen.add(identity)
        thread = row['pid'], row['tid']
        start, end = float(row['start']), float(row['end'])
        pos = bisect_right(starts.get(thread, []), start)-1
        if pos < 0 or end > intervals[thread][pos][1]:
            outside += 1
            continue
        entry = levels[intervals[thread][pos][2]]
        entry['mpi_ns'] += end-start
        if 'wait' in row['event'].lower(): entry['mpi_wait_ns'] += end-start
    for entry in levels.values():
        for field in ('host', 'kernel', 'projected_span', 'mpi', 'mpi_wait'):
            entry[field+'_ms_per_amg_cycle'] = entry.pop(field+'_ns') / 1e6 / cycles
    return dict(rank=rank, captured_amg_cycles=cycles, levels=levels, mpi_events_outside_levels=outside)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('log', type=Path)
    parser.add_argument('--reference', type=Path)
    parser.add_argument('--audits', type=Path)
    parser.add_argument('--ranks', type=int, default=1)
    parser.add_argument('--trace', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    answer = {'solve': result(args.log)}
    if args.reference:
        reference = result(args.reference)
        if answer['solve']['iterations'] != reference['iterations']:
            raise ValueError('Instrumentation changed iteration count')
        h1, h2 = answer['solve']['history'], reference['history']
        if len(h1) != len(h2) or any(a[0] != b[0] or not math.isclose(float(a[1]), float(b[1]), rel_tol=5e-6)
                                   for a, b in zip(h1, h2)):
            raise ValueError('Instrumentation changed the printed FCG history')
        answer['hook_to_stock_solve_ratio'] = answer['solve']['solve_ms']/reference['solve_ms']
    if args.audits: check_audits(args.audits, args.ranks)
    if args.trace:
        answer['ranks'] = [trace_summary(args.audits, rank, args.ranks > 1) for rank in range(args.ranks)]
        answer['interpretation'] = ('Per-rank measurements; never add CPU, GPU and MPI times. '
            'GPU projected spans include gaps and may overlap across levels. Kernel totals exclude copies. '
            'Level roots exclude their duplicate child summaries; MPI waits may reflect remote imbalance.')
    args.output.write_text(json.dumps(answer, indent=2)+'\n')
    print(json.dumps(answer, indent=2))


if __name__ == '__main__': main()
