"""Match host trace intervals to the unchanged benchmark's MPI timers."""
import collections
import json
from pathlib import Path

import performance


def intervals(events):
    timers = sorted((e for e in events if e['kind'] == 'timer'), key=lambda e: e['start'])
    # Each workload has mesh, DOF and twenty halo timer pairs, repeated after warmup.
    if len(timers) != 12 * 44:
        raise ValueError('Expected 528 benchmark timer calls, got ' + str(len(timers)))
    result = {}
    for i, case in enumerate(performance.CASES):
        markers = timers[(6 + i) * 44:(7 + i) * 44]
        row = {}
        for phase, offset in (('mesh_seconds', 0), ('dof_seconds', 2)):
            first, last = markers[offset:offset + 2]
            elapsed = last['value'] - first['value']
            if elapsed <= 0 or last['start'] <= first['start']:
                raise ValueError('Invalid setup timer interval')
            inside = [e for e in events if e['kind'] != 'timer' and
                      first['start'] <= e['start'] <= last['start']]
            totals = collections.defaultdict(lambda: {'calls': 0, 'seconds': 0.0})
            for e in inside:
                name = e['kind'] + ':' + e['name']
                totals[name]['calls'] += 1
                totals[name]['seconds'] += e['seconds']
            row[phase] = {'seconds': elapsed, 'totals': dict(totals),
                          'longest_calls': sorted(inside, key=lambda e: e['seconds'], reverse=True)[:20]}
        row['halo_seconds'] = sum(markers[j + 1]['value'] - markers[j]['value']
                                  for j in range(4, 44, 2)) / 20
        result[case] = row
    return result


def audit(prefix, ranks, benchmark, cuda=False):
    records = []
    for rank in range(ranks):
        data = json.loads(Path(str(prefix) + '.rank' + str(rank) + '.json').read_text())
        if data['global_fencing']:
            raise ValueError('Unfenced tracing added global fences')
        if cuda:
            calls = {e['name'] for e in data['events'] if e['kind'] == 'cuda'}
            if not {'cudaMalloc', 'cudaFree', 'cudaLaunchKernel', 'cudaStreamSynchronize'} <= calls:
                raise ValueError('CUDA runtime interception was not active: ' + str(calls))
            if any(e['kind'] == 'cuda' and e['error'] for e in data['events']):
                raise ValueError('CUDA returned an error during tracing')
            if 'cudaMallocAsync' in calls:
                raise ValueError('Unexpected pooled allocation in the corrected dependency')
        records.append({'rank': rank, 'cases': intervals(data['events'])})
    for case in performance.CASES:
        for metric in performance.METRICS:
            values = [r['cases'][case][metric] for r in records]
            if metric != 'halo_seconds':
                values = [v['seconds'] for v in values]
            if abs(max(values) - benchmark['cases'][case][metric]) > 1.e-12:
                raise ValueError('Trace timer differs from benchmark: ' + case + '/' + metric)
    return records
