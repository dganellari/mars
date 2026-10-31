#!/usr/bin/env python3
"""Keep probe overhead and per-rank trace attribution distinct from stock timings."""
import json
from pathlib import Path
import statistics
import sys

root = Path(sys.argv[1])
for trace in sorted(root.glob('p*_r*_trace/summary.json')):
    case = trace.parent.name.removesuffix('_trace')
    print(case)
    for mode in ('stock', 'hook'):
        times = [json.loads(p.read_text())['solve']['solve_ms'] for p in sorted(root.glob(f'{case}_{mode}*/summary.json'))]
        if len(times) != 3: raise ValueError(f'{case}: expected three {mode} rounds')
        print(f'  {mode} driver solve ms min/median/max: {min(times):.3f} / {statistics.median(times):.3f} / {max(times):.3f}')
    ratios = [json.loads(p.read_text())['hook_to_stock_solve_ratio'] for p in root.glob(f'{case}_hook*/summary.json')]
    print(f'  hook/stock median paired ratio: {statistics.median(ratios):.4f}')
    data = json.loads(trace.read_text())
    print('  Level: median / max across ranks, ms per AMG cycle (host, kernels, MPI waits)')
    levels = sorted({int(level) for r in data['ranks'] for level in r['levels']})
    for level in levels:
        parts = []
        for field in ('host_ms_per_amg_cycle', 'kernel_ms_per_amg_cycle', 'mpi_wait_ms_per_amg_cycle'):
            values = [r['levels'].get(str(level), {}).get(field, 0) for r in data['ranks']]
            parts.append(f'{statistics.median(values):.4f}/{max(values):.4f}')
        print(f'    L{level:02d}: '+', '.join(parts))
print('Host, GPU and MPI times overlap. Do not add them or sum different ranks\' maxima.')
