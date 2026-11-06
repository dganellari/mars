"""Read native hierarchy receipts without changing raw logs or rank audits."""
import json
import math
import re
from common import parse, save

def parse_h(out, ranks, ne, depth, checked, solves=3):
    text = (out / 'solver.log').read_text()
    geometry = re.findall(r'BP3 geometric levels: (\d+), global cells per axis ([\d ]+); bottom rows (\d+); smoother bounds ([\d.eE+ -]+)', text)
    if len(geometry) != 1:
        raise ValueError('Missing unique h hierarchy receipt')
    d, cells, rows, bounds = geometry[0]
    if (int(d) != depth or [int(c) for c in cells.split()] != [ne >> k for k in range(depth + 1)] or
            int(rows) != ((ne >> depth) - 1)**3):
        raise ValueError('Wrong geometric hierarchy or free row count')
    values = [float(x) for x in bounds.split()]
    if len(values) != depth or not all(math.isfinite(v) and v > 0 for v in values):
        raise ValueError('Invalid smoother bounds')
    checks = re.findall(r'PASS h validation: levels=(\d+), adjoint/Galerkin/symmetry/linearity/continuity/repeat worst (\S+)', text)
    if (len(checks) != int(checked) or checked and
            (int(checks[0][0]) != depth or not 0 <= float(checks[0][1]) < 5e-11)):
        raise ValueError('Missing or unexpected direct h validation')
    # Ordinary sizes use the unchanged strict parser. Tiny coarse grids can reach
    # Hypre's nine-row stopping rule before its consolidation callback is needed.
    audits = [json.loads((out / ('rank%d.json' % r)).read_text()) for r in range(ranks)]
    inactive = ranks > 1 and all(a['setups'] == 0 for a in audits)
    if inactive:
        tail = text.split('Operator Matrix Information:', 1)[1].split('Interpolation Matrix Information:', 1)[0]
        matrix_rows = [(int(a), int(b)) for a, b in re.findall(r'^\s+(\d+)\s+(\d+)\s+\d+\s+[\d.]+\s+', tail, re.M)]
        if int(rows) >= 64 or not matrix_rows or matrix_rows[0] != (0, int(rows)) or matrix_rows[-1][1] > 9:
            raise ValueError('Inactive adapter without a verified tiny terminal matrix')
        limit = min(50000, max(5000, 1000 * ranks))
        if any(a != dict(version=1, setups=0, destroyed=0, cycles=0, empty_setups=0, skipped_setups=0,
                         threshold_rows=limit, level=-1, rows=0, nnz=0) for a in audits):
            raise ValueError('Unexpected inactive adapter state')
        matches = re.findall(r'FCG \+ V-cycle: (\d+) iterations, sqrt\(<r,z>/<r0,z0>\) = (\S+), true assembled residual (\S+)', text)
        if len(matches) != 1 or not re.search(r'\bPASS\b', text) or re.search(r'\bFAIL\b', text):
            raise ValueError('Missing or failed small-grid solve')
        it, relative, residual = int(matches[0][0]), float(matches[0][1]), float(matches[0][2])
        if not (0 < it < 2000 and 0 <= relative <= 1e-13 and 0 <= residual <= 1e-10):
            raise ValueError('Small-grid accuracy gate failed')
        wall = re.search(r'host wall solves (.+) ms \(max over ranks; last is reported\)', text)
        event = re.search(r'FCG iteration\s+[\d.]+ ms\s+\(solve ([\d.]+) ms, \+ setup ([\d.]+) ms', text)
        coarse = re.search(r'coarse solve\s+([\d.]+) ms per V-cycle', text)
        times = [float(x) for x in wall[1].split()] if wall else []
        if len(times) != solves or not event or not coarse or not all(math.isfinite(x) and x > 0 for x in times):
            raise ValueError('Missing small-grid timings')
        if not float(event[1]) > 0 or not float(coarse[1]) > 0 or (times[-1] > 10 and abs(times[-1] / float(event[1]) - 1) > .03):
            raise ValueError('Small-grid timings disagree')
        if 'Hypre: SpMV hypre, SpGEMM hypre, GPU-aware MPI on' not in text:
            raise ValueError('Wrong small-grid Hypre policy')
        if len(list(out.glob('rank*.json'))) != ranks:
            raise ValueError('Unexpected small-grid rank receipt count')
        result = dict(iterations=it, relative=relative, true_residual=residual, wall_ms=times[-1],
                      wall_solves_ms=times, setup_ms=float(event[2]), coarse_ms=float(coarse[1]),
                      hierarchy=[-1, 0, 0, 0], tiny_terminal_without_consolidation=True)
    else:
        result = parse(out, ranks, -1, solves=solves)
    result.update(h_levels=depth, h_bottom_rows=int(rows), h_bounds=values,
                  h_check_worst=float(checks[0][1]) if checked else None)
    save(out / 'summary.json', result)
    return result
