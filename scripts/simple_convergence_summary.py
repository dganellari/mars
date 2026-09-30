#!/usr/bin/env python3
"""Summarize a completed SIMPLE run locally without exporting private values."""

import csv
import json
import math
from pathlib import Path
import re
import statistics
import sys

from simple_public_diagnostics import DiagnosticState, SafeParser


COLUMNS = ('iteration', 'momentum', 'continuity', 'mass_balance', 'du', 'dp',
           'dflux', 'cancellation', 'inlet_kg_s', 'outlet_kg_s', 'umax_m_s',
           'closed_faces', 'changed_faces')
INTEGER_COLUMNS = ('iteration', 'closed_faces', 'changed_faces')
LOG_COLUMNS = dict(iteration='iteration', momentum='momentum', continuity='continuity',
                   balance='mass_balance', du='du', dp='dp', dflux='dflux',
                   umax='umax_m_s', closed='closed_faces', changed='changed_faces')
HEADER = re.compile(r'SIMPLE Tet4, ([1-9][0-9]*) ranks \(ElementDomain/cstone\), '
                    r'(upwind|high-resolution), laminar')
COMPLETION = re.compile(r'(CONVERGED|NOT CONVERGED: iteration limit) '
                        r'iterations=([0-9]+) ranks=([1-9][0-9]*) exchange_rounds=[0-9]+')


def require(condition):
    if not condition:
        raise ValueError('Incomplete or inconsistent convergence evidence')


def trend(rows, key, limit):
    # Three equal blocks cover the last 30% of the completed iterations.
    size = (len(rows) - 1) // 10
    if size < 5:
        return 'insufficient_data'
    blocks = [rows[-3*size:-2*size], rows[-2*size:-size], rows[-size:]]
    if all(row[key] <= limit for row in blocks[-1]):
        return 'below_target'
    medians = [statistics.median_low(row[key] for row in block) for block in blocks]
    if medians[0] > medians[2] and all(b <= .9*a for a, b in zip(medians, medians[1:])):
        return 'decreasing'
    if medians[2] > medians[0] and all(b >= 1.1*a for a, b in zip(medians, medians[1:])):
        return 'increasing'
    if max(medians) <= 1.1*min(medians):
        return 'flat'
    return 'mixed'


def summarize(log, metrics, exit_text, residual, mass, change):
    require(all(math.isfinite(value) and value > 0 for value in (residual, mass, change)))
    state = DiagnosticState()
    headers, completions, reports = [], [], []
    for raw in log:
        state.feed(raw)
        line = raw.strip()
        if line.startswith('SIMPLE Tet4,'):
            match = HEADER.fullmatch(line)
            require(match is not None)
            headers.append(int(match.group(1)))
        if line.startswith(('CONVERGED', 'NOT CONVERGED')):
            match = COMPLETION.fullmatch(line)
            require(match is not None)
            completions.append((int(match.group(2)), int(match.group(3))))
        if line.startswith('[simple]'):
            tokens = line.split()[1:]
            pairs = [token.split('=', 1) for token in tokens]
            require(all(len(pair) == 2 for pair in pairs))
            record = dict(pairs)
            require(len(record) == len(pairs) and set(record) == set(LOG_COLUMNS))
            report = {target: (int(record[source]) if target in INTEGER_COLUMNS else float(record[source]))
                      for source, target in LOG_COLUMNS.items()}
            require(not reports or report['iteration'] > reports[-1]['iteration'])
            reports.append(report)

    status = state.result(exit_text)['run_status']
    require(status in ('converged', 'iteration_limit'))
    require(len(headers) == len(completions) == 1 and bool(reports))
    final_iteration, ranks = completions[0]
    require(headers[0] == ranks and reports[-1]['iteration'] == final_iteration)

    reader = csv.DictReader(metrics)
    require(reader.fieldnames == list(COLUMNS))
    rows = []
    for record in reader:
        require(set(record) == set(COLUMNS) and all(value is not None for value in record.values()))
        row = {key: (int(record[key]) if key in INTEGER_COLUMNS else float(record[key])) for key in COLUMNS}
        require(row['iteration'] == len(rows))
        require(all(math.isfinite(value) for value in row.values()))
        require(all(value >= 0 for key, value in row.items() if key not in ('inlet_kg_s', 'outlet_kg_s')))
        rows.append(row)
    require(bool(rows) and rows[-1]['iteration'] == final_iteration)
    # Both streams print the same double values with 17 significant digits.
    for report in reports:
        require(0 <= report['iteration'] < len(rows))
        require(all(rows[report['iteration']][key] == value for key, value in report.items()))

    last = rows[-1]
    limits = dict(momentum=residual, continuity=residual, mass_balance=mass,
                  du=change, dp=change, dflux=change)
    checks = {key: last[key] <= limit for key, limit in limits.items()}
    checks.update(minimum_iterations=final_iteration >= 2, finite=True,
                  conservation=last['cancellation'] <= 1e-10,
                  stable_outlet_flags=last['changed_faces'] == 0)
    result = dict(schema='mars-simple-convergence-v1', run_status=status,
                  targets_source='command_line', all_final_checks_passed=all(checks.values()),
                  completion_matches_supplied_targets=(status == 'converged') == all(checks.values()),
                  closed_outlet_faces_present=last['closed_faces'] > 0,
                  small_updates_with_residual_failure=(all(checks[key] for key in ('du', 'dp', 'dflux'))
                                                      and not all(checks[key] for key in ('momentum', 'continuity'))))
    result.update(('final_' + key + '_passed', value) for key, value in checks.items())
    result.update(('tail_' + key + '_trend', trend(rows, key, limit)) for key, limit in limits.items())
    size = (len(rows) - 1) // 10
    result['outlet_flags_changed_in_tail'] = any(row['changed_faces'] > 0 for row in rows[-3*size:]) if size >= 5 else None
    return result


def main():
    parser = SafeParser(description=__doc__)
    parser.add_argument('--log', required=True)
    parser.add_argument('--metrics', required=True)
    parser.add_argument('--exit-file', required=True)
    parser.add_argument('--residual-tol', required=True, type=float)
    parser.add_argument('--mass-tol', required=True, type=float)
    parser.add_argument('--change-tol', required=True, type=float)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    try:
        with open(args.log, encoding='utf-8', errors='replace') as log, open(args.metrics, newline='') as metrics:
            result = summarize(log, metrics, Path(args.exit_file).read_text(),
                               args.residual_tol, args.mass_tol, args.change_tol)
        with open(args.output, 'x') as output:
            json.dump(result, output, sort_keys=True, indent=2, allow_nan=False)
            output.write('\n')
    except Exception:
        print('ERROR: missing, invalid or conflicting convergence evidence; check files privately.', file=sys.stderr)
        return 1
    print('Convergence summary written. Share only that JSON file; export success is not convergence.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
