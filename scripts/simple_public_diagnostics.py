#!/usr/bin/env python3
"""Run beside a private SIMPLE log; export fixed labels and booleans only."""

import argparse
import json
import math
import re
import sys


TAGS = ('[simple-linear]', '[HypreGMRES] rejected:', '[hypre-spmv]')


def fields(text):
    result = {}
    for key, value in re.findall(r'\b([a-z_]+)=([^\s]+)', text):
        # A merged MPI line must not silently replace an earlier value.
        result[key] = None if key in result else value
    return result


def number(record, key):
    try:
        return float(record[key])
    except (KeyError, TypeError, ValueError, OverflowError):
        return None


def finite(record, key):
    value = number(record, key)
    return None if value is None else math.isfinite(value)


def passed(record, residual_key, limit_key):
    residual = number(record, residual_key)
    limit = number(record, limit_key)
    if residual is None or limit is None:
        return None
    return (math.isfinite(residual) and math.isfinite(limit)
            and 0 <= residual <= limit)


def enum_value(record, key, allowed):
    value = record.get(key)
    return value if value in allowed else 'unknown'


def flag(record, key):
    return {'0': False, '1': True}.get(record.get(key))


def cap_reached(record):
    match = re.fullmatch(r'([0-9]+)/([0-9]+)', record.get('iterations') or '')
    if not match:
        return None
    done, cap = map(int, match.groups())
    return done >= cap if cap > 0 else None


def error_nonzero(record):
    value = record.get('solve_error') or ''
    return int(value) != 0 if re.fullmatch(r'[0-9]+', value) else None


def consensus(records, extract, unknown=None):
    values = [extract(record) for record in records]
    return values[0] if values and all(v == values[0] for v in values) else unknown


def summarize(lines, exit_text):
    records = {tag: [] for tag in TAGS}
    completions = []
    false_convergence = False
    for line in lines:
        line = line.strip()
        for tag in TAGS:
            if line.startswith(tag):
                record = fields(line[len(tag):])
                if record not in records[tag]:
                    records[tag].append(record)
        if re.fullmatch(r'CONVERGED iterations=[0-9]+ ranks=[0-9]+ exchange_rounds=[0-9]+', line):
            completions.append('converged')
        if re.fullmatch(r'NOT CONVERGED: iteration limit iterations=[0-9]+ ranks=[0-9]+ exchange_rounds=[0-9]+', line):
            completions.append('iteration_limit')
        if re.match(r'^false convergence [12](?:,|$)', line):
            false_convergence = True

    exit_text = exit_text.strip()
    exit_code = int(exit_text) if re.fullmatch(r'[0-9]{1,3}', exit_text) else None
    status = 'incomplete'
    if exit_code is not None and exit_code not in (0, 2):
        status = 'failed'
    elif len(completions) == 1:
        if exit_code == 0 and completions[0] == 'converged':
            status = 'converged'
        elif exit_code == 2 and completions[0] == 'iteration_limit':
            status = 'iteration_limit'

    simple = records['[simple-linear]']
    hypre = records['[HypreGMRES] rejected:']
    spmv = records['[hypre-spmv]']
    if simple or hypre:
        # A concatenated successful run cannot hide a rejection.
        status = 'failed'
    return {
        'schema': 'mars-simple-public-diagnostics-v1',
        'run_status': status,
        'multiple_completions': len(completions) > 1,
        'simple_diagnostic_present': bool(simple),
        'hypre_rejection_present': bool(hypre),
        'linear_stage': consensus(simple, lambda r: enum_value(r, 'stage', ('momentum', 'pressure')), 'unknown'),
        'backend': consensus(hypre, lambda r: enum_value(r, 'backend', ('GMRES', 'FlexGMRES')), 'unknown'),
        'solver_accepted': consensus(simple, lambda r: flag(r, 'solver_accepted')),
        'mars_passed': consensus(simple, lambda r: flag(r, 'mars_passed')),
        'iteration_cap_reached': consensus(hypre, cap_reached),
        'hypre_error_nonzero': consensus(hypre, error_nonzero),
        'reported_residual_finite': consensus(hypre, lambda r: finite(r, 'reported_relative')),
        'reported_residual_passed': consensus(hypre, lambda r: passed(r, 'reported_relative', 'tolerance')),
        'hypre_explicit_residual_finite': consensus(hypre, lambda r: finite(r, 'absolute_residual')),
        'hypre_explicit_residual_passed': consensus(hypre, lambda r: passed(r, 'absolute_residual', 'acceptance_limit')),
        'mars_explicit_residual_finite': consensus(simple, lambda r: finite(r, 'mars_absolute_residual')),
        'mars_explicit_residual_passed': consensus(simple, lambda r: passed(r, 'mars_absolute_residual', 'acceptance_limit')),
        'spmv_backend': consensus(spmv, lambda r: {'0': 'native', '1': 'vendor'}.get(r.get('vendor_requested'), 'unknown'), 'unknown'),
        'false_convergence_message_seen': false_convergence,
    }


class SafeParser(argparse.ArgumentParser):
    def error(self, message):
        # argparse normally echoes unrecognized arguments and their private values.
        self.exit(2, 'ERROR: invalid diagnostic arguments; use --help.\n')


def main():
    parser = SafeParser(description=__doc__)
    parser.add_argument('--log', required=True)
    parser.add_argument('--exit-file', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    try:
        with open(args.log, encoding='utf-8', errors='replace') as log:
            with open(args.exit_file, encoding='ascii') as exit_file:
                result = summarize(log, exit_file.read(64))
        # An exclusive create prevents accidental replacement of an input or old summary.
        with open(args.output, 'x', encoding='ascii') as output:
            json.dump(result, output, indent=2, sort_keys=True, allow_nan=False)
            output.write('\n')
    except Exception:
        print('ERROR: diagnostic export failed; check files privately.', file=sys.stderr)
        return 1
    print('Public diagnostic summary written. Share only that JSON file.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
