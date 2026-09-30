#!/usr/bin/env python3
"""Prepare an early matched-state run privately; never launch a solver."""

import json
import os
from pathlib import Path
import re
import sys

from simple_convergence_summary import COMPLETION, summarize
from simple_public_diagnostics import SafeParser
from simple_snapshot_compare import controls, digest, finite, options, require, result_files


def first_snapshot(paths, maximum, before):
    import numpy as np
    from netCDF4 import Dataset
    times = None
    for path in paths:
        with Dataset(str(path)) as ds:
            values = finite(ds.variables['time_whole'][:])
        require(values.ndim == 1 and len(values) > 0 and np.all(np.diff(values) > 0))
        require(np.all(values >= 0))
        if times is None:
            times = values
        else:
            require(np.array_equal(times, values))
    positive = times[times > 0] if times is not None else []
    require(len(positive) > 0)
    first = float(positive[0])
    require(first.is_integer() and first <= maximum and first < before)
    return int(first)


def prepare(args, public):
    public['failed_check'] = 'baseline_completion'
    log = (args.baseline / 'run.log').read_text(errors='replace').splitlines()
    with (args.baseline / 'flow-metrics.csv').open() as metrics:
        status = summarize(log, metrics, (args.baseline / 'run.exit').read_text(),
                           args.residual_tol, args.mass_tol, args.change_tol)
    require(status['completion_matches_supplied_targets'])
    endings = [COMPLETION.fullmatch(line.strip()) for line in log
               if line.startswith(('CONVERGED', 'NOT CONVERGED'))]
    require(len(endings) == 1 and endings[0] is not None)
    baseline_iteration, ranks = int(endings[0].group(2)), int(endings[0].group(3))
    require(1 <= ranks <= 4)  # The documented launch uses one four-GPU node.
    public['baseline_completion_checked'] = True

    public['failed_check'] = 'binary_changed_or_unavailable'
    recorded = (args.baseline / 'executable.sha256').read_text().splitlines()
    require(len(recorded) == 1)
    match = re.fullmatch(r'([0-9a-f]{64})[ \t]+.+', recorded[0])
    require(match is not None and args.executable.is_file() and os.access(str(args.executable), os.X_OK))
    executable_hash = digest(args.executable)
    require(match.group(1) == executable_hash)
    public['binary_matches_baseline'] = True

    public['failed_check'] = 'saved_arguments_or_controls'
    prepared, _, reference_status, same_deck = controls(args.case, args.reference_dir, None, log)
    require(same_deck and reference_status == 'mapped_controls_match')
    require(Path(prepared['--mesh']).is_file())
    saved = json.loads(args.case.read_text())['arguments']
    raw = (args.case.parent / 'args.nul').read_bytes()
    require(raw.endswith(b'\0') and raw[:-1].split(b'\0') == [value.encode('utf-8') for value in saved])
    require(prepared == options(saved))
    public['prepared_and_reference_controls_match'] = True

    public['failed_check'] = 'baseline_runtime_options'
    runtime = [line.strip() for line in log if line.startswith('linear_cache=')]
    require(len(runtime) == 1)
    flags = re.fullmatch(r'linear_cache=([01]) halo_overlap=([01]) field_output=(gathered|distributed) profile=0', runtime[0])
    require(flags is not None)
    public['runtime_options_and_ranks_preserved'] = True

    public['failed_check'] = 'early_reference_state'
    require(0 < args.max_iteration <= 100)
    paths = result_files(args.reference_dir)
    iteration = first_snapshot(paths, args.max_iteration, baseline_iteration)
    public['early_reference_state_available'] = True
    extra = ['--iterations', str(iteration), '--report-every', '1',
             '--residual-tol', str(args.residual_tol), '--mass-tol', str(args.mass_tol),
             '--change-tol', str(args.change_tol), '--linear-cache', flags.group(1),
             '--halo-overlap', flags.group(2), '--field-output', flags.group(3), '--profile', '0']
    require(not (set(saved[::2]) & set(extra[::2])))
    arguments = saved + extra

    public['failed_check'] = 'private_output'
    args.output_dir.mkdir(mode=0o700)
    (args.output_dir / 'args.nul').write_bytes(b'\0'.join(value.encode('utf-8') for value in arguments) + b'\0')
    (args.output_dir / 'iteration.txt').write_text(str(iteration) + '\n')
    (args.output_dir / 'ranks.txt').write_text(str(ranks) + '\n')
    (args.output_dir / 'probe.json').write_text(json.dumps(dict(
        format='mars-simple-snapshot-probe-v1', iteration=iteration, ranks=ranks,
        baseline=str(args.baseline.resolve()), baseline_iteration=baseline_iteration,
        executable=str(args.executable.resolve()), executable_sha256=executable_hash,
        case_sha256=digest(args.case), arguments=arguments,
        reference_files=[str(path.resolve()) for path in paths],
        note='Same binary and recorded controls; shared-library and unrecorded environment identity not attested.'
    ), indent=2, sort_keys=True) + '\n')
    public.update(preparation_status='ready', failed_check='none')


def main(argv=None):
    parser = SafeParser(description=__doc__)
    for name in ('baseline', 'case', 'reference-dir', 'executable', 'output-dir', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--max-iteration', type=int, default=100)
    for name in ('residual-tol', 'mass-tol', 'change-tol'):
        parser.add_argument('--' + name, type=float, default=1e-6)
    args = parser.parse_args(argv)
    os.umask(0o077)
    public = dict(schema='mars-simple-snapshot-probe-v1', preparation_status='rejected',
                  failed_check='output_preflight', baseline_completion_checked=False,
                  binary_matches_baseline=False, prepared_and_reference_controls_match=False,
                  runtime_options_and_ranks_preserved=False, early_reference_state_available=False)
    try:
        with args.output.open('x') as output:
            try:
                require(not args.output_dir.exists())
                prepare(args, public)
            except Exception:
                pass  # File, YAML and netCDF exceptions can contain private data.
            json.dump(public, output, indent=2, sort_keys=True)
            output.write('\n')
    except Exception:
        print('ERROR: cannot create a fresh public preparation summary.', file=sys.stderr)
        return 1
    if public['preparation_status'] != 'ready':
        print('Probe preparation rejected. Share only the public preparation JSON.')
        return 1
    print('Probe preparation complete. No solver was launched; private launch files are ready.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
