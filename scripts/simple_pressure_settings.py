#!/usr/bin/env python3
"""Compare saved pressure-solver controls locally; export fixed labels, never private values."""

import math
from pathlib import Path
import re
import sys

from prepare_simple_deck import load_deck
from simple_pressure_probe import resolved_solver
from simple_public_diagnostics import SafeParser
from simple_snapshot_compare import digest, options, require
import simple_startup_probe as startup


SCHEMA = 'mars-simple-pressure-settings-v1'
UNKNOWN = object()
# These are wrapper calls, not inferred defaults of the linked Hypre library.
AMG = {
    'coarsen_type': ('coarsentype', 'MARS_AMG_COARSEN', 8),
    'interp_type': ('interptype', 'MARS_AMG_INTERP', 6),
    'relax_type': ('relaxtype', 'MARS_AMG_RELAX', 18),
    'relax_order': ('relaxorder', 'MARS_AMG_RELAXORDER', 0),
    'strong_threshold': ('strongthreshold', 'MARS_AMG_STRONG', 0.25),
    'p_max_elmts': (None, 'MARS_AMG_PMAX', 4),
    'agg_num_levels': ('aggnumlevels', 'MARS_AMG_AGG', 0),
    'num_sweeps': ('numsweeps', 'MARS_AMG_SWEEPS', 2),
    'max_levels': ('maxlevels', None, 25),
    'min_coarse_size': (None, None, 32),
    'max_coarse_size': ('maxcoarsesize', None, 128),
    'coarse_relax_type': (None, None, 18),
    'keep_transpose': (None, None, 1),
}
AMG_LIBRARY = ('numpaths', 'agginterptype', 'nodal', 'nodaldiag', 'numfunctions',
               'coarsencutfactor', 'cycletype', 'fcycle', 'truncfactor',
               'aggtruncfactor', 'jacobitruncthreshold', 'maxrowsum')


def number(value, integer=False):
    if isinstance(value, bool) or not isinstance(value, (str, int, float)):
        return UNKNOWN
    text = str(value).strip()
    if integer and not re.fullmatch(r'[+-]?[0-9]+', text):
        return UNKNOWN
    try:
        result = int(text) if integer else float(text)
        if not math.isfinite(result) or (integer and not -2**31 <= result < 2**31):
            return UNKNOWN
        return result
    except (ValueError, OverflowError):
        return UNKNOWN


def environment_number(environment, name, default):
    if name is None or name not in environment or environment[name] == '':
        return default
    # The C++ parser accepts numeric prefixes. Do not certify malformed overrides.
    return number(environment[name], integer=isinstance(default, int))


def lower_options(value):
    require(isinstance(value, dict))
    result = {}
    for key, item in value.items():
        require(isinstance(key, str) and key.lower() not in result)
        result[key.lower()] = item
    return result


def compare_settings(config, arguments, environment):
    require(isinstance(config, dict) and isinstance(environment, dict))
    require(all(isinstance(k, str) and isinstance(v, str) for k, v in environment.items()))
    require(str(config.get('family', '')).lower() == 'hypre')
    reference_options = config.get('options', {})
    require(isinstance(reference_options, dict))
    solver = lower_options(reference_options.get('solver', {}))
    precond = lower_options(reference_options.get('precond', {}))
    method = str(solver.get('type', 'gmres')).lower()
    require(method in ('gmres', 'flexgmres'))
    require('options' not in config or 'type' in reference_options.get('solver', {}))
    preconditioner = str(precond.get('type', 'none')).lower()
    require(preconditioner in ('none', 'boomeramg', 'mgr'))
    require('precond' not in reference_options or
            ('type' in reference_options['precond'] and preconditioner != 'none'))
    flex = environment.get('MARS_HYPRE_FLEXGMRES')
    status = {}

    def compare(label, reference, mars):
        status[label] = ('unresolved' if reference is UNKNOWN or mars is UNKNOWN else
                         'equal' if reference == mars else 'different')

    compare('krylov_method', method, 'flexgmres' if flex is not None and flex != '0' else 'gmres')
    compare('preconditioner', preconditioner, 'boomeramg')
    compare('restart_dimension', number(solver.get('kdim'), integer=True), 100)
    compare('maximum_iterations', number(config.get('max_iterations', 20), integer=True), 2000)
    # ContextHYPRE never forwards the base class's min_iterations to Hypre.
    compare('minimum_iterations', UNKNOWN, environment_number(environment, 'MARS_HYPRE_MINITER', 3))
    for short, default in (('rtol', 1e-6), ('atol', 1e-16)):
        compare('pressure_' + short, number(config.get(short, default)),
                number(arguments.get('--pressure-linear-' + short)))
    for label in ('normalize_matrix', 'diagonal_scaling'):
        value = config.get(label, False)
        compare(label, value if type(value) is bool else UNKNOWN, False)
    compare('zero_initial_guess', True, True)

    for label, (key, variable, default) in AMG.items():
        reference = number(precond.get(key), integer=isinstance(default, int)) if key else UNKNOWN
        compare('amg_' + label, reference, environment_number(environment, variable, default))
    for key in AMG_LIBRARY:
        compare('amg_' + key, UNKNOWN, UNKNOWN)
    # The reference overwrites these two preconditioner options after parsing.
    compare('amg_max_iterations', 1, 1)
    compare('amg_tolerance', 0., 0.)
    if preconditioner != 'boomeramg':
        for label in status:
            if label.startswith('amg_'):
                status[label] = 'not_applicable'

    supported_amg = {key for key, _, _ in AMG.values() if key} | set(AMG_LIBRARY)
    supported_amg |= {'type', 'printlevel', 'logging', 'maxiter', 'tol'}
    ignored = bool(set(solver) - {'type', 'printlevel', 'logging', 'kdim'})
    ignored |= bool(set(reference_options) - {'solver', 'precond'})
    if preconditioner == 'boomeramg':
        ignored |= bool(set(precond) - supported_amg)
    return dict(
        known_matches=sorted(k for k, value in status.items() if value == 'equal'),
        known_differences=sorted(k for k, value in status.items() if value == 'different'),
        unresolved_settings=sorted(k for k, value in status.items() if value == 'unresolved'),
        amg_comparison_applicable=preconditioner == 'boomeramg',
        reference_ignored_solver_option_keys_present=ignored,
        reference_min_iterations_forwarded=False,
        reference_preconditioner_modelled=preconditioner != 'mgr')


def saved_configuration(pair):
    record = startup.read_json(pair / 'pair.json')
    require(record['schema'] == startup.SCHEMA)
    require(record['steps'] == (1 if record.get('first_step_audit', False) else startup.STEPS))
    require(record['case_sha256'] == digest(pair / 'case.json'))
    deck_hash = digest(pair / 'reference/input.i')
    require(record['deck_sha256'] == deck_hash)
    case = startup.read_json(pair / 'case.json')
    require(case['format'] == 'mars-simple-deck-v1' and case['deck_sha256'] == deck_hash)
    require(case.get('pressure_linear_policy') == 'reference')
    arguments = options(startup.solver_arguments(pair, 'mars'))
    require('--pressure-linear-rtol' in arguments and '--pressure-linear-atol' in arguments)
    deck = load_deck((pair / 'reference/input.i').read_bytes())
    return resolved_solver(deck, 'pressure_correction'), arguments


def saved_launch(pair, solver):
    directory = pair / ('reference' if solver == 'openaccel' else 'mars')
    start = startup.read_json(directory / 'launch-start.json')
    record = startup.read_json(directory / 'launch.json')
    require(record['schema'] == startup.SCHEMA and record['solver'] == solver)
    require(record['status'] in ('finished', 'failed') and start['status'] == 'started')
    keys = ('schema', 'solver', 'ranks', 'command', 'executable', 'executable_sha256',
            'pair_sha256', 'libraries', 'environment')
    require(all(record[key] == start[key] for key in keys))
    require(type(record['ranks']) is int and 1 <= record['ranks'] <= 4)
    require(record['pair_sha256'] == digest(pair / 'pair.json'))
    require(isinstance(record['executable'], str) and record['executable'])
    suffix = [record['executable']] + startup.solver_arguments(pair, solver)
    require(isinstance(record['command'], list) and len(record['command']) > len(suffix))
    require(all(isinstance(item, str) for item in record['command']))
    require(record['command'][-len(suffix):] == suffix)
    require(isinstance(record['libraries'], dict) and record['libraries'])
    require(all(isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value)
                for value in [record['executable_sha256']] + list(record['libraries'].values())))
    require(isinstance(record['environment'], dict) and
            all(isinstance(k, str) and isinstance(v, str) for k, v in record['environment'].items()))
    require(type(record['exit_code']) is int and
            record['exit_code'] == int((directory / 'run.exit').read_text()))
    accepted = record['exit_code'] in ((0,) if solver == 'openaccel' else (0, 2))
    require(record['status'] == ('finished' if accepted else 'failed'))
    for name in ('run.log', 'run.exit'):
        require(record['files'][name] == digest(directory / name))
    if solver == 'openaccel' and record['status'] == 'finished':
        require(record['files']['input.i'] == digest(directory / 'input.i'))
    return record


def inspect(pair):
    public = dict(schema=SCHEMA, comparison_status='invalid_evidence', failed_check='saved_configuration',
                  scope='source_calls_from_saved_configuration_and_launch_metadata',
                  binary_source_identity_verified=False, effective_library_defaults_verified=False,
                  launcher_child_environment_verified=False,
                  field_provenance_verified=False, identical_linear_solvers_verified=False)
    try:
        pair = pair.resolve()
        config, arguments = saved_configuration(pair)
        public['failed_check'] = 'saved_launch_metadata'
        launches = {solver: saved_launch(pair, solver) for solver in ('openaccel', 'mars')}
        public['failed_check'] = 'experimental_recovery'
        require(arguments.get('--pressure-refinement', '0') == '0')
        with (pair / 'mars/run.log').open(encoding='utf-8', errors='replace') as log:
            modes = []
            for line in log:
                modes.extend(re.findall(r'(?:^|\s)pressure_refinement=([^\s]*)', line))
                require('[simple-pressure-refinement]' not in line)
            require(modes in ([], ['0']))
        public['failed_check'] = 'pressure_configuration'
        result = compare_settings(config, arguments, launches['mars']['environment'])
        public.update(result)
        public.update(comparison_status='completed', failed_check='none',
                      assessment='source_settings_differ' if result['known_differences'] else
                                 'effective_settings_not_fully_resolved',
                      saved_launch_metadata_consistent=True,
                      same_recorded_rank_count=launches['mars']['ranks'] == launches['openaccel']['ranks'],
                      reference_runtime_convergence_verified=False,
                      acceptance_checks_identical=False)
    except Exception:
        # Exceptions can contain private paths, YAML fragments or option values.
        pass
    return public


def main(argv=None):
    parser = SafeParser(description=__doc__)
    parser.add_argument('--pair', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    result = inspect(args.pair)
    try:
        startup.write_json(args.output, result)
    except Exception:
        print('ERROR: cannot create a new public summary; existing files were not replaced.', file=sys.stderr)
        return 1
    print('Pressure settings inspected. No solver launched; share only the public JSON.')
    return 0 if result['comparison_status'] == 'completed' else 1


if __name__ == '__main__':
    sys.exit(main())
