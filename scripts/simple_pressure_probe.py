#!/usr/bin/env python3
"""Prepare a pressure-only accuracy experiment; keep case data and detailed results local."""

import copy
import json
import os
from pathlib import Path
import sys

import simple_startup_probe as startup
from prepare_simple_deck import load_deck, translate
from simple_public_diagnostics import SafeParser, summarize
from simple_snapshot_compare import EvidenceError, digest, options, require

SCHEMA = 'mars-simple-pressure-probe-v1'
ERRORS = frozenset(('baseline_identity', 'baseline_first_step_required', 'baseline_reference_policy_required',
    'pressure_configuration', 'pressure_target_unchanged', 'momentum_configuration_changed',
    'nonpressure_controls_changed', 'experiment_identity', 'executable_changed', 'libraries_changed',
    'solver_environment_changed', 'launcher_changed', 'baseline_stage_mismatch', 'launcher_exit',
    'environment_restore_requires_mars'))

# Only scalar solver/execution controls can be restored. Never restore paths or device binding.
RESTORABLE_ENVIRONMENT = frozenset(('MARS_AMG_AGG', 'MARS_AMG_COARSEN', 'MARS_AMG_INTERP',
    'MARS_AMG_PMAX', 'MARS_AMG_RELAX', 'MARS_AMG_RELAXORDER', 'MARS_AMG_STRONG', 'MARS_AMG_SWEEPS',
    'MARS_HYPRE_ABSTOL', 'MARS_HYPRE_FLEXGMRES', 'MARS_HYPRE_MAXX_RATIO', 'MARS_HYPRE_MINITER',
    'MARS_HYPRE_NULLX_RATIO', 'MARS_HYPRE_RESIDUAL_AUDIT', 'MARS_HYPRE_SPMV_VENDOR', 'MARS_HYPRE_VERBOSE',
    'MARS_HYPRE_SPGEMM_VENDOR', 'MARS_HYPRE_GPU_AWARE', 'MARS_HYPRE_POOL_MAX_MIB', 'MARS_HYPRE_POOL_CACHE_MIB',
    'MARS_OWNERSHIP', 'MARS_HALO_FACTOR', 'MARS_NODEHALO_HOST', 'MARS_NODEHALO_VALIDATE',
    'MARS_SIMPLE_PRESSURE_AUDIT', 'CUDA_LAUNCH_BLOCKING', 'CUDA_DEVICE_MAX_CONNECTIONS',
    'MPICH_GPU_SUPPORT_ENABLED', 'OMP_NUM_THREADS', 'OMP_PROC_BIND', 'OMP_PLACES'))
PUBLIC_ENVIRONMENT_NAMES = RESTORABLE_ENVIRONMENT | frozenset(('CUDA_VISIBLE_DEVICES',
    'CUDA_HOME', 'CUDA_PATH', 'CUDA_ROOT', 'CUDA_VERSION', 'CUDA_MODULE_LOADING',
    'MARS_SIGNAL_DIR', 'MARS_VERBOSE_MESH', 'MARS_HALO_DEBUG', 'MARS_SS_RESOLVE_DEBUG',
    'MARS_SYNC_TRACE', 'MARS_BLOCK_NODE_IDENTITY', 'PETSC_OPTIONS', 'PETSC_OPTIONS_YAML'))


def resolved_solver(deck, equation):
    solver = deck['simulation']['solver']
    settings = solver['solver_control']['advanced_options']['linear_solver_settings']
    config = next((settings[k] for k in (equation, 'segregated_flow', 'default') if k in settings), None)
    require(isinstance(config, dict), 'pressure_configuration')
    if 'lookup' in config:
        config = solver[config['lookup']]
    require(isinstance(config, dict) and 'lookup' not in config, 'pressure_configuration')
    return copy.deepcopy(config)


def tightened_deck(original):
    mapped, _ = translate(original, 'reference')
    old = options(mapped)
    modified = copy.deepcopy(original)
    config = resolved_solver(original, 'pressure_correction')
    rtol = min(float(old['--pressure-linear-rtol']), 1e-10)
    require(rtol < float(old['--pressure-linear-rtol']) or float(old['--pressure-linear-atol']) > 0,
            'pressure_target_unchanged')
    config.update(rtol=rtol, atol=0.)
    settings = modified['simulation']['solver']['solver_control']['advanced_options']['linear_solver_settings']
    # A fallback or named block may also configure momentum. Give pressure its own copy.
    settings['pressure_correction'] = config
    require(resolved_solver(original, 'coupled_navier_stokes') == resolved_solver(modified, 'coupled_navier_stokes'),
            'momentum_configuration_changed')
    translated, _ = translate(modified, 'reference')
    filtered = lambda values: {k: v for k, v in options(values).items()
                               if k not in ('--pressure-linear-rtol', '--pressure-linear-atol')}
    require(filtered(mapped) == filtered(translated), 'nonpressure_controls_changed')
    return modified, translated


def baseline_records(baseline):
    record = startup.pair_inputs(baseline)
    require(record.get('first_step_audit') is True and record['steps'] == 1, 'baseline_first_step_required')
    case = startup.read_json(baseline / 'case.json')
    require(case.get('pressure_linear_policy') == 'reference', 'baseline_reference_policy_required')
    launches = {solver: startup.verified_launch(baseline, solver) for solver in ('openaccel', 'mars')}
    return record, case, launches


def experiment_inputs(baseline):
    baseline = baseline.resolve()
    old_record, old_case, launches = baseline_records(baseline)
    deck, mapped = tightened_deck(load_deck((baseline / 'reference/input.i').read_bytes()))
    values = options(old_case['arguments'])
    args = ['--mesh', old_record['mesh'], '--mesh-format', 'exodus'] + mapped + [
        '--reference-length', values['--reference-length']]
    case = dict(old_case, arguments=args, pressure_linear_policy='reference')
    identity = dict(schema=SCHEMA, baseline_pair=str(baseline), baseline_pair_sha256=digest(baseline / 'pair.json'),
        baseline_launch_hashes={solver: digest(baseline / ('reference' if solver == 'openaccel' else 'mars') / 'launch.json')
                               for solver in launches})
    return old_record, case, deck, identity, launches


def prepare(baseline, output):
    import yaml
    record, case, deck, identity, _ = experiment_inputs(baseline)
    output.mkdir(mode=0o700)
    (output / 'reference').mkdir(mode=0o700)
    path = output / 'reference/input.i'
    path.write_text(yaml.safe_dump(deck, default_flow_style=False))
    case['deck_sha256'] = digest(path)
    startup.write_json(output / 'case.json', case)
    record = dict(record, deck_sha256=digest(path), case_sha256=digest(output / 'case.json'),
                  pressure_experiment=identity, changed_controls=['pressure_linear_rtol', 'pressure_linear_atol'])
    startup.write_json(output / 'pair.json', record)
    check_inputs(output)


def check_inputs(pair):
    inputs = startup.pair_inputs(pair)
    identity = inputs['pressure_experiment']
    require(identity['schema'] == SCHEMA, 'experiment_identity')
    baseline = Path(identity['baseline_pair'])
    old, case, deck, expected_identity, launches = experiment_inputs(baseline)
    require(identity == expected_identity, 'baseline_identity')
    require(load_deck((pair / 'reference/input.i').read_bytes()) == deck, 'nonpressure_controls_changed')
    case['deck_sha256'] = digest(pair / 'reference/input.i')
    require(startup.read_json(pair / 'case.json') == case, 'experiment_identity')
    expected = dict(old, deck_sha256=case['deck_sha256'], case_sha256=digest(pair / 'case.json'),
                    pressure_experiment=identity, changed_controls=['pressure_linear_rtol', 'pressure_linear_atol'])
    require(inputs == expected, 'experiment_identity')
    return baseline, launches


def solver_environment(environment):
    # Scheduler IDs and paths change between launches; solver-affecting overrides must not.
    return {k: v for k, v in environment.items() if k.startswith(('MARS_', 'HYPRE_', 'CUDA_', 'MPICH_', 'OMP_'))
            or k in ('PETSC_OPTIONS', 'PETSC_OPTIONS_YAML')}


def environment_changes(current, recorded):
    current, recorded = solver_environment(current), solver_environment(recorded)
    return {key for key in set(current) | set(recorded) if current.get(key) != recorded.get(key)}


def environment_report(changes):
    # Unknown names can carry private identifiers; export only literal labels from this module.
    return dict(changed_known_names=sorted(changes & PUBLIC_ENVIRONMENT_NAMES),
                other_changed_names_present=bool(changes - PUBLIC_ENVIRONMENT_NAMES),
                restorable_changes_present=bool(changes & RESTORABLE_ENVIRONMENT),
                protected_or_unknown_changes_present=bool(changes - RESTORABLE_ENVIRONMENT))


def restore_environment(current, recorded):
    result = dict(current)
    for key in RESTORABLE_ENVIRONMENT:
        if key in recorded:
            result[key] = recorded[key]
        else:
            result.pop(key, None)
    return result


def launcher(record, pair, solver):
    suffix = [record['executable']] + startup.solver_arguments(pair.resolve(), solver)
    require(record['command'][-len(suffix):] == suffix, 'launcher_changed')
    prefix = record['command'][:-len(suffix)]
    require(prefix, 'launcher_changed')
    return prefix


def runtime_matches(baseline, record, pair, solver):
    for key, label in (('executable', 'executable_changed'), ('executable_sha256', 'executable_changed'),
                       ('libraries', 'libraries_changed'), ('ranks', 'launcher_changed')):
        require(record[key] == baseline[key], label)
    require(solver_environment(record['environment']) == solver_environment(baseline['environment']),
            'solver_environment_changed')
    old_pair = Path(startup.read_json(pair / 'pair.json')['pressure_experiment']['baseline_pair'])
    require(launcher(record, pair, solver) == launcher(baseline, old_pair, solver), 'launcher_changed')


def run(pair, solver, restore_solver_environment=False, public=None):
    baseline, records = check_inputs(pair)
    old = records[solver]
    executable = Path(old['executable'])
    require(executable.is_file() and digest(executable) == old['executable_sha256'], 'executable_changed')
    require(startup.runtime_libraries(executable) == old['libraries'], 'libraries_changed')
    environment = dict(os.environ)
    for key in ('MARS_OPENACCEL_EXPORT_DIR', 'MARS_OPENACCEL_PUBLIC_FIXTURE'):
        environment.pop(key, None)
    if solver == 'openaccel':
        environment.update(OMP_NUM_THREADS='1', OMP_PROC_BIND='close', OMP_PLACES='cores')
    changed = environment_changes(environment, old['environment'])
    if public is not None:
        public['environment_check'] = environment_report(changed)
    if restore_solver_environment:
        require(solver == 'mars', 'environment_restore_requires_mars')
        environment = restore_environment(environment, old['environment'])
        if public is not None:
            public['restored_solver_controls'] = sorted(changed & RESTORABLE_ENVIRONMENT)
    require(solver_environment(environment) == solver_environment(old['environment']), 'solver_environment_changed')
    startup.launch(pair, solver, executable, old['ranks'], launcher(old, baseline, solver), environment=environment)
    runtime_matches(old, startup.verified_launch(pair, solver), pair, solver)


def compare(pair, public, details):
    baseline, records = check_inputs(pair)
    for solver, old in records.items():
        runtime_matches(old, startup.verified_launch(pair, solver), pair, solver)
    old_details, new_details = details / 'baseline', details / 'tightened'
    old_details.mkdir(mode=0o700)
    new_details.mkdir(mode=0o700)
    old_public, new_public = {}, {}
    startup.compare(baseline, old_public, old_details)
    startup.compare(pair, new_public, new_details)
    stages = ('momentum_matrix', 'momentum_rhs', 'momentum_predictor', 'momentum_influence',
              'pressure_matrix', 'pressure_rhs')
    require(all(old_public['first_step_stage_matches'][key] for key in stages), 'baseline_stage_mismatch')
    checks = new_public['pressure_solve_checks']
    targets = all(checks[key] for key in ('same_declared_tolerances', 'mars_referenced_copies_equal_owners',
        'mars_local_residual_meets_runtime_limit', 'mars_owner_residual_meets_runtime_limit',
        'mars_residual_below_declared_unpreconditioned_limit', 'reference_residual_below_declared_unpreconditioned_limit'))
    matched = all(new_public['first_step_stage_matches'].values())
    upstream = all(new_public['first_step_stage_matches'][key] for key in stages)
    outcome = ('upstream_stage_mismatch' if not upstream else
               'pressure_target_not_met' if not targets else
               'first_step_matches' if matched else 'first_step_still_differs')
    public.update(comparison_status='completed', failed_check='none', outcome=outcome,
        pressure_only_controls_verified=True, captured_runtime_settings_match=True,
        runtime_scope='recorded_binaries_libraries_launchers_and_solver_environment',
        baseline_first_differing_stage=old_public['first_differing_stage'],
        tighter_pressure_residuals_pass=targets, tightened_first_step=new_public,
        root_cause_proven=False, nonlinear_convergence_verified=False)


def failure_summary(pair, solver):
    directory = pair / ('reference' if solver == 'openaccel' else 'mars')
    if not (directory / 'run.log').is_file() or not (directory / 'run.exit').is_file():
        return None
    if solver == 'mars':
        with (directory / 'run.log').open(errors='replace') as stream:
            return summarize(stream, (directory / 'run.exit').read_text())
    executable = Path(startup.read_json(directory / 'launch-start.json')['executable'])
    return startup.inspect_launch(pair, solver, executable)


def main(argv=None):
    parser = SafeParser(description=__doc__)
    sub = parser.add_subparsers(dest='action')
    prep = sub.add_parser('prepare')
    prep.add_argument('--baseline-pair', type=Path, required=True)
    prep.add_argument('--output-dir', type=Path, required=True)
    launch = sub.add_parser('run')
    launch.add_argument('--pair', type=Path, required=True)
    launch.add_argument('--solver', choices=('openaccel', 'mars'), required=True)
    launch.add_argument('--output', type=Path, required=True, help='New public launch status JSON')
    launch.add_argument('--restore-solver-environment', action='store_true',
                        help='MARS only: restore allowlisted scalar controls from the verified baseline in the child process')
    check = sub.add_parser('compare')
    check.add_argument('--pair', type=Path, required=True)
    check.add_argument('--output', type=Path, required=True)
    check.add_argument('--detail-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    os.umask(0o077)
    public = dict(schema=SCHEMA, comparison_status='invalid_evidence', failed_check='private_preflight_or_capture')
    try:
        if args.action == 'prepare':
            prepare(args.baseline_pair.resolve(), args.output_dir.resolve())
            print('Pressure-only experiment prepared. Baseline preserved; no solver launched.')
            return 0
        if args.action not in ('run', 'compare'):
            parser.error('subcommand required')
        # Reserve the report first so a repeated command cannot launch another job.
        with args.output.open('x') as stream:
            try:
                pair = args.pair.resolve()
                if args.action == 'run':
                    run(pair, args.solver, args.restore_solver_environment, public)
                    public.update(comparison_status='capture_complete', failed_check='none')
                else:
                    args.detail_dir.mkdir(mode=0o700)
                    compare(pair, public, args.detail_dir)
            except Exception as error:
                public['failed_check'] = str(error) if isinstance(error, EvidenceError) and str(error) in ERRORS else 'private_preflight_or_capture'
                if args.action == 'run':
                    public['outcome'] = 'inconclusive'
                    try:
                        public['launch_diagnostics'] = failure_summary(pair, args.solver)
                    except Exception:
                        public['launch_diagnostics'] = None
            json.dump(public, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')
        print('Pressure experiment status written. Share only the public JSON.')
        return 0 if public['comparison_status'] in ('completed', 'capture_complete') else 1
    except Exception as error:
        label = str(error) if isinstance(error, EvidenceError) and str(error) in ERRORS else 'private_preflight_or_capture'
        print('ERROR: pressure experiment failed (' + label + ').', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
