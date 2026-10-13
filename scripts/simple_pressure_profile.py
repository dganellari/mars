#!/usr/bin/env python3
"""Compare original-control startup fields with an explicit GPU pressure configuration."""
import copy
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import sys

import simple_hypre_defaults as defaults
import simple_pressure_probe as pressure
import simple_pressure_settings as settings
from simple_pressure_replay import reference_controls
import simple_startup_probe as startup
from simple_snapshot_compare import EvidenceError, band, digest, options, require
from simple_public_diagnostics import SafeParser

SCHEMA = 'mars-simple-pressure-profile-v1'
HISTORY = 'original_twenty_steps'
DEFAULT_KEYS = dict(coarsentype='coarsen_type', interptype='interp_type', relaxorder='relax_order',
    pmax='p_max_elmts', maxlevels='max_levels', mincoarsesize='min_coarse_size', maxcoarsesize='max_coarse_size',
    numfunctions='num_functions', coarsencutfactor='coarsen_cut_factor', cycletype='cycle_type', fcycle='fcycle',
    strongthreshold='strong_threshold', truncfactor='trunc_factor', jacobitruncthreshold='jacobi_trunc_threshold',
    maxrowsum='max_row_sum', aggnumlevels='agg_num_levels', agginterptype='agg_interp_type',
    aggtruncfactor='agg_trunc_factor', numpaths='num_paths', keeptranspose='keep_transpose', nodal='nodal',
    nodaldiag='nodal_diag', coarserelax='relax_coarse', relax_down='relax_down', relax_up='relax_up',
    sweeps_1='sweeps_down', sweeps_2='sweeps_up', sweeps_3='sweeps_coarse')
REAL_KEYS = {'rtol', 'atol', 'strongthreshold', 'truncfactor', 'jacobitruncthreshold', 'maxrowsum', 'aggtruncfactor'}
KEYS = set(DEFAULT_KEYS) | {'method', 'kdim', 'miniter', 'maxiter', 'rtol', 'atol', 'relaxtype', 'numsweeps'}
DEFAULT_FILES = ('public.json', 'probe.log', 'private-provenance.json', 'private-probe-libraries.json',
                 'compile.exit', 'probe.exit')
ERRORS = {'original_baseline_required', 'defaults_evidence', 'unsupported_gpu_profile', 'profile_identity',
          'runtime_profile_mismatch', 'executable_changed', 'libraries_changed', 'solver_environment_changed',
          'launcher_changed', 'launcher_exit', 'runtime_recovery_mismatch', 'initial_fields_mismatch',
          'reference_pressure_log', 'history_required', 'unsupported_momentum_profile'}


def profile_violations(v):
    if set(v) != KEYS:
        return ['profile_keys']
    invalid = []
    for k, value in v.items():
        if (type(value) not in (int, float) or not math.isfinite(value)
                or (k not in REAL_KEYS and not (value == int(value) and 0 <= value < 2**31))):
            invalid.append(k)
    if invalid:
        return sorted(invalid)
    checks = dict(method=v['method'] in (0, 1), kdim=v['kdim'] > 0, maxiter=v['maxiter'] > 0,
        miniter=v['miniter'] <= v['maxiter'], rtol=0 < v['rtol'] < 1, atol=v['atol'] >= 0,
        coarsentype=v['coarsentype'] == 8, interptype=v['interptype'] == 6,
        numfunctions=v['numfunctions'] == 1, cycletype=v['cycletype'] == 1,
        keeptranspose=v['keeptranspose'] in (0, 1), maxlevels=v['maxlevels'] > 0,
        maxcoarsesize=v['maxcoarsesize'] > 0, mincoarsesize=v['mincoarsesize'] <= v['maxcoarsesize'],
        numpaths=v['numpaths'] > 0, coarserelax=v['coarserelax'] == 18)
    checks.update((k, v[k] == 0) for k in ('relaxorder', 'aggnumlevels', 'nodal', 'nodaldiag', 'fcycle'))
    checks.update((k, v[k] in (6, 13, 14, 18)) for k in ('relaxtype', 'relax_down', 'relax_up'))
    checks.update((k, v[k] > 0) for k in ('numsweeps', 'sweeps_1', 'sweeps_2', 'sweeps_3'))
    checks.update((k, 0 <= v[k] <= 1) for k in REAL_KEYS - {'rtol', 'atol'})
    return sorted(k for k, passed in checks.items() if not passed)


def validate(v):
    require(not profile_violations(v), 'unsupported_gpu_profile')


def read_defaults(directory, reference):
    # Recheck saved probe evidence without requiring both uenv mounts to be visible.
    public = startup.read_json(directory / 'public.json')
    raw = defaults.raw_output(directory / 'probe.log')
    projection = defaults.public_defaults(raw)
    metadata = startup.read_json(directory / 'private-provenance.json')
    libraries = startup.read_json(directory / 'private-probe-libraries.json')
    require(public.get('schema') == defaults.SCHEMA and public.get('comparison_status') == 'completed'
        and public.get('failed_check') == 'none' and all(public.get(k) == v for k, v in projection.items())
        and all(public.get(k) is True for k in ('saved_capture_identity_verified', 'captured_libraries_unchanged',
            'loaded_hypre_library_matches_capture', 'mpi_library_identity_verified', 'mpi_enabled'))
        and all((directory / k).read_text().strip() == '0' for k in ('compile.exit', 'probe.exit')),
        'defaults_evidence')
    captured = metadata['capture']
    require(all(captured[k] == reference[k] for k in ('executable', 'executable_sha256', 'libraries')),
            'defaults_evidence')
    try:
        defaults.probe_dependencies(libraries, reference['libraries'], True)
    except ValueError:
        raise EvidenceError('defaults_evidence')
    for kind, path in (('hypre', raw['library_path']), ('mpi', raw['mpi_library_path'])):
        require(Path(path).is_absolute() and defaults.library_kind(path) == kind, 'defaults_evidence')
        # dladdr and ldd may spell the same library through different symlinks.
        # The successful saved probe checked the actual loaded path by content.
        if path in libraries:
            require(libraries[path] in {h for p, h in reference['libraries'].items()
                                       if defaults.library_kind(p) == kind}, 'defaults_evidence')
    # The omitted one-level relaxation/sweep fallback below is inspected in this version.
    require(projection['version'] == '3.1.0' and projection['build_cuda'] is False, 'defaults_evidence')
    return projection['defaults']


def resolved_values(config, target, measured):
    explicit = reference_controls(config, target)
    v = {key: measured['amg_' + name] for key, name in DEFAULT_KEYS.items()}
    method = 'flexgmres' if explicit['method'] else 'gmres'
    v.update(kdim=measured[method + '_restart_dimension'], miniter=measured[method + '_minimum_iterations'],
             relaxtype=6, numsweeps=1)
    v.update(explicit)
    # Hypre generic setters also overwrite the cycle settings, in this order.
    if 'relaxtype' in explicit:
        v.update(relax_down=v['relaxtype'], relax_up=v['relaxtype'], coarserelax=9)
    if 'numsweeps' in explicit:
        v.update(sweeps_1=v['numsweeps'], sweeps_2=v['numsweeps'], sweeps_3=1)
    changes = []
    if v['coarsentype'] == 10:
        v['coarsentype'] = 8
        changes.append('hmis_to_device_pmis')
    if v['coarserelax'] == 9:
        v['coarserelax'] = 18
        changes.append('coarse_direct_to_device_l1_jacobi')
    return v, changes


def resolve(config, target, measured):
    values, changes = resolved_values(config, target, measured)
    validate(values)
    return values, changes


def inspect_configuration(config, measured, stage):
    from prepare_simple_deck import reference_linear_target

    def label(value, allowed):
        value = value.lower() if isinstance(value, str) else ''
        return value if value in allowed else 'other_or_invalid'

    report = dict(mapping_supported=False, failed_component='definition')
    if not isinstance(config, dict):
        return report
    outer = config.get('options', {})
    report.update(family=label(config.get('family'), {'hypre', 'petsc', 'trilinos', 'amgsolver', 'gmres'}),
        normalize_matrix_requested=config.get('normalize_matrix', False) is not False,
        diagonal_scaling_requested=config.get('diagonal_scaling', False) is not False)
    report['failed_component'] = 'options_structure'
    if not isinstance(outer, dict):
        return report
    try:
        solver = settings.lower_options(outer.get('solver', {}))
        precond = settings.lower_options(outer.get('precond', {}))
    except ValueError:
        return report
    supported_amg = {key for key, _, _ in settings.AMG.values() if key} | set(settings.AMG_LIBRARY)
    report.update(solver=label(solver.get('type', 'gmres'), {'gmres', 'flexgmres', 'boomeramg', 'mgr'}),
        preconditioner=label(precond.get('type', 'none'), {'none', 'boomeramg', 'mgr'}),
        unrecognized_option_keys_present=bool(set(outer) - {'solver', 'precond'}
            or set(solver) - {'type', 'kdim', 'printlevel', 'logging'}
            or set(precond) - supported_amg - {'type', 'printlevel', 'logging', 'maxiter', 'tol'}))
    for component, passed in (
            ('family', report['family'] == 'hypre'),
            ('normalize_matrix', not report['normalize_matrix_requested']),
            ('diagonal_scaling', not report['diagonal_scaling_requested']),
            ('solver', report['solver'] in ('gmres', 'flexgmres')),
            ('preconditioner', report['preconditioner'] == 'boomeramg'),
            ('option_keys', not report['unrecognized_option_keys_present'])):
        if not passed:
            report['failed_component'] = component
            return report
    report['failed_component'] = 'residual_target_or_solver_definition'
    try:
        arguments = options(reference_linear_target({'equation': config}, {}, 'equation', stage))
        target = tuple(float(arguments['--' + stage + '-linear-' + key]) for key in ('atol', 'rtol'))
    except ValueError:
        return report
    report['failed_component'] = 'reference_options'
    try:
        values, changes = resolved_values(config, target, measured)
    except ValueError:
        return report
    report['unsupported_gpu_controls'] = profile_violations(values)
    report['gpu_adaptations'] = changes
    report['mapping_supported'] = not report['unsupported_gpu_controls']
    report['failed_component'] = 'none' if report['mapping_supported'] else 'gpu_controls'
    return report


def inspect_profiles(baseline, probe, public):
    before = pressure.baseline_records(baseline)
    measured = read_defaults(probe, before[2]['openaccel'])
    deck = startup.load_deck((baseline / 'reference/input.i').read_bytes())
    reports = {}
    for stage, equation in (('momentum', 'coupled_navier_stokes'), ('pressure', 'pressure_correction')):
        try:
            config = pressure.resolved_solver(deck, equation)
        except (ValueError, KeyError, TypeError):
            reports[stage] = dict(mapping_supported=False, failed_component='definition')
        else:
            reports[stage] = inspect_configuration(config, measured, stage)
    require(pressure.baseline_records(baseline) == before
            and read_defaults(probe, before[2]['openaccel']) == measured, 'profile_identity')
    public.update(comparison_status='inspection_complete', failed_check='none', profiles=reports,
        both_profiles_supported=all(r['mapping_supported'] for r in reports.values()),
        solver_launched=False, scope='saved_configuration_mapping_not_runtime_or_field_parity')


def text_profile(values):
    return ''.join('{} {:.17g}\n'.format(k, values[k]) for k in sorted(values))


def baseline_data(baseline, probe):
    record, case, launches = pressure.baseline_records(baseline)
    require('pressure_experiment' not in record and 'pressure_profile' not in record
            and launches['mars']['ranks'] == launches['openaccel']['ranks'], 'original_baseline_required')
    require(launches['mars']['environment'].get('MARS_HYPRE_SPMV_VENDOR', '0') == '0', 'unsupported_gpu_profile')
    config, arguments = settings.saved_configuration(baseline)
    target = (float(arguments['--pressure-linear-atol']), float(arguments['--pressure-linear-rtol']))
    values, adaptations = resolve(config, target, read_defaults(probe, launches['openaccel']))
    return record, case, launches, values, adaptations


def momentum_data(baseline, probe, reference):
    try:
        deck = startup.load_deck((baseline / 'reference/input.i').read_bytes())
        config = pressure.resolved_solver(deck, 'coupled_navier_stokes')
        from prepare_simple_deck import reference_momentum_target
        solver = deck['simulation']['solver']
        arguments = reference_momentum_target(solver['solver_control']['advanced_options']['linear_solver_settings'], solver)
        target = options(arguments)
        values, adaptations = resolve(config, (float(target['--momentum-linear-atol']),
                                              float(target['--momentum-linear-rtol'])), read_defaults(probe, reference))
        return values, adaptations, arguments
    except Exception:
        # Never print a private solver name or approximate a different backend.
        raise EvidenceError('unsupported_momentum_profile')


def pair_content(baseline, pair, record, case, identity):
    case = copy.deepcopy(case)
    record = dict(record)
    if identity.get('kind') == HISTORY:
        original = startup.load_deck((baseline / 'reference/input.i').read_bytes())
        require(startup.load_deck((pair / 'reference/input.i').read_bytes()) == pressure.history_deck(original),
                'profile_identity')
        case['deck_sha256'] = digest(pair / 'reference/input.i')
        record.update(steps=startup.STEPS, first_step_audit=False, deck_sha256=case['deck_sha256'],
                      changed_controls=['iteration_limits', 'output_control', 'linear_system_output', 'gpu_pressure_profile'])
    case['arguments'] += ['--pressure-solver-profile', str(pair / 'pressure.profile')]
    case['pressure_solver_profile_sha256'] = identity['sha256']
    if 'momentum' in identity:
        _, _, arguments = momentum_data(baseline, pair / 'defaults', startup.verified_launch(baseline, 'openaccel'))
        case['arguments'] += arguments + ['--momentum-solver-profile', str(pair / 'momentum.profile')]
        case['momentum_linear_policy'] = 'reference'
        case['momentum_solver_profile_sha256'] = identity['momentum']['sha256']
        record['changed_controls'] = record['changed_controls'] + ['gpu_momentum_profile']
        case['notes'] = [note for note in case.get('notes', []) if not note.startswith((
            'MARS retains its momentum solver settings', 'Reference pressure residual target copied;'))] + [
            'Both linear targets are copied; explicit GPU profiles record backend adaptations.',
            'Nonlinear MARS convergence norms are unchanged; this compares saved iteration states.']
    return record, case


def prepare(baseline, probe, executable, output, history=False, both_linear_solvers=False):
    record, case, launches, values, adaptations = baseline_data(baseline, probe)
    momentum = momentum_data(baseline, probe, launches['openaccel']) if both_linear_solvers else None
    require(executable.is_file() and os.access(str(executable), os.X_OK), 'executable_changed')
    libraries = startup.runtime_libraries(executable)
    require(libraries == launches['mars']['libraries'], 'libraries_changed')
    output.mkdir(mode=0o700)
    (output / 'reference').mkdir(mode=0o700)
    if history:
        import yaml
        original = startup.load_deck((baseline / 'reference/input.i').read_bytes())
        (output / 'reference/input.i').write_text(yaml.safe_dump(pressure.history_deck(original), default_flow_style=False))
    else:
        shutil.copyfile(baseline / 'reference/input.i', output / 'reference/input.i')
    (output / 'defaults').mkdir(mode=0o700)
    for name in DEFAULT_FILES:
        shutil.copyfile(probe / name, output / 'defaults' / name)
    (output / 'pressure.profile').write_text(text_profile(values))
    identity = dict(schema=SCHEMA, reference_pair=str(baseline), reference_pair_sha256=digest(baseline / 'pair.json'),
        sha256=digest(output / 'pressure.profile'), adaptations=adaptations,
        defaults_hashes={name: digest(output / 'defaults' / name) for name in DEFAULT_FILES},
        baseline_launch_hashes={s: digest(baseline / ('reference' if s == 'openaccel' else 'mars') / 'launch.json') for s in launches},
        executable=str(executable), executable_sha256=digest(executable), libraries=libraries)
    if momentum is not None:
        (output / 'momentum.profile').write_text(text_profile(momentum[0]))
        identity['momentum'] = dict(sha256=digest(output / 'momentum.profile'), adaptations=momentum[1])
    if history:
        identity['kind'] = HISTORY
    record, case = pair_content(baseline, output, record, case, identity)
    startup.write_json(output / 'case.json', case)
    startup.write_json(output / 'pair.json', dict(record, case_sha256=digest(output / 'case.json'), pressure_profile=identity))
    check_inputs(output)


def check_inputs(pair):
    record = startup.pair_inputs(pair)
    identity = record['pressure_profile']
    require(identity['schema'] == SCHEMA and identity.get('kind') in (None, HISTORY), 'profile_identity')
    baseline = Path(identity['reference_pair'])
    require(digest(baseline / 'pair.json') == identity['reference_pair_sha256'], 'profile_identity')
    require(identity['defaults_hashes'] == {name: digest(pair / 'defaults' / name) for name in DEFAULT_FILES}, 'defaults_evidence')
    old, case, launches, values, adaptations = baseline_data(baseline, pair / 'defaults')
    require(identity['baseline_launch_hashes'] == {s: digest(baseline / ('reference' if s == 'openaccel' else 'mars') / 'launch.json')
                                                  for s in launches}, 'profile_identity')
    require(text_profile(values) == (pair / 'pressure.profile').read_text()
            and identity['adaptations'] == adaptations and identity['libraries'] == launches['mars']['libraries'], 'profile_identity')
    if 'momentum' in identity:
        momentum, changes, _ = momentum_data(baseline, pair / 'defaults', launches['openaccel'])
        require(text_profile(momentum) == (pair / 'momentum.profile').read_text()
                and identity['momentum'] == dict(sha256=digest(pair / 'momentum.profile'), adaptations=changes), 'profile_identity')
    old, expected_case = pair_content(baseline, pair, old, case, identity)
    require(startup.read_json(pair / 'case.json') == expected_case
            and record == dict(old, case_sha256=digest(pair / 'case.json'), pressure_profile=identity), 'profile_identity')
    return baseline, launches, values, identity


def check_launch(pair, baseline, old, identity, solver='mars'):
    record = startup.verified_launch(pair, solver)
    expected = identity if solver == 'mars' else old
    require(all(record[k] == expected[k] for k in ('executable', 'executable_sha256', 'libraries')), 'profile_identity')
    require(record['ranks'] == old['ranks'] and pressure.launcher(record, pair, solver) ==
            pressure.launcher(old, baseline, solver), 'launcher_changed')
    require(pressure.solver_environment(record['environment']) == pressure.solver_environment(old['environment']),
            'solver_environment_changed')
    return record


def run(pair, public, solver='mars'):
    baseline, launches, _, identity = check_inputs(pair)
    require(solver == 'mars' or (solver == 'openaccel' and identity.get('kind') == HISTORY), 'profile_identity')
    old = launches[solver]
    expected = identity if solver == 'mars' else old
    exe = Path(expected['executable'])
    require(exe.is_file() and os.access(str(exe), os.X_OK) and digest(exe) == expected['executable_sha256'], 'executable_changed')
    require(startup.runtime_libraries(exe) == expected['libraries'], 'libraries_changed')
    environment = dict(os.environ)
    for key in ('MARS_OPENACCEL_EXPORT_DIR', 'MARS_OPENACCEL_PUBLIC_FIXTURE'):
        environment.pop(key, None)
    if solver == 'openaccel':
        environment.update(OMP_NUM_THREADS='1', OMP_PROC_BIND='close', OMP_PLACES='cores')
    changes = pressure.environment_changes(environment, old['environment'])
    public['environment_check'] = pressure.environment_report(changes)
    if solver == 'mars':
        environment = pressure.restore_environment(environment, old['environment'])
    require(not pressure.environment_changes(environment, old['environment']), 'solver_environment_changed')
    launcher = pressure.launcher(old, baseline, solver)
    # Keep the complete direct command private; the wrapper also supplies the checked environment.
    with (pair / (solver + '-launch-command.txt')).open('x') as stream:
        stream.write(' '.join(shlex.quote(word) for word in launcher + [str(exe)] +
                              startup.solver_arguments(pair, solver)) + '\n')
    startup.launch(pair, solver, exe, old['ranks'], launcher, environment)
    check_launch(pair, baseline, old, identity, solver)


def check_runtime(pair, ranks, values, history=False, stage='pressure'):
    expected = dict(values)
    expected['effective_relax_1'] = expected.pop('relax_down')
    expected['effective_relax_2'] = expected.pop('relax_up')
    expected['effective_relax_3'] = expected['coarserelax']
    prefix = 'flow-' + stage + '-'
    paths = sorted((pair / 'mars').glob(prefix + 'rank-*.settings'))
    names = [(prefix + 'rank-{}.settings').format(rank) for rank in range(ranks)]
    if history:
        paths += sorted((pair / 'mars').glob(prefix + 'step-*-rank-*.settings'))
        names = [(prefix + 'step-{}-rank-{}.settings').format(step, rank)
                 for step in range(1, startup.STEPS + 1) for rank in range(ranks)]
    require({p.name for p in paths} == set(names), 'runtime_profile_mismatch')
    for name in names:
        lines = (pair / 'mars' / name).read_text().splitlines()
        parsed = [line.split() for line in lines]
        require(all(len(row) == 2 for row in parsed), 'runtime_profile_mismatch')
        actual = {k: float(v) for k, v in parsed}
        require(len(actual) == len(parsed) and set(actual) == set(expected) | {'hypre_release', 'effective_levels'}
                and all(math.isfinite(v) for v in actual.values()) and all(actual[k] == v for k, v in expected.items())
                and actual['effective_levels'] >= 1, 'runtime_profile_mismatch')


def compare(pair, public, details):
    baseline, launches, values, identity = check_inputs(pair)
    record = check_launch(pair, baseline, launches['mars'], identity)
    history = identity.get('kind') == HISTORY
    if history:
        check_launch(pair, baseline, launches['openaccel'], identity, 'openaccel')
        lines = (pair / 'mars/run.log').read_text(errors='replace').splitlines()
        require([line for line in lines if line.startswith('pressure_refinement=')] ==
                ['pressure_refinement=0 pressure_expansion=0'] and
                not any(line.startswith(('[simple-pressure-expansion]', '[simple-pressure-refinement]')) for line in lines),
                'runtime_recovery_mismatch')
    check_runtime(pair, record['ranks'], values, history)
    if 'momentum' in identity:
        momentum, _, _ = momentum_data(baseline, pair / 'defaults', launches['openaccel'])
        check_runtime(pair, record['ranks'], momentum, history, 'momentum')
    result = {}
    startup.compare(pair, result, details, reference_pair=None if history else baseline)
    check_inputs(pair)
    if history:
        require(result['initial_field_parity_verified'], 'initial_fields_mismatch')
        saved = startup.read_json(details / 'comparison-private.json')
        public.update(experiment=HISTORY, recovery_enabled=False,
            pressure_residuals_recomputed_at_every_step=False,
            outcome='twenty_step_fields_match' if result['all_twenty_snapshots_match'] else 'twenty_step_fields_differ',
            snapshot_agreement=[dict(iteration=row['iteration'], **{
                field + '_band': band(row['errors'][field]) for field in
                ('velocity_max_scaled', 'velocity_rms_scaled', 'pressure_max_scaled', 'pressure_rms_scaled')})
                for row in saved['snapshots']])
    public.update(comparison_status='completed', failed_check='none', comparison=result,
        original_pressure_targets_preserved=True, recorded_gpu_pressure_settings_verified=True,
        gpu_adaptations=identity['adaptations'], reference_solver_launched=history,
        identical_linear_solvers_verified=False, nonlinear_convergence_verified=False)
    if 'momentum' in identity:
        public.update(original_momentum_targets_preserved=True, recorded_gpu_momentum_settings_verified=True,
                      momentum_gpu_adaptations=identity['momentum']['adaptations'])


def inspect_linear(pair, public):
    from simple_pressure_log import inspect_pressure_log
    baseline, launches, values, identity = check_inputs(pair)
    require(identity.get('kind') == HISTORY, 'history_required')
    records = {solver: check_launch(pair, baseline, launches[solver], identity, solver)
               for solver in ('mars', 'openaccel')}
    check_runtime(pair, records['mars']['ranks'], values, history=True)
    checks = inspect_pressure_log((pair / 'reference/run.log').read_text(errors='replace'),
                                  startup.STEPS, values['rtol'], values['atol'], int(values['maxiter']))
    # Recheck hashed evidence after reading the log; this command never launches a solver.
    require(check_inputs(pair) == (baseline, launches, values, identity), 'profile_identity')
    for solver in records:
        require(check_launch(pair, baseline, launches[solver], identity, solver) == records[solver], 'profile_identity')
    public.update(comparison_status='completed', failed_check='none', experiment=HISTORY,
                  solver_launched=False, reference_pressure_log_checks=checks,
                  original_pressure_targets_preserved=True, nonlinear_convergence_verified=False,
                  openaccel_field_parity_verified=False)


def main(argv=None):
    parser = SafeParser(description=__doc__)
    sub = parser.add_subparsers(dest='action')
    prep = sub.add_parser('prepare')
    for name in ('baseline-pair', 'defaults-probe', 'executable', 'output-dir'):
        prep.add_argument('--' + name, type=Path, required=True)
    prep.add_argument('--history', action='store_true')
    prep.add_argument('--both-linear-solvers', action='store_true',
                      help='Also resolve momentum from the original deck; unsupported backends/settings fail before launch')
    inspect = sub.add_parser('inspect-profiles')
    for name in ('baseline-pair', 'defaults-probe', 'output'):
        inspect.add_argument('--' + name, type=Path, required=True)
    for action in ('run', 'compare', 'inspect-linear'):
        command = sub.add_parser(action)
        command.add_argument('--pair', type=Path, required=True)
        command.add_argument('--output', type=Path, required=True)
        if action == 'compare':
            command.add_argument('--detail-dir', type=Path, required=True)
        elif action == 'run':
            command.add_argument('--solver', choices=('mars', 'openaccel'), default='mars')
    args = parser.parse_args(argv)
    os.umask(0o077)
    public = dict(schema=SCHEMA, comparison_status='invalid_evidence', failed_check='private_preflight_or_capture')
    try:
        if args.action == 'prepare':
            prepare(args.baseline_pair.resolve(), args.defaults_probe.resolve(), args.executable.resolve(), args.output_dir.resolve(), args.history, args.both_linear_solvers)
            print('Original linear profiles prepared. No solver launched; baseline preserved.')
            return 0
        require(args.action in ('run', 'compare', 'inspect-linear', 'inspect-profiles'))
        with args.output.open('x') as stream:
            try:
                if args.action == 'inspect-profiles':
                    inspect_profiles(args.baseline_pair.resolve(), args.defaults_probe.resolve(), public)
                elif args.action == 'run':
                    pair = args.pair.resolve()
                    run(pair, public, args.solver)
                    public.update(comparison_status='capture_complete', failed_check='none')
                elif args.action == 'inspect-linear':
                    inspect_linear(args.pair.resolve(), public)
                else:
                    args.detail_dir.mkdir(mode=0o700)
                    compare(args.pair.resolve(), public, args.detail_dir)
            except Exception as error:
                public['failed_check'] = str(error) if isinstance(error, EvidenceError) and str(error) in ERRORS else 'private_preflight_or_capture'
                if args.action == 'run':
                    try:
                        public['launch_diagnostics'] = pressure.failure_summary(pair, args.solver)
                    except Exception:
                        public['launch_diagnostics'] = None
            json.dump(public, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')
        print('Pressure profile status written. Share only the public JSON.')
        return 0 if public['comparison_status'] in ('completed', 'capture_complete', 'inspection_complete') else 1
    except Exception as error:
        label = str(error) if isinstance(error, EvidenceError) and str(error) in ERRORS else 'private_preflight_or_capture'
        print('ERROR: pressure profile failed (' + label + ').', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
