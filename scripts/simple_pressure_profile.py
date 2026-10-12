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
          'launcher_changed', 'launcher_exit', 'runtime_recovery_mismatch', 'initial_fields_mismatch'}


def validate(v):
    require(set(v) == KEYS, 'unsupported_gpu_profile')
    for k, value in v.items():
        require(type(value) in (int, float) and math.isfinite(value), 'unsupported_gpu_profile')
        if k not in REAL_KEYS:
            require(value == int(value) and 0 <= value < 2**31, 'unsupported_gpu_profile')
    require(v['method'] in (0, 1) and v['kdim'] > 0 and v['maxiter'] > 0 and v['miniter'] <= v['maxiter']
        and 0 < v['rtol'] < 1 and v['atol'] >= 0 and v['coarsentype'] == 8 and v['interptype'] == 6
        and v['relaxorder'] == v['aggnumlevels'] == v['nodal'] == v['nodaldiag'] == v['fcycle'] == 0
        and v['numfunctions'] == v['cycletype'] == 1 and v['keeptranspose'] in (0, 1)
        and v['maxlevels'] > 0 and v['maxcoarsesize'] > 0 and v['mincoarsesize'] <= v['maxcoarsesize']
        and v['numpaths'] > 0 and v['coarserelax'] == 18
        and all(v[k] in (6, 13, 14, 18) for k in ('relaxtype', 'relax_down', 'relax_up'))
        and all(v[k] > 0 for k in ('numsweeps', 'sweeps_1', 'sweeps_2', 'sweeps_3'))
        and all(0 <= v[k] <= 1 for k in REAL_KEYS - {'rtol', 'atol'}), 'unsupported_gpu_profile')


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


def resolve(config, target, measured):
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
    validate(v)
    return v, changes


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
    return record, case


def prepare(baseline, probe, executable, output, history=False):
    record, case, launches, values, adaptations = baseline_data(baseline, probe)
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


def check_runtime(pair, ranks, values, history=False):
    expected = dict(values)
    expected['effective_relax_1'] = expected.pop('relax_down')
    expected['effective_relax_2'] = expected.pop('relax_up')
    expected['effective_relax_3'] = expected['coarserelax']
    paths = sorted((pair / 'mars').glob('flow-pressure-rank-*.settings'))
    names = ['flow-pressure-rank-{}.settings'.format(rank) for rank in range(ranks)]
    if history:
        paths += sorted((pair / 'mars').glob('flow-pressure-step-*-rank-*.settings'))
        names = ['flow-pressure-step-{}-rank-{}.settings'.format(step, rank)
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


def main(argv=None):
    parser = SafeParser(description=__doc__)
    sub = parser.add_subparsers(dest='action')
    prep = sub.add_parser('prepare')
    for name in ('baseline-pair', 'defaults-probe', 'executable', 'output-dir'):
        prep.add_argument('--' + name, type=Path, required=True)
    prep.add_argument('--history', action='store_true')
    for action in ('run', 'compare'):
        command = sub.add_parser(action)
        command.add_argument('--pair', type=Path, required=True)
        command.add_argument('--output', type=Path, required=True)
        if action == 'compare':
            command.add_argument('--detail-dir', type=Path, required=True)
        else:
            command.add_argument('--solver', choices=('mars', 'openaccel'), default='mars')
    args = parser.parse_args(argv)
    os.umask(0o077)
    public = dict(schema=SCHEMA, comparison_status='invalid_evidence', failed_check='private_preflight_or_capture')
    try:
        if args.action == 'prepare':
            prepare(args.baseline_pair.resolve(), args.defaults_probe.resolve(), args.executable.resolve(), args.output_dir.resolve(), args.history)
            print('Original pressure profile prepared. No solver launched; baseline preserved.')
            return 0
        require(args.action in ('run', 'compare'))
        with args.output.open('x') as stream:
            try:
                pair = args.pair.resolve()
                if args.action == 'run':
                    run(pair, public, args.solver)
                    public.update(comparison_status='capture_complete', failed_check='none')
                else:
                    args.detail_dir.mkdir(mode=0o700)
                    compare(pair, public, args.detail_dir)
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
        return 0 if public['comparison_status'] in ('completed', 'capture_complete') else 1
    except Exception as error:
        label = str(error) if isinstance(error, EvidenceError) and str(error) in ERRORS else 'private_preflight_or_capture'
        print('ERROR: pressure profile failed (' + label + ').', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
