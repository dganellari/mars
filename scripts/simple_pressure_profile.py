#!/usr/bin/env python3
"""Compare the original first step with an explicit GPU pressure configuration."""
import copy
import json
import math
import os
from pathlib import Path
import shutil
import sys

import simple_hypre_defaults as defaults
import simple_pressure_probe as pressure
import simple_pressure_settings as settings
from simple_pressure_replay import reference_controls
import simple_startup_probe as startup
from simple_snapshot_compare import EvidenceError, digest, options, require
from simple_public_diagnostics import SafeParser

SCHEMA = 'mars-simple-pressure-profile-v1'
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
          'launcher_changed', 'launcher_exit'}


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


def prepare(baseline, probe, executable, output):
    record, case, launches, values, adaptations = baseline_data(baseline, probe)
    require(executable.is_file() and os.access(str(executable), os.X_OK), 'executable_changed')
    libraries = startup.runtime_libraries(executable)
    require(libraries == launches['mars']['libraries'], 'libraries_changed')
    output.mkdir(mode=0o700)
    (output / 'reference').mkdir(mode=0o700)
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
    case = copy.deepcopy(case)
    case['arguments'] += ['--pressure-solver-profile', str(output / 'pressure.profile')]
    case['pressure_solver_profile_sha256'] = identity['sha256']
    startup.write_json(output / 'case.json', case)
    startup.write_json(output / 'pair.json', dict(record, case_sha256=digest(output / 'case.json'), pressure_profile=identity))
    check_inputs(output)


def check_inputs(pair):
    record = startup.pair_inputs(pair)
    identity = record['pressure_profile']
    require(identity['schema'] == SCHEMA, 'profile_identity')
    baseline = Path(identity['reference_pair'])
    require(digest(baseline / 'pair.json') == identity['reference_pair_sha256'], 'profile_identity')
    require(identity['defaults_hashes'] == {name: digest(pair / 'defaults' / name) for name in DEFAULT_FILES}, 'defaults_evidence')
    old, case, launches, values, adaptations = baseline_data(baseline, pair / 'defaults')
    require(identity['baseline_launch_hashes'] == {s: digest(baseline / ('reference' if s == 'openaccel' else 'mars') / 'launch.json')
                                                  for s in launches}, 'profile_identity')
    require(text_profile(values) == (pair / 'pressure.profile').read_text()
            and identity['adaptations'] == adaptations and identity['libraries'] == launches['mars']['libraries'], 'profile_identity')
    expected_case = copy.deepcopy(case)
    expected_case['arguments'] += ['--pressure-solver-profile', str(pair / 'pressure.profile')]
    expected_case['pressure_solver_profile_sha256'] = identity['sha256']
    require(startup.read_json(pair / 'case.json') == expected_case
            and record == dict(old, case_sha256=digest(pair / 'case.json'), pressure_profile=identity), 'profile_identity')
    return baseline, launches, values, identity


def check_launch(pair, baseline, old, identity):
    record = startup.verified_launch(pair, 'mars')
    require(all(record[k] == identity[k] for k in ('executable', 'executable_sha256', 'libraries')), 'profile_identity')
    require(record['ranks'] == old['ranks'] and pressure.launcher(record, pair, 'mars') ==
            pressure.launcher(old, baseline, 'mars'), 'launcher_changed')
    require(pressure.solver_environment(record['environment']) == pressure.solver_environment(old['environment']),
            'solver_environment_changed')
    return record


def run(pair, public):
    baseline, launches, _, identity = check_inputs(pair)
    old = launches['mars']
    exe = Path(identity['executable'])
    require(digest(exe) == identity['executable_sha256'], 'executable_changed')
    require(startup.runtime_libraries(exe) == identity['libraries'], 'libraries_changed')
    environment = dict(os.environ)
    for key in ('MARS_OPENACCEL_EXPORT_DIR', 'MARS_OPENACCEL_PUBLIC_FIXTURE'):
        environment.pop(key, None)
    changes = pressure.environment_changes(environment, old['environment'])
    public['environment_check'] = pressure.environment_report(changes)
    environment = pressure.restore_environment(environment, old['environment'])
    require(not pressure.environment_changes(environment, old['environment']), 'solver_environment_changed')
    startup.launch(pair, 'mars', exe, old['ranks'], pressure.launcher(old, baseline, 'mars'), environment)
    check_launch(pair, baseline, old, identity)


def check_runtime(pair, ranks, values):
    expected = dict(values)
    expected['effective_relax_1'] = expected.pop('relax_down')
    expected['effective_relax_2'] = expected.pop('relax_up')
    expected['effective_relax_3'] = expected['coarserelax']
    paths = sorted((pair / 'mars').glob('flow-pressure-rank-*.settings'))
    require(len(paths) == ranks, 'runtime_profile_mismatch')
    for rank in range(ranks):
        lines = (pair / 'mars' / ('flow-pressure-rank-{}.settings'.format(rank))).read_text().splitlines()
        parsed = [line.split() for line in lines]
        require(all(len(row) == 2 for row in parsed), 'runtime_profile_mismatch')
        actual = {k: float(v) for k, v in parsed}
        require(len(actual) == len(parsed) and set(actual) == set(expected) | {'hypre_release', 'effective_levels'}
                and all(math.isfinite(v) for v in actual.values()) and all(actual[k] == v for k, v in expected.items())
                and actual['effective_levels'] >= 1, 'runtime_profile_mismatch')


def compare(pair, public, details):
    baseline, launches, values, identity = check_inputs(pair)
    record = check_launch(pair, baseline, launches['mars'], identity)
    check_runtime(pair, record['ranks'], values)
    result = {}
    startup.compare(pair, result, details, reference_pair=baseline)
    check_inputs(pair)
    public.update(comparison_status='completed', failed_check='none', comparison=result,
        original_pressure_targets_preserved=True, recorded_gpu_pressure_settings_verified=True,
        gpu_adaptations=identity['adaptations'], reference_solver_launched=False,
        identical_linear_solvers_verified=False, nonlinear_convergence_verified=False)


def main(argv=None):
    parser = SafeParser(description=__doc__)
    sub = parser.add_subparsers(dest='action')
    prep = sub.add_parser('prepare')
    for name in ('baseline-pair', 'defaults-probe', 'executable', 'output-dir'):
        prep.add_argument('--' + name, type=Path, required=True)
    for action in ('run', 'compare'):
        command = sub.add_parser(action)
        command.add_argument('--pair', type=Path, required=True)
        command.add_argument('--output', type=Path, required=True)
        if action == 'compare':
            command.add_argument('--detail-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    os.umask(0o077)
    public = dict(schema=SCHEMA, comparison_status='invalid_evidence', failed_check='private_preflight_or_capture')
    try:
        if args.action == 'prepare':
            prepare(args.baseline_pair.resolve(), args.defaults_probe.resolve(), args.executable.resolve(), args.output_dir.resolve())
            print('Original pressure profile prepared. No solver launched; baseline preserved.')
            return 0
        require(args.action in ('run', 'compare'))
        with args.output.open('x') as stream:
            try:
                pair = args.pair.resolve()
                if args.action == 'run':
                    run(pair, public)
                    public.update(comparison_status='capture_complete', failed_check='none')
                else:
                    args.detail_dir.mkdir(mode=0o700)
                    compare(pair, public, args.detail_dir)
            except Exception as error:
                public['failed_check'] = str(error) if isinstance(error, EvidenceError) and str(error) in ERRORS else 'private_preflight_or_capture'
                if args.action == 'run':
                    try:
                        public['launch_diagnostics'] = pressure.failure_summary(pair, 'mars')
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
