#!/usr/bin/env python3
"""User-local capture and replay of a rejected pressure system; only fixed summaries are public."""
import json
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys

import simple_hypre_defaults as defaults
import simple_pressure_probe as probe
import simple_pressure_settings as settings
import simple_startup_probe as startup
from simple_snapshot_compare import digest, options, require
from simple_public_diagnostics import DiagnosticState, SafeParser

SCHEMA = 'mars-simple-pressure-replay-v1'
ARCHIVE_SCHEMA = 'mars-pressure-replay-inputs-v1'
CAPTURE_MARKER = '[simple-pressure-capture] complete; files are private; original rejection retained'
ROOT = Path(__file__).resolve().parent.parent
CPP = ROOT / 'tests/reference/openaccel/simple_performance/pressure_replay.cpp'
HEADERS = [ROOT / 'backend/distributed/unstructured/solvers/mars_hypre_pressure_settings.hpp',
           ROOT / 'backend/distributed/unstructured/fem/segregated/mars_segregated_pressure_capture.hpp']


def hashes(paths):
    return {str(Path(p).resolve()): digest(p) for p in paths}


def verify(record):
    defaults.unchanged(record)


def file_identity_checks(files, category, replacements=None):
    issues = {name: set() for name in ('missing', 'changed', 'unreadable')}
    valid = isinstance(files, dict)
    if valid:
        for path, expected in files.items():
            if (not isinstance(path, str) or not Path(path).is_absolute()
                    or not isinstance(expected, str) or not re.fullmatch('[0-9a-f]{64}', expected)):
                valid = False
                continue
            kind = category(path)
            source = (replacements or {}).get(path, path)
            try:
                if not Path(source).is_file():
                    issues['missing'].add(kind)
                elif digest(source) != expected:
                    issues['changed'].add(kind)
            except OSError:
                issues['unreadable'].add(kind)
    return dict(record_valid=valid, matched=valid and not any(issues.values()),
                **{name: sorted(values) for name, values in issues.items()})


def replay_input_kind(path, record, capture_dir, capture):
    command = record.get('command', [])
    if len(command) >= 4 and path == command[-4]:
        return 'executable'
    if path == str(capture_dir / 'capture.json') or path in capture['files']:
        return 'capture_input'
    kind = defaults.library_kind(path)
    if kind is not None:
        return kind + '_library'
    if re.search(r'\.so(?:\.[0-9]+)*$|\.dylib$', Path(path).name):
        return 'other_runtime_library'
    return 'source_or_build_input'


def external_inputs(record, capture_dir, capture):
    live = set(capture['files']) | {str(capture_dir / 'capture.json')}
    return {p: value for p, value in record['inputs'].items() if p not in live}


def replay_input_checks(path, record, capture_dir, capture, archive=None):
    replacements = {}
    if archive is not None:
        saved = startup.read_json(archive / 'archive.json')
        require(saved['schema'] == ARCHIVE_SCHEMA and saved['backend'] == record['backend']
                and saved['capture_sha256'] == record['capture_sha256']
                and saved['replay_sha256'] == digest(path / 'replay.json')
                and saved['inputs'] == external_inputs(record, capture_dir, capture))
        # Capture and result files stay live. Only dependencies from another
        # environment may be supplied by an exact, hash-checked byte copy.
        for source, value in saved['inputs'].items():
            require(isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value))
            replacements[source] = archive / 'objects' / value
    checks = file_identity_checks(record['inputs'],
        lambda p: replay_input_kind(p, record, capture_dir, capture), replacements)
    checks['scope'] = 'archived_dependencies_and_live_capture_inputs' if archive is not None else 'live_inputs'
    return checks


def checked_replay_record(path, name, capture_dir, capture, check, archive=None):
    check.update(verified=False, failed_check='replay_record')
    record = startup.read_json(path / 'replay.json')
    check['failed_check'] = 'replay_binding'
    require(record['schema'] == SCHEMA and record['backend'] == name and record['exit_code'] == 0
            and record['capture_sha256'] == digest(capture_dir / 'capture.json'))
    require(record['profile'] in ('captured', 'reference') and (name != 'mars' or record['profile'] == 'captured'))
    check['failed_check'] = 'input_archive' if archive is not None else 'replay_inputs'
    check['inputs'] = replay_input_checks(path, record, capture_dir, capture, archive)
    check['failed_check'] = 'replay_outputs'
    check['outputs'] = file_identity_checks(record['files'], lambda p: 'replay_output')
    check['failed_check'] = 'replay_inputs' if not check['inputs']['matched'] else 'replay_outputs'
    require(check['inputs']['matched'] and check['outputs']['matched'])
    check.update(verified=True, failed_check='none')
    return record


def archive_inputs(args, output, public):
    public.update(solver_launched=False, convergence_verified=False, failed_check='capture_identity')
    directory = args.capture_run.resolve(); capture = captured_record(directory)
    path = args.replay_run.resolve()
    public['failed_check'] = 'replay_record'
    record_hash = digest(path / 'replay.json')
    name = startup.read_json(path / 'replay.json')['backend']
    require(name in ('mars', 'reference'))
    check = public['replay_evidence_checks'] = {}
    try:
        record = checked_replay_record(path, name, directory, capture, check)
    except Exception:
        public['failed_check'] = check['failed_check']
        raise
    public['failed_check'] = 'loaded_library_identity'
    expected = capture['libraries'] if name == 'mars' else capture['reference']['libraries']
    loaded_libraries_match(path / 'result', capture['ranks'], expected, Path(record['command'][-4]), public)
    public['failed_check'] = 'input_archive_copy'
    inputs = external_inputs(record, directory, capture)
    objects = output / 'objects'; objects.mkdir(mode=0o700)
    for source, value in inputs.items():
        target = objects / value
        if not target.exists():
            with Path(source).open('rb') as src, target.open('xb') as dst:
                shutil.copyfileobj(src, dst)
        require(digest(target) == value)
    public['failed_check'] = 'evidence_changed'
    require(digest(path / 'replay.json') == record_hash)
    verify(record['inputs']); verify(record['files']); verify(capture['files'])
    # This records reproducible bytes, not a signed attestation or a solver verdict.
    startup.write_json(output / 'archive.json', dict(schema=ARCHIVE_SCHEMA, backend=name,
        capture_sha256=record['capture_sha256'], replay_sha256=record_hash, inputs=inputs))
    public.update(comparison_status='input_archive_complete', failed_check='none', backend=name,
                  original_inputs_verified=True, archived_dependencies_verified=True)


def marker(directory, schema):
    lines = (directory / 'complete').read_text().splitlines()
    require(len(lines) == 2 and lines[0] == schema and lines[1].isdigit())
    ranks = int(lines[1]); require(1 <= ranks <= 1000000)
    return ranks


def rank_file(directory, rank, suffix):
    return directory / ('rank-{:06d}'.format(rank) + suffix)


def numeric_file(path):
    result = {}
    for line in path.read_text().splitlines():
        key, value = line.split()
        require(key not in result)
        result[key] = float(value)
    return result


def stopping_checks(reports, residual):
    require(bool(reports))
    integer_keys = ('method', 'maxiter', 'result_iterations', 'result_converged',
                    'result_solve_error', 'result_global_error', 'result_fatal_error')
    for report in reports:
        for key in integer_keys:
            value = report[key]
            require(math.isfinite(value) and 0 <= value <= 2**31-1 and value == int(value))
        require(report['method'] in (0, 1) and report['result_converged'] in (0, 1)
                and report['result_fatal_error'] in (0, 1))
        require(math.isfinite(report['rtol']) and 0 < report['rtol'] < 1
                and math.isfinite(report['atol']) and report['atol'] >= 0)
        require(not math.isfinite(report['result_reported']) or report['result_reported'] >= 0)
    relations = {('below' if r['result_iterations'] < r['maxiter'] else
                  'at' if r['result_iterations'] == r['maxiter'] else 'above') for r in reports}
    relation = next(iter(relations)) if len(relations) == 1 else 'mixed'
    # Compare exit metadata, not rank-dependent AMG hierarchy sizes.
    keys = integer_keys + ('rtol', 'atol')
    consistent = all(tuple(r[k] for k in keys) == tuple(reports[0][k] for k in keys) for r in reports)
    finite = all(math.isfinite(r['result_reported']) for r in reports)
    zero = all(r['result_iterations'] == 0 for r in reports)
    fatal = any(r['result_fatal_error'] != 0 for r in reports)
    if not consistent:
        assessment = 'rank_reports_disagree'
    elif fatal:
        assessment = 'fatal_backend_error'
    elif not finite:
        assessment = 'nonfinite_reported_residual'
    elif not residual['finite']:
        assessment = 'nonfinite_candidate_residual'
    elif residual['residual_passed']:
        assessment = 'independent_residual_passed'
    elif residual['residual_inconclusive']:
        assessment = 'independent_residual_inconclusive'
    else:
        require(residual['residual_failed'])
        assessment = ('failed_without_iterations' if zero else
                      'failed_before_iteration_limit' if relation == 'below' else
                      'failed_at_or_above_iteration_limit')
    methods = {r['method'] for r in reports}
    return dict(scope='saved_exit_metadata_not_internal_branch_trace', assessment=assessment,
        backend='mixed' if len(methods) != 1 else 'FlexGMRES' if 1 in methods else 'GMRES',
        exit_metadata_agrees_across_ranks=consistent, iteration_limit_relation=relation,
        iterations_zero_on_all_ranks=zero,
        solve_return_nonzero=any(r['result_solve_error'] != 0 for r in reports),
        global_error_nonzero=any(r['result_global_error'] != 0 for r in reports),
        fatal_backend_error_seen=fatal,
        reported_relative_residual_finite=finite,
        reported_relative_residual_below_rtol=(all(r['result_reported'] <= r['rtol'] for r in reports)
                                               if finite else None),
        absolute_tolerance_enabled=any(r['atol'] != 0 for r in reports),
        convergence_claim_contradicted=(residual['residual_failed'] and any(r['result_converged'] == 1 for r in reports)))


def capture_inputs(directory):
    ranks = marker(directory, 'mars-pressure-capture-v1')
    files = [directory / 'complete']
    previous, first = 0, None
    for rank in range(ranks):
        path = rank_file(directory, rank, '.bin')
        with path.open('rb') as stream:
            h = struct.unpack('<13Q', stream.read(104))
            target = struct.unpack('<2d', stream.read(16))
        require(h[:2] == (0x4d53505245535331, 1) and h[2:4] == (rank, ranks))
        begin, end, total, nodes, nnz, iteration = h[4:10]
        require(begin == previous and begin < end <= total <= 2**31-1 and nodes <= 2**31-1 and end-begin <= nnz <= 2**31-1)
        require(iteration > 0 and h[10] in (0, 1) and h[11] in (0, 1) and not (h[10] and h[11]) and h[12] == 1)
        require(all(math.isfinite(x) for x in target) and target[0] >= 0 and 0 < target[1] < 1)
        require(path.stat().st_size == 120 + 4*(end-begin+1) + 12*nnz + 16*nodes + 8*(end-begin))
        identity = (total, iteration, target)
        require(first is None or identity == first)
        first, previous = identity, end
        config = rank_file(directory, rank, '.settings')
        values = numeric_file(config)
        require(all(math.isfinite(x) for x in values.values()))
        require(values['rtol'] == target[1] and values['atol'] == target[0])
        files += [path, config]
    require(previous == first[0])
    require(set(directory.iterdir()) == set(files))
    return ranks, first[2], hashes(files)


def reference_controls(config, target):
    require(str(config.get('family', '')).lower() == 'hypre')
    require(config.get('normalize_matrix', False) is False and config.get('diagonal_scaling', False) is False)
    outer = config.get('options', {})
    require(isinstance(outer, dict) and not set(outer) - {'solver', 'precond'})
    solver = settings.lower_options(outer.get('solver', {}))
    amg = settings.lower_options(outer.get('precond', {}))
    method = str(solver.get('type', 'gmres')).lower()
    require(method in ('gmres', 'flexgmres') and str(amg.get('type', 'none')).lower() == 'boomeramg')
    require(not set(solver) - {'type', 'kdim', 'printlevel', 'logging'})
    supported = {key for key, _, _ in settings.AMG.values() if key} | set(settings.AMG_LIBRARY)
    require(not set(amg) - supported - {'type', 'printlevel', 'logging', 'maxiter', 'tol'})
    def number(value):
        result = settings.number(value)
        require(result is not settings.UNKNOWN)
        return result
    result = {key: number(value) for key, value in amg.items() if key in supported}
    result.update(method=int(method == 'flexgmres'), maxiter=number(config.get('max_iterations', 20)),
                  rtol=number(config.get('rtol', 1e-6)), atol=number(config.get('atol', 1e-16)))
    if 'kdim' in solver:
        result['kdim'] = number(solver['kdim'])
    require(result['rtol'] == target[1] and result['atol'] == target[0])
    require(all(math.isfinite(x) for x in result.values()))
    # Omitted options remain omitted, so the captured reference library supplies them.
    return result


def capture(pair, executable, output, public):
    public['failed_check'] = 'saved_pair'
    pair = pair.resolve(); config, _ = settings.saved_configuration(pair)
    require(not startup.read_json(pair / 'pair.json').get('first_step_audit', False))
    baseline = settings.saved_launch(pair, 'mars'); reference = settings.saved_launch(pair, 'openaccel')
    require(reference['status'] == 'finished' and baseline['ranks'] == reference['ranks'])
    arguments = options(startup.solver_arguments(pair, 'mars'))
    require(arguments.get('--pressure-refinement', '0') == '0' and '--first-step-audit' not in arguments)
    public['failed_check'] = 'reference_configuration'
    requested_target = (float(arguments['--pressure-linear-atol']), float(arguments['--pressure-linear-rtol']))
    controls = reference_controls(config, requested_target)
    executable = executable.resolve(strict=True)
    arguments.update({'--output-prefix': str(output / 'flow'), '--field-output': 'none',
                      '--snapshot-iterations': '0', '--pressure-failure-capture': str(output / 'system')})
    environment = probe.restore_environment(os.environ, baseline['environment'])
    public['failed_check'] = 'solver_environment_changed'
    require(not probe.environment_changes(environment, baseline['environment']))
    require(environment.get('MARS_HYPRE_SPMV_VENDOR', '0') == '0')
    launch = probe.launcher(baseline, pair, 'mars')
    command = launch + [str(executable)] + [word for key, value in arguments.items() for word in (key, value)]
    public['failed_check'] = 'runtime_libraries'
    libraries = startup.runtime_libraries(executable)
    defaults.probe_dependencies(libraries, libraries, True)
    inputs = hashes([executable, pair / 'pair.json', pair / 'case.json', pair / 'reference/input.i'] +
                    [pair / name / file for name in ('reference', 'mars')
                     for file in ('launch-start.json', 'launch.json', 'run.log', 'run.exit')])
    inputs.update(libraries)
    record = dict(schema=SCHEMA, pair=str(pair), command=command, launcher=launch,
                  ranks=baseline['ranks'], environment=probe.solver_environment(environment),
                  libraries=libraries, reference=reference, inputs=inputs)
    startup.write_json(output / 'launch-start.json', record)
    public['failed_check'] = 'capture_launch'
    with (output / 'run.log').open('xb') as log:
        process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=environment)
    (output / 'run.exit').write_text(str(process.returncode) + '\n')
    public['failed_check'] = 'capture_missing'
    require(process.returncode not in (0, 2))
    require((output / 'run.log').read_text(errors='replace').count(CAPTURE_MARKER) == 1)
    ranks, target, files = capture_inputs(output / 'system')
    require(ranks == record['ranks'] and target == requested_target)
    with (output / 'reference.settings').open('x') as file:
        for key, value in sorted(controls.items()):
            file.write('{} {:.17g}\n'.format(key, value))
    public['failed_check'] = 'inputs_changed'
    verify(inputs)
    record.update(exit_code=process.returncode, files=dict(files, **hashes(
        [output / 'reference.settings', output / 'run.log', output / 'run.exit', output / 'launch-start.json'])))
    startup.write_json(output / 'capture.json', record)
    public.update(comparison_status='capture_complete', failed_check='none', original_rejection_preserved=True)


def captured_record(directory):
    record = startup.read_json(directory / 'capture.json')
    require(record['schema'] == SCHEMA and record['exit_code'] not in (0, 2))
    verify(record['files'])
    ranks, _, files = capture_inputs(directory / 'system')
    require(ranks == record['ranks'] and all(record['files'][key] == value for key, value in files.items()))
    return record


def loaded_library_checks(directory, ranks, expected, executable):
    accepted = {kind: {value for name, value in expected.items() if defaults.library_kind(name) == kind}
                for kind in ('hypre', 'mpi')}
    failures = {kind: set() for kind in accepted}
    for rank in range(ranks):
        try:
            identity = rank_file(directory, rank, '.libraries').read_text().splitlines()
        except (OSError, UnicodeError):
            for issues in failures.values(): issues.add('record_unreadable')
            continue
        if len(identity) != 2:
            for issues in failures.values(): issues.add('record_malformed')
            continue
        for path, kind in zip(identity, ('hypre', 'mpi')):
            if not Path(path).is_absolute():
                failures[kind].add('non_absolute_path')
                continue
            try:
                value = digest(path)
                if value not in accepted[kind]:
                    failures[kind].add('executable_instead_of_library' if Path(path).resolve() == executable.resolve()
                                       else 'hash_mismatch')
            except (OSError, ValueError):
                failures[kind].add('file_unreadable')
    return {kind: dict(matched=not issues, failures=sorted(issues)) for kind, issues in failures.items()}


def loaded_libraries_match(directory, ranks, expected, executable, public):
    checks = loaded_library_checks(directory, ranks, expected, executable)
    public['library_identity_checks'] = checks
    public['loaded_library_identity_verified'] = all(check['matched'] for check in checks.values())
    require(public['loaded_library_identity_verified'])


def inspect_replay(directory, public):
    public.update(comparison_status='invalid_evidence', failed_check='saved_launch_record',
                  convergence_verified=False, solver_launched_by_inspection=False)
    record = startup.read_json(directory / 'launch-start.json')
    require(record['schema'] == SCHEMA and record['backend'] in ('mars', 'reference'))
    require(record['profile'] in ('captured', 'reference'))
    command = record['command']
    require(isinstance(command, list) and len(command) >= 4)
    capture_dir = Path(command[-3]).parent
    require(Path(command[-3]).name == 'system' and Path(command[-1]).resolve() == (directory / 'result').resolve())
    require(digest(capture_dir / 'capture.json') == record['capture_sha256'])
    capture = startup.read_json(capture_dir / 'capture.json')
    require(capture['schema'] == SCHEMA and type(capture['ranks']) is int and 1 <= capture['ranks'] <= 1000000)
    ranks = capture['ranks']
    expected = capture['libraries'] if record['backend'] == 'mars' else capture['reference']['libraries']
    public.update(backend=record['backend'], profile=record['profile'])
    checks = {}
    def check(label, action):
        try:
            action()
            checks[label] = True
        except Exception:
            checks[label] = False
    check('launch_inputs_unchanged', lambda: verify(record['inputs']))
    status_file = directory / 'run.exit'
    public['exit_file_present'] = status_file.is_file()
    public['process_exit_code'] = None
    if status_file.is_file():
        text = status_file.read_text().strip()
        if re.fullmatch(r'-?[0-9]{1,3}', text) and -127 <= int(text) <= 255:
            public['process_exit_code'] = int(text)
    state = DiagnosticState()
    public['log_present'] = (directory / 'run.log').is_file()
    public['library_load_error_seen'] = False
    public['replay_error_seen'] = False
    if public['log_present']:
        with (directory / 'run.log').open(errors='replace') as log:
            for line in log:
                state.feed(line)
                public['library_load_error_seen'] |= any(message in line for message in (
                    'error while loading shared libraries:', 'symbol lookup error:',
                    'Library not loaded:', 'Symbol not found:')) or ('version ' in line and ' not found (required by ' in line)
                public['replay_error_seen'] |= line.strip() == 'ERROR: private pressure replay failed'
    for key in ('mpi_abort_seen', 'scheduler_time_limit_seen', 'scheduler_out_of_memory_seen', 'scheduler_signal_seen'):
        public[key] = key in state.failures
    public['scheduler_messages'] = list(state.scheduler_messages)
    result = directory / 'result'
    public['result_directory_present'] = result.is_dir()
    for name in ('libraries', 'solution', 'report'):
        present = [rank_file(result, rank, '.' + name).is_file() for rank in range(ranks)]
        public['any_' + name + '_parts_present'] = any(present)
        public['all_' + name + '_parts_present'] = all(present)
    public['completion_marker_present'] = (result / 'complete').is_file()
    check('completion_marker_valid', lambda: require(marker(result, 'mars-pressure-replay-v1') == ranks))
    check('loaded_library_identity_verified', lambda: loaded_libraries_match(
        result, ranks, expected, Path(command[-4]), public))
    public.update(checks)
    if not checks['launch_inputs_unchanged']:
        failure = 'inputs_changed'
    elif public['process_exit_code'] is None:
        failure = 'launcher_exit_missing_or_invalid'
    elif public['process_exit_code'] != 0:
        failure = 'launcher_exit'
    elif not checks['completion_marker_valid']:
        failure = 'completion_marker'
    elif not checks['loaded_library_identity_verified']:
        failure = 'loaded_library_identity'
    elif not public['all_solution_parts_present'] or not public['all_report_parts_present']:
        failure = 'replay_parts_missing'
    else:
        failure = 'none'
    public.update(comparison_status='inspection_complete', failed_check=failure,
                  scope='saved_launch_log_and_file_checks_not_residual_or_convergence')


def replay(args, output, public):
    public['failed_check'] = 'capture_identity'
    capture_dir = args.capture_run.resolve(); record = captured_record(capture_dir)
    profile = args.profile or ('reference' if args.backend == 'reference' else 'captured')
    require(not (args.backend == 'mars' and profile != 'captured'))
    expected = record['libraries'] if args.backend == 'mars' else record['reference']['libraries']
    public['failed_check'] = 'runtime_libraries'
    defaults.unchanged(expected)
    if args.backend == 'reference':
        public['failed_check'] = 'reference_compile'
        require(args.build_cache is not None and args.executable is None)
        library, include = defaults.matching_install(expected)
        compiler, compile_flags, link_flags, cache_inputs = defaults.cached_toolchain(args.build_cache.resolve(), True)
        require(defaults.find_compiler(compiler) is not None)
        executable = output / 'pressure-replay'
        # dladdr needs the shared-library address, not an executable PLT stub.
        command = [compiler, '-std=c++17', '-O2', '-I' + str(include)] + compile_flags + [
            str(CPP), str(library)] + link_flags + ['-fPIC', '-Wl,-rpath,' + str(library.parent), '-ldl', '-o', str(executable)]
        with (output / 'compile.log').open('xb') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        require(result.returncode == 0)
        inputs = dict(cache_inputs, **hashes([CPP] + HEADERS + list(include.glob('*.h'))))
    else:
        require(args.executable is not None and args.build_cache is None)
        executable = args.executable.resolve(strict=True); inputs = {}
    libraries = startup.runtime_libraries(executable)
    defaults.probe_dependencies(libraries, expected, True)
    inputs.update(hashes([executable, capture_dir / 'capture.json'])); inputs.update(libraries); inputs.update(record['files'])
    configuration = capture_dir / ('reference.settings' if profile == 'reference' else 'system')
    command = record['launcher'] + [str(executable), str(capture_dir / 'system'), str(configuration), str(output / 'result')]
    launch = dict(schema=SCHEMA, backend=args.backend, profile=profile, command=command, inputs=inputs,
                  environment=probe.solver_environment(os.environ), capture_sha256=digest(capture_dir / 'capture.json'))
    startup.write_json(output / 'launch-start.json', launch)
    public['failed_check'] = 'replay_launch'
    with (output / 'run.log').open('xb') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    (output / 'run.exit').write_text(str(result.returncode) + '\n')
    public['failed_check'] = 'launcher_exit'
    require(result.returncode == 0)
    public['failed_check'] = 'completion_marker'
    require(marker(output / 'result', 'mars-pressure-replay-v1') == record['ranks'])
    public['failed_check'] = 'loaded_library_identity'
    loaded_libraries_match(output / 'result', record['ranks'], expected, executable, public)
    public['failed_check'] = 'inputs_changed'
    verify(inputs)
    launch.update(exit_code=0, files=hashes(list((output / 'result').iterdir()) +
                  [output / 'run.exit', output / 'run.log', output / 'launch-start.json']))
    startup.write_json(output / 'replay.json', launch)
    public.update(comparison_status='replay_complete', failed_check='none', loaded_library_identity_verified=True,
                  backend=args.backend, profile=profile, convergence_verified=False)


def compare(args, public):
    public['failed_check'] = 'capture_identity'
    directory = args.capture_run.resolve(); capture = captured_record(directory)
    public['failed_check'] = 'checker_executable'
    checker = args.checker.resolve(strict=True); checker_hash = digest(checker)
    checks = {}
    stopping = {}
    profiles = {}
    records = {}
    archives = {name: getattr(args, name + '_input_archive') for name in ('mars', 'reference')}
    evidence = public['replay_evidence_checks'] = {}
    # Check both environments before launching the checker, so one missing mount
    # does not hide a second problem in the other replay.
    for name, path in (('mars', args.mars_run), ('reference', args.reference_run)):
        evidence[name] = {}
        try:
            records[name] = checked_replay_record(path, name, directory, capture, evidence[name], archives[name])
        except Exception:
            pass
    for name in ('mars', 'reference'):
        if not evidence[name]['verified']:
            public['failed_candidate'] = name
            public['failed_check'] = name + '_' + evidence[name]['failed_check']
            require(False)
    for name, path in (('original', None), ('mars', args.mars_run), ('reference', args.reference_run)):
        public['failed_candidate'] = name
        if path is not None:
            record = records[name]
            profiles[name] = record['profile']
            public['failed_check'] = name + '_replay_inputs'
            require(replay_input_checks(path, record, directory, capture, archives[name])['matched'])
            public['failed_check'] = name + '_replay_outputs'
            verify(record['files'])
        command = [str(checker), str(directory / 'system'), '-' if path is None else str(path / 'result')]
        public['failed_check'] = name + '_checker_launch'
        try:
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except OSError as error:
            public['checker_launch_error'] = ('not_found' if isinstance(error, FileNotFoundError) else
                                             'permission_denied' if isinstance(error, PermissionError) else 'os_error')
            raise
        public['failed_check'] = name + '_checker_exit'
        if result.returncode != 0:
            stderr = result.stderr.decode('utf-8', errors='replace')
            public['checker_diagnostics'] = dict(process_exit_code=result.returncode,
                library_load_error_seen=any(message in stderr for message in (
                    'error while loading shared libraries:', 'symbol lookup error:',
                    'Library not loaded:', 'Symbol not found:')) or ('version ' in stderr and ' not found (required by ' in stderr),
                residual_check_error_seen='ERROR: private pressure residual check failed' in stderr)
        require(result.returncode == 0)
        public['failed_check'] = name + '_checker_output'
        value = json.loads(result.stdout.decode('utf-8'))
        keys = {'capture_valid', 'original_referenced_copies_equal_owners', 'finite', 'residual_passed', 'residual_failed', 'residual_inconclusive'}
        require(value['schema'] == 'mars-pressure-residual-v1' and all(type(value[key]) is bool for key in keys))
        checks[name] = {key: value[key] for key in sorted(keys)}
        if path is not None:
            public['failed_check'] = name + '_replay_report'
            reports = [numeric_file(rank_file(path / 'result', rank, '.report')) for rank in range(capture['ranks'])]
            checks[name]['backend_fatal_error_seen'] = any(r['result_fatal_error'] != 0 for r in reports)
            checks[name]['backend_converged_flag'] = all(r['result_converged'] == 1 for r in reports)
            control_matches = True
            for rank, report in enumerate(reports):
                config = (rank_file(directory / 'system', rank, '.settings') if record['profile'] == 'captured'
                          else directory / 'reference.settings')
                expected = numeric_file(config)
                control_matches = control_matches and all(report.get(k) == v for k, v in expected.items()
                                                          if k != 'hypre_release' and not k.startswith('effective_'))
            checks[name]['recorded_controls_match_requested_profile'] = control_matches
            public['failed_check'] = name + '_replay_stopping_report'
            stopping[name] = stopping_checks(reports, checks[name])
            public['failed_check'] = name + '_replay_inputs_changed'
            require(replay_input_checks(path, record, directory, capture, archives[name])['matched'])
            public['failed_check'] = name + '_replay_outputs_changed'
            verify(record['files'])
    public.pop('failed_candidate', None)
    public['failed_check'] = 'checker_changed'
    require(digest(checker) == checker_hash)
    public['failed_check'] = 'capture_changed'
    verify(capture['files'])
    public.update(comparison_status='completed', failed_check='none', residual_checks=checks, stopping_checks=stopping,
                  same_frozen_system_verified=True, profiles=profiles, actual_exit_branch_verified=False,
                  original_amg_hierarchy_reused=False, nonlinear_convergence_verified=False)


def main(argv=None):
    parser = SafeParser(description=__doc__)
    sub = parser.add_subparsers(dest='action')
    c = sub.add_parser('capture'); c.add_argument('--pair', type=Path, required=True); c.add_argument('--executable', type=Path, required=True)
    c.add_argument('--output-dir', type=Path, required=True)
    r = sub.add_parser('replay'); r.add_argument('--capture-run', type=Path, required=True)
    r.add_argument('--backend', choices=('mars', 'reference'), required=True); r.add_argument('--profile', choices=('captured', 'reference'))
    r.add_argument('--executable', type=Path); r.add_argument('--build-cache', type=Path); r.add_argument('--output-dir', type=Path, required=True)
    c = sub.add_parser('compare'); c.add_argument('--capture-run', type=Path, required=True)
    c.add_argument('--mars-run', type=Path, required=True); c.add_argument('--reference-run', type=Path, required=True)
    c.add_argument('--mars-input-archive', type=Path); c.add_argument('--reference-input-archive', type=Path)
    c.add_argument('--checker', type=Path, required=True); c.add_argument('--output', type=Path, required=True)
    c = sub.add_parser('archive-inputs'); c.add_argument('--capture-run', type=Path, required=True)
    c.add_argument('--replay-run', type=Path, required=True); c.add_argument('--output-dir', type=Path, required=True)
    c = sub.add_parser('inspect'); c.add_argument('--replay-run', type=Path, required=True)
    c.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv); os.umask(0o077)
    public = dict(schema=SCHEMA, comparison_status='invalid_evidence', failed_check='preflight')
    destination = None
    try:
        if args.action in ('compare', 'inspect'):
            destination = args.output
            if args.action == 'compare':
                compare(args, public)
            else:
                inspect_replay(args.replay_run.resolve(), public)
        else:
            require(args.action in ('capture', 'replay', 'archive-inputs'))
            args.output_dir.mkdir(mode=0o700); directory = args.output_dir.resolve(); destination = directory / 'public.json'
            if args.action == 'capture':
                capture(args.pair, args.executable, directory, public)
            elif args.action == 'replay':
                replay(args, directory, public)
            else:
                archive_inputs(args, directory, public)
    except Exception:
        if args.action in ('compare', 'inspect'):
            destination = args.output
        elif destination is None:
            print('ERROR: new private output directory required.', file=sys.stderr); return 1
    try:
        startup.write_json(destination, public)
    except Exception:
        print('ERROR: cannot create new public summary.', file=sys.stderr); return 1
    print('Pressure replay status written. Share only public.json; completion alone is not convergence.')
    return 1 if public['comparison_status'] == 'invalid_evidence' else 0


if __name__ == '__main__':
    sys.exit(main())
