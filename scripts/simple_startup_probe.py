#!/usr/bin/env python3
"""Capture fresh SIMPLE startup histories locally; export only fixed diagnostic labels."""

import argparse
import copy
import csv
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys

from prepare_simple_deck import load_deck, translate
from simple_public_diagnostics import SafeParser
from simple_snapshot_compare import (COMPLETION, HEADER, EvidenceError, controls, coordinates, digest,
                                     field_errors, mars_fields, node_ids, options,
                                     reference_fields, require, result_files)

STEPS = 20
SCHEMA = 'mars-simple-startup-v1'
PREPARATION_ERRORS = {'saved_case_format', 'saved_deck_identity', 'saved_control_identity',
                      'saved_mesh_path', 'reference_length', 'modified_controls', 'launcher_exit',
                      'runtime_library_probe', 'runtime_libraries_unresolved',
                      'runtime_library_paths', 'runtime_library_unreadable', 'reference_outputs'}


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def prepare(case_path, reference_dir, output):
    import yaml
    case = read_json(case_path)
    require(case['format'] == 'mars-simple-deck-v1', 'saved_case_format')
    saved = options(case['arguments'])
    decks = [p for p in reference_dir.iterdir() if p.is_file()
             and p.suffix.lower() in ('.i', '.yaml', '.yml') and digest(p) == case['deck_sha256']]
    require(len(decks) == 1, 'saved_deck_identity')
    original = load_deck(decks[0].read_bytes())
    policy = case.get('pressure_linear_policy', 'mars')
    mapped, mesh_name = translate(original, policy)  # Reject restarts and nonzero initialization.
    expected = options(mapped)
    require(set(saved) == set(expected) | {'--mesh', '--mesh-format', '--reference-length'}, 'saved_control_identity')
    require(all(saved[k] == v for k, v in expected.items()) and saved['--mesh-format'] == 'exodus', 'saved_control_identity')
    mesh = Path(saved['--mesh']).resolve()
    require(mesh.is_file() and (decks[0].parent / mesh_name).resolve() == mesh, 'saved_mesh_path')
    require(math.isfinite(float(saved['--reference-length'])) and float(saved['--reference-length']) > 0, 'reference_length')
    modified = copy.deepcopy(original)
    modified['mesh']['file_path'] = str(mesh)
    solver = modified['simulation']['solver']
    convergence = solver['solver_control']['basic_settings']['convergence_controls']
    convergence['min_iterations'] = STEPS
    convergence['max_iterations'] = STEPS
    solver['output_control'] = dict(file_path='results.e', output_frequency=1,
                                  output_fields=['velocity', 'pressure'], corrected_boundary_values=False)
    require(translate(modified, policy)[0] == mapped, 'modified_controls')
    output.mkdir(mode=0o700)
    reference = output / 'reference'
    reference.mkdir()
    deck = reference / 'input.i'
    deck.write_text(yaml.safe_dump(modified, default_flow_style=False))
    case = dict(case, deck_sha256=digest(deck))
    args = list(case['arguments'])
    args[args.index('--mesh') + 1] = str(mesh)
    case['arguments'] = args
    write_json(output / 'case.json', case)
    write_json(output / 'pair.json', dict(schema=SCHEMA, steps=STEPS, mesh=str(mesh), mesh_sha256=digest(mesh),
        original_deck_sha256=digest(decks[0]), original_case_sha256=digest(case_path),
        deck_sha256=digest(deck), case_sha256=digest(output / 'case.json'),
        initialization='declared_zero_fields_no_restart',
        changed_controls=['mesh_path_spelling', 'iteration_limits', 'output_control']))


def pair_inputs(pair):
    record = read_json(pair / 'pair.json')
    require(record['schema'] == SCHEMA and record['steps'] == STEPS)
    require(record['case_sha256'] == digest(pair / 'case.json'))
    require(record['deck_sha256'] == digest(pair / 'reference/input.i'))
    require(record['mesh_sha256'] == digest(Path(record['mesh'])))
    case = read_json(pair / 'case.json')
    mapped, mesh = translate(load_deck((pair / 'reference/input.i').read_bytes()),
                             case.get('pressure_linear_policy', 'mars'))
    values = options(case['arguments'])
    require(Path(mesh).resolve() == Path(record['mesh']) == Path(values['--mesh']))
    require(values == dict(options(mapped), **{'--mesh': mesh, '--mesh-format': 'exodus',
                                             '--reference-length': values['--reference-length']}))
    return record


def solver_arguments(pair, solver):
    if solver == 'openaccel':
        return ['-i', 'input.i']
    return read_json(pair / 'case.json')['arguments'] + [
        '--iterations', str(STEPS), '--snapshot-iterations', str(STEPS), '--report-every', '1',
        '--residual-tol', '1e-6', '--mass-tol', '1e-6', '--change-tol', '1e-6',
        '--linear-cache', '1', '--halo-overlap', '1', '--field-output', 'distributed',
        '--profile', '0', '--output-prefix', str(pair / 'mars/flow')]


def runtime_libraries(executable):
    try:
        result = subprocess.run(['ldd', str(executable)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except OSError:
        raise EvidenceError('runtime_library_probe')
    text = result.stdout.decode('utf-8', errors='replace')
    require('not found' not in text, 'runtime_libraries_unresolved')
    require(result.returncode == 0, 'runtime_library_probe')
    paths = re.findall(r'(?:=>\s+|^\s*)(/\S+)\s+\(', text, re.M)
    require(paths, 'runtime_library_paths')
    try:
        return dict((p, digest(Path(p))) for p in paths)
    except OSError:
        raise EvidenceError('runtime_library_unreadable')


def run_files(directory, solver):
    files = [directory / 'run.log', directory / 'run.exit']
    if solver == 'openaccel':
        try:
            files += [directory / 'input.i'] + result_files(directory)
        except Exception:
            raise EvidenceError('reference_outputs')
    else:
        files += sorted(directory.glob('flow-*.csv')) + sorted(directory.glob('flow-*.json'))
    require(all(p.is_file() for p in files))
    return dict((p.name, digest(p)) for p in files)


class ReferenceLogState:
    # Ordered milestones from the pinned OpenAccel source; never export captured text.
    stages = (
        ('banner', r'^.*\bOpenAccel 3D\b'),
        ('controls_read', r'^Reading controls \.\.$'),
        ('controls_ready', r'^Finished reading controls \.\.$'),
        ('mesh_read', r'^Reading mesh \.\.$'),
        ('mesh_validation', r'^Validating YAML input against Exodus file$'),
        ('mesh_validated', r'^Finished validating YAML input$'),
        ('zone_registration', r'^Registering zones$'),
        ('part_registration', r'^Registering mesh parts$'),
        ('mesh_ready', r'^Finished reading mesh \.\.$'),
        ('domain_setup', r'^Setting up simulation domains:$'),
        ('equation_initialization', r'^Initializing equation `'),
        ('iteration', r'^Iter = [0-9]+$'),
        ('complete', r'^.*\bSimulation is complete\b'),
    )
    categories = (
        ('yaml', r'yaml-cpp|YAML::'),
        ('allocation', r'std::bad_alloc|std::bad_array_new_length|out of memory|cannot allocate memory'),
        ('mpi_initialization', r'MPI_Init(?:_thread)?|MPIR_Init|PMPI_Init'),
        ('mpi_thread_support', r'Provided MPI thread-level support is not sufficient'),
        ('kokkos', r'Kokkos'),
        ('mesh_io', r'Ioss::|IOSS ERROR|Exodus|ex_open|netCDF|netcdf'),
        ('decomposition', r'decomposition|decompos|Zoltan|METIS|ParMETIS'),
        ('file_access', r'No such file or directory|Permission denied|could not open|cannot open|unable to open|failed to open|does not exist'),
        ('disk_space', r'No space left on device|Disk quota exceeded'),
        ('input_validation', r'invalid boundary part|invalid side[12] part|Mesh dimension mismatch|not provided in the yaml input file'),
        ('assertion', r'Assertion .* failed|assertion .* failed|Requirement\('),
        ('linear_solver', r'Belos::|Tpetra::|Amesos2::|Ifpack2::|MueLu::|HYPRE ERROR'),
    )

    def __init__(self):
        self.seen_stages = set()
        self.seen_categories = set()
        self.exception_seen = False
        self.what_seen = False
        self.abort_seen = False
        self.in_exception = False

    def feed(self, line):
        line = re.sub(r'^\[[0-9]+\]\s*', '', line.strip())
        progress = False
        for label, pattern in self.stages:
            if re.search(pattern, line):
                self.seen_stages.add(label)
                progress = True
        exception = bool(re.match(r'^(?:terminate called|libc\+\+abi: terminating)', line))
        what = bool(re.match(r'^what\(\)\s*:', line))
        self.exception_seen |= exception
        self.what_seen |= what
        self.in_exception |= exception or what
        self.abort_seen |= bool(re.search(r'\b(?:SIGABRT|Aborted)\b', line))
        # Do not classify routine mesh/library banners as failures. Multiline what()
        # messages stay local; only matches to these fixed categories leave the log.
        error_line = bool(re.match(r'^(?:ERROR\b|Error\b|IOSS ERROR\b|Kokkos.*(?:Error|error)|Assertion\b)', line))
        if (self.in_exception or error_line) and not progress:
            for label, pattern in self.categories:
                if re.search(pattern, line, re.I):
                    self.seen_categories.add(label)

    def result(self):
        stages = [label for label, _ in self.stages if label in self.seen_stages]
        return dict(reference_stages_seen=stages,
                    reference_progress_scope='any_logged_rank',
                    reference_cpp_termination_seen=self.exception_seen,
                    reference_exception_message_seen=self.what_seen,
                    reference_abort_seen=self.abort_seen,
                    reference_error_categories=sorted(self.seen_categories))


def launch(pair, solver, executable, ranks, launcher):
    pair = pair.resolve()
    pair_inputs(pair)
    require(1 <= ranks <= 4 and launcher)
    executable = executable.resolve()
    require(executable.is_file() and os.access(str(executable), os.X_OK))
    directory = pair / ('reference' if solver == 'openaccel' else 'mars')
    if solver == 'mars':
        directory.mkdir()
    require(not (directory / 'run.log').exists() and not (directory / 'launch.json').exists())
    environment = dict(os.environ)
    # Old public instrumentation must never capture this potentially private case.
    for key in ('MARS_OPENACCEL_EXPORT_DIR', 'MARS_OPENACCEL_PUBLIC_FIXTURE'):
        environment.pop(key, None)
    if solver == 'openaccel':
        environment.update(OMP_NUM_THREADS='1', OMP_PROC_BIND='close', OMP_PLACES='cores')
    command = list(launcher) + [str(executable)] + solver_arguments(pair, solver)
    record = dict(schema=SCHEMA, solver=solver, ranks=ranks, command=command,
                  executable=str(executable), executable_sha256=digest(executable),
                  pair_sha256=digest(pair / 'pair.json'), libraries=runtime_libraries(executable),
                  environment={k: v for k, v in environment.items() if k.startswith(
                      ('MARS_', 'HYPRE_', 'CUDA_', 'MPICH_', 'OMP_', 'SLURM_'))
                      or k in ('PATH', 'LD_LIBRARY_PATH', 'UENV_VIEW')},
                  status='started')
    write_json(directory / 'launch-start.json', record)
    with (directory / 'run.log').open('xb') as log:
        child = subprocess.Popen(command, cwd=str(directory), env=environment,
                                 stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
        try:
            code = child.wait()
        except BaseException:
            child.terminate()
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait()
            raise
    code = code if code >= 0 else 128 - code
    (directory / 'run.exit').write_text(str(code) + '\n')
    record.update(status='finished', exit_code=code)
    if code not in ((0,) if solver == 'openaccel' else (0, 2)):
        # A failed process may produce no Exodus files. Preserve its exit first.
        record.update(status='failed', files={name: digest(directory / name) for name in ('run.log', 'run.exit')})
        write_json(directory / 'launch.json', record)
        raise EvidenceError('launcher_exit')
    require(digest(executable) == record['executable_sha256'])
    require(runtime_libraries(executable) == record['libraries'])
    pair_inputs(pair)
    record['files'] = run_files(directory, solver)
    write_json(directory / 'launch.json', record)
    require(code in ((0,) if solver == 'openaccel' else (0, 2)), 'launcher_exit')


def inspect_launch(pair, solver, executable):
    """Inspect a saved attempt without launching or changing its solver files."""
    directory = pair / ('reference' if solver == 'openaccel' else 'mars')
    result = dict(schema='mars-simple-startup-inspection-v1',
                  launch_start_present=(directory / 'launch-start.json').is_file(),
                  log_present=(directory / 'run.log').is_file(),
                  exit_present=(directory / 'run.exit').is_file(),
                  launch_record_present=(directory / 'launch.json').is_file(),
                  process_exit_code=None, input_check='not_checked',
                  runtime_check='not_checked', runtime_check_scope='current_environment',
                  outputs_check='not_checked')
    def check(label, function):
        try:
            function()
            result[label] = 'passed'
        except Exception as error:
            result[label] = (str(error) if isinstance(error, EvidenceError)
                             and str(error) in PREPARATION_ERRORS else 'rejected')
            result[label + '_exception'] = next((name for kind, name in (
                (ImportError, 'import_error'), (FileNotFoundError, 'file_missing'),
                (PermissionError, 'permission_denied'), (OSError, 'os_error'),
                (ValueError, 'value_error'), (KeyError, 'key_error'), (TypeError, 'type_error'))
                if isinstance(error, kind)), 'other')
    if result['exit_present']:
        def exit_code():
            raw = (directory / 'run.exit').read_text().strip()
            require(re.fullmatch(r'[0-9]{1,3}', raw) and 0 <= int(raw) <= 255)
            result['process_exit_code'] = int(raw)
        check('exit_check', exit_code)
    if result['log_present']:
        from simple_public_diagnostics import DiagnosticState
        def log_flags():
            state = DiagnosticState()
            reference = ReferenceLogState() if solver == 'openaccel' else None
            result['library_load_error_seen'] = False
            with (directory / 'run.log').open(errors='replace') as stream:
                for line in stream:
                    state.feed(line)
                    if reference is not None:
                        reference.feed(line)
                    result['library_load_error_seen'] |= 'error while loading shared libraries:' in line
            diagnostics = state.result(str(result['process_exit_code']) if result['process_exit_code'] is not None else '')
            for key in ('mpi_abort_seen', 'scheduler_time_limit_seen', 'scheduler_out_of_memory_seen',
                        'scheduler_signal_seen', 'application_error_seen'):
                result[key] = diagnostics[key]
            if reference is not None:
                result.update(reference.result())
        check('log_check', log_flags)
    check('input_check', lambda: pair_inputs(pair))
    result['executable_available'] = executable.is_file() and os.access(str(executable), os.X_OK)
    if result['executable_available']:
        check('runtime_check', lambda: runtime_libraries(executable))
    if result['exit_present']:
        check('outputs_check', lambda: run_files(directory, solver))
    # Absence of an exit file is not proof that no job started or is still running.
    return result


def verified_launch(pair, solver):
    directory = pair / ('reference' if solver == 'openaccel' else 'mars')
    record = read_json(directory / 'launch.json')
    require(record['schema'] == SCHEMA and record['solver'] == solver and record['status'] == 'finished')
    require(type(record['ranks']) is int and 1 <= record['ranks'] <= 4)
    require(record['pair_sha256'] == digest(pair / 'pair.json'))
    require(record['files'] == run_files(directory, solver))
    require(record['exit_code'] == int((directory / 'run.exit').read_text()))
    require(record['exit_code'] in ((0,) if solver == 'openaccel' else (0, 2)))
    start = read_json(directory / 'launch-start.json')
    require(start['status'] == 'started')
    require(all(record[key] == value for key, value in start.items() if key != 'status'))
    suffix = [record['executable']] + solver_arguments(pair.resolve(), solver)
    require(record['command'][-len(suffix):] == suffix)
    require(bool(record['libraries']) and re.fullmatch('[0-9a-f]{64}', record['executable_sha256']))
    return record


def compare(pair, public):
    import numpy as np
    from netCDF4 import Dataset
    pair = pair.resolve()
    public['failed_check'] = 'input_identity'
    inputs = pair_inputs(pair)
    public['failed_check'] = 'launch_records'
    mars_record = verified_launch(pair, 'mars')
    reference_record = verified_launch(pair, 'openaccel')
    public['fresh_launch_records_verified'] = True
    public['declared_zero_initialization_verified'] = True
    public['failed_check'] = 'completion_and_controls'
    mars_log = (pair / 'mars/run.log').read_text(errors='replace').splitlines()
    reference_log = (pair / 'reference/run.log').read_text(errors='replace')
    endings = [COMPLETION.fullmatch(line.strip()) for line in mars_log
               if line.startswith(('CONVERGED', 'NOT CONVERGED'))]
    require(len(endings) == 1 and endings[0] and int(endings[0].group(2)) == STEPS
            and int(endings[0].group(3)) == mars_record['ranks'])
    require(mars_record['exit_code'] == (0 if endings[0].group(1) == 'CONVERGED' else 2))
    headers = [HEADER.fullmatch(line.strip()) for line in mars_log if line.startswith('SIMPLE Tet4,')]
    require(len(headers) == 1 and headers[0] and int(headers[0].group(1)) == mars_record['ranks'])
    prepared, _, matched, exact = controls(pair / 'case.json', pair / 'reference', Path(inputs['mesh']), mars_log)
    require(exact and matched == 'mapped_controls_match')
    require([int(x) for x in re.findall(r'^Iter = (\d+)\s*$', reference_log, re.M)] == list(range(1, STEPS + 1)))
    require('Simulation is complete' in reference_log)
    reference_paths = result_files(pair / 'reference')
    require(len(reference_paths) == reference_record['ranks'])
    with (pair / 'mars/flow-metrics.csv').open() as stream:
        metrics = list(csv.DictReader(stream))
    require([int(m['iteration']) for m in metrics] == list(range(STEPS + 1)))
    for row in metrics:
        require(all(math.isfinite(float(value)) for value in row.values()))
    public['mapped_controls_verified'] = True
    public['failed_check'] = 'source_nodes'
    with Dataset(inputs['mesh']) as ds:
        count = len(ds.dimensions['num_nodes'])
        require(count > 0)
        ids, xyz = node_ids(ds, count), coordinates(ds)
    require(xyz.shape == (count, 3))
    tol = 64 * np.finfo(float).eps * max(1., float(np.max(np.abs(xyz))))
    u, rho = float(prepared['--inlet-velocity']), float(prepared['--rho'])
    scales = np.array([u, u, u, rho*u*u])
    require(np.all(np.isfinite(scales)) and np.all(scales > 0))
    rows = []
    first = None
    first_fields = []
    for iteration in range(STEPS + 1):
        public['failed_check'] = 'snapshot_coverage_or_mapping'
        mars, _ = mars_fields(pair / ('mars/flow-step-' + str(iteration)), xyz, tol, mars_record['ranks'])
        peak = float(np.max(np.sqrt(np.sum(mars[:, :3]**2, axis=1))))
        require(math.isclose(peak, float(metrics[iteration]['umax_m_s']), rel_tol=1e-12, abs_tol=1e-12*u))
        if iteration == STEPS:
            final, _ = mars_fields(pair / 'mars/flow', xyz, tol, mars_record['ranks'])
            require(np.array_equal(mars, final))
        reference = reference_fields(reference_paths, ids, xyz, iteration, tol, scales)
        public['failed_check'] = 'field_arithmetic'
        errors = field_errors((mars - reference) / scales)
        require(all(math.isfinite(value) for value in errors.values()))
        failed = [name for name in ('velocity', 'pressure') if errors[name + '_max_scaled'] > 1e-5]
        if iteration == 0:
            public['initial_field_parity_verified'] = not failed
        if failed and first is None:
            first, first_fields = iteration, failed
        rows.append(dict(iteration=iteration, errors=errors))
    write_json(pair / 'comparison-private.json', dict(schema=SCHEMA, snapshots=rows,
        first_differing_iteration=first, first_differing_fields=first_fields,
        field_tolerance=1e-5, pressure_gauge_shift_applied=False,
        mars_launch_sha256=digest(pair / 'mars/launch.json'),
        reference_launch_sha256=digest(pair / 'reference/launch.json')))
    public.update(comparison_status='completed', failed_check='none',
                  first_differing_iteration=first, first_differing_fields=first_fields,
                  all_twenty_snapshots_match=first is None)


def main(argv=None):
    parser = SafeParser(description=__doc__)
    sub = parser.add_subparsers(dest='action')
    prep = sub.add_parser('prepare')
    for name in ('case', 'reference-dir', 'output-dir'):
        prep.add_argument('--' + name, type=Path, required=True)
    run = sub.add_parser('run')
    run.add_argument('--pair', type=Path, required=True)
    run.add_argument('--solver', choices=('openaccel', 'mars'), required=True)
    run.add_argument('--executable', type=Path, required=True)
    run.add_argument('--ranks', type=int, required=True)
    run.add_argument('launcher', nargs=argparse.REMAINDER)
    check = sub.add_parser('compare')
    check.add_argument('--pair', type=Path, required=True)
    check.add_argument('--output', type=Path, required=True)
    inspect = sub.add_parser('inspect')
    inspect.add_argument('--pair', type=Path, required=True)
    inspect.add_argument('--solver', choices=('openaccel', 'mars'), required=True)
    inspect.add_argument('--executable', type=Path, required=True)
    inspect.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    os.umask(0o077)
    try:
        if args.action == 'prepare':
            prepare(args.case, args.reference_dir, args.output_dir.resolve())
            print('Startup preparation complete. All launch files are private; no solver launched.')
        elif args.action == 'run':
            launcher = args.launcher[1:] if args.launcher[:1] == ['--'] else args.launcher
            launch(args.pair, args.solver, args.executable, args.ranks, launcher)
            print('Startup capture complete. Detailed logs and fields stay private.')
        elif args.action == 'compare':
            public = dict(schema=SCHEMA, comparison_status='invalid_evidence', failed_check='inputs',
                          fresh_launch_records_verified=False, mapped_controls_verified=False,
                          declared_zero_initialization_verified=False, initial_field_parity_verified=False,
                          identical_linear_solvers_verified=False, nonlinear_convergence_required=False)
            with args.output.open('x') as stream:
                try:
                    compare(args.pair, public)
                except Exception:
                    pass  # Paths, geometry and numerical values never enter the public report.
                json.dump(public, stream, indent=2, sort_keys=True, allow_nan=False)
                stream.write('\n')
            print('Startup comparison written. Share only the public JSON.')
            return 0 if public['comparison_status'] == 'completed' else 1
        elif args.action == 'inspect':
            with args.output.open('x') as stream:
                result = inspect_launch(args.pair.resolve(), args.solver, args.executable.resolve())
                json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
                stream.write('\n')
            print('Startup inspection written. No solver launched; share only the public JSON.')
        else:
            parser.error('subcommand required')
    except Exception as error:
        label = str(error) if isinstance(error, EvidenceError) and str(error) in PREPARATION_ERRORS else 'private_preflight_or_capture'
        print('ERROR: startup probe failed (' + label + '); inspect the private launch files locally.', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
