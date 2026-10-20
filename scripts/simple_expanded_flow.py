#!/usr/bin/env python3
"""Run a saved SIMPLE case with retained pressure; keep case data private."""
import math
import os
from pathlib import Path
import re
import shlex
import sys

import run_simple_diagnostics as diagnostics
import simple_convergence_summary as convergence
import simple_pressure_replay as replay
import simple_pressure_profile as profile
from simple_public_diagnostics import SafeParser, summarize
from simple_snapshot_compare import options, require

SCHEMA = 'mars-simple-expanded-flow-v1'
TARGETS = ('--residual-tol', '--mass-tol', '--change-tol')


def extended_launcher(saved, time_limit):
    result = list(saved)
    if time_limit is None:
        return result
    require(result and Path(result[0]).name == 'srun')
    require(re.fullmatch(r'[0-9]{2}:[0-5][0-9]:[0-5][0-9]', time_limit)
            and time_limit != '00:00:00')
    # Stop at the binding program: its arguments are not scheduler options.
    index, found = 1, False
    valued = {'--account', '-A', '--time', '-t', '--nodes', '-N', '--ntasks', '-n',
              '--ntasks-per-node', '--cpus-per-task', '-c', '--cpu-bind', '--export',
              '--kill-on-bad-exit', '--gpus-per-node', '--gpus-per-task', '--partition', '-p'}
    while index < len(result) and result[index].startswith('-'):
        token = result[index]
        if token in ('--time', '-t'):
            require(not found and index + 1 < len(result))
            result[index + 1] = time_limit
            found = True
            index += 2
        elif token.startswith('--time=') or (token.startswith('-t') and not token.startswith('--')):
            require(not found)
            result[index] = '--time=' + time_limit
            found = True
            index += 1
        elif '=' in token:
            index += 1
        elif token in valued:
            require(index + 1 < len(result))
            index += 2
        else:
            # Do not guess whether an unknown flag consumes a following argument.
            require(token in ('--exclusive', '--overlap', '--unbuffered', '-u'))
            index += 1
    if not found:
        result.insert(1, '--time=' + time_limit)
    return result


def flow_arguments(saved, values, output, iterations=None):
    target = (float(saved['--pressure-linear-atol']), float(saved['--pressure-linear-rtol']))
    require(target == (values['atol'], values['rtol']))
    require(saved.get('--pressure-expansion', '0') == '0'
            and saved.get('--pressure-refinement', '0') == '0'
            and saved.get('--first-step-audit', '0') == '0'
            and saved.get('--setup-only', '0') == '0')
    result = dict(saved)
    result.pop('--pressure-failure-capture', None)
    result.update({'--pressure-expansion': '1', '--snapshot-iterations': '0',
                   '--field-output': 'none', '--output-prefix': str(output / 'flow'),
                   '--pressure-solver-profile': str(output / 'pressure.profile')})
    if iterations is not None:
        require(type(iterations) is int and int(saved['--iterations']) <= iterations < 2**31)
        result['--iterations'] = str(iterations)
    return result


def prepare(capture_dir, profile_pair, executable, output, public, iterations=None, time_limit=None):
    public['failed_check'] = 'saved_capture'
    old = replay.captured_record(capture_dir)
    pair = Path(old['pair'])
    public['failed_check'] = 'saved_case'
    case = replay.startup.pair_inputs(pair)
    require(not case.get('first_step_audit', False))
    saved = options(replay.startup.solver_arguments(pair, 'mars'))
    expected = dict(saved, **{'--output-prefix': str(capture_dir / 'flow'),
                             '--field-output': 'none', '--snapshot-iterations': '0',
                             '--pressure-failure-capture': str(capture_dir / 'system')})
    command = old['command']; launcher = old['launcher']
    require(launcher and command[:len(launcher)] == launcher
            and options(command[len(launcher)+1:]) == expected)
    public['failed_check'] = 'saved_capture_inputs'
    old_executable = str(Path(command[len(launcher)]).resolve())
    require(old_executable in old['inputs'])
    # Rebuilding the executable is intended; every other captured input must match.
    inputs = {key: value for key, value in old['inputs'].items() if key != old_executable}
    public['saved_input_checks'] = replay.file_identity_checks(inputs,
        lambda key: 'runtime_library' if key in old['libraries'] else 'saved_case_or_launch')
    require(public['saved_input_checks']['matched'])
    require(int(saved['--iterations']) == case['steps'] == replay.startup.STEPS)
    public['failed_check'] = 'gpu_profile'
    values, dependencies = replay.gpu_configuration(profile_pair, capture_dir, old)
    arguments = flow_arguments(saved, values, output, iterations)
    targets = {key: float(arguments[key]) for key in TARGETS}
    require(all(math.isfinite(value) and value > 0 for value in targets.values()))
    public['failed_check'] = 'launcher_controls'
    new_launcher = extended_launcher(launcher, time_limit)
    public['failed_check'] = 'solver_environment_changed'
    public['environment_check'] = replay.probe.environment_report(
        replay.probe.environment_changes(os.environ, old['environment']))
    environment = replay.probe.restore_environment(os.environ, old['environment'])
    require(not replay.probe.environment_changes(environment, old['environment']))
    require(environment.get('MARS_HYPRE_SPMV_VENDOR', '0') == '0')
    public['failed_check'] = 'runtime_libraries'
    require(executable.is_file() and os.access(str(executable), os.X_OK))
    libraries = replay.startup.runtime_libraries(executable)
    require(libraries == old['libraries'])
    (output / 'pressure.profile').write_text(profile.text_profile(values))
    inputs.update(replay.hashes([executable, capture_dir / 'capture.json', pair / 'pair.json',
                           pair / 'case.json', pair / 'reference/input.i', Path(case['mesh']),
                           output / 'pressure.profile']))
    inputs.update(dependencies); inputs.update(libraries)
    launch = new_launcher + [str(executable)] + [word for item in arguments.items() for word in item]
    record = dict(schema=SCHEMA, command=launch, inputs=inputs, ranks=old['ranks'],
                  steps=int(arguments['--iterations']), source_steps=case['steps'], record_version=2,
                  nonlinear_targets=targets, launcher=new_launcher, source_launcher=launcher,
                  time_limit_override=time_limit,
                  environment=replay.probe.solver_environment(environment),
                  pressure_target_preserved=True, saved_gpu_profile_verified=True)
    replay.startup.write_json(output / 'launch-start.json', record)
    # Keep the exact direct launch visible to the user without publishing case arguments.
    (output / 'launch-command.txt').write_text(' '.join(shlex.quote(word) for word in launch) + '\n')
    public.update(pressure_target_preserved=True, saved_gpu_profile_verified=True,
                  saved_case_verified=True, launcher_preserved=new_launcher == launcher,
                  launcher_resources_preserved_except_time=True,
                  iteration_budget_extended=record['steps'] > case['steps'],
                  loaded_libraries_preflight_verified=True, fresh_zero_start=True)
    return record, environment


def run(capture_dir, profile_pair, executable, output, public, iterations=None, time_limit=None):
    record, environment = prepare(capture_dir, profile_pair, executable, output, public, iterations, time_limit)
    public['failed_check'] = 'launch'
    code = diagnostics.capture(record['command'], output / 'run.log', output / 'run.exit',
                               output / 'diagnostics.json', environment=environment)
    public['diagnostics'] = replay.startup.read_json(output / 'diagnostics.json')
    public['failed_check'] = 'inputs_changed'
    replay.verify(record['inputs'])
    public['inputs_unchanged'] = True
    files = [output / name for name in ('run.log', 'run.exit', 'diagnostics.json', 'launch-start.json')]
    if (output / 'flow-metrics.csv').is_file():
        files.append(output / 'flow-metrics.csv')
    record.update(exit_code=code, outputs=replay.hashes(files))
    replay.startup.write_json(output / 'launch.json', record)
    return completion(record, output, code, public)


def completion(record, output, code, public):
    public['failed_check'] = 'completion'
    status = public['diagnostics']['run_status']
    completions = re.findall(r'^(CONVERGED|NOT CONVERGED: iteration limit) iterations=(\d+) ranks=(\d+)\b',
                             (output / 'run.log').read_text(errors='replace'), re.M)
    valid = len(completions) == 1
    if valid:
        label, steps, ranks = completions[0]
        valid = int(ranks) == record['ranks'] and (
            (code == 2 and status == 'iteration_limit' and label.startswith('NOT ')
             and int(steps) == record['steps']) or
            (code == 0 and status == 'converged' and label == 'CONVERGED'
             and 0 < int(steps) <= record['steps']))
    public.update(run_completed=valid, short_run_completed=valid and record['steps'] == replay.startup.STEPS,
                  nonlinear_convergence_reported=status == 'converged',
                  nonlinear_convergence_verified=False,
                  comparison_status='completed' if valid else 'run_failed',
                  failed_check='none' if valid else 'completion')
    if valid and record.get('record_version') == 2:
        public['comparison_status'] = 'invalid_evidence'
        public['failed_check'] = 'convergence_evidence'
        metrics = output / 'flow-metrics.csv'
        require(str(metrics.resolve()) in record['outputs'])
        with (output / 'run.log').open(errors='replace') as log, metrics.open(newline='') as rows:
            result = convergence.summarize(log, rows, str(code),
                                           *(record['nonlinear_targets'][key] for key in TARGETS))
        result['targets_source'] = 'verified_launch_arguments'
        public['convergence'] = result
        require(result['completion_matches_supplied_targets'])
        public['nonlinear_convergence_verified'] = status == 'converged' and result['all_final_checks_passed']
        public['comparison_status'] = 'completed'
        public['failed_check'] = 'none'
    return 0 if valid else 1


def inspect_run(directory, public):
    public['failed_check'] = 'saved_launch'
    record = replay.startup.read_json(directory / 'launch.json')
    require(record['schema'] == SCHEMA and record['pressure_target_preserved'] and record['saved_gpu_profile_verified'])
    if record.get('record_version') == 2:
        require(record['source_steps'] == replay.startup.STEPS
                and type(record['steps']) is int and replay.startup.STEPS <= record['steps'] < 2**31)
        require(record['launcher'] == extended_launcher(record['source_launcher'], record['time_limit_override']))
        require(record['command'][:len(record['launcher'])] == record['launcher'])
        arguments = options(record['command'][len(record['launcher'])+1:])
        require(int(arguments['--iterations']) == record['steps']
                and arguments['--pressure-expansion'] == '1'
                and arguments['--output-prefix'] == str(directory / 'flow')
                and record['nonlinear_targets'] == {key: float(arguments[key]) for key in TARGETS})
    else:
        require('record_version' not in record and record['steps'] == replay.startup.STEPS)
    public['failed_check'] = 'saved_inputs_or_outputs'
    required = {str((directory / name).resolve()) for name in
                ('run.log', 'run.exit', 'diagnostics.json', 'launch-start.json')}
    require(required <= set(record['outputs']))
    replay.verify(record['inputs']); replay.verify(record['outputs'])
    start = replay.startup.read_json(directory / 'launch-start.json')
    require(start == {key: value for key, value in record.items() if key not in ('exit_code', 'outputs')})
    code = int((directory / 'run.exit').read_text().strip())
    require(code == record['exit_code'])
    with (directory / 'run.log').open(errors='replace') as log:
        public['diagnostics'] = summarize(log, str(code))
    public.update(inputs_unchanged=True, saved_outputs_verified=True, solver_launched=False)
    return completion(record, directory, code, public)


def main(argv=None):
    parser = SafeParser(description=__doc__)
    parser.add_argument('--capture-run', type=Path)
    parser.add_argument('--gpu-profile-pair', type=Path)
    parser.add_argument('--executable', type=Path)
    parser.add_argument('--inspect-run', type=Path)
    parser.add_argument('--iterations', type=int, help='Optional larger iteration budget; starts from the saved zero state')
    parser.add_argument('--time-limit', help='Optional srun wall-time override as HH:MM:SS')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    launch_args = (args.capture_run, args.gpu_profile_pair, args.executable)
    if (args.inspect_run and (any(launch_args) or args.iterations is not None or args.time_limit is not None)) or (not args.inspect_run and not all(launch_args)):
        parser.error('Choose inspection or launch arguments')
    old_umask = os.umask(0o077)
    public = dict(schema=SCHEMA, comparison_status='invalid_evidence', failed_check='output_directory',
                  run_completed=False, short_run_completed=False, nonlinear_convergence_reported=False,
                  nonlinear_convergence_verified=False,
                  openaccel_field_parity_verified=False)
    try:
        output = args.output_dir.resolve()
        output.mkdir(mode=0o700)
    except Exception:
        os.umask(old_umask)
        print('ERROR: a fresh private output directory is required.', file=sys.stderr)
        return 1
    try:
        if args.inspect_run:
            code = inspect_run(args.inspect_run.resolve(), public)
        else:
            code = run(args.capture_run.resolve(), args.gpu_profile_pair.resolve(),
                       args.executable.resolve(), output, public, args.iterations, args.time_limit)
    except KeyboardInterrupt:
        public.update(comparison_status='interrupted', failed_check='interrupted')
        code = 130
    except Exception:
        # Exception text can contain private boundary names, fields or paths.
        code = 1
    finally:
        diagnostics.publish(output / 'public.json', public)
        os.umask(old_umask)
    print('Flow status written. Share only public.json; completion is not field parity.')
    return code


if __name__ == '__main__':
    sys.exit(main())
