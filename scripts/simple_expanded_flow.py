#!/usr/bin/env python3
"""Resume a saved short SIMPLE case with retained pressure; keep case data private."""
import os
from pathlib import Path
import re
import sys

import run_simple_diagnostics as diagnostics
import simple_pressure_replay as replay
import simple_pressure_profile as profile
from simple_public_diagnostics import SafeParser
from simple_snapshot_compare import options, require

SCHEMA = 'mars-simple-expanded-flow-v1'


def flow_arguments(saved, values, output):
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
    return result


def prepare(capture_dir, profile_pair, executable, output, public):
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
    arguments = flow_arguments(saved, values, output)
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
    launch = list(launcher) + [str(executable)] + [word for item in arguments.items() for word in item]
    record = dict(schema=SCHEMA, command=launch, inputs=inputs, ranks=old['ranks'],
                  steps=case['steps'], environment=replay.probe.solver_environment(environment),
                  pressure_target_preserved=True, saved_gpu_profile_verified=True)
    replay.startup.write_json(output / 'launch-start.json', record)
    public.update(pressure_target_preserved=True, saved_gpu_profile_verified=True,
                  saved_case_verified=True, launcher_preserved=True,
                  loaded_libraries_preflight_verified=True)
    return record, environment


def run(capture_dir, profile_pair, executable, output, public):
    record, environment = prepare(capture_dir, profile_pair, executable, output, public)
    public['failed_check'] = 'launch'
    code = diagnostics.capture(record['command'], output / 'run.log', output / 'run.exit',
                               output / 'diagnostics.json', environment=environment)
    public['diagnostics'] = replay.startup.read_json(output / 'diagnostics.json')
    public['failed_check'] = 'inputs_changed'
    replay.verify(record['inputs'])
    public['inputs_unchanged'] = True
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
    record.update(exit_code=code, outputs=replay.hashes(
        [output / name for name in ('run.log', 'run.exit', 'diagnostics.json', 'launch-start.json')]))
    replay.startup.write_json(output / 'launch.json', record)
    public.update(short_run_completed=valid, nonlinear_convergence_reported=status == 'converged',
                  comparison_status='completed' if valid else 'run_failed',
                  failed_check='none' if valid else 'completion')
    return 0 if valid else 1


def main(argv=None):
    parser = SafeParser(description=__doc__)
    parser.add_argument('--capture-run', type=Path, required=True)
    parser.add_argument('--gpu-profile-pair', type=Path, required=True)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args(argv)
    old_umask = os.umask(0o077)
    public = dict(schema=SCHEMA, comparison_status='invalid_evidence', failed_check='output_directory',
                  short_run_completed=False, nonlinear_convergence_reported=False,
                  openaccel_field_parity_verified=False)
    try:
        output = args.output_dir.resolve()
        output.mkdir(mode=0o700)
    except Exception:
        os.umask(old_umask)
        print('ERROR: a fresh private output directory is required.', file=sys.stderr)
        return 1
    try:
        code = run(args.capture_run.resolve(), args.gpu_profile_pair.resolve(),
                   args.executable.resolve(), output, public)
    except KeyboardInterrupt:
        public.update(comparison_status='interrupted', failed_check='interrupted')
        code = 130
    except Exception:
        # Exception text can contain private boundary names, fields or paths.
        code = 1
    finally:
        diagnostics.publish(output / 'public.json', public)
        os.umask(old_umask)
    print('Short-flow status written. Share only public.json; completion is not field parity.')
    return code


if __name__ == '__main__':
    sys.exit(main())
