#!/usr/bin/env python3
"""User-run, bounded setup diagnostic. Detailed debugger output stays private."""

import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import time

from simple_setup_stack_gdb import safe_report, write_json

SCHEMA = 'mars-simple-setup-stack-v1'
HERE = Path(__file__).resolve()
GDB_SCRIPT = HERE.with_name('simple_setup_stack_gdb.py')


def read_json(path):
    return json.loads(path.read_text())


def debugger(gdb, command, directory, rank, ranks, environment, self_test=False):
    config = directory / ('gdb-{}.json'.format(rank))
    write_json(config, dict(directory=str(directory), rank=rank, ranks=ranks, self_test=self_test))
    env = dict(environment, SIMPLE_SETUP_STACK_CONFIG=str(config))
    argv = [gdb, '-nx', '-nh', '-batch', '-iex', 'set auto-load off',
            '-x', str(GDB_SCRIPT), '--args'] + command
    with (directory / ('debugger-{}.log'.format(rank))).open('xb') as log:
        child = subprocess.Popen(argv, cwd=str(directory), env=env, stdin=subprocess.DEVNULL,
                                 stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        begin = time.monotonic()
        deadline, hard_deadline = begin + 30, begin + (20 if self_test else 300)
        running = entered = interrupted = captured = False
        try:
            while child.poll() is None:
                now = time.monotonic()
                if not running and (directory / ('running-{}'.format(rank))).exists():
                    running = True
                    deadline = now + (1 if self_test else 120)
                if not entered and (directory / ('entered-{}'.format(rank))).exists():
                    entered = True
                    deadline = now + 30
                captured |= (directory / ('stack-{}.json'.format(rank))).is_file()
                if now >= hard_deadline or (interrupted and not captured and now >= deadline + 20):
                    break
                if not interrupted and not captured and now >= deadline:
                    child.send_signal(signal.SIGINT)
                    interrupted = True
                time.sleep(.1)
        finally:
            if child.poll() is None:
                # Include the inferior if GDB cannot stop or exits unexpectedly.
                os.killpg(child.pid, signal.SIGKILL)
            child.wait()
    report = directory / ('stack-{}.json'.format(rank))
    return safe_report(read_json(report)) if report.is_file() else None


def rank_run(config_path):
    config = read_json(config_path)
    rank = int(os.environ['SLURM_PROCID'])
    if not 0 <= rank < config['ranks'] <= 4:
        raise ValueError('rank outside diagnostic scope')
    report = debugger(config['gdb'], config['solver_command'], Path(config['directory']),
                      rank, config['ranks'], os.environ)
    # A deliberately killed inferior is not a failed CFD convergence test.
    return 0 if report is not None and not report['debugger_error'] else 1


def saved_command(pair, directory, public):
    import simple_pressure_profile as profile
    import simple_pressure_probe as pressure
    import simple_startup_probe as startup
    from simple_snapshot_compare import digest, require

    baseline, launches, _, identity = profile.check_inputs(pair)
    require(identity.get('pressure_amg_mode') == 'reference' and
            startup.read_json(pair / 'pair.json').get('first_step_audit') is True and
            startup.read_json(pair / 'case.json').get('momentum_backend') == 'petsc',
            'reference_amg_petsc_first_step_required')
    old = launches['mars']
    require(1 <= old['ranks'] <= 4, 'unsupported_rank_count')
    executable = Path(identity['executable'])
    require(executable.is_file() and os.access(str(executable), os.X_OK)
            and digest(executable) == identity['executable_sha256'], 'executable_changed')
    public['library_checks'] = profile.profile_library_checks(
        startup.runtime_libraries(executable), identity['libraries'], False)
    require(public['library_checks']['matched'], 'libraries_changed')
    environment = dict(os.environ)
    for key in ('MARS_OPENACCEL_EXPORT_DIR', 'MARS_OPENACCEL_PUBLIC_FIXTURE'):
        environment.pop(key, None)
    public['environment_check'] = pressure.environment_report(
        pressure.environment_changes(environment, old['environment']))
    environment = pressure.restore_environment(environment, old['environment'])
    require(not pressure.environment_changes(environment, old['environment']), 'solver_environment_changed')
    arguments = startup.solver_arguments(pair, 'mars')
    arguments[arguments.index('--output-prefix') + 1] = str(directory / 'flow')
    launcher = pressure.launcher(old, baseline, 'mars')
    require('srun' == Path(launcher[0]).name, 'saved_srun_launcher_required')
    return launcher, [str(executable)] + arguments, environment, old['ranks'], identity


def wait_with_progress(child, directory, ranks):
    seen = set()
    while True:
        try:
            code = child.wait(timeout=1)
        except subprocess.TimeoutExpired:
            code = None
        for rank in range(ranks):
            for name, label in (('running-{}', 'debugger running'), ('entered-{}', 'Hypre setup entered'),
                                ('stack-{}.json', 'stack saved')):
                path = directory / name.format(rank)
                if path not in seen and path.is_file():
                    print('MARS rank {}: {}'.format(rank, label), flush=True)
                    seen.add(path)
        if code is not None:
            return code


def run(pair, output, public):
    output.mkdir(mode=0o700)
    directory = output / 'private'
    directory.mkdir(mode=0o700)
    public.update(failed_check='saved_case_preflight')
    launcher, command, environment, ranks, identity = saved_command(pair, directory, public)
    public.update(failed_check='debugger_unavailable')
    gdb = shutil.which('gdb')
    if gdb is None:
        return 1
    test = directory / 'self-test'
    test.mkdir()
    public.update(failed_check='debugger_self_test')
    # Require readable stopped frames as well as interruption before requesting GPUs.
    check = debugger(gdb, ['/bin/sleep', '60'], test, 0, 1, environment, self_test=True)
    public['debugger_self_test_passed'] = (check is not None and check['interrupt_seen']
                                          and check['inferior_alive'] and not check['debugger_error']
                                          and not check['unexpected_signal_seen']
                                          and check['main_thread_stopped']
                                          and check['main_thread_frame_readable']
                                          and not check['unwind_incomplete'])
    if not public['debugger_self_test_passed']:
        public['debugger_self_test_stack'] = check
        return 1
    config = directory / 'launch-input.json'
    write_json(config, dict(gdb=gdb, solver_command=command, directory=str(directory), ranks=ranks,
                           pair=str(pair), executable_identity=identity))
    argv = launcher + [sys.executable, str(HERE), '_rank', str(config)]
    with (directory / 'launch-command.txt').open('x') as stream:
        stream.write(' '.join(shlex.quote(word) for word in argv) + '\n')
    public.update(failed_check='diagnostic_launch', solver_launch_requested=True,
                  solver_controls_preserved=True, saved_launcher_prefix_preserved=True,
                  output_prefix_changed=True, ranks_expected=ranks)
    print('GDB self-test passed. Launching the saved MARS case; setup sampling stops it deliberately.', flush=True)
    with (directory / 'launch.log').open('xb') as log:
        child = subprocess.Popen(argv, cwd=str(directory), env=environment, stdin=subprocess.DEVNULL,
                                 stdout=log, stderr=subprocess.STDOUT)
        try:
            code = wait_with_progress(child, directory, ranks)
        except BaseException:
            child.terminate()
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait()
            raise
    (directory / 'launch.exit').write_text(str(code) + '\n')
    public['launcher_exit_code'] = code
    public['solver_launched'] = any((directory / ('running-{}'.format(i))).exists() for i in range(ranks))
    public['rank_stacks'] = {}
    for rank in range(ranks):
        path = directory / ('stack-{}.json'.format(rank))
        if path.is_file():
            public['rank_stacks'][str(rank)] = safe_report(read_json(path))
    public['all_rank_stacks_present'] = len(public['rank_stacks']) == ranks
    public['failed_check'] = 'saved_case_recheck'
    _, _, _, _, after = saved_command(pair, directory, {})
    if after != identity:
        return 1
    reports = list(public['rank_stacks'].values())
    public['all_rank_stacks_usable'] = public['all_rank_stacks_present'] and all(
        r['inferior_alive'] and not r['debugger_error'] and not r['unexpected_signal_seen']
        and r['main_thread_stopped'] and r['main_thread_frame_readable']
        and (bool(r['main_thread_functions']) or r['setup_return_seen']) for r in reports)
    if code != 0 or not public['all_rank_stacks_present'] or any(r['debugger_error'] for r in reports):
        public.update(comparison_status='diagnostic_incomplete', failed_check='rank_stack_capture')
    elif not public['all_rank_stacks_usable']:
        public.update(comparison_status='diagnostic_incomplete', failed_check='stack_unreadable')
    elif any(r['unwind_incomplete'] for r in reports):
        public.update(comparison_status='diagnostic_partial', failed_check='stack_unwind')
    else:
        public.update(comparison_status='diagnostic_complete', failed_check='none')
        return 0
    return 1


def main(argv=None):
    os.umask(0o077)
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == '_rank':
        try:
            return rank_run(Path(argv[1]))
        except Exception:
            return 1
    from simple_public_diagnostics import SafeParser
    parser = SafeParser(description=__doc__)
    parser.add_argument('--pair', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args(argv)
    output = args.output_dir.resolve()
    public = dict(schema=SCHEMA, comparison_status='diagnostic_incomplete',
                  failed_check='output_directory', solver_launched=False,
                  scope='allowlisted_host_stack_names_not_gpu_state_or_deadlock_proof',
                  intentional_stop=True, convergence_verified=False, field_parity_verified=False)
    status = 1
    if output.exists():
        print('STOP: use a new diagnostic output directory. Existing evidence was not changed.')
        return 1
    try:
        status = run(args.pair.resolve(), output, public)
    except KeyboardInterrupt:
        public.update(failed_check='interrupted')
    except Exception:
        # Exception text can contain private command arguments and input paths.
        pass
    if output.is_dir() and not (output / 'public.json').exists():
        write_json(output / 'public.json', public)
    print('Setup diagnostic finished. Share only public.json; raw debugger output stays private.')
    return status


if __name__ == '__main__':
    sys.exit(main())
