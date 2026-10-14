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


def debugger(gdb, command, directory, rank, ranks, environment, self_test=False, communication_progress=False):
    config = directory / ('gdb-{}.json'.format(rank))
    write_json(config, dict(directory=str(directory), rank=rank, ranks=ranks, self_test=self_test,
                           communication_progress=communication_progress))
    env = dict(environment, SIMPLE_SETUP_STACK_CONFIG=str(config))
    argv = [gdb, '-nx', '-nh', '-batch', '-iex', 'set auto-load off',
            '-x', str(GDB_SCRIPT), '--args'] + command
    with (directory / ('debugger-{}.log'.format(rank))).open('xb') as log:
        child = subprocess.Popen(argv, cwd=str(directory), env=env, stdin=subprocess.DEVNULL,
                                 stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        begin = time.monotonic()
        deadline, hard_deadline = begin + 30, begin + (20 if self_test else 300)
        running = entered = interrupted = captured = resumed = False
        try:
            while child.poll() is None:
                now = time.monotonic()
                if not running and (directory / ('running-{}'.format(rank))).exists():
                    running = True
                    deadline = now + (1 if self_test else 120)
                if not entered and (directory / ('entered-{}'.format(rank))).exists():
                    entered = True
                    deadline = now + 30
                if communication_progress and not resumed and (directory / ('resumed-{}'.format(rank))).exists():
                    resumed = True
                    interrupted = False
                    deadline = now + 15
                first_saved = communication_progress and (directory / ('first-{}.json'.format(rank))).is_file()
                captured |= (directory / ('stack-{}.json'.format(rank))).is_file()
                if now >= hard_deadline or (interrupted and not captured and (not first_saved or resumed)
                                           and now >= deadline + 20):
                    break
                if not interrupted and not captured and not (first_saved and not resumed) and now >= deadline:
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
                      rank, config['ranks'], os.environ,
                      communication_progress=config.get('communication_progress', False))
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
                                ('first-{}.json', 'first communication sample saved'),
                                ('resumed-{}', 'resumed for second communication sample'),
                                ('stack-{}.json', 'stack saved')):
                path = directory / name.format(rank)
                if path not in seen and path.is_file():
                    print('MARS rank {}: {}'.format(rank, label), flush=True)
                    seen.add(path)
        if code is not None:
            return code


def communication_report(path, final_stack):
    value = read_json(path)
    if (not isinstance(value, dict) or set(value) != {'first', 'final', 'two_samples_taken', 'peer_samples_present'}
            or any(type(value[k]) is not bool for k in ('two_samples_taken', 'peer_samples_present'))
            or not isinstance(value['first'], dict) or set(value['first']) != {'stack', 'exchange'}):
        raise ValueError('invalid communication report')
    first_stack = safe_report(value['first']['stack'])
    first, last = value['first']['exchange'], value['final']
    for counts in (first, last):
        if (not isinstance(counts, dict) or set(counts) != {'entered', 'returned', 'tracking_complete'}
                or type(counts['tracking_complete']) is not bool
                or any(type(counts[k]) is not int or counts[k] < 0 for k in ('entered', 'returned'))
                or counts['returned'] > counts['entered']):
            raise ValueError('invalid communication counters')
    if any(last[k] < first[k] for k in ('entered', 'returned')):
        raise ValueError('communication counters moved backwards')
    usable = value['two_samples_taken'] and all(
        s['inferior_alive'] and s['main_thread_stopped'] and s['main_thread_frame_readable']
        and not s['debugger_error'] and not s['unexpected_signal_seen'] and not s['unwind_incomplete']
        and (s['interrupt_seen'] or s['setup_return_seen']) for s in (first_stack, final_stack))
    tracked = first['tracking_complete'] and last['tracking_complete']
    return dict(first_main_thread_functions=first_stack['main_thread_functions'],
                two_samples_taken=value['two_samples_taken'], two_samples_usable=usable,
                peer_samples_present=value['peer_samples_present'], exchange_tracking_complete=tracked,
                exchange_entry_seen=last['entered'] > 0, exchange_return_seen=last['returned'] > 0,
                exchange_advanced_between_samples=(any(last[k] > first[k] for k in ('entered', 'returned'))
                                                   if tracked and usable else None),
                exchange_active_at_first_sample=(first['entered'] > first['returned'] if tracked else None),
                exchange_active_at_final_sample=(last['entered'] > last['returned'] if tracked else None))


def run(pair, output, public, communication_progress=False):
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
    # Reading stopped frames does not require unwinding the whole program stack.
    check = debugger(gdb, ['/bin/sleep', '60'], test, 0, 1, environment, self_test=True)
    public['debugger_self_test_stack'] = check
    public['debugger_self_test_scope'] = 'stopped_named_frames_not_complete_unwind'
    public['debugger_self_test_passed'] = (check is not None and check['interrupt_seen']
                                          and check['inferior_alive'] and not check['debugger_error']
                                          and not check['unexpected_signal_seen']
                                          and check['main_thread_stopped']
                                          and check['main_thread_frame_readable']
                                          and not (set(check['stack_errors']) - {'unwind_stopped'})
                                          and (not check['unwind_incomplete'] or bool(check['stack_errors'])))
    if not public['debugger_self_test_passed']:
        return 1
    if check['unwind_incomplete']:
        print('GDB can read stopped frames; the self-test unwind is partial.', flush=True)
    config = directory / 'launch-input.json'
    write_json(config, dict(gdb=gdb, solver_command=command, directory=str(directory), ranks=ranks,
                           pair=str(pair), executable_identity=identity,
                           communication_progress=communication_progress))
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
    if communication_progress:
        public.update(failed_check='communication_reports', communication_progress={})
        for rank, stack in public['rank_stacks'].items():
            public['communication_progress'][rank] = communication_report(
                directory / ('exchange-{}.json'.format(rank)), stack)
        public['communication_scope'] = 'host_exchange_entry_return_between_samples_not_deadlock_proof'
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
    elif communication_progress and any(
            not c['exchange_tracking_complete'] or not c['peer_samples_present']
            or (not c['two_samples_usable'] and not any(r['setup_return_seen'] for r in reports))
            for c in public['communication_progress'].values()):
        public.update(comparison_status='diagnostic_partial', failed_check='communication_tracking')
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
    parser.add_argument('--communication-progress', action='store_true',
                        help='Track host exchange entry/return and take a second sample after all ranks resume.')
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
        status = run(args.pair.resolve(), output, public, args.communication_progress)
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
