"""Synthetic debugger and saved-case checks; no private geometry or cluster launch."""

import contextlib
import io
import json
import os
from pathlib import Path
import signal
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import simple_setup_stack as trace
import simple_setup_stack_gdb as sampler
import simple_startup_probe as startup
import test_simple_pressure_profile as fixtures


def report(**changes):
    value = dict.fromkeys(sampler.BOOL_KEYS, False)
    value.update({key: [] for key in sampler.LIST_VALUES})
    value.update(changes)
    return value


def readable_report(**changes):
    value = report(interrupt_seen=True, inferior_alive=True, main_thread_stopped=True,
                   main_thread_frame_readable=True)
    value.update(changes)
    return value


class Frame:
    def __init__(self, name, older=None, reason=0):
        self._name, self._older, self.reason = name, older, reason

    def name(self):
        return self._name

    def older(self):
        return self._older

    def unwind_stop_reason(self):
        return self.reason


class FakeGdb:
    # These are the normal Python exports in GDB 13.2; FIRST_ERROR is not exported.
    FRAME_UNWIND_NO_REASON = 0
    FRAME_UNWIND_NULL_ID = 1
    FRAME_UNWIND_OUTERMOST = 2
    FRAME_UNWIND_UNAVAILABLE = 3
    FRAME_UNWIND_INNER_ID = 4
    FRAME_UNWIND_SAME_ID = 5
    FRAME_UNWIND_NO_SAVED_PC = 6
    FRAME_UNWIND_MEMORY_ERROR = 7

    def __init__(self, mode='sample'):
        self.mode, self.commands, self.points = mode, [], []
        self.alive, self.selected = True, None
        self.events = SimpleNamespace(cont=SimpleNamespace(connect=Mock()),
                                      stop=SimpleNamespace(connect=Mock()),
                                      exited=SimpleNamespace(connect=Mock()))
        owner = self

        class Breakpoint:
            def __init__(self, *args, **kwargs):
                owner.points.append(self)

            def delete(self):
                pass

        self.Breakpoint = Breakpoint
        self.FinishBreakpoint = Breakpoint
        self.threads = [SimpleNamespace(num=i, switch=lambda i=i: setattr(self, 'selected', i),
                                       is_stopped=lambda: True)
                        for i in (1, 2)]
        self.frames = {
            1: Frame('MPI_Waitall', Frame('SECRET_PUMP_FUNCTION',
                     Frame('hypre_BoomerAMGCoarsenRuge', Frame('hypre_BoomerAMGSetup')))),
            2: Frame('SECRET_OTHER_THREAD', Frame('cudaEventSynchronize'))}

    def selected_inferior(self):
        return SimpleNamespace(pid=123 if self.alive else 0,
                               threads=lambda: self.threads if self.alive else [])

    def newest_frame(self):
        return self.frames[self.selected]

    def execute(self, command):
        self.commands.append(command)
        if command == 'run':
            self.events.cont.connect.call_args[0][0](None)
            if self.mode == 'before_setup':
                raise KeyboardInterrupt()
            if self.mode == 'error':
                raise RuntimeError('PRIVATE_PATH_AND_VALUES')
            if self.mode == 'exit':
                self.alive = False
                self.events.exited.connect.call_args[0][0](SimpleNamespace(exit_code=127))
            else:
                self.points[0].stop()
        elif command == 'continue':
            if self.mode == 'return':
                self.points[-1].stop()
            else:
                raise KeyboardInterrupt()
        elif command == 'kill':
            self.alive = False


class StackTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name).resolve()

    def test_only_literal_function_names_are_exported(self):
        gdb = FakeGdb()
        self.assertFalse(hasattr(gdb, 'FRAME_UNWIND_FIRST_ERROR'))
        gdb.frames[1]._older._older._older.reason = gdb.FRAME_UNWIND_OUTERMOST
        value = sampler.stack_functions(gdb)
        self.assertEqual(value['main_thread_functions'],
                         ['MPI_Waitall', 'hypre_BoomerAMGCoarsenRuge', 'hypre_BoomerAMGSetup'])
        self.assertEqual(value['other_thread_functions'], ['cudaEventSynchronize'])
        self.assertFalse(value['unwind_incomplete'])
        self.assertTrue(value['main_thread_stopped'])
        self.assertTrue(value['main_thread_frame_readable'])
        self.assertNotIn('SECRET', json.dumps(value))

    def test_unknown_names_and_raw_extra_keys_cannot_cross_export(self):
        for value in (report(main_thread_functions=['PRIVATE']), dict(report(), raw='PRIVATE'),
                      report(setup_entry_seen='PRIVATE'), report(other_thread_functions=['MPI_Wait PRIVATE']),
                      report(main_thread_functions=['MPI_Wait'] * 257), report(stack_errors=['PRIVATE']),
                      report(stack_errors=['MPI_Wait']), report(main_thread_stopped=1),
                      report(unwind_stop_reasons=['PRIVATE']), report(unwind_stop_reasons=[3])):
            with self.assertRaisesRegex(ValueError, 'invalid stack report'):
                sampler.safe_report(value)

    def test_broken_unwind_is_explicit(self):
        gdb = FakeGdb()
        gdb.frames[1].reason = 3
        self.assertEqual(sampler.stack_functions(gdb)['stack_errors'], ['unwind_stopped'])
        gdb.frames[1].name = Mock(side_effect=RuntimeError('SECRET'))
        value = sampler.stack_functions(gdb)
        self.assertTrue(value['unwind_incomplete'])
        self.assertIn('frame_name', value['stack_errors'])
        self.assertIn('hypre_BoomerAMGSetup', value['main_thread_functions'])
        gdb.frames[1].name.assert_called_once()
        self.assertNotIn('SECRET', json.dumps(value))

    def test_unwind_reasons_are_fixed_labels_and_legacy_end_is_normal(self):
        cases = [(name, getattr(FakeGdb, 'FRAME_UNWIND_' + name)) for name in sampler.UNWIND_ERRORS]
        cases += [('unknown', 1000), ('normal', FakeGdb.FRAME_UNWIND_NULL_ID)]
        for name, reason in cases:
            with self.subTest(name=name):
                gdb = FakeGdb()
                gdb.frames[1]._older._older._older.reason = reason
                value = sampler.stack_functions(gdb)
                self.assertEqual(value['unwind_stop_reasons'], [] if name == 'normal' else [name.lower()])
                self.assertEqual(value['unwind_incomplete'], name != 'normal')
                sampler.safe_report(report(**value))

    def test_partial_unwind_keeps_names_collected_before_exception(self):
        gdb = FakeGdb()
        gdb.frames[1]._older._older.older = Mock(side_effect=RuntimeError('PRIVATE'))
        gdb.frames[2]._older.older = Mock(side_effect=RuntimeError('PRIVATE'))
        value = sampler.stack_functions(gdb)
        self.assertEqual(value['main_thread_functions'], ['MPI_Waitall', 'hypre_BoomerAMGCoarsenRuge'])
        self.assertEqual(value['other_thread_functions'], ['cudaEventSynchronize'])
        self.assertEqual(value['stack_errors'], ['older_frame'])
        self.assertTrue(value['unwind_incomplete'])
        self.assertNotIn('PRIVATE', json.dumps(value))

    def test_reason_failure_does_not_stop_frame_name_collection(self):
        gdb = FakeGdb()
        gdb.frames[1].unwind_stop_reason = Mock(side_effect=AttributeError('PRIVATE'))
        value = sampler.stack_functions(gdb)
        self.assertEqual(value['main_thread_functions'],
                         ['MPI_Waitall', 'hypre_BoomerAMGCoarsenRuge', 'hypre_BoomerAMGSetup'])
        self.assertEqual(value['stack_errors'], ['unwind_reason'])
        self.assertTrue(value['unwind_incomplete'])

    def test_stop_reason_is_checked_after_attempting_unwind(self):
        gdb = FakeGdb()
        frame = Frame('MPI_Wait')
        def older():
            frame.reason = 3
            return None
        frame.older = older
        gdb.frames[1] = frame
        self.assertEqual(sampler.stack_functions(gdb)['stack_errors'], ['unwind_stopped'])

    def test_running_main_thread_is_not_read(self):
        gdb = FakeGdb()
        gdb.threads[0].is_stopped = lambda: False
        gdb.threads[0].switch = Mock()
        value = sampler.stack_functions(gdb)
        gdb.threads[0].switch.assert_not_called()
        self.assertEqual(value['main_thread_functions'], [])
        self.assertFalse(value['main_thread_stopped'])
        self.assertFalse(value['main_thread_frame_readable'])
        self.assertEqual(value['stack_errors'], ['thread_not_stopped'])

    def test_missing_frames_and_missing_main_thread_are_distinct(self):
        gdb = FakeGdb()
        gdb.frames[1] = None
        self.assertEqual(sampler.stack_functions(gdb)['stack_errors'], ['no_frames'])
        gdb.threads = gdb.threads[1:]
        self.assertEqual(sampler.stack_functions(gdb)['stack_errors'], ['main_thread_missing'])

    def test_frame_limit_keeps_bounded_partial_stack(self):
        gdb = FakeGdb()
        gdb.frames[1]._older = gdb.frames[1]
        value = sampler.stack_functions(gdb)
        self.assertEqual(len(value['main_thread_functions']), 256)
        self.assertEqual(value['stack_errors'], ['frame_limit'])
        sampler.safe_report(report(**value))

    def test_setup_sample_and_normal_return_are_distinct(self):
        for mode in ('sample', 'return', 'before_setup', 'exit', 'error'):
            with self.subTest(mode=mode):
                directory = self.root / mode
                directory.mkdir()
                gdb = FakeGdb(mode)
                sampler.capture(gdb, dict(directory=str(directory), rank=0, ranks=1, self_test=False))
                value = sampler.safe_report(trace.read_json(directory / 'stack-0.json'))
                self.assertEqual(value['setup_entry_seen'], mode in ('sample', 'return'))
                self.assertEqual(value['setup_return_seen'], mode == 'return')
                self.assertEqual(value['interrupt_seen'], mode in ('sample', 'before_setup'))
                self.assertEqual(value['debugger_error'], mode == 'error')
                self.assertEqual(value['inferior_exit_seen'], mode == 'exit')
                self.assertEqual(value['inferior_exit_nonzero'], mode == 'exit')
                self.assertNotIn('PRIVATE', json.dumps(value))
                self.assertNotIn('SECRET', json.dumps(value))
                self.assertFalse(gdb.alive)
                self.assertEqual(gdb.commands[-1], 'quit')
                self.assertIn('unset environment SIMPLE_SETUP_STACK_CONFIG', gdb.commands)

    def test_no_rank_is_killed_before_other_stack_is_saved(self):
        gdb = FakeGdb()

        def other_rank(_):
            self.assertNotIn('kill', gdb.commands)
            self.assertTrue((self.root / 'stack-0.json').exists())
            sampler.write_json(self.root / 'stack-1.json', report())

        with patch.object(sampler.time, 'sleep', side_effect=other_rank):
            sampler.capture(gdb, dict(directory=str(self.root), rank=0, ranks=2, self_test=False))
        self.assertIn('kill', gdb.commands)

    def test_watchdog_interrupts_gdb_and_keeps_raw_output_private(self):
        child = Mock(pid=123)
        child.poll.side_effect = [None, None, None, 0, 0]

        def interrupted(sig):
            self.assertEqual(sig, signal.SIGINT)
            sampler.write_json(self.root / 'stack-0.json', report(interrupt_seen=True, inferior_alive=True))

        child.send_signal.side_effect = interrupted
        with patch.object(trace.subprocess, 'Popen', return_value=child) as popen, \
                patch.object(trace.time, 'monotonic', side_effect=[0, 0, 31, 32]), \
                patch.object(trace.time, 'sleep'), patch.object(trace.os, 'killpg') as kill:
            value = trace.debugger('/gdb', ['/solver', '--secret', 'PRIVATE'], self.root, 0, 1, {})
        self.assertTrue(value['interrupt_seen'])
        self.assertEqual(popen.call_args[0][0][-3:], ['/solver', '--secret', 'PRIVATE'])
        self.assertEqual(Path(popen.call_args[1]['stdout'].name), self.root / 'debugger-0.log')
        self.assertEqual(popen.call_args[1]['stderr'], trace.subprocess.STDOUT)
        self.assertNotIn('PRIVATE', json.dumps(value))
        child.send_signal.assert_called_once()
        kill.assert_not_called()

    def test_unresponsive_debugger_is_killed_at_hard_deadline(self):
        child = Mock(pid=123)
        child.poll.return_value = None
        with patch.object(trace.subprocess, 'Popen', return_value=child), \
                patch.object(trace.time, 'monotonic', side_effect=[0, 301]), \
                patch.object(trace.os, 'killpg') as kill:
            self.assertIsNone(trace.debugger('/gdb', ['/solver'], self.root, 0, 1, {}))
        kill.assert_called_once_with(123, signal.SIGKILL)

    def test_real_child_receives_only_a_timed_interrupt_and_private_output(self):
        child = self.root / 'fake-gdb'
        child.write_text('#!' + sys.executable + '\n' + '''
import json, os, signal, time
from pathlib import Path
with open(os.environ['SIMPLE_SETUP_STACK_CONFIG']) as stream:
    config = json.load(stream)
root = Path(config['directory'])
value = ''' + repr(report(interrupt_seen=True, inferior_alive=True)) + '''
def interrupted(sig, frame):
    (root / 'stack-0.json').write_text(json.dumps(value))
    raise SystemExit(0)
signal.signal(signal.SIGINT, interrupted)
(root / 'running-0').touch()
print('PRIVATE_RAW_DEBUGGER_TEXT', flush=True)
while True:
    time.sleep(.1)
''')
        child.chmod(0o700)
        value = trace.debugger(str(child), ['/unused'], self.root, 0, 1, os.environ, self_test=True)
        self.assertTrue(value['interrupt_seen'])
        self.assertNotIn('PRIVATE', json.dumps(value))
        self.assertIn('PRIVATE_RAW_DEBUGGER_TEXT', (self.root / 'debugger-0.log').read_text())

    def test_missing_finish_breakpoint_still_allows_stack_sampling(self):
        gdb = FakeGdb()
        class Unavailable(gdb.Breakpoint):
            def __init__(self, *args, **kwargs):
                raise RuntimeError('PRIVATE')
        gdb.FinishBreakpoint = Unavailable
        sampler.capture(gdb, dict(directory=str(self.root), rank=0, ranks=1, self_test=False))
        value = trace.read_json(self.root / 'stack-0.json')
        self.assertTrue(value['setup_entry_seen'])
        self.assertTrue(value['interrupt_seen'])
        self.assertFalse(value['return_watch_available'])
        self.assertFalse(value['debugger_error'])

    def test_existing_output_is_not_modified(self):
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(trace.main(['--pair', '/unused', '--output-dir', str(self.root)]), 1)
        self.assertEqual(list(self.root.iterdir()), [])

    def test_live_progress_uses_fixed_markers_only(self):
        (self.root / 'running-0').touch()
        (self.root / 'entered-0').write_text('PRIVATE 123')
        child = Mock()
        child.wait.side_effect = [trace.subprocess.TimeoutExpired('PRIVATE', 1), 0]
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream):
            self.assertEqual(trace.wait_with_progress(child, self.root, 1), 0)
        self.assertEqual(stream.getvalue(), 'MARS rank 0: debugger running\nMARS rank 0: Hypre setup entered\n')


class SavedCaseTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixtures.ProfileTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.fixture.reference_amg = True
        self.fixture.enable_petsc_momentum()
        self.fixture.prepare()
        self.root, self.pair = self.fixture.root, self.fixture.pair
        baseline = self.fixture.baseline / 'mars'
        for name in ('launch.json', 'launch-start.json'):
            path = baseline / name
            value = startup.read_json(path)
            value['command'][0] = 'srun'
            path.write_text(json.dumps(value))
        # Reprepare after updating the synthetic saved launcher identity.
        import shutil
        shutil.rmtree(self.pair)
        self.fixture.prepare()
        self.libs = patch.object(startup, 'runtime_libraries', return_value=self.fixture.rebuilt_libs)
        self.libs.start()
        self.addCleanup(self.libs.stop)
        self.environment = patch.dict(os.environ, {}, clear=True)
        self.environment.start()
        self.addCleanup(self.environment.stop)

    def test_saved_controls_and_binding_are_preserved_except_output_prefix(self):
        directory = self.root / 'new'
        public = {}
        launcher, command, environment, ranks, identity = trace.saved_command(self.pair, directory, public)
        expected = startup.solver_arguments(self.pair, 'mars')
        expected[expected.index('--output-prefix') + 1] = str(directory / 'flow')
        self.assertEqual(command, [identity['executable']] + expected)
        self.assertEqual(launcher[0], 'srun')
        self.assertEqual(environment, {})
        self.assertTrue(public['library_checks']['matched'])
        self.assertNotIn(str(self.root), json.dumps(public))

    def test_missing_debugger_never_launches_solver(self):
        public = {}
        with patch.object(trace.shutil, 'which', return_value=None), \
                patch.object(trace.subprocess, 'Popen') as launch:
            self.assertEqual(trace.run(self.pair, self.root / 'trace', public), 1)
        launch.assert_not_called()
        self.assertEqual(public['failed_check'], 'debugger_unavailable')

    def test_failed_debugger_self_test_never_allocates_gpus(self):
        cases = (None, report(debugger_error=True),
                 report(interrupt_seen=True, inferior_alive=True, unwind_incomplete=True),
                 readable_report(unwind_incomplete=True, stack_errors=['older_frame']),
                 readable_report(unwind_incomplete=True, stack_errors=['unwind_reason']),
                 readable_report(unwind_incomplete=True, stack_errors=['unwind_stopped', 'frame_name']),
                 readable_report(unwind_incomplete=True),
                 readable_report(main_thread_stopped=False),
                 readable_report(main_thread_frame_readable=False),
                 readable_report(unexpected_signal_seen=True), readable_report(interrupt_seen=False))
        for i, check in enumerate(cases):
            public = {}
            with self.subTest(i=i), patch.object(trace.shutil, 'which', return_value='/gdb'), \
                    patch.object(trace, 'debugger', return_value=check), \
                    patch.object(trace.subprocess, 'Popen') as launch:
                self.assertEqual(trace.run(self.pair, self.root / ('trace-' + str(i)), public), 1)
            launch.assert_not_called()
            self.assertFalse(public['debugger_self_test_passed'])
            self.assertEqual(public['debugger_self_test_stack'], check)

    def test_readable_partial_self_test_allows_capture_without_hiding_incomplete_unwind(self):
        gdb = FakeGdb('before_setup')
        gdb.threads = gdb.threads[:1]
        gdb.frames[1] = Frame('clock_nanosleep', reason=gdb.FRAME_UNWIND_UNAVAILABLE)
        control = self.root / 'control'
        control.mkdir()
        sampler.capture(gdb, dict(directory=str(control), rank=0, ranks=1, self_test=True))
        check = sampler.safe_report(trace.read_json(control / 'stack-0.json'))
        self.assertTrue(check['main_thread_frame_readable'])
        self.assertTrue(check['unwind_incomplete'])
        self.assertEqual(check['main_thread_functions'], [])
        self.assertEqual(check['stack_errors'], ['unwind_stopped'])

        def launch(command, **kwargs):
            directory = Path(kwargs['cwd'])
            ranks = trace.read_json(directory / 'launch-input.json')['ranks']
            for rank in range(ranks):
                sampler.write_json(directory / ('stack-{}.json'.format(rank)),
                                   readable_report(setup_entry_seen=True,
                                                   main_thread_functions=['hypre_BoomerAMGCoarsenRuge'],
                                                   unwind_incomplete=True, stack_errors=['unwind_stopped'],
                                                   unwind_stop_reasons=['unavailable']))
            return Mock(wait=Mock(return_value=0))

        public, stream = {}, io.StringIO()
        with patch.object(trace.shutil, 'which', return_value='/gdb'), \
                patch.object(trace, 'debugger', return_value=check), \
                patch.object(trace.subprocess, 'Popen', side_effect=launch) as launched, \
                contextlib.redirect_stdout(stream):
            self.assertEqual(trace.run(self.pair, self.root / 'trace', public), 1)
        launched.assert_called_once()
        self.assertTrue(public['debugger_self_test_passed'])
        self.assertEqual(public['debugger_self_test_stack'], check)
        self.assertTrue(public['solver_launch_requested'])
        self.assertTrue(public['all_rank_stacks_usable'])
        self.assertEqual(public['comparison_status'], 'diagnostic_partial')
        self.assertEqual(public['failed_check'], 'stack_unwind')
        self.assertIn('the self-test unwind is partial', stream.getvalue())
        self.assertNotIn('clock_nanosleep', json.dumps(public))

    def test_success_aggregates_only_validated_rank_stacks(self):
        output, public = self.root / 'trace', {}
        before = {str(p): startup.digest(p) for p in self.pair.rglob('*') if p.is_file()}

        def launch(command, **kwargs):
            directory = Path(kwargs['cwd'])
            config = trace.read_json(directory / 'launch-input.json')
            for rank in range(config['ranks']):
                sampler.write_json(directory / ('stack-{}.json'.format(rank)),
                                   readable_report(setup_entry_seen=True,
                                                   main_thread_functions=['hypre_BoomerAMGCoarsenRuge']))
            return Mock(wait=Mock(return_value=0))

        with patch.object(trace.shutil, 'which', return_value='/gdb'), \
                patch.object(trace, 'debugger', return_value=readable_report()), \
                patch.object(trace.subprocess, 'Popen', side_effect=launch):
            self.assertEqual(trace.run(self.pair, output, public), 0)
        self.assertTrue(public['all_rank_stacks_present'])
        self.assertTrue(public['all_rank_stacks_usable'])
        self.assertEqual(public['comparison_status'], 'diagnostic_complete')
        self.assertEqual(before, {str(p): startup.digest(p) for p in self.pair.rglob('*') if p.is_file()})
        self.assertNotIn(str(self.root), json.dumps(public))

    def test_unreadable_partial_and_missing_reports_never_claim_complete(self):
        cases = (
            (report(setup_entry_seen=True, inferior_alive=True, interrupt_seen=True, unwind_incomplete=True),
             0, 'diagnostic_incomplete', 'stack_unreadable'),
            (readable_report(), 0, 'diagnostic_incomplete', 'stack_unreadable'),
            (readable_report(main_thread_functions=['MPI_Wait'], unwind_incomplete=True,
                             stack_errors=['older_frame']), 0, 'diagnostic_partial', 'stack_unwind'),
            (readable_report(main_thread_functions=['MPI_Wait'], unexpected_signal_seen=True),
             0, 'diagnostic_incomplete', 'stack_unreadable'),
            (None, 0, 'diagnostic_incomplete', 'rank_stack_capture'),
            (readable_report(main_thread_functions=['MPI_Wait']),
             255, 'diagnostic_incomplete', 'rank_stack_capture'),
        )
        for i, (value, code, status, failure) in enumerate(cases):
            def launch(command, **kwargs):
                directory = Path(kwargs['cwd'])
                ranks = trace.read_json(directory / 'launch-input.json')['ranks']
                if value is not None:
                    for rank in range(ranks):
                        sampler.write_json(directory / ('stack-{}.json'.format(rank)), value)
                return Mock(wait=Mock(return_value=code))
            public = {}
            with self.subTest(i=i), patch.object(trace.shutil, 'which', return_value='/gdb'), \
                    patch.object(trace, 'debugger', return_value=readable_report()), \
                    patch.object(trace.subprocess, 'Popen', side_effect=launch), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(trace.run(self.pair, self.root / ('trace-' + str(i)), public), 1)
            self.assertEqual(public['comparison_status'], status)
            self.assertEqual(public['failed_check'], failure)

    def test_changed_library_or_protected_environment_refuses_before_launch(self):
        for mode in ('library', 'environment'):
            with self.subTest(mode=mode), patch.object(trace.subprocess, 'Popen') as launch:
                context = (patch.object(startup, 'runtime_libraries', return_value={}) if mode == 'library'
                           else patch.dict(os.environ, {'MARS_PRIVATE_UNKNOWN': 'PRIVATE'}))
                with context, self.assertRaises(startup.EvidenceError):
                    trace.saved_command(self.pair, self.root / mode, {})
                launch.assert_not_called()


if __name__ == '__main__':
    unittest.main()
