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
                if getattr(owner, 'in_stop', False):
                    raise AssertionError('breakpoint changed inside stop callback')
                self.spec = args[0] if args else None
                self.pending = False
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


class CommunicationGdb(FakeGdb):
    def __init__(self, advance=False, setup_returns=False):
        super().__init__()
        self.continues, self.samples = 0, 0
        self.actions = iter(['entry', 'sample'] + (['return_entry'] if advance else [])
                            + ['setup_return' if setup_returns else 'sample'])
        self.frames[1] = Frame('PMPI_Iprobe', Frame('hypre_DataExchangeList',
            Frame('hypre_ParCSRCommPkgCreateApart_core', Frame('hypre_MatvecCommPkgCreate',
            Frame('hypre_ParCSRMatrixGenerateFFFCDevice_core', Frame('hypre_BoomerAMGBuildExtPIInterpDevice'))))))

    def execute(self, command):
        if command != 'continue':
            return super().execute(command)
        self.commands.append(command)
        self.continues += 1
        entry = next(p for p in self.points if p.spec == 'hypre_DataExchangeList')
        action = next(self.actions)
        if self.samples == 1:
            getattr(self, 'assert_continues', lambda: None)()
        self.in_stop = True
        try:
            if action in ('entry', 'return_entry'):
                if action == 'return_entry':
                    self.points[-1].stop()
                entry.stop()
                return
            if action == 'setup_return':
                next(p for p in self.points if type(p).__name__ == 'Returned').stop()
                return
        finally:
            self.in_stop = False
        self.samples += 1
        raise KeyboardInterrupt()


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

    def test_communication_helpers_are_literal_names_not_argument_text(self):
        gdb = CommunicationGdb()
        value = sampler.stack_functions(gdb)
        self.assertEqual(value['main_thread_functions'], [
            'PMPI_Iprobe', 'hypre_DataExchangeList', 'hypre_ParCSRCommPkgCreateApart_core',
            'hypre_MatvecCommPkgCreate', 'hypre_ParCSRMatrixGenerateFFFCDevice_core',
            'hypre_BoomerAMGBuildExtPIInterpDevice'])
        for name in ('hypre_DataExchangeList PRIVATE', 'hypre_PrivateMesh', 'MPI_Test PRIVATE'):
            with self.assertRaises(ValueError):
                sampler.safe_report(report(main_thread_functions=[name]))

    def test_progress_distinguishes_same_named_stack_with_completed_exchanges(self):
        for advance in (False, True):
            directory = self.root / str(advance)
            directory.mkdir()
            gdb = CommunicationGdb(advance=advance)
            sampler.capture(gdb, dict(directory=str(directory), rank=0, ranks=1,
                                     self_test=False, communication_progress=True))
            final = trace.read_json(directory / 'stack-0.json')
            value = trace.communication_report(directory / 'exchange-0.json', final)
            self.assertTrue(value['two_samples_usable'])
            self.assertTrue(value['exchange_tracking_complete'])
            self.assertEqual(value['exchange_advanced_between_samples'], advance)
            self.assertEqual(value['first_main_thread_functions'], final['main_thread_functions'])
            self.assertTrue(value['exchange_active_at_first_sample'])
            self.assertTrue(value['exchange_active_at_final_sample'])
            self.assertNotIn('entered', value)
            self.assertNotIn('returned', value)
            self.assertNotIn('SECRET', json.dumps(value))

    def test_second_sample_can_observe_setup_return(self):
        gdb = CommunicationGdb(advance=True, setup_returns=True)
        sampler.capture(gdb, dict(directory=str(self.root), rank=0, ranks=1,
                                 self_test=False, communication_progress=True))
        final = trace.read_json(self.root / 'stack-0.json')
        self.assertTrue(final['setup_return_seen'])
        self.assertFalse(final['interrupt_seen'])
        self.assertTrue(trace.communication_report(self.root / 'exchange-0.json', final)['two_samples_usable'])

    def test_all_ranks_must_save_first_sample_before_any_resume(self):
        gdb = CommunicationGdb()

        def peer(_):
            self.assertEqual(gdb.samples, 1)
            self.assertNotIn('kill', gdb.commands)
            first = trace.read_json(self.root / 'first-0.json')
            sampler.write_json(self.root / 'first-1.json', first)
            # The peer's final report also lets cleanup finish without waiting.
            sampler.write_json(self.root / 'stack-1.json', first['stack'])
            gdb.assert_continues = lambda: self.assertTrue((self.root / 'first-1.json').exists())

        with patch.object(sampler.time, 'sleep', side_effect=peer):
            sampler.capture(gdb, dict(directory=str(self.root), rank=0, ranks=2,
                                     self_test=False, communication_progress=True))
        self.assertEqual(gdb.samples, 2)

    def test_returned_peer_prevents_resuming_into_a_flow_run(self):
        first = dict(stack=readable_report(setup_entry_seen=True, setup_return_seen=True), exchange={})
        sampler.write_json(self.root / 'first-1.json', first)
        sampler.write_json(self.root / 'stack-1.json', first['stack'])
        gdb = CommunicationGdb()
        sampler.capture(gdb, dict(directory=str(self.root), rank=0, ranks=2,
                                 self_test=False, communication_progress=True))
        value = trace.communication_report(self.root / 'exchange-0.json', trace.read_json(self.root / 'stack-0.json'))
        self.assertEqual(gdb.samples, 1)
        self.assertFalse(value['two_samples_taken'])
        self.assertIsNone(value['exchange_advanced_between_samples'])

    def test_missing_peer_times_out_without_resuming(self):
        gdb = CommunicationGdb()
        with patch.object(sampler.time, 'monotonic', side_effect=(0, 61, 61, 122)):
            sampler.capture(gdb, dict(directory=str(self.root), rank=0, ranks=2,
                                     self_test=False, communication_progress=True))
        value = trace.communication_report(self.root / 'exchange-0.json', trace.read_json(self.root / 'stack-0.json'))
        self.assertEqual(gdb.samples, 1)
        self.assertFalse(value['peer_samples_present'])
        self.assertFalse(value['two_samples_taken'])
        self.assertIsNone(value['exchange_advanced_between_samples'])
        self.assertIn('kill', gdb.commands)

    def test_unavailable_return_tracking_is_inconclusive(self):
        gdb = CommunicationGdb()
        old = gdb.FinishBreakpoint
        class Unavailable(old):
            def __init__(self, *args, **kwargs):
                if type(self).__name__ == 'ExchangeReturned':
                    raise RuntimeError('PRIVATE')
                super().__init__(*args, **kwargs)
        gdb.FinishBreakpoint = Unavailable
        sampler.capture(gdb, dict(directory=str(self.root), rank=0, ranks=1,
                                 self_test=False, communication_progress=True))
        value = trace.communication_report(self.root / 'exchange-0.json', trace.read_json(self.root / 'stack-0.json'))
        self.assertFalse(value['exchange_tracking_complete'])
        self.assertIsNone(value['exchange_advanced_between_samples'])
        self.assertNotIn('PRIVATE', json.dumps(value))

    def test_communication_report_refuses_untrusted_counters_and_partial_samples(self):
        counts = dict(entered=2, returned=1, tracking_complete=True)
        value = dict(first=dict(stack=readable_report(), exchange=counts), final=dict(counts),
                     two_samples_taken=True, peer_samples_present=True)
        path = self.root / 'exchange.json'
        for final in (dict(counts, entered='PRIVATE'), dict(counts, entered=True),
                      dict(counts, returned=3), dict(counts, entered=1, returned=0),
                      dict(counts, raw='PRIVATE')):
            path.write_text(json.dumps(dict(value, final=final)))
            with self.assertRaises(ValueError):
                trace.communication_report(path, readable_report())
        path.write_text(json.dumps(value))
        result = trace.communication_report(path, readable_report(unwind_incomplete=True))
        self.assertFalse(result['two_samples_usable'])
        self.assertIsNone(result['exchange_advanced_between_samples'])

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

    def test_watchdog_waits_for_peers_then_times_second_interrupt_from_resume(self):
        child = Mock(pid=123)
        child.poll.return_value = None
        clock = iter((0, 0, 31, 60, 61, 76))
        signals = []

        def now():
            value = next(clock)
            if value == 61:
                (self.root / 'resumed-0').touch()
            return value

        def interrupted(sig):
            signals.append(sig)
            if len(signals) == 1:
                (self.root / 'first-0.json').write_text('{}')
            else:
                self.assertTrue((self.root / 'resumed-0').exists())
                sampler.write_json(self.root / 'stack-0.json', readable_report())
                child.poll.return_value = 0

        child.send_signal.side_effect = interrupted
        with patch.object(trace.subprocess, 'Popen', return_value=child), \
                patch.object(trace.time, 'monotonic', side_effect=now), \
                patch.object(trace.time, 'sleep'), patch.object(trace.os, 'killpg') as kill:
            value = trace.debugger('/gdb', ['/solver'], self.root, 0, 2, {}, communication_progress=True)
        self.assertTrue(value['interrupt_seen'])
        self.assertEqual(signals, [signal.SIGINT, signal.SIGINT])
        kill.assert_not_called()

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

    def test_communication_capture_exports_progress_without_counters_or_paths(self):
        for tracking in (True, False):
            def launch(command, **kwargs):
                directory = Path(kwargs['cwd'])
                config = trace.read_json(directory / 'launch-input.json')
                self.assertTrue(config['communication_progress'])
                for rank in range(config['ranks']):
                    stack = readable_report(setup_entry_seen=True, main_thread_functions=['hypre_DataExchangeList'])
                    sampler.write_json(directory / ('stack-{}.json'.format(rank)), stack)
                    sampler.write_json(directory / ('exchange-{}.json'.format(rank)), dict(
                        first=dict(stack=stack, exchange=dict(entered=123, returned=122, tracking_complete=tracking)),
                        final=dict(entered=124, returned=123, tracking_complete=tracking),
                        two_samples_taken=True, peer_samples_present=True))
                return Mock(wait=Mock(return_value=0))

            public = {}
            with patch.object(trace.shutil, 'which', return_value='/gdb'), \
                    patch.object(trace, 'debugger', return_value=readable_report()), \
                    patch.object(trace.subprocess, 'Popen', side_effect=launch), \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(trace.run(self.pair, self.root / str(tracking), public, True), 0 if tracking else 1)
            self.assertEqual(public['failed_check'], 'none' if tracking else 'communication_tracking')
            self.assertEqual(public['comparison_status'], 'diagnostic_complete' if tracking else 'diagnostic_partial')
            for progress in public['communication_progress'].values():
                self.assertEqual(progress['exchange_advanced_between_samples'], True if tracking else None)
            self.assertNotIn('123', json.dumps(public))
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


class SavedCommunicationTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.source = self.root / 'capture'
        self.private = self.source / 'private'
        self.private.mkdir(parents=True)

    def fixture(self, entries=(12345, 12345, 12345, 12345), change=None):
        saved = dict(schema=trace.SCHEMA, comparison_status='diagnostic_complete', failed_check='none',
                     ranks_expected=len(entries), rank_stacks={}, communication_progress={})
        for rank, entered in enumerate(entries):
            stack = readable_report(setup_entry_seen=True, main_thread_functions=['hypre_DataExchangeList'])
            counts = dict(entered=entered, returned=entered - 1, tracking_complete=True)
            value = dict(first=dict(stack=stack, exchange=dict(counts)), final=dict(counts),
                         two_samples_taken=True, peer_samples_present=True)
            if change:
                change(rank, value, stack)
            trace.write_json(self.private / ('first-{}.json'.format(rank)), value['first'])
            trace.write_json(self.private / ('exchange-{}.json'.format(rank)), value)
            trace.write_json(self.private / ('stack-{}.json'.format(rank)), stack)
            saved['rank_stacks'][str(rank)] = stack
            saved['communication_progress'][str(rank)] = trace.communication_value(value, stack)
        # An unrelated private string must never be copied into the output.
        saved['private_extra'] = 'SECRET_CASE_NAME'
        trace.write_json(self.source / 'public.json', saved)

    def inspect(self):
        output = self.root / 'inspection'
        before = {p: p.read_bytes() for p in self.source.rglob('*') if p.is_file()}
        with patch.object(trace, 'run', side_effect=AssertionError('solver launch')), \
                patch.object(trace, 'debugger', side_effect=AssertionError('debugger launch')), \
                contextlib.redirect_stdout(io.StringIO()) as console:
            code = trace.main(['--inspect-dir', str(self.source), '--output-dir', str(output)])
        self.assertEqual(before, {p: p.read_bytes() for p in self.source.rglob('*') if p.is_file()})
        result = trace.read_json(output / 'public.json')
        text = json.dumps(result) + console.getvalue()
        for private in ('12345', '12346', 'SECRET', str(self.source)):
            self.assertNotIn(private, text)
        self.assertFalse(result['solver_launched'])
        return code, result

    def test_matching_ordinals_do_not_claim_protocol_identity(self):
        self.fixture()
        code, result = self.inspect()
        self.assertEqual(code, 0)
        self.assertEqual(result['assessment'], 'same_observed_exchange_ordinal')
        self.assertTrue(result['same_observed_exchange_ordinal'])
        self.assertFalse(result['protocol_phase_verified'])
        self.assertFalse(result['communicator_identity_verified'])

    def test_different_ordinals_are_reported_without_counts(self):
        self.fixture((12345, 12345, 12346, 12346))
        code, result = self.inspect()
        self.assertEqual(code, 0)
        self.assertEqual(result['assessment'], 'different_observed_exchange_ordinals')
        self.assertFalse(result['same_observed_exchange_ordinal'])

    def test_tracking_loss_is_inconclusive(self):
        self.fixture(change=lambda rank, value, stack: value['final'].update(tracking_complete=False))
        code, result = self.inspect()
        self.assertEqual(code, 1)
        self.assertIsNone(result['same_observed_exchange_ordinal'])

    def test_advancing_exchange_is_not_compared_as_stationary(self):
        self.fixture(change=lambda rank, value, stack: value['final'].update(entered=12346, returned=12345))
        code, result = self.inspect()
        self.assertEqual(code, 1)
        self.assertFalse(result['no_rank_completed_or_entered_exchange_between_samples'])

    def test_nested_exchanges_are_inconclusive(self):
        def nested(rank, value, stack):
            value['first']['exchange']['returned'] -= 1
            value['final']['returned'] -= 1
        self.fixture(change=nested)
        code, result = self.inspect()
        self.assertEqual(code, 1)
        self.assertFalse(result['single_active_exchange_per_rank'])

    def test_missing_rank_record_is_invalid(self):
        self.fixture()
        (self.private / 'exchange-3.json').unlink()
        code, result = self.inspect()
        self.assertEqual(code, 1)
        self.assertFalse(result['saved_reports_consistent'])

    def test_changed_public_report_is_invalid(self):
        self.fixture()
        path = self.source / 'public.json'
        value = trace.read_json(path)
        value['communication_progress']['0']['exchange_advanced_between_samples'] = True
        path.write_text(json.dumps(value))
        code, result = self.inspect()
        self.assertEqual(code, 1)
        self.assertFalse(result['saved_reports_consistent'])

    def test_duplicate_first_record_must_match(self):
        self.fixture()
        path = self.private / 'first-0.json'
        value = trace.read_json(path)
        value['exchange']['entered'] += 1
        path.write_text(json.dumps(value))
        code, result = self.inspect()
        self.assertEqual(code, 1)
        self.assertFalse(result['saved_reports_consistent'])


if __name__ == '__main__':
    unittest.main()
