"""GDB-side setup sampler. Only fixed library symbols enter the exported report."""

import json
import os
from pathlib import Path
import time


SETUP_FUNCTIONS = ('HYPRE_ParCSRGMRESSetup', 'HYPRE_ParCSRFlexGMRESSetup')
FUNCTIONS = frozenset(SETUP_FUNCTIONS + (
    'hypre_GMRESSetup', 'hypre_FlexGMRESSetup', 'HYPRE_BoomerAMGSetup',
    'hypre_BoomerAMGSetup', 'hypre_BoomerAMGCreateS', 'hypre_BoomerAMGCreateSHost',
    'hypre_BoomerAMGCreateSDevice', 'hypre_BoomerAMGCoarsenHMIS',
    'hypre_BoomerAMGCoarsenRuge', 'hypre_BoomerAMGCoarsenPMIS',
    'hypre_BoomerAMGCoarsenPMISHost', 'hypre_BoomerAMGCoarsenPMISDevice',
    'hypre_BoomerAMGBuildExtInterp', 'hypre_BoomerAMGBuildExtInterpDevice',
    'hypre_BoomerAMGBuildExtPIInterp', 'hypre_BoomerAMGBuildExtPIInterpHost',
    'hypre_BoomerAMGBuildExtPIInterpDevice', 'hypre_BoomerAMGBuildModExtInterp',
    'hypre_BoomerAMGBuildModExtPIInterp', 'hypre_BoomerAMGBuildCoarseOperator',
    'hypre_BoomerAMGBuildCoarseOperatorKT', 'hypre_ParCSRMatrixRAP',
    'hypre_ParCSRMatrixRAPKT', 'hypre_ParCSRMatrixRAPKTDevice',
    'hypre_ParCSRMatrixClone', 'hypre_ParCSRMatrixClone_v2',
    'hypre_GaussElimSetup', 'hypre_GaussElimSetupDevice', 'hypre_GenerateSubComm',
    'hypre_ParKrylovCreateVectorArray', 'hypre_ParVectorInitialize_v2',
    'hypre_MAlloc', 'hypre_CAlloc', 'hypre_Memcpy', 'hypre_SyncComputeStream',
    'cudaMalloc', 'cudaMallocManaged', 'cudaMemcpy', 'cudaMemcpyAsync',
    'cudaStreamSynchronize', 'cudaDeviceSynchronize', 'cudaEventSynchronize',
    'cuStreamSynchronize', 'cuCtxSynchronize', 'cuEventSynchronize',
    'KSPSetUp', 'KSPSolve', 'PCSetUp', 'PetscInitialize', 'MPI_Init_thread',
) + tuple(prefix + operation for prefix in ('MPI_', 'PMPI_', 'hypre_MPI_')
          for operation in ('Barrier', 'Wait', 'Waitall', 'Waitany', 'Test', 'Testall',
                            'Allreduce', 'Allgather', 'Allgatherv', 'Bcast', 'Recv',
                            'Send', 'Sendrecv', 'Comm_split')))
BOOL_KEYS = ('setup_entry_seen', 'setup_return_seen', 'return_watch_available',
             'interrupt_seen', 'unexpected_signal_seen', 'inferior_alive',
             'inferior_exit_seen', 'inferior_exit_nonzero', 'debugger_error', 'unwind_incomplete',
             'main_thread_stopped', 'main_thread_frame_readable')
STACK_ERRORS = frozenset(('thread_list', 'thread_state', 'thread_not_stopped', 'thread_switch',
                          'newest_frame', 'no_frames', 'frame_name', 'older_frame', 'unwind_reason',
                          'unwind_stopped', 'frame_limit', 'main_thread_missing'))
UNWIND_ERRORS = ('UNAVAILABLE', 'INNER_ID', 'SAME_ID', 'NO_SAVED_PC', 'MEMORY_ERROR')
LIST_VALUES = dict(main_thread_functions=FUNCTIONS, other_thread_functions=FUNCTIONS,
                   stack_errors=STACK_ERRORS,
                   unwind_stop_reasons=frozenset(name.lower() for name in UNWIND_ERRORS) | {'unknown'})


def write_json(path, value):
    # Readers must never mistake a partially written report for a stopped rank.
    temporary = path.with_suffix('.tmp')
    with temporary.open('x') as stream:
        json.dump(value, stream, sort_keys=True, allow_nan=False)
        stream.write('\n')
    temporary.replace(path)


def safe_report(value):
    if (not isinstance(value, dict) or set(value) != set(BOOL_KEYS) | set(LIST_VALUES)
            or any(type(value[k]) is not bool for k in BOOL_KEYS)
            or any(not isinstance(value[k], list) or len(value[k]) > 256
                   or any(not isinstance(s, str) or s not in allowed for s in value[k])
                   for k, allowed in LIST_VALUES.items())):
        raise ValueError('invalid stack report')
    return {k: value[k] for k in BOOL_KEYS + tuple(LIST_VALUES)}


def stack_functions(gdb):
    result = dict(main_thread_functions=[], other_thread_functions=[], main_thread_stopped=False,
                  main_thread_frame_readable=False)
    others, errors, stops, main_seen = set(), set(), set(), False
    # GDB 13.2 exports these reasons but not FRAME_UNWIND_FIRST_ERROR.
    normal_reasons = {getattr(gdb, name) for name in
                      ('FRAME_UNWIND_NO_REASON', 'FRAME_UNWIND_NULL_ID', 'FRAME_UNWIND_OUTERMOST')
                      if hasattr(gdb, name)}
    stop_names = {getattr(gdb, 'FRAME_UNWIND_' + name): name.lower() for name in UNWIND_ERRORS
                  if hasattr(gdb, 'FRAME_UNWIND_' + name)}
    try:
        threads = gdb.selected_inferior().threads()
    except Exception:
        threads = ()
        errors.add('thread_list')
    for thread in threads:
        names, is_main, stage = [], False, 'thread_state'
        try:
            is_main = thread.num == 1
            main_seen |= is_main
            stopped = bool(thread.is_stopped())
            if is_main:
                result['main_thread_stopped'] = stopped
            if not stopped:
                errors.add('thread_not_stopped')
                continue
            stage = 'thread_switch'
            thread.switch()
            stage = 'newest_frame'
            frame = gdb.newest_frame()
            if frame is None:
                errors.add('no_frames')
            for _ in range(256):
                if frame is None:
                    break
                try:
                    name = frame.name()
                    if is_main and isinstance(name, str) and name:
                        result['main_thread_frame_readable'] = True
                    if name in FUNCTIONS:
                        names.append(name)
                except Exception:
                    errors.add('frame_name')
                stage = 'older_frame'
                older = frame.older()
                try:
                    # Asking for the older frame first lets GDB establish the stop reason.
                    reason = frame.unwind_stop_reason()
                    if reason not in normal_reasons:
                        errors.add('unwind_stopped')
                        stops.add(stop_names.get(reason, 'unknown'))
                except Exception:
                    errors.add('unwind_reason')
                frame = older
            if frame is not None:
                errors.add('frame_limit')
        except Exception:
            errors.add(stage)
        finally:
            # A later unreadable frame must not discard the frames already captured.
            if is_main:
                result['main_thread_functions'] = names
            else:
                others.update(names)
    if not main_seen:
        errors.add('main_thread_missing')
    result.update(other_thread_functions=sorted(others), stack_errors=sorted(errors),
                  unwind_stop_reasons=sorted(stops), unwind_incomplete=bool(errors))
    return result


def capture(gdb, config):
    output = Path(config['directory'])
    rank, ranks = config['rank'], config['ranks']
    report = dict.fromkeys(BOOL_KEYS, False)
    report.update({key: [] for key in LIST_VALUES})

    class Entry(gdb.Breakpoint):
        def stop(self):
            report['setup_entry_seen'] = True
            return True

    class Returned(gdb.FinishBreakpoint):
        def stop(self):
            report['setup_return_seen'] = True
            return True

    def running(event):
        (output / ('running-{}'.format(rank))).touch()

    def stopped(event):
        sig = getattr(event, 'stop_signal', None)
        report['interrupt_seen'] |= sig == 'SIGINT'
        report['unexpected_signal_seen'] |= sig is not None and sig != 'SIGINT'

    def exited(event):
        report['inferior_exit_seen'] = True
        code = getattr(event, 'exit_code', None)
        report['inferior_exit_nonzero'] = code is not None and code != 0

    def execute(command):
        try:
            gdb.execute(command)
        except KeyboardInterrupt:
            report['interrupt_seen'] = True

    try:
        for command in ('set confirm off', 'set pagination off', 'set breakpoint pending on',
                        'set print frame-arguments none', 'set startup-with-shell off',
                        'set disable-randomization off', 'unset environment SIMPLE_SETUP_STACK_CONFIG'):
            gdb.execute(command)
        gdb.events.cont.connect(running)
        gdb.events.stop.connect(stopped)
        gdb.events.exited.connect(exited)
        points = [] if config['self_test'] else [Entry(name, internal=True) for name in SETUP_FUNCTIONS]
        execute('run')
        if report['setup_entry_seen']:
            for point in points:
                point.delete()
            try:
                Returned(internal=True)
                report['return_watch_available'] = True
            except Exception:
                # Sampling still works when a stripped frame cannot be watched on return.
                pass
            (output / ('entered-{}'.format(rank))).touch()
            execute('continue')
        report['inferior_alive'] = gdb.selected_inferior().pid != 0
        if report['inferior_alive']:
            report.update(stack_functions(gdb))
    except Exception:
        report['debugger_error'] = True
    finally:
        write_json(output / ('stack-{}.json'.format(rank)), safe_report(report))
        # Do not kill one MPI rank until the others have saved their stacks.
        deadline = time.monotonic() + 130
        while not all((output / ('stack-{}.json'.format(i))).is_file() for i in range(ranks)):
            if time.monotonic() >= deadline:
                break
            time.sleep(.2)
        try:
            if gdb.selected_inferior().pid:
                gdb.execute('kill')
        finally:
            gdb.execute('quit')


if __name__ == '__main__':
    import gdb
    with open(os.environ['SIMPLE_SETUP_STACK_CONFIG']) as stream:
        capture(gdb, json.load(stream))
