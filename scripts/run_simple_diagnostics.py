#!/usr/bin/env python3
"""Capture a user-launched command privately and maintain a shareable SIMPLE report."""

import argparse
import contextlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from simple_public_diagnostics import DiagnosticState, SafeParser


def create_private(path, mode):
    fd = os.open(str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    return os.fdopen(fd, mode)


def publish(path, result):
    # Readers see either the previous complete JSON or the new one.
    fd, temporary = tempfile.mkstemp(prefix='.simple-public-', dir=str(path.parent))
    try:
        with os.fdopen(fd, 'w') as output:
            json.dump(result, output, indent=2, sort_keys=True, allow_nan=False)
            output.write('\n')
        os.replace(temporary, str(path))
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def stop(child):
    if child is not None and child.poll() is None:
        child.terminate()
        try:
            child.wait(timeout=10)
        except subprocess.TimeoutExpired:
            child.kill()
            child.wait()


def capture(command, log_path, exit_path, output_path, environment=None):
    paths = [Path(p).resolve() for p in (log_path, exit_path, output_path)]
    if len(set(paths)) != 3:
        raise ValueError('Capture paths must differ')
    state = DiagnosticState()
    child = None
    with contextlib.ExitStack() as stack:
        log = stack.enter_context(create_private(paths[0], 'wb'))
        status = stack.enter_context(create_private(paths[1], 'w'))
        with create_private(paths[2], 'w'):
            pass
        previous = state.result('')
        previous['capture_status'] = 'running'
        publish(paths[2], previous)
        code, capture_status = 1, 'failed'
        try:
            child = subprocess.Popen(command, stdin=subprocess.DEVNULL,
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                     env=environment)
            for raw in iter(child.stdout.readline, b''):
                log.write(raw)
                log.flush()
                state.feed(raw.decode('utf-8', errors='replace'))
                current = state.result('')
                current['capture_status'] = 'running'
                if current != previous:
                    publish(paths[2], current)
                    previous = current
            code = child.wait()
            code = code if code >= 0 else 128 - code
            capture_status = 'finished'
        except KeyboardInterrupt:
            code, capture_status = 130, 'interrupted'
        except Exception:
            # Never print an exception containing the command, path or private output.
            code, capture_status = (127, 'launch_failed') if child is None else (1, 'failed')
        finally:
            stop(child)
            if child is not None:
                child.stdout.close()
        status.write(str(code) + '\n')
        status.flush()
        result = state.result(str(code))
        result['capture_status'] = capture_status
        publish(paths[2], result)
    return 1 if code == 0 and result['run_status'] == 'failed' else code


def main():
    parser = SafeParser(description=__doc__)
    parser.add_argument('--log', required=True)
    parser.add_argument('--exit-file', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command:
        parser.error('A command is required')
    try:
        code = capture(command, args.log, args.exit_file, args.output)
    except KeyboardInterrupt:
        return 130
    except Exception:
        print('ERROR: diagnostic capture failed; check files privately.', file=sys.stderr)
        return 1
    print('Diagnostic capture finished. Share only the public JSON file.')
    return code


if __name__ == '__main__':
    sys.exit(main())
