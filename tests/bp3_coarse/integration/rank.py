"""Keep optional audit writes off the shared filesystem while the solver runs."""
import argparse
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
from common import digest, save, unblock_sigchld


def run(out, audited, args, scratch_base=Path('/dev/shm')):
    rank = int(os.environ['SLURM_PROCID'])
    out = out.resolve()
    with tempfile.TemporaryDirectory(prefix='bp3-integration-', dir=str(scratch_base)) as scratch:
        audit = Path(scratch) / 'audit.json'
        env = dict(os.environ, BP3_CONSOLIDATION_AUDIT=str(audit))
        receipt = dict(rank=rank, hostname=os.uname().nodename, status='starting',
                       parent_blocked_signals=sorted(int(s) for s in signal.pthread_sigmask(signal.SIG_BLOCK, [])),
                       binary_sha256=digest(args[0]), audit_expected=audited)
        save(out / ('process%d.json' % rank), receipt)
        # Slurm can pass MPI bootstrap sockets through inherited descriptors.
        result = subprocess.run(args, env=env, stdin=subprocess.DEVNULL, close_fds=False,
                                preexec_fn=unblock_sigchld)
        receipt.update(returncode=result.returncode, status='complete')
        if audit.is_file():
            target = out / ('rank%d.json' % rank)
            shutil.copyfile(str(audit), str(target))
            receipt['audit_sha256'] = digest(target)
        save(out / ('process%d.json' % rank), receipt)
        if result.returncode:
            return result.returncode
        if audit.exists() != audited:
            raise ValueError('Wrong audit presence for this build')
    return 0


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('out', type=Path)
    p.add_argument('audited', choices=('yes', 'no'))
    p.add_argument('command', nargs=argparse.REMAINDER)
    a = p.parse_args()
    raise SystemExit(run(a.out, a.audited == 'yes', a.command))
