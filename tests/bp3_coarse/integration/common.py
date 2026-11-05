"""Pinned inputs and checks for the standard-build integration gate."""
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess

IMAGE = '/ritom/scratch/cscs/gandanie/.uenv-images/images/fe67719844db4b44e7a05cd88bc0dd1cfd5e3aeea9ba1a0734d98184c1e3e2d2/store.squashfs'
EXCLUDE = 'nid005323,nid005617,nid005669,nid005679,nid005742,nid005820,nid005854,nid006377,nid005927'
HYPRE = Path('/user-environment/linux-neoverse_v2/hypre-2.33.0-vnuxfgup65dzsv7rsjf7fnxly3h3cozo/lib/libHYPRE-2.33.0.so')
HYPRE_SHA = '0c2c23c7a5d4537c842ccadfb06a2c2d45f65c3de47c80fda3ffba9e46c58dc6'
CONFIG_SHA = 'c81ae5ae99c632c154a46cbe0aae3419c2a15c59d94c415062a0a8bdf03d2b09'
ACCEPTED = Path('/capstor/scratch/cscs/gandanie/bp3-coarse-finish/5253439/native')
ACCEPTED_SHA = '685a2166cc0ce340c830961ddbd5bab5667ea3ba960f488cacf1743263a93b53'
CALLBACKS = ('hypre_seqAMGSetup', 'hypre_seqAMGCycle', 'mars_bp3_amg_setup', 'mars_bp3_amg_destroy')


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def verify_files(root, files):
    for name, expected in files.items():
        path = Path(name)
        if path.is_absolute() or '..' in path.parts:
            raise ValueError('Unsafe manifest name: ' + name)
        path = root / path
        if path.is_symlink() or not path.is_file() or digest(path) != expected:
            raise ValueError('Changed or missing input: ' + str(path))


def verify_bundle(bundle):
    record = json.loads((bundle / 'manifest.json').read_text())
    actual = set()
    for path in bundle.rglob('*'):
        if path.is_symlink():
            raise ValueError('Bundle contains a symbolic link')
        if path.is_file():
            actual.add(str(path.relative_to(bundle)))
    if actual != set(record['files']) | {'manifest.json'}:
        raise ValueError('Bundle inventory changed')
    verify_files(bundle, record['files'])
    return record


def unblock_sigchld():
    # The login environment can otherwise leave CMake waiting for completed children.
    signal.pthread_sigmask(signal.SIG_UNBLOCK, {signal.SIGCHLD})


def command(args, log, env=None, expected=0, timeout=None):
    args = [str(arg) for arg in args]
    save(Path(str(log) + '.command.json'), args)
    with Path(log).open('w') as stream:
        result = subprocess.run(args, stdin=subprocess.DEVNULL, stdout=stream, stderr=subprocess.STDOUT,
                                env=env, preexec_fn=unblock_sigchld, timeout=timeout)
    if result.returncode != expected:
        raise RuntimeError('Exit {} (expected {}): {}'.format(result.returncode, expected, log))
    return Path(log).read_text()


def step(ranks, minutes=20):
    return ['srun', '--uenv-passthrough=use', '--uenv=' + IMAGE + ':/user-environment', '--view=default',
            '--nodes=' + str((ranks + 3) // 4), '-n' + str(ranks),
            '--ntasks-per-node=' + str(min(4, ranks)), '--cpus-per-task=8', '--export=ALL',
            '--kill-on-bad-exit=1', '--exclude=' + EXCLUDE, '--time={:02d}:{:02d}:00'.format(minutes // 60, minutes % 60)]


def solver_args(generated, p, ne, cutoff):
    return ['--ptx-dir', str(generated), '--spec', 'bp3', '--p', str(p), '--ne', str(ne),
            '--deform', '0.05', '--mg', '--amg', '--amg-cycles', '3', '--amg-consolidate', str(cutoff),
            '--pre', '2', '--post', '2', '--cheb-interval', 'mars', '--rhs', 'discrete',
            '--tol', '1e-13', '--true-tol', '1e-10', '--maxit', '2000', '--solves', '3', '--reps', '1', '--amg-print']


def parse(directory, ranks, cutoff, audited=True):
    log = (directory / 'solver.log').read_text()
    matches = re.findall(r'FCG \+ V-cycle: (\d+) iterations, sqrt\(<r,z>/<r0,z0>\) = (\S+), true assembled residual (\S+)', log)
    if len(matches) != 1 or not re.search(r'\bPASS\b', log) or re.search(r'\bFAIL\b', log):
        raise ValueError('Missing or failed solve gate: ' + str(directory))
    it, relative, residual = int(matches[0][0]), float(matches[0][1]), float(matches[0][2])
    if not (0 < it < 2000 and 0 <= relative <= 1e-13 and 0 <= residual <= 1e-10):
        raise ValueError('Accuracy gate failed')
    wall = re.search(r'host wall solves (.+) ms \(max over ranks; last is reported\)', log)
    times = [float(x) for x in wall[1].split()] if wall else []
    event = re.search(r'FCG iteration\s+[\d.]+ ms\s+\(solve ([\d.]+) ms, \+ setup ([\d.]+) ms', log)
    coarse = re.search(r'coarse solve\s+([\d.]+) ms per V-cycle', log)
    if len(times) != 3 or not event or not coarse:
        raise ValueError('Missing solve timings')
    if not all(math.isfinite(x) and x > 0 for x in times + [float(event[1]), float(coarse[1])]):
        raise ValueError('Invalid solve timings')
    # Rounded event printing is too coarse to compare with sub-millisecond smoke runs.
    if times[-1] > 10 and abs(times[-1] / float(event[1]) - 1) > 0.03:
        raise ValueError('Event and host solve times disagree')
    if 'Hypre: SpMV hypre, SpGEMM hypre, GPU-aware MPI on' not in log:
        raise ValueError('Wrong Hypre execution policy')
    count = len(list(directory.glob('rank*.json')))
    if count != (ranks if audited else 0):
        raise ValueError('Unexpected audit count')
    audits = [json.loads((directory / ('rank%d.json' % r)).read_text()) for r in range(ranks)] if audited else []
    limit = 0 if ranks == 1 or cutoff == 0 else min(50000, max(5000, ranks * 1000)) if cutoff == -1 else cutoff
    for audit in audits:
        if audit['version'] != 1 or audit['threshold_rows'] != limit or audit['setups'] != audit['destroyed'] or audit['skipped_setups'] != 0:
            raise ValueError('Coarse lifecycle or cutoff failed')
        if not limit:
            if audit['setups'] != 0 or audit['cycles'] != 0:
                raise ValueError('Disabled adapter ran')
        elif not (audit['setups'] == 1 and audit['cycles'] >= 3 * (it + 1) and 9 < audit['rows'] <= limit and audit['level'] > 0):
            raise ValueError('Expected coarse path did not run')
    hierarchy = {(a['level'], a['rows'], a['nnz'], a['cycles']) for a in audits}
    if len(hierarchy) > 1:
        raise ValueError('Rank audits disagree')
    result = dict(iterations=it, relative=relative, true_residual=residual, wall_ms=times[-1],
                  wall_solves_ms=times, setup_ms=float(event[2]), coarse_ms=float(coarse[1]),
                  hierarchy=list(next(iter(hierarchy))) if hierarchy else None)
    save(directory / 'summary.json', result)
    return result


def clean_environment():
    if os.environ.get('LD_PRELOAD'):
        raise ValueError('Unset LD_PRELOAD')
    for name in list(os.environ):
        if name.startswith('BP3_'):
            os.environ.pop(name)
    os.environ.update(MPICH_GPU_SUPPORT_ENABLED='1', MARS_HYPRE_SPMV_VENDOR='0',
                      MARS_HYPRE_SPGEMM_VENDOR='0', MARS_HYPRE_GPU_AWARE='1',
                      MARS_HYPRE_POOL_MAX_MIB='128', MARS_HYPRE_POOL_CACHE_MIB='2048',
                      CUDA_CACHE_PATH='/dev/shm/gandanie/cuda_cache', CUDA_MODULE_LOADING='LAZY',
                      PYTHONDONTWRITEBYTECODE='1')
