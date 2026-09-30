#!/usr/bin/env python3
"""Run four short, public-only SIMPLE references with one controlled change per axis."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from netCDF4 import Dataset

from openaccel_simple_convergence import PIN, convergence_deck, digest, read_reference, require
from prepare_openaccel_geometry import DECK_SHA256, MESH_SHA256

CASES = ('straight-average', 'straight-static', 'warped-average', 'warped-static')
STEPS = 20


def probe_deck(original, case):
    require(case in CASES, 'unknown public probe')
    text = convergence_deck(original, 'high-resolution', 'linear-linear')
    changes = {'min_iterations: 2\n': 'min_iterations: 20\n',
               'max_iterations: 5000\n': 'max_iterations: 20\n'}
    if case.endswith('-static'):
        changes['option: average_static_pressure'] = 'option: static_pressure'
        changes['                pressure_profile_blend: 0.05\n'] = ''
    for old, new in changes.items():
        require(text.count(old) == 1, 'unexpected public deck control')
        text = text.replace(old, new)
    # Historical input comments are not publication attribution for a new fixture.
    return '\n'.join(line for line in text.splitlines()
                     if not line.startswith('# Author:')) + '\n'


def warp(xyz):
    # A polynomial warp and proper rotation exercise nonplanar, oblique boundaries.
    points = np.array(xyz, dtype=float, copy=True)
    points[:, 0] += .2 * xyz[:, 1] ** 2 + .1 * xyz[:, 2] ** 2
    points[:, 1] += .05 * xyz[:, 0] ** 2
    rotation = np.array([[.8, -.6, 0], [.48, .64, -.6], [.36, .48, .8]])
    return points.dot(rotation.T)


def prepare_mesh(source, destination, case):
    require(case in CASES and digest(source) == MESH_SHA256,
            'only the pinned public mesh is accepted')
    require(not destination.exists(), 'mesh output exists')
    shutil.copyfile(str(source), str(destination))
    with Dataset(str(destination), 'r+') as ds:
        xyz = np.column_stack([ds.variables['coord' + a][:] for a in 'xyz'])
        cells = np.asarray(ds.variables['connect1'][:], dtype=np.int64) - 1
        require(xyz.shape == (425, 3) and cells.shape == (1536, 4), 'public mesh topology')
        if case.startswith('warped-'):
            transformed = warp(xyz)
            before = np.linalg.det(xyz[cells[:, 1:]] - xyz[cells[:, :1]])
            after = np.linalg.det(transformed[cells[:, 1:]] - transformed[cells[:, :1]])
            require(np.all(np.isfinite(after)) and np.all(after * before > 0),
                    'warp inverted or collapsed an element')
            for j, a in enumerate('xyz'):
                ds.variables['coord' + a][:] = transformed[:, j]


def capture_identity(capture, executable):
    record = json.loads((capture / 'run.json').read_text())
    require(record['fixture'] == 'public_channel' and record['returncode'] == 0
            and record['status'].endswith('_capture_completed'), 'incomplete public capture')
    require(record['bundle']['sha256']['input.i'] == DECK_SHA256
            and record['bundle']['sha256']['channel.exo'] == MESH_SHA256,
            'capture is not the pinned public fixture')
    require(digest(capture / 'input.i') == DECK_SHA256
            and digest(capture / 'channel.exo') == MESH_SHA256, 'public inputs changed')
    require(executable.is_file() and os.access(str(executable), os.X_OK), 'executable unavailable')
    require(digest(executable) == record['binary_sha256'], 'capture executable changed')
    return record['binary_sha256']


def check_result(directory):
    log = (directory / 'run.log').read_text()
    require(re.findall(r'^Iter = (\d+)\s*$', log, re.M) == [str(i) for i in range(1, STEPS + 1)],
            'reference did not complete the requested iterations')
    require('Git hash: ' + PIN in log, 'unexpected reference source revision')
    files = sorted(directory.glob('results.e*'))
    require(len(files) == 1 and files[0].is_file(), 'expected one serial Exodus output')
    with Dataset(str(files[0])) as ds:
        times = ds.variables['time_whole'][:]
        require(not np.any(np.ma.getmaskarray(times)) and np.all(np.isfinite(times)),
                'missing reference output times')
        require(np.array_equal(times, np.arange(1, STEPS + 1))
                or np.array_equal(times, np.arange(STEPS + 1)), 'missing or repeated output states')
    read_reference(files[0], STEPS)
    return files[0]


def run(capture, executable, output):
    capture, executable, output = (p.resolve() for p in (capture, executable, output))
    binary_hash = capture_identity(capture, executable)
    require(not output.exists(), 'choose a fresh output directory')
    for key in ('SLURM_NTASKS', 'OMPI_COMM_WORLD_SIZE', 'PMI_SIZE', 'PMIX_SIZE'):
        require(int(os.environ.get(key, '1')) == 1, 'reference probe requires one MPI rank')
    output.mkdir(parents=True)
    original = (capture / 'input.i').read_text()
    environment = os.environ.copy()
    for key in ('MARS_OPENACCEL_EXPORT_DIR', 'MARS_OPENACCEL_PUBLIC_FIXTURE'):
        environment.pop(key, None)
    environment['OMP_NUM_THREADS'] = '1'
    for case in CASES:
        directory = output / case
        directory.mkdir()
        prepare_mesh(capture / 'channel.exo', directory / 'channel.exo', case)
        (directory / 'input.i').write_text(probe_deck(original, case))
        command = [str(executable), '-i', 'input.i']
        record = dict(schema='mars-public-simple-geometry-probe-v1', case=case,
                      steps=STEPS, reference_revision=PIN, binary_sha256=binary_hash,
                      source_mesh_sha256=MESH_SHA256, source_deck_sha256=DECK_SHA256,
                      capture_record_sha256=digest(capture / 'run.json'), command=command,
                      mesh_sha256=digest(directory / 'channel.exo'),
                      deck_sha256=digest(directory / 'input.i'), status='started',
                      scope='public serial snapshots; no convergence, GPU or pump claim')
        manifest = directory / 'probe.json'
        manifest.write_text(json.dumps(record, indent=2) + '\n')
        with (directory / 'run.log').open('w') as log:
            result = subprocess.run(command, cwd=str(directory), env=environment,
                                    stdout=log, stderr=subprocess.STDOUT)
        record.update(returncode=result.returncode, status='failed', log_sha256=digest(directory / 'run.log'))
        (directory / 'run.exit').write_text(str(result.returncode) + '\n')
        try:
            require(result.returncode == 0, 'OpenAccel failed; inspect public run.log')
            result_file = check_result(directory)
            record.update(status='snapshots_completed', result_file=result_file.name,
                          result_sha256=digest(result_file))
        finally:
            manifest.write_text(json.dumps(record, indent=2) + '\n')
        print(case + ': 20 public snapshots saved', flush=True)
    print('Reference probes complete: ' + str(output), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('capture', 'executable', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    run(args.capture, args.executable, args.output)


if __name__ == '__main__':
    main()
