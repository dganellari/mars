#!/usr/bin/env python3
"""Run and compare the pinned public SIMPLE channel; never use pump results."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

from prepare_openaccel_geometry import DECK_SHA256, MESH_SHA256

PIN = '0d69041ba1afda63e9e4328d9e0d9834bba37756'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def convergence_deck(text):
    changes = {'max_iterations: 2\n': 'max_iterations: 5000\n',
               'residual_target: 1.0e-6': 'residual_target: 1.0e-10',
               'rtol: 1.0e-8': 'rtol: 1.0e-12',
               'atol: 1.0e-12': 'atol: 1.0e-14'}
    for old, new in changes.items():
        require(text.count(old) == 1, 'unexpected pinned deck control: ' + old)
        text = text.replace(old, new)
    return text


def last_iteration(log):
    iterations = [int(i) for i in re.findall(r'^Iter = (\d+)\s*$', log, re.M)]
    require(iterations and iterations == list(range(1, iterations[-1] + 1)),
            'missing/duplicated reference iterations')
    require(re.search(r'^Converged\.\s*$', log, re.M), 'OpenAccel did not report convergence')
    require('Git hash: ' + PIN in log, 'reference executable is not at the pinned revision')
    return iterations[-1]


def run(capture, executable, output):
    capture, executable, output = [p.resolve() for p in (capture, executable, output)]
    require(not output.exists(), 'output exists; choose a fresh directory')
    record = json.loads((capture / 'run.json').read_text())
    require(record['fixture'] == 'public_channel' and record['returncode'] == 0
            and record['status'].endswith('_capture_completed'), 'not a completed public capture')
    require(record['bundle']['sha256']['input.i'] == DECK_SHA256
            and record['bundle']['sha256']['channel.exo'] == MESH_SHA256,
            'capture has another mesh/deck')
    require(digest(capture / 'input.i') == DECK_SHA256
            and digest(capture / 'channel.exo') == MESH_SHA256, 'public input checksum mismatch')
    require(executable.is_file() and os.access(str(executable), os.X_OK), 'executable missing')
    require(digest(executable) == record['binary_sha256'], 'use the executable from this completed capture')
    for key in ('SLURM_NTASKS', 'OMPI_COMM_WORLD_SIZE', 'PMI_SIZE', 'PMIX_SIZE'):
        require(int(os.environ.get(key, '1')) == 1, 'reference runner requires one rank')
    output.mkdir(parents=True)
    (output / 'input.i').write_text(convergence_deck((capture / 'input.i').read_text()))
    shutil.copyfile(str(capture / 'channel.exo'), str(output / 'channel.exo'))
    env = os.environ.copy()
    for key in ('MARS_OPENACCEL_EXPORT_DIR', 'MARS_OPENACCEL_PUBLIC_FIXTURE'):
        env.pop(key, None)
    env['OMP_NUM_THREADS'] = '1'
    result = dict(fixture='public_simple_convergence_v1', source_capture=str(capture),
                  source_deck_sha256=DECK_SHA256, mesh_sha256=MESH_SHA256,
                  deck_sha256=digest(output / 'input.i'), binary_sha256=digest(executable),
                  executable=str(executable), status='started')
    manifest = output / 'comparison-run.json'
    manifest.write_text(json.dumps(result, indent=2) + '\n')
    with (output / 'run.log').open('w') as log:
        process = subprocess.Popen([str(executable), '-i', 'input.i'], cwd=str(output),
                                   env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   universal_newlines=True, errors='replace')
        with process.stdout:
            for line in process.stdout:
                log.write(line)
                log.flush()
                print(line, end='', flush=True)
        result['returncode'] = process.wait()
    result['log_sha256'] = digest(output / 'run.log')
    result['status'] = 'failed'
    try:
        require(result['returncode'] == 0, 'OpenAccel failed; see run.log')
        result['iteration'] = last_iteration((output / 'run.log').read_text())
        files = list(output.glob('results.e*'))
        require(len(files) == 1 and files[0].is_file(), 'expected exactly one serial Exodus result')
        result['result_file'] = files[0].name
        result['result_sha256'] = digest(files[0])
        result['status'] = 'native_convergence_reported'
    finally:
        manifest.write_text(json.dumps(result, indent=2) + '\n')
    print('Reference complete. Run compare; native convergence alone is not field parity.')


def finite(np, data):
    require(not np.any(np.ma.getmaskarray(data)), 'missing/unwritten Exodus values')
    data = np.asarray(data, dtype=float)
    require(np.all(np.isfinite(data)), 'nonfinite field values')
    return data


def read_reference(path, iteration):
    import numpy as np
    from netCDF4 import Dataset, chartostring
    with Dataset(str(path)) as ds:
        require(len(ds.dimensions['num_nodes']) == 425, 'not the public channel node count')
        ids = finite(np, ds.variables['node_num_map'][:])
        require(np.all(ids == ids.astype('int64')) and len(set(ids)) == 425, 'invalid global IDs')
        if 'coord' in ds.variables:
            xyz = finite(np, ds.variables['coord'][:]).T
        else:
            xyz = np.column_stack([finite(np, ds.variables['coord' + a][:]) for a in 'xyz'])
        times = finite(np, ds.variables['time_whole'][:])
        require(times.ndim == 1 and len(times) >= 2 and np.all(np.diff(times) > 0)
                and times[-1] == iteration and times[-2] == iteration - 1,
                'final two output states do not match the converged iteration')
        names = [str(s).strip('\x00 ').lower() for s in
                 chartostring(np.ma.filled(ds.variables['name_nod_var'][:], b'\0'))]
        fields = []
        for name in ('velocity_x', 'velocity_y', 'velocity_z', 'pressure'):
            require(names.count(name) == 1, 'missing/ambiguous Exodus field: ' + name)
            k = names.index(name)
            if 'vals_nod_var' in ds.variables:
                v = ds.variables['vals_nod_var']
                require(v.dimensions == ('time_step', 'num_nod_var', 'num_nodes'), 'nodal dimensions')
                data = v[-2:, k, :]
            else:
                v = ds.variables['vals_nod_var' + str(k + 1)]
                require(v.dimensions == ('time_step', 'num_nodes'), 'nodal dimensions')
                data = v[-2:, :]
            fields.append(finite(np, data))
    require(xyz.shape == (425, 3), 'coordinate dimensions')
    return ids.astype('int64'), xyz, np.stack(fields, axis=2)


def compare(reference, mars, output):
    import numpy as np
    require(not output.exists(), 'comparison output exists')
    r = json.loads((reference / 'comparison-run.json').read_text())
    require(r['fixture'] == 'public_simple_convergence_v1' and r['status'] == 'native_convergence_reported'
            and r['returncode'] == 0 and r['mesh_sha256'] == MESH_SHA256
            and r['source_deck_sha256'] == DECK_SHA256, 'reference run is incomplete or unsupported')
    require(Path(r['result_file']).name == r['result_file'] and r['result_file'].startswith('results.e'),
            'invalid reference result filename')
    for filename, key in [('input.i', 'deck_sha256'), ('channel.exo', 'mesh_sha256'),
                          ('run.log', 'log_sha256'), (r['result_file'], 'result_sha256')]:
        require(digest(reference / filename) == r[key], 'reference file changed: ' + filename)
    require(last_iteration((reference / 'run.log').read_text()) == r['iteration'], 'reference iteration mismatch')
    meta = json.loads((mars / 'channel.json').read_text())
    require(digest(mars / 'channel.txt') == meta['packed_sha256'], 'MARS mesh checksum mismatch')
    require((mars / 'channel.txt').read_text().splitlines()[0] == 'MARS_PUBLIC_SIMPLE_MESH_V1 425 1536 576',
            'not public mesh-only input')
    ids = np.asarray(meta['node_global_ids'])
    require(ids.shape == (425,) and len(set(ids)) == 425 and np.all(ids == ids.astype('int64')),
            'MARS global ID map invalid')
    with (mars / 'channel-fields.csv').open() as f:
        rows = list(csv.DictReader(f))
    require(len(rows) == 425 and {int(row['node']) for row in rows} == set(range(425)), 'MARS row coverage')
    rows.sort(key=lambda row: int(row['node']))
    xyz = finite(np, [[row[k] for k in ('x', 'y', 'z')] for row in rows])
    values = finite(np, [[row[k] for k in ('u', 'v', 'w', 'p')] for row in rows])
    rid, rxyz, states = read_reference(reference / r['result_file'], r['iteration'])
    require(set(rid) == set(ids), 'global node ID sets differ')
    index = {n:i for i,n in enumerate(rid)}
    order = [index[n] for n in ids]
    require(np.max(np.abs(xyz - rxyz[order])) <= 1e-12, 'coordinates differ after global ID mapping')
    with (mars / 'channel-metrics.csv').open() as f:
        metrics = list(csv.DictReader(f))[-1]
    end = re.findall(r'^CONVERGED iterations=(\d+)\s*$', (mars / 'run.log').read_text(), re.M)
    require(len(end) == 1 and int(end[0]) == int(metrics['iteration']), 'MARS completion/metrics mismatch')
    for k in ('momentum', 'continuity', 'mass_balance', 'du', 'dp', 'dflux', 'cancellation'):
        x = float(metrics[k])
        require(np.isfinite(x) and 0 <= x <= (1e-10 if k == 'cancellation' else 1e-6), 'MARS gate failed: ' + k)
    require(int(metrics['changed_faces']) == 0, 'MARS outlet flags not stable')
    # Fixed physical scales make this an absolute test, including the pressure level.
    scale = np.array([0.1, 0.1, 0.1, 0.01])
    delta = (values - states[-1, order]) / scale
    drift = (states[-1, order] - states[-2, order]) / scale
    uerr = np.linalg.norm(delta[:, :3], axis=1)
    udrift = np.linalg.norm(drift[:, :3], axis=1)
    report = dict(scope='converged public single-rank nodal fields; no MPI or pump claim',
                  reference_iteration=r['iteration'], mars_iteration=int(end[0]),
                  velocity_max_scaled=float(uerr.max()), velocity_rms_scaled=float(np.sqrt(np.mean(uerr**2))),
                  pressure_max_scaled=float(np.abs(delta[:, 3]).max()),
                  pressure_rms_scaled=float(np.sqrt(np.mean(delta[:, 3]**2))),
                  reference_velocity_change_scaled=float(udrift.max()),
                  reference_pressure_change_scaled=float(np.abs(drift[:, 3]).max()),
                  pressure_mean_shift_pa=float(np.mean(values[:, 3] - states[-1, order, 3])),
                  field_tolerance=1e-5, reference_change_tolerance=1e-6,
                  hashes={str(p):digest(p) for p in [reference / 'comparison-run.json', mars / 'channel.json',
                          mars / 'channel-fields.csv', mars / 'channel-metrics.csv', mars / 'run.log']})
    passed = all(report[k] <= 1e-5 for k in ('velocity_max_scaled', 'pressure_max_scaled'))
    passed = passed and all(report[k] <= 1e-6 for k in ('reference_velocity_change_scaled', 'reference_pressure_change_scaled'))
    report['passed'] = bool(passed)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps(report, indent=2))
    require(passed, 'field parity or final reference change failed; inspect report, do not relax tolerances blindly')
    print('PASS: converged public SIMPLE velocity and absolute pressure agree by global node ID')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command')
    r = sub.add_parser('run')
    for name in ('capture', 'executable', 'output'):
        r.add_argument('--' + name, type=Path, required=True)
    c = sub.add_parser('compare')
    for name in ('reference', 'mars', 'output'):
        c.add_argument('--' + name, type=Path, required=True)
    args = p.parse_args()
    try:
        if args.command == 'run':
            run(args.capture, args.executable, args.output)
        elif args.command == 'compare':
            compare(args.reference, args.mars, args.output)
        else:
            p.error('choose run or compare')
    except (OSError, ValueError, KeyError, TypeError, IndexError, ImportError) as error:
        p.exit(1, 'ERROR: ' + str(error) + '\n')


if __name__ == '__main__':
    main()
