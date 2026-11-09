"""Record every node's driver and toolchain before the timed solves start."""
import csv
import json
import os
from pathlib import Path
import shutil
import socket
import sys

HERE = Path(__file__).resolve().parent
KIT = HERE.parent / 'kit' if (HERE.parent / 'kit').is_dir() else HERE.parents[1] / 'bp3_coarse/integration'
sys.path.insert(0, str(KIT))
from common import CONFIG_SHA, HYPRE, HYPRE_SHA, IMAGE, command, digest, save, step


def node(directory):
    rank = int(os.environ['SLURM_PROCID'])
    prefix = directory / ('node%d' % rank)
    hypre = digest(HYPRE)
    config = digest(HYPRE.parent.parent / 'include/HYPRE_config.h')
    if (hypre, config) != (HYPRE_SHA, CONFIG_SHA):
        raise ValueError('Runtime Hypre changed')
    text = command(['nvidia-smi', '--query-gpu=name,uuid,driver_version', '--format=csv,noheader'],
                   prefix.with_suffix('.gpus.txt'), timeout=30)
    gpus = []
    for row in csv.reader(text.splitlines(), skipinitialspace=True):
        if len(row) != 3 or not all(v.strip() for v in row):
            raise ValueError('Incomplete GPU identity')
        gpus.append(dict(zip(('name', 'uuid', 'driver'), (v.strip() for v in row))))
    tools = {}
    for name in ('nvcc', 'mpicxx', 'cmake'):
        executable = shutil.which(name)
        if executable is None:
            raise ValueError('Missing build tool: ' + name)
        version = command([executable, '--version'], prefix.with_suffix('.' + name + '.txt'), timeout=30)
        tools[name] = dict(path=executable, version=version.strip())
    save(prefix.with_suffix('.json'), dict(rank=rank, hostname=socket.gethostname(), gpus=gpus,
         tools=tools, image=IMAGE, hypre_sha256=hypre, config_sha256=config,
         visibility={k: os.environ.get(k) for k in ('CUDA_VISIBLE_DEVICES', 'SLURM_LOCALID')}))


def validate(rows, count, expected=None):
    if len(rows) != count or [r['rank'] for r in rows] != list(range(count)):
        raise ValueError('Missing node runtime receipt')
    if len({r['hostname'] for r in rows}) != count:
        raise ValueError('Runtime probe did not visit every node')
    signatures = []
    uuids = []
    for row in rows:
        gpus = row['gpus']
        if len(gpus) != 4 or len({(g['name'], g['driver']) for g in gpus}) != 1:
            raise ValueError('Expected four matching GPUs per node')
        uuids += [g['uuid'] for g in gpus]
        if (row['hypre_sha256'], row['config_sha256'], row['image']) != (HYPRE_SHA, CONFIG_SHA, IMAGE):
            raise ValueError('Runtime dependency differs from the pinned input')
        signatures.append(dict(gpu=gpus[0]['name'], driver=gpus[0]['driver'], tools=row['tools'],
                               image=row['image'], hypre_sha256=row['hypre_sha256'],
                               config_sha256=row['config_sha256']))
    if len(set(uuids)) != 4 * count:
        raise ValueError('Duplicate GPU identities')
    signature = signatures[0]
    if any(s != signature for s in signatures) or (expected is not None and signature != expected):
        raise ValueError('Driver or toolchain differs across the paper allocations')
    return signature


def capture(bundle, root, nodes, prepared=None):
    if int(os.environ.get('SLURM_JOB_NUM_NODES', '0')) != nodes:
        raise ValueError('Wrong allocation size for this phase')
    directory = root / 'runtime'
    directory.mkdir()
    # One probe per node avoids four identical nvidia-smi requests on each node.
    args = step(nodes * 4, 5)
    args[args.index('-n' + str(nodes * 4))] = '-n' + str(nodes)
    args[args.index('--ntasks-per-node=4')] = '--ntasks-per-node=1'
    command(args + ['/usr/bin/python3', HERE / 'runtime.py', directory], root / 'runtime.log', timeout=360)
    rows = [json.loads((directory / ('node%d.json' % r)).read_text()) for r in range(nodes)]
    expected = json.loads((prepared / 'runtime.json').read_text())['signature'] if prepared else None
    signature = validate(rows, nodes, expected)
    save(root / 'runtime.json', dict(signature=signature, nodes=rows, job_id=os.environ['SLURM_JOB_ID'],
                                   manifest_sha256=digest(bundle / 'manifest.json')))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: runtime.py receipt-directory')
    node(Path(sys.argv[1]))
