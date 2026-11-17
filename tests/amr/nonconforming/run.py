#!/usr/bin/env python3
"""Run the fixed diagnostic matrix and preserve every failure."""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys


def cases():
    result = []
    for mode in ('mesh', 'full'):
        for repeat in range(1, 4):
            result.append((f'tet_{mode}_r8_repeat{repeat}', 8, 'tet', mode))
        result.append((f'tet_{mode}_r4', 4, 'tet', mode))
    result.append(('hex_full_r8', 8, 'hex', 'full'))
    return result


def mesh_hash(path):
    files = [path] if path.is_file() else sorted(p for p in path.rglob('*') if p.is_file())
    if not files or any(p.is_symlink() for p in files):
        raise ValueError(f'Missing or linked mesh input: {path}')
    inventory = {p.relative_to(path.parent).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    return hashlib.sha256(json.dumps(inventory, sort_keys=True).encode()).hexdigest()


def assess(text, code, element, mode):
    if code or re.search(r'\bFAIL\b|\b(?:nan|[-+]?inf)\b', text, re.I):
        return 'FAIL'
    if '[amr-galerkin] PASS' not in text:
        return 'INCOMPLETE'
    prefix = 'amr-mesh-only' if mode == 'mesh' else 'amr-galerkin'
    levels = re.findall(r'\[' + prefix + r'\] level=(\d+) elements=', text)
    if levels != ['0', '1', '2', '3']:
        return 'INCOMPLETE'
    if element == 'tet':
        if text.count('stage=connectivity-rebuilt ') != 3 or text.count(' expected=') != 3:
            return 'INCOMPLETE'
    if mode == 'full':
        for level in range(4):
            for field in ('adjoint_relative=', 'symmetry_relative=', 'true_residual_assembled='):
                if f'[amr-algebra-audit] level={level} {field}' not in text:
                    return 'INCOMPLETE'
            for trial in range(3):
                if f'level={level} trial={trial} operator_relative=' not in text:
                    return 'INCOMPLETE'
            for basis in range(4):
                if f'level={level} basis={basis} max=' not in text:
                    return 'INCOMPLETE'
    return 'PASS'


def main():
    kit, results = map(Path, sys.argv[1:])
    executable = results / 'build/examples/distributed/unstructured/mars_amr_galerkin'
    affinity = Path.home() / 'affinity/bind_numa.sh'
    if not os.access(executable, os.X_OK) or not os.access(affinity, os.X_OK):
        raise SystemExit('Missing executable or existing affinity launcher')
    receipts = []
    for name, ranks, element, mode in cases():
        case = results / name
        case.mkdir()
        mesh = kit / 'meshes' / ('cube16_tet.mesh' if element == 'tet' else 'cube16.mesh')
        command = ['srun', '--export=ALL', '--kill-on-bad-exit=1', '--time=10',
                   f'--nodes={(ranks + 3) // 4}', f'--ntasks={ranks}', '--ntasks-per-node=4',
                   '--cpus-per-task=8', str(affinity), str(executable), f'--mesh={mesh}',
                   f'--element={element}', '--problem=linear', '--levels=3']
        if element == 'tet':
            command.append('--audit-mesh')
        command += ['--mesh-only'] if mode == 'mesh' else ['--audit-algebra', '--diagnose']
        env = os.environ.copy()
        env.pop('MARS_AMR_AUDIT_TET', None)
        env['MARS_AMR_AUDIT_DIR'] = str(case / 'cells')
        env['MPICH_GPU_SUPPORT_ENABLED'] = '1'
        receipt = {'name': name, 'ranks': ranks, 'element': element, 'mode': mode, 'command': command,
                   'binary_sha256': hashlib.sha256(executable.read_bytes()).hexdigest(),
                   'mesh_inventory_sha256': mesh_hash(mesh)}
        with (case / 'run.log').open('w') as log:
            completed = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
        text = (case / 'run.log').read_text()
        receipt.update(exit_code=completed.returncode, status=assess(text, completed.returncode, element, mode))
        (case / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        receipts.append(receipt)
        print(name, receipt['status'], 'exit', completed.returncode, flush=True)
    (results / 'cases.json').write_text(json.dumps(receipts, indent=2) + '\n')
    return 0 if all(r['status'] == 'PASS' for r in receipts) else 1


if __name__ == '__main__':
    sys.exit(main())
