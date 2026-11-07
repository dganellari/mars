"""Build and measure a frozen MARS/MFEM pair; failed samples remain in the report."""
import json
import math
from pathlib import Path
import re
import shutil
import statistics
import subprocess
import sys

HERE = Path(__file__).resolve().parent
HCOARSE = HERE.parent / 'hcoarse' if (HERE.parent / 'hcoarse').is_dir() else HERE.parent / 'bp3_hcoarse'
sys.path.insert(0, str(HCOARSE))
import compat
import integrate
from common import (CONFIG_SHA, HYPRE, HYPRE_SHA, clean_environment, command, digest, parse, save,
                    step, verify_bundle, verify_files)
from receipts import parse_h
from plan import VARIANTS, cases, depth, mars_args, mfem_args, validate_case


def mfem_parse(out):
    text = (out / 'solver.log').read_text()
    found = re.findall(r'CG \+ V-cycle: (\d+) iterations, sqrt\(<r,z>/<r0,z0>\) = (\S+), true residual (\S+)', text)
    if len(found) != 1 or not re.search(r'\bPASS\b', text) or re.search(r'\bFAIL\b|did not reach tol|above --true-tol', text):
        raise ValueError('Missing or failed MFEM convergence gate')
    it, relative, true = int(found[0][0]), float(found[0][1]), float(found[0][2])
    if not (0 < it < 2000 and 0 <= relative <= 1e-13 and 0 <= true <= 1e-10):
        raise ValueError('MFEM accuracy gate failed')
    if "Hypre: MARS's settings" not in text or 'GPU-aware MPI on' not in text:
        raise ValueError('Wrong MFEM Hypre policy')
    found = re.findall(r'^  solves (.+) ms \(the last is reported\)', text, re.M)
    times = [float(v) for v in found[0].split()] if len(found) == 1 else []
    setup = re.findall(r'\(solve [\d.]+ ms, \+ setup ([\d.]+) ms\)', text)
    if len(times) != 5 or len(setup) != 1 or not all(math.isfinite(t) and t > 0 for t in times):
        raise ValueError('Missing MFEM timings')
    result = dict(iterations=it, relative=relative, true_residual=true, wall_ms=times[-1],
                  wall_solves_ms=times, setup_ms=float(setup[0]))
    save(out / 'summary.json', result)
    return result


def executable(prepared, variant, checked=False):
    if variant == 'mfem':
        return prepared / 'mfem/mfem_bp3'
    return integrate.binary(prepared, 'check' if checked else 'on')


def solve(bundle, prepared, root, name, variant, ranks, p, ne, checked=False):
    out = root / name
    out.mkdir()
    exe = executable(prepared, variant, checked)
    args = mfem_args(p, ne) if variant == 'mfem' else mars_args(bundle, variant, p, ne, checked)
    command(step(ranks, 20) + ['/usr/bin/python3', bundle / 'kit/rank.py', out,
            'no' if variant == 'mfem' else 'yes', exe] + args, out / 'solver.log', timeout=1250)
    if variant == 'mfem':
        row = mfem_parse(out)
    elif depth(variant, p):
        row = parse_h(out, ranks, ne, depth(variant, p), checked, solves=5)
    else:
        row = parse(out, ranks, 0 if variant == 'matched' else -1, solves=5)
    text = (out / 'solver.log').read_text()
    if variant != 'mfem':
        if re.findall(r'BP3 h selection: requested=(\w+) levels=(\d+)', text) != [('explicit', str(depth(variant, p)))]:
            raise ValueError('Wrong h selection')
        if not checked and ('PASS h validation:' in text or 'BP3 h memory after setup:' in text):
            raise ValueError('Diagnostics entered a timing binary')
    build = json.loads((prepared / 'paper-build.json').read_text())
    hosts = []
    for rank in range(ranks):
        process = json.loads((out / ('process%d.json' % rank)).read_text())
        if process['rank'] != rank or process['returncode'] or process['status'] != 'complete' or process['binary_sha256'] != build['binaries'][str(exe.relative_to(prepared))]:
            raise ValueError('Wrong process receipt')
        if variant != 'mfem' and process['audit_sha256'] != digest(out / ('rank%d.json' % rank)):
            raise ValueError('Changed audit')
        hosts.append(process['hostname'])
    if len(set(hosts)) != (ranks + 3) // 4:
        raise ValueError('Unexpected node placement')
    row.update(hosts=hosts, variant=variant, ranks=ranks, p=p, ne=ne, h_depth=depth(variant, p), checked=checked)
    save(out / 'summary.json', row)
    print('PASS ' + name, flush=True)
    return row


def build_mfem(bundle, root):
    out = root / 'mfem'
    out.mkdir()
    prefix = out / 'install'
    command(['cmake', '-S', bundle / 'dependencies/mfem', '-B', out / 'build',
             '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_INSTALL_PREFIX=' + str(prefix),
             '-DCMAKE_CXX_COMPILER=' + str(shutil.which('mpicxx')),
             '-DCMAKE_CUDA_HOST_COMPILER=' + str(shutil.which('mpicxx')),
             '-DMFEM_USE_MPI=YES', '-DMFEM_USE_METIS=NO', '-DMFEM_USE_CUDA=YES', '-DCUDA_ARCH=sm_90',
             '-DHYPRE_DIR=' + str(HYPRE.parent.parent), '-DMFEM_ENABLE_EXAMPLES=NO',
             '-DMFEM_ENABLE_MINIAPPS=NO', '-DMFEM_ENABLE_TESTING=NO'], out / 'configure.log', timeout=600)
    command(['cmake', '--build', out / 'build', '--parallel', '8'], out / 'compile.log', timeout=6000)
    command(['cmake', '--install', out / 'build'], out / 'install.log', timeout=600)
    command(['cmake', '-E', 'chdir', out, 'make', '-f', bundle / 'source/tests/performance/mfem_bp3/Makefile',
             'MFEM_DIR=' + str(prefix)], out / 'driver.log', timeout=600)
    link = command(['ldd', '-r', out / 'mfem_bp3'], out / 'linkage.txt')
    libs = re.findall(r'libHYPRE[^\s]*\s+=>\s+(\S+)', link)
    if 'not found' in link or 'undefined symbol:' in link or len(libs) != 1 or Path(libs[0]).resolve() != HYPRE.resolve():
        raise ValueError('MFEM resolves a different Hypre')
    save(out / 'binary.json', dict(sha256=digest(out / 'mfem_bp3')))


def seal(root, bundle, name, passed=True):
    files = {str(p.relative_to(root)): digest(p) for p in root.rglob('*') if p.is_file() and
             p.suffix in ('.json', '.log', '.txt') and not {'CMakeFiles', '_deps', 'mfem'}.intersection(p.relative_to(root).parts)}
    # MFEM's build files are numerous; retain its explicit logs and binary receipt.
    if (root / 'mfem').exists():
        files.update({str(p.relative_to(root)): digest(p) for p in (root / 'mfem').iterdir()
                      if p.is_file() and p.suffix in ('.json', '.log', '.txt')})
    save(root / name, dict(passed=passed, manifest_sha256=digest(bundle / 'manifest.json'), evidence=files))


def verify_prepared(prepared, bundle):
    record = json.loads((prepared / 'prepared.json').read_text())
    if not record['passed'] or record['manifest_sha256'] != digest(bundle / 'manifest.json'):
        raise ValueError('Wrong preparation')
    verify_files(prepared, record['evidence'])
    verify_files(prepared, json.loads((prepared / 'paper-build.json').read_text())['binaries'])


def prepare(bundle, root):
    compat.check(bundle, root)
    command(step(1, 115) + ['/usr/bin/python3', HERE / 'driver.py', 'build-mfem', bundle, root], root / 'mfem-build.log')
    built = json.loads((root / 'build/build.json').read_text())
    binaries = {'build/' + n: sha for n, sha in built['binaries'].items()}
    binaries['mfem/mfem_bp3'] = digest(root / 'mfem/mfem_bp3')
    save(root / 'paper-build.json', dict(binaries=binaries, manifest_sha256=digest(bundle / 'manifest.json'),
                                      hypre_sha256=HYPRE_SHA, config_sha256=CONFIG_SHA))
    # Small checks include every new baseline setting before any large allocation.
    for ranks in (1, 4):
        for p in (3, 7):
            for variant in VARIANTS:
                solve(bundle, root, root, 'baseline_%s_p%d_r%d' % (variant, p, ranks), variant, ranks, p, 32)
    verify_bundle(bundle)
    seal(root, bundle, 'prepared.json')
    print('PASS preparation: combined operators and matched MFEM baseline', flush=True)


def sample_set(bundle, prepared, root, workload, ranks, p, ne):
    rows = {v: [] for v in VARIANTS}
    failures = []
    metadata = validate_case(ranks, p, ne)
    for rd in range(3):
        order = VARIANTS[rd:] + VARIANTS[:rd]
        for variant in order:
            if any(f.get('variant') == variant for f in failures):
                continue
            name = '%s_r%d_p%d_%s_%d' % (workload, ranks, p, variant, rd + 1)
            try:
                rows[variant].append(solve(bundle, prepared, root, name, variant, ranks, p, ne))
            except (ValueError, RuntimeError, OSError, subprocess.TimeoutExpired) as error:
                failures.append(dict(name=name, variant=variant, error=str(error)))
                print('FAIL ' + name + ': ' + str(error), flush=True)
    completed = [r for values in rows.values() for r in values]
    placement_changed = bool(completed) and any(r['hosts'] != completed[0]['hosts'] for r in completed)
    if placement_changed:
        failures.append(dict(name=workload, error='Rank placement changed within the comparison'))
    summary = {}
    for variant, values in rows.items():
        if len(values) == 3 and not placement_changed:
            times = [r['wall_ms'] for r in values]
            summary[variant] = dict(median_ms=statistics.median(times), min_ms=min(times), max_ms=max(times),
                                    iterations=[r['iterations'] for r in values],
                                    cold_ms=[r['setup_ms'] + r['wall_solves_ms'][0] for r in values],
                                    true_residuals=[r['true_residual'] for r in values],
                                    iteration_ms=[r['wall_ms'] / r['iterations'] for r in values],
                                    coarse_ms=[r.get('coarse_ms') for r in values])
    return dict(workload=workload, ranks=ranks, p=p, ne=ne, metadata=metadata, samples=rows,
                summary=summary, failures=failures)


def benchmark(bundle, prepared, root, phase):
    verify_prepared(prepared, bundle)
    if phase == 'large':
        # No timed 256-rank experiment runs unless all pilot solves and direct checks pass.
        for p in (3, 7):
            for variant in VARIANTS:
                solve(bundle, prepared, root, 'pilot_%s_p%d' % (variant, p), variant, 256, p, 64)
            checked = solve(bundle, prepared, root, 'pilot_check_p%d' % p, 'best_h', 256, p, 64, checked=True)
            normal = json.loads((root / ('pilot_best_h_p%d' % p) / 'summary.json').read_text())
            if not integrate.checked_equivalent(normal, checked):
                raise ValueError('256-rank direct check changes the hierarchy')
        seal(root, bundle, 'pilot-passed.json')
    reports = [sample_set(bundle, prepared, root, *case) for case in cases(phase)]
    verify_bundle(bundle)
    verify_prepared(prepared, bundle)
    failures = [f for c in reports for f in c['failures']]
    save(root / 'results.json', dict(passed=not failures, phase=phase, cases=reports, failures=failures,
                                    manifest_sha256=digest(bundle / 'manifest.json')))
    seal(root, bundle, 'evidence.json', passed=not failures)
    if failures:
        raise RuntimeError('%d failed samples retained; do not use them as timings' % len(failures))
    print('PASS ' + phase, flush=True)


if __name__ == '__main__':
    clean_environment()
    if len(sys.argv) not in (4, 5) or sys.argv[1] not in ('prepare', 'standard', 'large', 'build-mfem'):
        raise SystemExit('usage: driver.py prepare|standard|large|build-mfem bundle results [prepared]')
    phase, bundle, root = sys.argv[1], Path(sys.argv[2]).resolve(), Path(sys.argv[3]).resolve()
    verify_bundle(bundle)
    if phase == 'build-mfem':
        build_mfem(bundle, root)
    else:
        root.mkdir(parents=True, exist_ok=False)
        save(root / 'inputs.json', dict(phase=phase, manifest_sha256=digest(bundle / 'manifest.json'),
                                       prepared=sys.argv[4] if len(sys.argv) == 5 else None))
        if phase == 'prepare':
            prepare(bundle, root)
        else:
            if len(sys.argv) != 5:
                raise SystemExit('benchmark phase needs a prepared result directory')
            benchmark(bundle, Path(sys.argv[4]).resolve(), root, phase)
