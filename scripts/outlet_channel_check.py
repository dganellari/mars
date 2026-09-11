#!/usr/bin/env python3
"""Compare public outlet-channel integration logs from 1/2/4 MPI ranks.

GPT/Codex, 2026-09-10. Scalars check partition agreement; they are not a pointwise
field comparison or an accuracy/convergence-order test. Uses the standard library.
"""
import argparse
import math
from pathlib import Path
import re


def fields(line):
    return {key: float(value) for key, value in re.findall(r'(\w+)=([^\s]+)', line)}


def read_run(text):
    configs = [fields(line) for line in text.splitlines() if line.startswith('[outlet-channel-config]')]
    if len(configs) != 1:
        raise ValueError('expected exactly one channel configuration')
    config = configs[0]
    config.setdefault('opening_area', 1.)
    config.setdefault('cut_check', 0.)
    required = {'ranks', 'steps', 'dt', 'nu', 'rho', 'inlet', 'ramp_steps', 'beta', 'p_ref', 'opening_area',
                'empty_opening_ranks'}
    if not required <= config.keys() or not all(math.isfinite(v) for v in config.values()):
        raise ValueError('incomplete or nonfinite configuration')
    ranks, steps = int(config['ranks']), int(config['steps'])
    if ranks not in (1, 2, 4) or steps < 2 or config['steps'] != steps or config['ranks'] != ranks:
        raise ValueError('expected 1/2/4 ranks and at least two complete steps')
    if (config['dt'] <= 0 or config['inlet'] <= 0
            or config['opening_area'] not in (1., .0625) or config['cut_check'] not in (0., 1.)):
        raise ValueError('invalid timestep or inlet speed')
    if not 0 <= config['empty_opening_ranks'] < ranks:
        raise ValueError('invalid empty-opening rank count')
    if text.count(f'PASS: public outlet channel integration steps={steps}\n') != 1:
        raise ValueError('missing or duplicated completion marker')
    if re.search(r'\bFAIL\b|\bERROR\b|\[diverged\]|Average-pressure outlet:', text):
        raise ValueError('failure reported in log')
    records = [fields(line) for line in text.splitlines() if line.startswith('[outlet-channel]')]
    if len(records) != steps:
        raise ValueError('missing or duplicated step records')
    required = {'step', 'ranks', 'bdf', 'dt_eff', 'rms', 'max', 'q_in', 'q_out', 'residual_sum',
                'u_mean', 'v_mean', 'w_mean', 'p_mean', 'u_rms', 'p_rms', 'trace_mean',
                'history_error', 'halo_error'}
    for step, record in enumerate(records, 1):
        if not required <= record.keys() or not all(math.isfinite(v) for v in record.values()):
            raise ValueError(f'step {step}: incomplete or nonfinite record')
        bdf = 1 if step == 1 else 2
        dt_eff = config['dt']*(1 if bdf == 1 else 2/3)
        ramp = min(1., step/config['ramp_steps']) if config['ramp_steps'] > 0 else 1.
        if record['step'] != step or record['ranks'] != ranks or record['bdf'] != bdf:
            raise ValueError(f'step {step}: wrong sequence, partition, or BDF stage')
        if not math.isclose(record['dt_eff'], dt_eff, rel_tol=1e-14, abs_tol=0.):
            raise ValueError(f'step {step}: wrong effective timestep')
        if not (0 <= record['rms'] <= 1e-7 and 0 <= record['max'] <= 1e-7):
            raise ValueError(f'step {step}: continuity tolerance failed')
        if (abs(record['q_in']+record['q_out']) > 1e-8 or record['q_out'] <= 0
                or abs(record['q_in']+config['inlet']*ramp*config['opening_area']) > 1e-10
                or abs(record['residual_sum']-record['q_in']-record['q_out']) > 1e-10):
            raise ValueError(f'step {step}: source or boundary balance failed')
        if record['history_error'] != 0 or record['halo_error'] != 0:
            raise ValueError(f'step {step}: history or halo mismatch')
        if abs(record['trace_mean']-config['p_ref']) > 1e-10:
            raise ValueError(f'step {step}: frozen trace mean mismatch')
    cuts = [fields(line) for line in text.splitlines() if line.startswith('[outlet-cut]')]
    if config['cut_check'] == 1 or cuts:
        if len(cuts) != steps:
            raise ValueError('missing or duplicated cut checks')
        for record, cut in zip(records, cuts):
            required_cut = {'step', 'q25', 'q50', 'q75', 'identity_error'}
            if (not required_cut <= cut.keys() or not all(math.isfinite(v) for v in cut.values())
                    or cut['step'] != record['step'] or not 0 <= cut['identity_error'] <= 1e-10
                    or any(abs(cut[k]+record['q_in']) > 1e-8 for k in ('q25', 'q50', 'q75'))):
                raise ValueError('cut flux or subset continuity identity failed')
    return config, records


def compare(texts, require_empty=False, coverage_text=None):
    runs = [read_run(text) for text in texts]
    if sorted(config['ranks'] for config, _ in runs) != [1, 2, 4]:
        raise ValueError('supply exactly one log for each of 1, 2, and 4 ranks')
    runs.sort(key=lambda run: run[0]['ranks'])
    base_config, base = runs[0]
    for config, records in runs[1:]:
        for key in ('steps', 'dt', 'nu', 'rho', 'inlet', 'ramp_steps', 'beta', 'p_ref', 'opening_area'):
            if config[key] != base_config[key]:
                raise ValueError(f'configuration differs between ranks: {key}')
        for a, b in zip(base, records):
            for key in ('q_in', 'q_out', 'u_mean', 'v_mean', 'w_mean', 'p_mean', 'u_rms', 'p_rms'):
                # Residuals have absolute gates; compare physical integrals at a looser solver-scale tolerance.
                if not math.isclose(a[key], b[key], rel_tol=5e-6, abs_tol=1e-8):
                    raise ValueError(f"ranks={config['ranks']:g} step={b['step']:g}: {key} partition mismatch")
    empty_covered = any(config['empty_opening_ranks'] > 0 for config, _ in runs)
    if coverage_text is not None:
        coverage, _ = read_run(coverage_text)
        if coverage['ranks'] != 4 or coverage['cut_check'] != 1 or coverage['empty_opening_ranks'] <= 0:
            raise ValueError('additional coverage run must exercise cuts and an empty-opening rank on four ranks')
        empty_covered = True
    if require_empty and not empty_covered:
        raise ValueError('no empty-opening rank exercised; additional partition coverage needed')
    return int(base_config['steps']), empty_covered


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('logs', nargs=3, type=Path)
    parser.add_argument('--require-empty', action='store_true')
    parser.add_argument('--coverage-log', type=Path, help='additional public four-rank run for empty-opening coverage')
    args = parser.parse_args()
    try:
        steps, empty = compare([path.read_text() for path in args.logs], args.require_empty,
                               args.coverage_log.read_text() if args.coverage_log else None)
    except (ValueError, OSError) as error:
        parser.exit(1, f'FAIL: {error}\n')
    print(f'PASS: {steps} steps, BDF startup, full continuity, halo/history, and 1/2/4-rank scalar agreement')
    if args.coverage_log:
        print('Additional public coverage run: cut flux, BDF startup, continuity and empty-opening rank passed')
    print('Empty-opening rank coverage: '+('exercised' if empty else 'NOT EXERCISED; do not claim this coverage'))
