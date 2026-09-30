import csv
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import simple_convergence_summary as summary


def fixture(count=100, value=lambda i: .01, status='iteration_limit'):
    rows = []
    for i in range(count + 1):
        row = dict.fromkeys(summary.COLUMNS, 0)
        row.update(iteration=i, momentum=value(i), inlet_kg_s=-987.654321,
                   outlet_kg_s=987.654321, umax_m_s=123.456789)
        rows.append(row)
    return rows, status


def serialize(rows, status='iteration_limit'):
    metrics = io.StringIO()
    writer = csv.DictWriter(metrics, fieldnames=summary.COLUMNS)
    writer.writeheader()
    writer.writerows(rows)
    log = ['SIMPLE Tet4, 4 ranks (ElementDomain/cstone), high-resolution, laminar',
           'PRIVATE boundary name and path /do/not/export']
    for row in rows:
        log.append('[simple] ' + ' '.join('{}={}'.format(key, row[source])
                                         for key, source in summary.LOG_COLUMNS.items()))
    ending = 'CONVERGED' if status == 'converged' else 'NOT CONVERGED: iteration limit'
    log.append('{} iterations={} ranks=4 exchange_rounds=401'.format(ending, rows[-1]['iteration']))
    return '\n'.join(log), metrics.getvalue(), '0' if status == 'converged' else '2'


def analyze(rows, status='iteration_limit', residual=1e-6, mass=1e-6, change=1e-6):
    log, metrics, code = serialize(rows, status)
    return summary.summarize(io.StringIO(log), io.StringIO(metrics), code, residual, mass, change)


class ConvergenceSummaryTests(unittest.TestCase):
    def test_flat_residual_with_small_updates(self):
        rows, status = fixture()
        result = analyze(rows, status)
        self.assertEqual(result['run_status'], 'iteration_limit')
        self.assertFalse(result['final_momentum_passed'])
        self.assertTrue(result['final_continuity_passed'])
        self.assertTrue(result['completion_matches_supplied_targets'])
        self.assertTrue(result['small_updates_with_residual_failure'])
        self.assertEqual(result['tail_momentum_trend'], 'flat')
        self.assertEqual(result['tail_dp_trend'], 'below_target')

    def test_decreasing_increasing_and_mixed_trends(self):
        for function, expected in ((lambda i: .97**i, 'decreasing'),
                                   (lambda i: 1.03**i, 'increasing'),
                                   (lambda i: 2 if i <= 80 or i > 90 else 1, 'mixed')):
            rows, status = fixture(value=function)
            self.assertEqual(analyze(rows, status)['tail_momentum_trend'], expected)

    def test_below_target_requires_whole_last_block(self):
        rows, status = fixture(value=lambda i: 1e-8)
        self.assertEqual(analyze(rows, status)['tail_momentum_trend'], 'below_target')
        rows[95]['momentum'] = 1e-4
        result = analyze(rows, status)
        self.assertNotEqual(result['tail_momentum_trend'], 'below_target')
        self.assertTrue(result['final_momentum_passed'])

    def test_large_finite_values_do_not_overflow_block_medians(self):
        rows, status = fixture(value=lambda i: 1.7e308)
        self.assertEqual(analyze(rows, status)['tail_momentum_trend'], 'flat')

    def test_exact_targets_and_every_stopping_criterion(self):
        rows, _ = fixture(value=lambda i: 1e-6)
        for key in ('continuity', 'mass_balance', 'du', 'dp', 'dflux'):
            rows[-1][key] = 1e-6
        rows[-1]['cancellation'] = 1e-10
        baseline = analyze(rows, 'converged')
        self.assertTrue(baseline['all_final_checks_passed'])
        self.assertTrue(baseline['completion_matches_supplied_targets'])
        for key, check in (('momentum', 'momentum'), ('continuity', 'continuity'),
                           ('mass_balance', 'mass_balance'), ('du', 'du'), ('dp', 'dp'),
                           ('dflux', 'dflux'), ('cancellation', 'conservation'),
                           ('changed_faces', 'stable_outlet_flags')):
            altered = [dict(row) for row in rows]
            altered[-1][key] = 1
            result = analyze(altered)
            self.assertFalse(result['final_' + check + '_passed'])
            self.assertFalse(result['all_final_checks_passed'])

    def test_targets_are_supplied_not_inferred_from_completion(self):
        rows, status = fixture()
        result = analyze(rows, status, residual=1)
        self.assertTrue(result['all_final_checks_passed'])
        self.assertFalse(result['completion_matches_supplied_targets'])
        self.assertEqual(result['run_status'], 'iteration_limit')
        self.assertEqual(result['targets_source'], 'command_line')

    def test_reversal_and_short_history(self):
        rows, status = fixture()
        rows[-10]['changed_faces'] = 5
        rows[-1]['closed_faces'] = 9
        result = analyze(rows, status)
        self.assertTrue(result['outlet_flags_changed_in_tail'])
        self.assertTrue(result['closed_outlet_faces_present'])
        self.assertTrue(result['final_stable_outlet_flags_passed'])
        rows, status = fixture(count=1)
        result = analyze(rows, status)
        self.assertFalse(result['final_minimum_iterations_passed'])
        self.assertEqual(result['tail_momentum_trend'], 'insufficient_data')
        self.assertIsNone(result['outlet_flags_changed_in_tail'])

    def test_nonfinite_negative_and_bad_thresholds_fail(self):
        for value in (float('nan'), float('inf'), -1):
            for key in ('momentum', 'mass_balance', 'du', 'umax_m_s'):
                rows, status = fixture()
                rows[-1][key] = value
                with self.subTest(value=value, key=key), self.assertRaises(ValueError):
                    analyze(rows, status)
        rows, status = fixture()
        for key in ('residual', 'mass', 'change'):
            for value in (0, -1, float('nan'), float('inf')):
                with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                    analyze(rows, status, **{key: value})

    def test_missing_duplicate_or_reordered_evidence_fails(self):
        rows, status = fixture()
        log, metrics, code = serialize(rows, status)
        cases = [('', metrics, code), (log, '', code), (log, metrics, '143'),
                 (log, metrics, '0'), (log + '\n' + log, metrics, code),
                 (log.replace('ranks=4', 'ranks=2'), metrics, code),
                 (log.replace('iterations=100', 'iterations=99'), metrics, code),
                 (log, metrics.replace('iteration,', 'secret,'), code),
                 (log, '\n'.join(metrics.splitlines()[:-1]), code),
                 (log, '\n'.join(metrics.splitlines()[:2] + metrics.splitlines()[3:]), code),
                 (log + '\nERROR: prepared Hypre: vector update failed', metrics, code),
                 (log.replace('momentum=0.01', 'momentum=0.01 momentum=0.02'), metrics, code),
                 (log, metrics.replace('100,0.01,', '100,0.01,SECRET,'), code)]
        for candidate_log, candidate_metrics, candidate_code in cases:
            with self.assertRaises(ValueError):
                summary.summarize(io.StringIO(candidate_log), io.StringIO(candidate_metrics),
                                  candidate_code, 1e-6, 1e-6, 1e-6)

    def test_final_and_earlier_logged_rows_must_match_csv(self):
        rows, status = fixture()
        log, metrics, code = serialize(rows, status)
        for position in (10, 100):
            changed = [dict(row) for row in rows]
            changed[position]['momentum'] = .125
            wrong_metrics = serialize(changed, status)[1]
            with self.assertRaises(ValueError):
                summary.summarize(io.StringIO(log), io.StringIO(wrong_metrics), code, 1e-6, 1e-6, 1e-6)

    def test_export_has_only_public_labels_and_booleans(self):
        rows, status = fixture()
        result = analyze(rows, status)
        allowed = {'mars-simple-convergence-v1', 'iteration_limit', 'command_line',
                   'flat', 'below_target'}
        for value in result.values():
            self.assertTrue(type(value) is bool or value is None or value in allowed)
        output = json.dumps(result)
        for secret in ('PRIVATE', '987.654321', '123.456789', '/do/not/export'):
            self.assertNotIn(secret, output)

    def test_cli_roundtrip_and_failure_redaction(self):
        rows, status = fixture()
        log, metrics, code = serialize(rows, status)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name, text in (('PRIVATE.log', log), ('PRIVATE.csv', metrics), ('PRIVATE.exit', code)):
                (root / name).write_text(text)
            command = [sys.executable, str(Path(summary.__file__).resolve()),
                       '--log', str(root / 'PRIVATE.log'), '--metrics', str(root / 'PRIVATE.csv'),
                       '--exit-file', str(root / 'PRIVATE.exit'), '--residual-tol', '1e-6',
                       '--mass-tol', '1e-6', '--change-tol', '1e-6', '--output', str(root / 'out.json')]
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(json.loads((root / 'out.json').read_text()), analyze(rows, status))
            original = (root / 'out.json').read_bytes()
            for altered in (command, command + ['--PRIVATE', 'SECRET'],
                            command[:-1] + [str(root / 'new.json')] + ['--residual-tol', 'SECRET']):
                result = subprocess.run(altered, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn('PRIVATE', result.stderr)
                self.assertNotIn('SECRET', result.stderr)
            self.assertEqual((root / 'out.json').read_bytes(), original)
            self.assertFalse((root / 'new.json').exists())


if __name__ == '__main__':
    unittest.main()
