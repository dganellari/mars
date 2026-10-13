"""Synthetic OpenAccel reports only; no case files or solver launches."""
from fractions import Fraction
import json
import unittest

from simple_pressure_log import inspect_pressure_log, printed_interval, residual_assessment
from simple_snapshot_compare import EvidenceError


def report(step=1, start='1.23e-02', end='1.23e-08', count=5, cap=20):
    return '''Iter = {step}
Equation name:  Coupled Navier-Stokes
Solver context: PETSc
Iterations:     2/200
| U-Mom | 1.0e+00 | 1.0e-02 | 1.00e-02/1.00e-08 | 1.0e-06 | CONV |
Equation name:  Pressure Correction
Solver context: HYPRE
Iterations:     {count}/{cap}
+----------------+---------+---------+-------------------+---------+------+
|    Equation    |  Rate   |   RMS   | LinSol: Start/End |  Drop   | Conv.|
| P-Mass         | 1.0e+00 | 9.9e-01 | {start}/{end} | 9.9e+99 | CONV |
'''.format(step=step, start=start, end=end, count=count, cap=cap)


class PressureLogTests(unittest.TestCase):
    def check(self, text, steps=1, rtol=1e-4, atol=0, maxiter=20):
        return inspect_pressure_log(text, steps, rtol, atol, maxiter)

    def test_relative_pass_failure_and_rounding_ambiguity(self):
        for end, label in (('1.23e-08', 'printed_relative_target_met'),
                           ('1.23e-04', 'printed_target_exceeded'),
                           ('1.23e-06', 'inconclusive_absolute_tolerance_or_print_rounding')):
            with self.subTest(end=end):
                result = self.check(report(end=end))
                self.assertEqual(result['steps'][0]['residual_assessment'], label)
                self.assertEqual(result['any_inconclusive'], label.startswith('inconclusive'))

    def test_absolute_tolerance_without_row_count_is_conservative(self):
        self.assertEqual(residual_assessment('1.00e+00', '1.00e-04', 1e-6, 1e-3),
                         'inconclusive_absolute_tolerance_or_print_rounding')
        self.assertEqual(residual_assessment('1.00e+00', '1.00e-02', 1e-6, 1e-3),
                         'printed_target_exceeded')

    def test_decimal_rounding_intervals_are_exact_and_include_transition(self):
        self.assertEqual(printed_interval('1.23e-06'), (Fraction(1225, 10**9), Fraction(1235, 10**9)))
        low, high = printed_interval('1.00e+00')
        self.assertLess(low, Fraction('0.9999'))
        self.assertGreater(high, Fraction('1.004'))
        low, high = printed_interval('4.94e-324')
        self.assertLess(low, Fraction.from_float(float('4.94e-324')))
        self.assertGreater(high, Fraction.from_float(float('4.94e-324')))

    def test_zero_and_iteration_cap_do_not_use_backend_hint(self):
        for count in (0, 20, 21):
            result = self.check(report(start='0.00e+00', end='0.00e+00', count=count))
            self.assertTrue(result['all_printed_relative_targets_met'])
            self.assertEqual(result['any_iteration_cap_reached'], count >= 20)
            self.assertFalse(result['backend_convergence_verified'])
        result = self.check(report(end='1.00e-02', count=5))
        self.assertTrue(result['any_printed_target_exceeded'])
        self.assertFalse(result['any_iteration_cap_reached'])

    def test_incomplete_duplicate_reordered_or_wrong_reports_fail(self):
        first, second = report(), report(step=2)
        bad = (first + first, second + first, first + second.replace('Pressure Correction', 'Other'),
               first + second.replace('Solver context: HYPRE', 'Solver context: Trilinos'),
               first + second.replace('5/20', '5/200'),
               first + second.replace('| P-Mass', '| Other'),
               first + second.replace('Iter = 2', 'Iter = 3'),
               first + second + 'Equation name:  Pressure Correction\n',
               first + second.replace('Iterations:     5/20', 'Iterations:     5/20\nIterations:     5/20'),
               first + second.replace('Solver context: HYPRE', 'Solver context: HYPRE\nSolver context: HYPRE'))
        for text in bad:
            with self.subTest(text=text), self.assertRaisesRegex(EvidenceError, '^reference_pressure_log$'):
                self.check(text, steps=2)

    def test_invalid_residual_tokens_fail_with_fixed_label(self):
        for token in ('nan', 'inf', '-1.00e-04', '1.234e-04', '0.01e-02', '0.00e-04',
                      '1.00e+309', '1.00e-325', '9.99e+308', 'PRIVATE_NAME'):
            with self.subTest(token=token), self.assertRaisesRegex(EvidenceError, '^reference_pressure_log$'):
                self.check(report(end=token))

    def test_invalid_controls_fail(self):
        for key, value in (('steps', 0), ('maxiter', 0), ('rtol', 0), ('rtol', float('nan')),
                           ('atol', -1), ('atol', float('inf'))):
            with self.subTest(key=key, value=value), self.assertRaises(EvidenceError):
                self.check(report(), **{key: value})

    def test_full_history_has_only_classifications_and_no_exact_metrics(self):
        text = 'PRIVATE_CASE_PATH\n' + ''.join(report(step=i) for i in range(1, 21))
        result = self.check(text, steps=20)
        self.assertTrue(result['all_steps_present'])
        self.assertTrue(result['all_printed_relative_targets_met'])
        self.assertFalse(result['residuals_recomputed'])
        self.assertEqual(len(result['steps']), 20)
        for forbidden in ('PRIVATE_CASE_PATH', '1.23', '5/20', '9.9e'):
            self.assertNotIn(forbidden, json.dumps(result))


if __name__ == '__main__':
    unittest.main()
