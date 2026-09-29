import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import simple_public_diagnostics as diagnostic


HYPRE = ('[HypreGMRES] rejected: backend=GMRES iterations=7/2000 restart=100 '
         'reported_relative=2e-5 tolerance=1e-12 solve_error=0 '
         'absolute_residual=3e-7 acceptance_limit=1e-11')
SIMPLE = ('[simple-linear] stage=pressure iteration=27 solver_accepted=0 '
          'mars_absolute_residual=3e-7 rhs_norm=0.02 acceptance_limit=1e-11 mars_passed=0')


def summarize(text, exit_text='255'):
    return diagnostic.summarize(io.StringIO(text), exit_text)


class PublicDiagnosticsTests(unittest.TestCase):
    def test_pressure_failure(self):
        result = summarize(HYPRE + '\n' + SIMPLE)
        self.assertEqual(result['linear_stage'], 'pressure')
        self.assertEqual(result['backend'], 'GMRES')
        self.assertEqual(result['run_status'], 'failed')
        for key in ('solver_accepted', 'mars_passed', 'iteration_cap_reached',
                    'hypre_error_nonzero', 'reported_residual_passed',
                    'hypre_explicit_residual_passed', 'mars_explicit_residual_passed'):
            self.assertIs(result[key], False)
        self.assertIs(result['reported_residual_finite'], True)

    def test_old_vendor_mismatch_stays_distinct(self):
        simple = SIMPLE.replace('stage=pressure', 'stage=momentum').replace('mars_passed=0', 'mars_passed=1').replace('mars_absolute_residual=3e-7', 'mars_absolute_residual=1e-16')
        result = summarize(HYPRE + '\n' + simple)
        self.assertIs(result['hypre_explicit_residual_passed'], False)
        self.assertIs(result['mars_explicit_residual_passed'], True)

    def test_unknowns_and_conflicts_do_not_pass(self):
        self.assertIsNone(summarize('')['solver_accepted'])
        result = summarize(HYPRE + '\n' + HYPRE.replace('solve_error=0', 'solve_error=256'))
        self.assertIsNone(result['hypre_error_nonzero'])
        result = summarize(HYPRE + ' iterations=2000/2000')
        self.assertIsNone(result['iteration_cap_reached'])
        result = summarize(SIMPLE + '\n' + SIMPLE.replace('stage=pressure', 'stage=momentum'))
        self.assertEqual(result['linear_stage'], 'unknown')

    def test_nonfinite_invalid_and_missing_values(self):
        for value in ('nan', 'inf', '-inf', '-1', 'secret', '1e999'):
            result = summarize(HYPRE.replace('absolute_residual=3e-7', 'absolute_residual=' + value))
            self.assertIsNot(result['hypre_explicit_residual_passed'], True)
        for text in ('iterations=bad', 'iterations=7/0'):
            self.assertIsNone(summarize(HYPRE.replace('iterations=7/2000', text))['iteration_cap_reached'])

    def test_cap_error_and_backend(self):
        text = HYPRE.replace('7/2000', '2000/2000').replace('solve_error=0', 'solve_error=256').replace('backend=GMRES', 'backend=FlexGMRES')
        result = summarize(text + '\n[hypre-spmv] vendor_requested=0\nfalse convergence 2, L2 norm of residual: 123')
        self.assertIs(result['iteration_cap_reached'], True)
        self.assertIs(result['hypre_error_nonzero'], True)
        self.assertEqual(result['backend'], 'FlexGMRES')
        self.assertEqual(result['spmv_backend'], 'native')
        self.assertIs(result['false_convergence_message_seen'], True)

    def test_pressure_audit_exports_only_flags(self):
        text = '[simple-pressure-audit] finite=1 zero_row=0 nonpositive_diagonal=0 positive_offdiagonal=1 constant_mode_detected=0 residual_within_roundoff_bound=1 roundoff_bound_exceeds_limit=1 compensated_residual_finite=1 compensated_residual_passed=0 secret=private'
        result = summarize(text)
        self.assertTrue(result['pressure_audit_present'])
        self.assertTrue(result['pressure_audit_finite'])
        self.assertFalse(result['pressure_audit_zero_row'])
        self.assertTrue(result['pressure_audit_roundoff_bound_exceeds_limit'])
        self.assertTrue(result['pressure_audit_compensated_residual_finite'])
        self.assertFalse(result['pressure_audit_compensated_residual_passed'])
        self.assertIsNone(summarize('')['pressure_audit_compensated_residual_passed'])
        self.assertNotIn('private', json.dumps(result))
        self.assertIsNone(summarize('')['pressure_audit_finite'])
        self.assertIsNone(summarize(text + '\n' + text.replace('finite=1', 'finite=0'))['pressure_audit_finite'])

    def test_completion_requires_matching_exit(self):
        converged = 'CONVERGED iterations=10 ranks=4 exchange_rounds=41\n'
        limited = 'NOT CONVERGED: iteration limit iterations=50 ranks=4 exchange_rounds=201\n'
        self.assertEqual(summarize(converged, '0')['run_status'], 'converged')
        self.assertEqual(summarize(limited, '2')['run_status'], 'iteration_limit')
        for text, code in ((converged, '2'), (limited, '0'), ('', '0'), (converged, 'bad')):
            self.assertEqual(summarize(text, code)['run_status'], 'incomplete')
        self.assertEqual(summarize(converged, '137')['run_status'], 'failed')
        self.assertEqual(summarize(HYPRE + '\n' + converged, '0')['run_status'], 'failed')
        self.assertTrue(summarize(converged * 2, '0')['multiple_completions'])
        self.assertEqual(summarize(converged * 2, '0')['run_status'], 'incomplete')

    def test_output_vocabulary_never_contains_input_text_or_numbers(self):
        allowed = {'mars-simple-public-diagnostics-v1', 'incomplete', 'failed',
                   'converged', 'iteration_limit', 'unknown', 'pressure', 'momentum',
                   'GMRES', 'FlexGMRES', 'native', 'vendor'}
        private = 'SECRET-path-boundary-value'
        text = (private + '\n' + HYPRE + ' ' + private + '=123456789\n' + SIMPLE
                + '\n[hypre-spmv] vendor_requested=' + private)
        for candidate in (text, text.replace('pressure', private).replace('GMRES', private)):
            result = summarize(candidate)
            self.assertNotIn(private, json.dumps(result))
            for value in result.values():
                self.assertTrue(value is None or type(value) is bool or value in allowed)

    def test_cli_private_errors_and_existing_output(self):
        script = str(Path(diagnostic.__file__).resolve())
        with tempfile.TemporaryDirectory() as work:
            root = Path(work)
            log, status, output = (root / name for name in ('SECRET.log', 'SECRET.exit', 'SECRET.json'))
            command = [sys.executable, script, '--log', str(log), '--exit-file', str(status), '--output', str(output)]
            log.write_text(HYPRE + '\n' + SIMPLE)
            status.write_text('255\n')
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            self.assertEqual(result.returncode, 0)
            original = output.read_bytes()
            for args in (command, command + ['--SECRET=secret'], command[:-4]):
                result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn('SECRET', result.stdout + result.stderr)
                self.assertNotIn(work, result.stdout + result.stderr)
                self.assertEqual(output.read_bytes(), original)
            log.unlink()
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertNotIn(work, result.stderr)


if __name__ == '__main__':
    unittest.main()
