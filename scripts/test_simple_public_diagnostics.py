import copy
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

    def reference_doc(self, config):
        return {'simulation': {'solver': {'solver_control': {'advanced_options': {
            'linear_solver_settings': {'default': config}}}}}}

    def test_reference_pressure_target_and_resolution(self):
        config = {'family': 'HYPRE', 'rtol': 1e-4, 'atol': 1e-16}
        doc = self.reference_doc(config)
        text = HYPRE + ' rhs_norm=0.02\n' + SIMPLE
        result = diagnostic.summarize(io.StringIO(text), '255', doc)
        self.assertTrue(result['reference_pressure_target_looser'])
        self.assertTrue(result['reference_pressure_mars_residual_passed'])
        self.assertTrue(result['reference_pressure_hypre_residual_passed'])
        self.assertEqual(result['run_status'], 'failed')
        self.assertFalse(result['mars_passed'])
        library = doc['simulation']['solver']
        settings = library['solver_control']['advanced_options']['linear_solver_settings']
        settings['segregated_flow'] = {'family': 'Hypre', 'rtol': 1e-8}
        result = diagnostic.summarize(io.StringIO(text), '255', doc)
        self.assertFalse(result['reference_pressure_mars_residual_passed'])
        settings['pressure_correction'] = {'lookup': 'SECRET pressure settings'}
        library['SECRET pressure settings'] = copy.deepcopy(config)
        result = diagnostic.summarize(io.StringIO(text), '255', doc)
        self.assertTrue(result['reference_pressure_mars_residual_passed'])
        self.assertNotIn('SECRET', json.dumps(result))
        settings['pressure_correction']['lookup'] = 'secret pressure settings'
        result = diagnostic.summarize(io.StringIO(text), '255', doc)
        self.assertFalse(result['reference_pressure_resolved'])
        self.assertIsNone(result['reference_pressure_mars_residual_passed'])

    def test_reference_backend_and_scaling_guards(self):
        config = {'family': 'HYPRE', 'rtol': 1e-4}
        text = HYPRE + ' rhs_norm=0.02\n' + SIMPLE
        for changes in ({'normalize_matrix': True}, {'diagonal_scaling': True},
                        {'normalize_matrix': 'false'}, {'family': 'PETSc'}, {'family': 'Trilinos'},
                        {'family': 'SECRET'}, {'rtol': True}, {'rtol': 'nan'}, {'rtol': -1},
                        {'options': {}}, {'options': {'solver': {'type': 'SECRET'}}},
                        {'options': {'solver': {'type': 'MGR'}}}):
            candidate = dict(config, **changes)
            result = diagnostic.summarize(io.StringIO(text), '255', self.reference_doc(candidate))
            self.assertFalse(result['reference_pressure_original_norm_comparable'])
            self.assertIsNone(result['reference_pressure_mars_residual_passed'])
            self.assertNotIn('SECRET', json.dumps(result))

    def test_reference_boomer_ignores_absolute_tolerance(self):
        config = {'family': 'Hypre', 'rtol': 1e-8, 'atol': 1.0,
                  'options': {'solver': {'type': 'BoomerAMG', 'tol': 1.0}}}
        text = HYPRE + ' rhs_norm=0.02\n' + SIMPLE
        result = diagnostic.summarize(io.StringIO(text), '255', self.reference_doc(config))
        self.assertEqual(result['reference_pressure_solver'], 'BoomerAMG')
        self.assertFalse(result['reference_pressure_mars_residual_passed'])
        config['options']['solver']['type'] = 'FlexGMRES'
        result = diagnostic.summarize(io.StringIO(text), '255', self.reference_doc(config))
        self.assertTrue(result['reference_pressure_mars_residual_passed'])

    def test_reference_missing_conflicting_and_rounded_evidence(self):
        doc = self.reference_doc({'family': 'Hypre', 'rtol': 1e-4})
        for text in ('', HYPRE, SIMPLE.replace('pressure', 'momentum'),
                     SIMPLE.replace('rhs_norm=0.02', 'rhs_norm=0'),
                     SIMPLE.replace('rhs_norm=0.02', 'rhs_norm=nan'),
                     SIMPLE.replace('mars_absolute_residual=3e-7', 'mars_absolute_residual=2e-6'),
                     SIMPLE + '\n' + SIMPLE.replace('mars_absolute_residual=3e-7', 'mars_absolute_residual=1')):
            result = diagnostic.summarize(io.StringIO(text), '255', doc)
            self.assertIsNone(result['reference_pressure_mars_residual_passed'])
        # The source's smaller atol can make its mixed limit tighter despite a larger rtol.
        text = SIMPLE.replace('rhs_norm=0.02', 'rhs_norm=1e-20')
        result = diagnostic.summarize(io.StringIO(text), '255', doc)
        self.assertFalse(result['reference_pressure_target_looser'])

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

    def test_application_failures_without_linear_rejection(self):
        cases = {
            'outlet_anchor_or_moment_error_seen': 'all outlet faces closed: no open pressure anchor or nonfinite outlet moments',
            'nonfinite_diagnostics_error_seen': 'nonfinite nonlinear diagnostics',
            'continuity_consistency_error_seen': 'assembled continuity does not match boundary mass flux',
            'momentum_assembly_error_seen': 'momentum assembly failed',
            'flux_update_error_seen': 'flux update failed',
            'output_error_seen': 'field manifest output failed',
            'cuda_error_seen': 'an illegal memory access was encountered',
            'hypre_wrapper_error_seen': 'prepared Hypre: GMRES/AMG refresh failed',
            'halo_error_seen': 'node-field halo: MPI exchange failed',
            'linear_candidate_missing_seen': 'pressure correction at SIMPLE iteration 99: linear solver returned no usable candidate',
        }
        for key, message in cases.items():
            for suffix in ('', ' (on this rank)', ' (on another rank)'):
                with self.subTest(key=key, suffix=suffix):
                    result = summarize('ERROR: ' + message + suffix)
                    self.assertTrue(result[key])
                    self.assertTrue(result['application_error_seen'])
                    self.assertFalse(result['unclassified_application_error_seen'])
                    self.assertFalse(result['hypre_rejection_present'])
                    self.assertIsNone(result['solver_accepted'])
                    self.assertEqual(result['run_status'], 'failed')
        result = summarize('ERROR: all outlet faces closed: no open pressure anchor; cannot solve this prescribed-inflow case')
        self.assertTrue(result['outlet_anchor_or_moment_error_seen'])

    def test_unknown_application_errors_are_redacted(self):
        result = summarize('Rank 2 ERROR: SECRET-name /private/SECRET values=1234')
        self.assertTrue(result['application_error_seen'])
        self.assertTrue(result['unclassified_application_error_seen'])
        self.assertNotIn('SECRET', json.dumps(result))
        self.assertNotIn('1234', json.dumps(result))

    def test_exact_public_messages_and_first_error_survive(self):
        first = 'prepared Hypre: nonfinite result, extraction error, or CUDA failure'
        second = 'prepared Hypre: true residual evaluation failed'
        result = summarize('ERROR: ' + first + '\nERROR: ' + second + '\nERROR: ' + first)
        self.assertEqual(result['software_errors'], [first, second])
        self.assertEqual(result['first_software_error'], first)
        self.assertFalse(result['unclassified_application_error_seen'])
        # Even a public prefix followed by private text must not enter the report.
        result = summarize('ERROR: ' + first + ' SECRET=123\nERROR: prepared Hypre: SECRET')
        self.assertEqual(result['software_errors'], [])
        self.assertEqual(result['first_software_error'], 'unknown')
        self.assertNotIn('SECRET', json.dumps(result))
        self.assertNotIn('123', json.dumps(result))

    def test_incremental_report_matches_saved_log(self):
        lines = ['SIMPLE Tet4, 4 ranks (ElementDomain/cstone), upwind, laminar',
                 '[simple] iteration=50 momentum=SECRET continuity=SECRET', HYPRE, SIMPLE,
                 'ERROR: prepared Hypre: vector update failed', 'ERROR: SECRET']
        state = diagnostic.DiagnosticState()
        for line in lines:
            state.feed(line)
        self.assertEqual(state.result('255'), summarize('\n'.join(lines)))
        self.assertTrue(state.result('255')['solver_started'])
        self.assertTrue(state.result('255')['iteration_report_seen'])
        self.assertNotIn('SECRET', json.dumps(state.result('255')))

    def test_scheduler_and_mpi_markers(self):
        cases = {
            'scheduler_time_limit_seen': '[2026-09-30T12:00:00.001] error: *** STEP 123.0 ON SECRET CANCELLED DUE TO TIME LIMIT ***',
            'scheduler_out_of_memory_seen': 'slurmstepd: error: Detected 1 oom_kill event in StepId=123.0',
            'scheduler_signal_seen': 'srun: error: SECRET: task 0: Segmentation fault',
            'mpi_abort_seen': 'MPICH ERROR [Rank 0] [job id 123.0] [SECRET] - Abort(1): application called MPI_Abort(MPI_COMM_WORLD, 1)',
        }
        for key, message in cases.items():
            result = summarize(message)
            self.assertTrue(result[key])
            self.assertFalse(result['application_error_seen'])
            self.assertNotIn('SECRET', json.dumps(result))
        combined = summarize('\n'.join(cases.values()))
        for key in cases:
            self.assertTrue(combined[key])
        self.assertTrue(summarize('srun: error: SECRET: task 0: Out Of Memory')['scheduler_out_of_memory_seen'])
        self.assertTrue(summarize('application called MPI_Abort(MPI_COMM_WORLD, 1) - process 0')['mpi_abort_seen'])

    def test_scheduler_labels_and_exit_code_are_software_metadata(self):
        text = ('srun: error: SECRET: task 0: Segmentation fault\n'
                'srun: error: SECRET: tasks 1-3: Terminated\n'
                'srun: error: SECRET: task 0: Killed\n'
                'srun: error: SECRET: task 0: Killed\n')
        result = summarize(text, '139')
        self.assertEqual(result['scheduler_messages'], ['segmentation_fault', 'terminated', 'killed'])
        self.assertEqual(result['process_exit_code'], 139)
        self.assertNotIn('SECRET', json.dumps(result))
        for invalid in ('999', '-9', 'SECRET', ''):
            self.assertIsNone(summarize('', invalid)['process_exit_code'])

    def test_failure_markers_cannot_be_hidden_by_completion(self):
        converged = 'CONVERGED iterations=10 ranks=4 exchange_rounds=41\n'
        for text in ('ERROR: nonfinite nonlinear diagnostics', 'ERROR: SECRET unknown error',
                     'srun: error: *** STEP 123.0 CANCELLED DUE TO TIME LIMIT ***',
                     'application called MPI_Abort(MPI_COMM_WORLD, 1) - process 0'):
            self.assertEqual(summarize(text + '\n' + converged, '0')['run_status'], 'failed')

    def test_failure_markers_require_error_context(self):
        text = ('boundary: all outlet faces closed: no open pressure anchor or nonfinite outlet moments\n'
                'path=/private/nonfinite nonlinear diagnostics\n'
                'saved path: ERROR: an illegal memory access was encountered\n'
                'file named DUE TO TIME LIMIT or oom_kill or MPI_Abort\n'
                'slurmstepd: info: no oom_kill events\n'
                'srun: job 123 queued and waiting for resources\n')
        result = summarize(text)
        for key in diagnostic.FAILURE_FLAGS:
            self.assertFalse(result[key])
        result = summarize('ERROR: /private/nonfinite nonlinear diagnostics')
        self.assertTrue(result['unclassified_application_error_seen'])
        self.assertFalse(result['nonfinite_diagnostics_error_seen'])
        limited = ('NOT CONVERGED: iteration limit iterations=50 ranks=4 exchange_rounds=201\n'
                   'srun: error: SECRET: task 0: Exited with exit code 2\n')
        self.assertEqual(summarize(limited, '2')['run_status'], 'iteration_limit')

    def test_output_vocabulary_never_contains_input_text_or_numbers(self):
        allowed = {'mars-simple-public-diagnostics-v1', 'incomplete', 'failed',
                   'converged', 'iteration_limit', 'unknown', 'pressure', 'momentum',
                   'GMRES', 'FlexGMRES', 'native', 'vendor'}
        private = 'SECRET-path-boundary-value'
        text = (private + '\n' + HYPRE + ' ' + private + '=123456789\n' + SIMPLE
                + '\n[hypre-spmv] vendor_requested=' + private
                + '\nERROR: prepared Hypre: ' + private
                + '\nsrun: error: ' + private + ': Out Of Memory')
        for candidate in (text, text.replace('pressure', private).replace('GMRES', private)):
            result = summarize(candidate)
            self.assertNotIn(private, json.dumps(result))
            for key, value in result.items():
                if key == 'software_errors':
                    self.assertTrue(all(v in diagnostic.SOFTWARE_ERRORS for v in value))
                elif key == 'scheduler_messages':
                    self.assertTrue(all(v in diagnostic.SCHEDULER_MESSAGES.values() for v in value))
                elif key == 'process_exit_code':
                    self.assertTrue(value is None or type(value) is int and 0 <= value <= 255)
                else:
                    self.assertTrue(value is None or type(value) is bool
                                    or value in allowed or value in diagnostic.SOFTWARE_ERRORS)

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
            deck = root / 'SECRET-deck.json'
            deck.write_text(json.dumps(self.reference_doc({'family': 'Hypre', 'rtol': 1e-4})))
            reference_output = root / 'reference.json'
            reference_command = command[:-1] + [str(reference_output), '--reference-deck', str(deck)]
            checked = subprocess.run(reference_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            self.assertEqual(checked.returncode, 0, checked.stderr)
            self.assertTrue(json.loads(reference_output.read_text())['reference_pressure_target_looser'])
            self.assertNotIn('SECRET', reference_output.read_text())
            reference_output.unlink()
            deck.write_text('SECRET: [bad')
            checked = subprocess.run(reference_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            self.assertNotEqual(checked.returncode, 0)
            self.assertFalse(reference_output.exists())
            self.assertNotIn('SECRET', checked.stdout + checked.stderr)
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
