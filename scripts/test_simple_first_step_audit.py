"""Synthetic matrices and fields only; no private case is used by these tests."""
import copy
import csv
import contextlib
import io
import json
import unittest

import numpy as np
from netCDF4 import Dataset

import simple_first_step_audit as audit
import simple_startup_probe as probe
import test_simple_snapshot_compare as fixtures
import test_simple_startup_probe as startup_tests


class FirstStepTests(unittest.TestCase):
    record = startup_tests.StartupTests.record
    compare = startup_tests.StartupTests.compare

    def setUp(self):
        self.fixture = fixtures.SnapshotTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.tearDown)
        f = self.fixture
        self.pair = f.root / 'first-step'
        probe.prepare(f.case, f.reference, self.pair, first_step_audit=True)
        self.mars = self.pair / 'mars'
        self.mars.mkdir()
        self.public = f.root / 'public-audit.json'
        self.numbering = np.array([5, 2, 4, 0, 1, 3])
        self.predictor = np.arange(18).reshape(6, 3)*.002 + .03
        self.influence = np.tile([.1, .12, .14], (6, 1))
        self.gradient = np.arange(18).reshape(6, 3)*.001
        self.phi = (np.arange(6)*.01).reshape(6, 1)
        self.final = np.column_stack([self.predictor-self.influence*self.gradient, .3*self.phi])
        self.matrices = {'momentum': np.eye(18)*3 + .01*np.ones((18, 18)),
                         'pressure': np.eye(6)*2 + .02*np.ones((6, 6))}
        self.solutions = {'momentum': self.predictor, 'pressure': self.phi}
        self.rhs = {name: matrix @ self.solutions[name].ravel() for name, matrix in self.matrices.items()}
        for step in (0, 1):
            self.write_fields(self.mars / ('flow-step-{}-fields.csv'.format(step)), np.zeros((6, 4)) if step == 0 else self.final)
        self.write_fields(self.mars / 'flow-fields.csv', self.final)
        log = f.log.read_text().replace('iterations=100', 'iterations=1')
        (self.mars / 'run.log').write_text(log)
        (self.mars / 'run.exit').write_text('2\n')
        with f.metrics.open() as stream:
            metrics = list(csv.DictReader(stream))[:2]
        for step, row in enumerate(metrics):
            row['umax_m_s'] = str(np.max(np.linalg.norm(self.final[:, :3], axis=1)) if step else 0.)
        with (self.mars / 'flow-metrics.csv').open('w') as stream:
            writer = csv.DictWriter(stream, list(metrics[0])); writer.writeheader(); writer.writerows(metrics)
        reference = self.pair / 'reference'
        (reference / 'run.log').write_text('Iter = 1\nSimulation is complete\n')
        (reference / 'run.exit').write_text('0\n')
        self.write_reference_fields()
        self.write_reference_matrices()
        self.write_mars_parts()
        self.refresh()

    def refresh(self):
        for solver in ('mars', 'openaccel'):
            self.record(solver)
        for path in (self.public, self.pair / 'comparison-private.json', self.pair / 'first-step-private.json'):
            if path.exists(): path.unlink()

    def write_fields(self, path, values):
        with path.open('w') as stream:
            writer = csv.writer(stream); writer.writerow(['node', 'x', 'y', 'z', 'u', 'v', 'w', 'p'])
            for n in range(6): writer.writerow([n] + list(self.fixture.xyz[n]) + list(values[n]))

    def write_reference_fields(self):
        names = ('velocity_x', 'velocity_y', 'velocity_z', 'pressure', 'aux', 'du_x', 'du_y', 'du_z',
                 'pressure_correction', 'pressure_correction_gradient_x', 'pressure_correction_gradient_y', 'pressure_correction_gradient_z')
        values = np.column_stack([self.final, self.numbering, self.influence, self.phi, self.gradient])
        for rank, nodes in enumerate(([4, 0, 2, 1], [3, 1, 5, 2])):
            path = self.pair / ('reference/results.e.2.' + str(rank))
            with Dataset(str(path), 'w') as ds:
                for name, count in (('num_nodes', len(nodes)), ('num_dim', 3), ('time_step', 2), ('num_nod_var', len(names)), ('len_name', 64)):
                    ds.createDimension(name, count)
                ds.createVariable('node_num_map', 'i8', ('num_nodes',))[:] = self.fixture.ids[nodes]
                ds.createVariable('coord', 'f8', ('num_dim', 'num_nodes'))[:] = self.fixture.xyz[nodes].T
                ds.createVariable('time_whole', 'f8', ('time_step',))[:] = [0., 1.]
                ds.createVariable('name_nod_var', 'S1', ('num_nod_var', 'len_name'))[:] = np.asarray(names, dtype='S64').view('S1').reshape(len(names), 64)
                field = ds.createVariable('vals_nod_var', 'f8', ('time_step', 'num_nod_var', 'num_nodes'))
                field[0] = 0.; field[1] = values[nodes].T

    def write_reference_matrices(self, width=4, row_factor=None):
        mapping = np.argsort(self.numbering)
        for stage, system, c in (('momentum', 'coupled_navier_stokes', 3), ('pressure', 'pressure_correction', 1)):
            dofs = (mapping[:, None]*c + np.arange(c)).ravel()
            # Different positive row factors exercise normalization and diagonal scaling.
            factors = np.arange(len(dofs))*.07 + .4
            if row_factor is not None:
                factors[:] = row_factor
            matrix = self.matrices[stage][dofs][:, dofs]*factors[:, None]
            prefix = self.pair / ('reference/' + system + '_0000')
            data = {'rows': (np.arange(len(dofs)+1)*len(dofs), '<i'+str(width)),
                    'cols': (np.tile(np.arange(len(dofs)), len(dofs)), '<i'+str(width)),
                    'vals': (matrix.ravel(), '<f8'), 'b': (self.rhs[stage][dofs]*factors, '<f8')}
            for name, (values, dtype) in data.items():
                with open(str(prefix) + '_' + name + '.bin', 'wb') as stream: np.asarray(values, dtype=dtype).tofile(stream)

    def write_mars_parts(self, fault=None):
        for rank, source in enumerate((np.array([5, 3, 1, 4, 2, 0]), np.array([2, 0, 4, 1, 5, 3]))):
            owned = np.flatnonzero(source % 2 == rank)
            if fault == 'duplicate_owner' and rank: owned = np.array([0, 1, 2])
            arrays = [(np.array([0x4d53415544495431, 1, rank, 2, 6, 3, 36]), '<u8'),
                      (source, '<i4'), (owned, '<i4'), (np.arange(7)*6, '<i4'), (np.tile(np.arange(6), 6), '<i4')]
            for stage, c in (('momentum', 3), ('pressure', 1)):
                dofs = (source[:, None]*c + np.arange(c)).ravel()
                matrix = self.matrices[stage][dofs][:, dofs].reshape(6, c, 6, c).transpose(0, 2, 1, 3).copy()
                rhs = self.rhs[stage][dofs].reshape(6, c).copy()
                # Ghost rows may be partial or poisoned and must not enter the comparison.
                ghosts = np.setdiff1d(np.arange(6), owned)
                matrix[ghosts] = np.nan; rhs[ghosts] = np.nan
                if fault == stage + '_matrix' and rank == 0: matrix[owned[0], 0, 0, 0] += .02
                if fault == stage + '_rhs' and rank == 0: rhs[owned[0], 0] += .02
                arrays += [(matrix, '<f8'), (rhs, '<f8'), (self.solutions[stage][source], '<f8')]
                if c == 3: arrays += [(self.predictor[source], '<f8'), (self.influence[source], '<f8')]
            arrays.append((self.gradient[source], '<f8'))
            with (self.mars / 'flow-audit-rank{:06d}.bin'.format(rank)).open('wb') as stream:
                for values, dtype in arrays: np.asarray(values, dtype=dtype).tofile(stream)

    def test_complete_pair_permutations_scaling_and_ghost_rows(self):
        code, result = self.compare()
        self.assertEqual(code, 0, result)
        self.assertEqual(result['first_differing_stage'], 'none')
        self.assertTrue(all(result['first_step_stage_matches'].values()))
        self.assertTrue(all(result['first_step_linear_accuracy'].values()))
        self.assertNotIn('all_twenty_snapshots_match', result)
        self.assertNotIn('matrix_max_row_scaled', self.public.read_text())

    def test_64_bit_reference_indices(self):
        self.write_reference_matrices(8); self.refresh()
        self.assertEqual(self.compare()[0], 0)

    def test_common_pressure_system_with_nonuniform_row_scaling(self):
        code, result = self.compare()
        self.assertEqual(code, 0, result)
        checks = result['pressure_common_system_checks']
        self.assertEqual(checks['row_scaling'], 'nonuniform')
        self.assertFalse(checks['unscaled_matrix_matches'])
        self.assertFalse(checks['unscaled_rhs_matches'])
        self.assertTrue(checks['both_solutions_meet_reference_target'])
        self.assertTrue(checks['reference_solution_meets_mars_runtime_target'])
        self.assertEqual(checks['assessment'], 'increments_match')

    def test_common_pressure_system_unit_and_uniform_scaling(self):
        for factor, expected in ((1., 'unit'), (7., 'uniform_nonunit')):
            with self.subTest(factor=factor):
                self.write_reference_matrices(row_factor=factor); self.refresh()
                code, result = self.compare()
                self.assertEqual(code, 0, result)
                checks = result['pressure_common_system_checks']
                self.assertEqual(checks['row_scaling'], expected)
                self.assertEqual(checks['unscaled_matrix_matches'], factor == 1.)
                self.assertEqual(checks['unscaled_rhs_matches'], factor == 1.)
                self.assertTrue(checks['both_solutions_meet_reference_target'])

    def test_distinct_increments_can_both_meet_same_reference_target(self):
        self.phi += 1e-6
        self.change_mars_final()
        code, result = self.compare()
        self.assertEqual(code, 0, result)
        self.assertFalse(result['first_step_stage_matches']['pressure_increment'])
        checks = result['pressure_common_system_checks']
        self.assertTrue(checks['both_solutions_meet_reference_target'])
        self.assertEqual(checks['assessment'], 'distinct_increments_pass_common_fp64_check')

    def test_own_residuals_can_pass_while_common_absolute_target_fails(self):
        reference_phi = self.phi.copy()
        self.write_reference_matrices(row_factor=1000.)
        self.phi += 1e-6
        self.write_mars_parts()
        parts = [audit.MarsPart(self.mars / 'flow-audit-rank{:06d}.bin'.format(rank), rank, 2, 6)
                 for rank in (0, 1)]
        reference = audit.ReferenceMatrix(self.pair / 'reference', 'pressure_correction', 1,
                                          np.argsort(self.numbering))
        system = audit.compare_system(parts, reference, 'pressure', self.phi, reference_phi, 1.)
        targets = dict(limits=dict(reference=1e-5, runtime=1e-5))
        self.assertLess(system['mars']['absolute_residual'], 1e-5)
        self.assertLess(system['reference']['absolute_residual'], 1e-5)
        checks = audit.pressure_common_system_checks(system, targets, False)
        self.assertEqual(checks['row_scaling'], 'uniform_nonunit')
        self.assertEqual(checks['assessment'], 'mars_candidate_fails_reference_target')

    def test_unscaled_rhs_overflow_does_not_report_zero_difference(self):
        columns = np.array([0])
        with np.errstate(over='ignore'):
            with self.assertRaisesRegex(ValueError, 'audit_residual_arithmetic'):
                audit.unscaled_row_difference((columns, np.array([1e-308]), 1.),
                                              (columns, np.array([1e-308]), .5), 1e308)

    def test_common_system_detects_candidate_outside_reference_target(self):
        self.phi += .001
        self.change_mars_final()
        code, result = self.compare()
        self.assertEqual(code, 0, result)
        checks = result['pressure_common_system_checks']
        self.assertTrue(checks['reference_solution_meets_reference_target'])
        self.assertFalse(checks['mars_solution_meets_reference_target'])
        self.assertEqual(checks['assessment'], 'mars_candidate_fails_reference_target')

    def test_reference_target_failure_is_not_diagnosed_as_tolerance_ambiguity(self):
        self.phi += .001
        self.final = np.column_stack([self.predictor-self.influence*self.gradient, .3*self.phi])
        self.write_reference_fields()
        self.phi += .001
        self.change_mars_final()
        code, result = self.compare()
        self.assertEqual(code, 0, result)
        checks = result['pressure_common_system_checks']
        self.assertFalse(checks['reference_solution_meets_reference_target'])
        self.assertEqual(checks['assessment'], 'reference_target_not_met')

    def test_common_checks_do_not_hide_bad_equations_or_ghost_copies(self):
        for fault in ('pressure_matrix', 'pressure_rhs'):
            with self.subTest(fault=fault):
                self.write_mars_parts(fault); self.refresh()
                code, result = self.compare()
                self.assertEqual(code, 0, result)
                checks = result['pressure_common_system_checks']
                self.assertEqual(checks['row_scaling'], 'not_equivalent')
                self.assertEqual(checks['assessment'], 'row_scaled_equivalence_failed')
        self.write_mars_parts(); self.change_pressure_ghost()
        code, result = self.compare()
        self.assertEqual(code, 0, result)
        self.assertEqual(result['pressure_common_system_checks']['assessment'], 'referenced_copies_differ')

    def test_each_assembly_error_is_localized(self):
        for fault in ('momentum_matrix', 'momentum_rhs', 'pressure_matrix', 'pressure_rhs'):
            with self.subTest(fault=fault):
                self.write_mars_parts(fault); self.refresh()
                code, result = self.compare()
                self.assertEqual(code, 0, result)
                self.assertEqual(result['first_differing_stage'], fault)

    def test_matching_fields_do_not_hide_bad_linear_equation(self):
        self.rhs['pressure'] += .1
        self.write_reference_matrices(); self.write_mars_parts(); self.refresh()
        code, result = self.compare()
        self.assertEqual(code, 0, result)
        self.assertTrue(all(result['first_step_stage_matches'].values()))
        self.assertFalse(result['first_step_linear_accuracy']['pressure_reference_relative_residual_below_1e_8'])
        self.assertFalse(result['first_step_linear_accuracy']['pressure_mars_relative_residual_below_1e_8'])

    def change_mars_final(self):
        values = np.column_stack([self.predictor-self.influence*self.gradient, .3*self.phi])
        for name in ('flow-step-1-fields.csv', 'flow-fields.csv'):
            self.write_fields(self.mars / name, values)
        path = self.mars / 'flow-metrics.csv'
        with path.open() as stream: rows = list(csv.DictReader(stream))
        rows[-1]['umax_m_s'] = str(np.max(np.linalg.norm(values[:, :3], axis=1)))
        with path.open('w') as stream:
            writer = csv.DictWriter(stream, list(rows[0])); writer.writeheader(); writer.writerows(rows)
        self.write_mars_parts(); self.refresh()

    def test_predictor_solve_difference_with_equal_assembly(self):
        self.predictor += .001
        self.change_mars_final()
        code, result = self.compare()
        self.assertEqual(code, 0, result)
        self.assertEqual(result['first_differing_stage'], 'momentum_predictor')
        self.assertTrue(result['first_step_stage_matches']['momentum_matrix'])
        self.assertTrue(result['first_step_stage_matches']['momentum_rhs'])
        self.assertFalse(result['first_step_linear_accuracy']['momentum_mars_relative_residual_below_1e_8'])

    def test_pressure_solve_difference_with_equal_assembly(self):
        self.phi += .001
        self.change_mars_final()
        code, result = self.compare()
        self.assertEqual(code, 0, result)
        self.assertEqual(result['first_differing_stage'], 'pressure_increment')
        self.assertTrue(result['first_step_stage_matches']['pressure_matrix'])
        self.assertTrue(result['first_step_stage_matches']['pressure_rhs'])
        self.assertFalse(result['first_step_linear_accuracy']['pressure_mars_relative_residual_below_1e_8'])

    def test_pressure_can_meet_declared_target_and_fail_common_threshold(self):
        self.phi += 1e-6
        self.final = np.column_stack([self.predictor-self.influence*self.gradient, .3*self.phi])
        self.write_reference_fields()
        self.change_mars_final()
        code, result = self.compare()
        self.assertEqual(code, 0, result)
        checks = result['pressure_solve_checks']
        self.assertTrue(checks['same_declared_tolerances'])
        self.assertTrue(checks['mars_owner_residual_meets_runtime_limit'])
        self.assertTrue(checks['mars_local_residual_meets_runtime_limit'])
        for name in ('mars', 'reference'):
            self.assertTrue(checks[name + '_residual_below_declared_unpreconditioned_limit'])
            self.assertTrue(checks[name + '_declared_limit_looser_than_common_relative_check'])
            self.assertFalse(result['first_step_linear_accuracy']['pressure_' + name + '_relative_residual_below_1e_8'])
        self.assertFalse(checks['reference_backend_convergence_verified'])

    def change_pressure_ghost(self, match_local_rhs=False, value=.05):
        path = self.mars / 'flow-audit-rank000000.bin'
        part = audit.MarsPart(path, 0, 2, 6)
        phi = part.phi.copy()
        ghost = np.setdiff1d(np.arange(6), part.owned)[0]
        phi[ghost] += value
        with path.open('r+b') as stream:
            stream.seek(part.phi.ctypes.data - part.data.ctypes.data)
            phi.tofile(stream)
            if match_local_rhs:
                rhs = part.pressure_rhs.copy()
                for local in part.owned:
                    start, end = part.offsets[local:local+2]
                    rhs[local] = part.pressure[start:end, 0, 0] @ phi[part.columns[start:end], 0]
                stream.seek(part.pressure_rhs.ctypes.data - part.data.ctypes.data)
                rhs.tofile(stream)
        self.refresh()

    def test_bad_ghost_is_visible_even_when_owner_residual_passes(self):
        self.change_pressure_ghost()
        code, result = self.compare()
        self.assertEqual(code, 0, result)
        checks = result['pressure_solve_checks']
        self.assertFalse(checks['mars_referenced_copies_equal_owners'])
        self.assertTrue(checks['mars_owner_residual_meets_runtime_limit'])
        self.assertFalse(checks['mars_local_residual_meets_runtime_limit'])

    def test_local_residual_can_pass_with_bad_owner_reconstruction(self):
        self.change_pressure_ghost(match_local_rhs=True)
        code, result = self.compare()
        self.assertEqual(code, 0, result)
        checks = result['pressure_solve_checks']
        self.assertFalse(checks['mars_referenced_copies_equal_owners'])
        self.assertFalse(checks['mars_owner_residual_meets_runtime_limit'])
        self.assertTrue(checks['mars_local_residual_meets_runtime_limit'])

    def test_nonfinite_referenced_ghost_rejects_evidence(self):
        self.change_pressure_ghost(value=float('nan'))
        code, result = self.compare()
        self.assertEqual(code, 1, result)
        self.assertEqual(result['failed_check'], 'audit_residual_arithmetic')

    def test_default_absolute_floor_is_distinct_from_stopping_tolerance(self):
        norms = dict(absolute_residual=5e-14, rhs_norm=1e-12)
        system = dict(mars=norms, mars_local=norms, reference=norms, mars_referenced_copies_equal_owners=True)
        checks, _ = audit.pressure_accuracy(self.pair, {}, system)
        self.assertEqual(checks['mars_runtime_target_source'], 'default_true_residual_check')
        self.assertTrue(checks['mars_owner_residual_meets_runtime_limit'])
        self.assertFalse(checks['mars_residual_below_declared_unpreconditioned_limit'])
        self.assertTrue(checks['mars_runtime_limit_looser_than_common_relative_check'])

    def test_reference_controls_are_not_reported_as_backend_convergence(self):
        import yaml
        path = self.pair / 'reference/input.i'
        doc = probe.load_deck(path.read_bytes())
        config = doc['simulation']['solver']['solver_control']['advanced_options']['linear_solver_settings']['pressure_correction']
        norms = dict(absolute_residual=1e-5, rhs_norm=1.)
        system = dict(mars=norms, mars_local=norms, reference=norms, mars_referenced_copies_equal_owners=True)
        for family, expected in (('PETSc', 'petsc'), ('PRIVATE', 'other')):
            config.update(family=family, diagonal_scaling=True)
            path.write_text(yaml.safe_dump(doc))
            checks, _ = audit.pressure_accuracy(self.pair, {}, system)
            self.assertEqual(checks['reference_family'], expected)
            self.assertTrue(checks['reference_diagonal_scaling'])
            self.assertFalse(checks['reference_backend_convergence_verified'])
            self.assertNotIn('PRIVATE', json.dumps(checks))

    def test_reanalysis_preserves_capture_and_previous_reports(self):
        self.assertEqual(self.compare()[0], 0)
        before = {path: probe.digest(path) for path in self.pair.rglob('*') if path.is_file()}
        detail = self.fixture.root / 'new-details'
        public = self.fixture.root / 'new-public.json'
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            code = probe.main(['compare', '--pair', str(self.pair), '--detail-dir', str(detail), '--output', str(public)])
        self.assertEqual(code, 0, public.read_text())
        self.assertEqual(before, {path: probe.digest(path) for path in self.pair.rglob('*') if path.is_file()})
        self.assertTrue((detail / 'first-step-private.json').is_file())
        text = public.read_text()
        self.assertNotIn(str(self.fixture.root), text)
        self.assertNotIn('absolute_residual', text)
        self.assertNotIn('pressure_targets', text)

    def test_reanalysis_still_rejects_modified_capture(self):
        self.assertEqual(self.compare()[0], 0)
        path = self.mars / 'flow-audit-rank000000.bin'
        with path.open('ab') as stream:
            stream.write(b'PRIVATE')
        public = self.fixture.root / 'tampered-public.json'
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            code = probe.main(['compare', '--pair', str(self.pair), '--detail-dir', str(self.fixture.root / 'details'),
                               '--output', str(public)])
        self.assertEqual(code, 1)
        self.assertEqual(json.loads(public.read_text())['failed_check'], 'launch_records')
        self.assertNotIn('PRIVATE', public.read_text())

    def test_missing_or_truncated_data_never_passes(self):
        for name in ('mars/flow-audit-rank000001.bin', 'reference/pressure_correction_0000_vals.bin'):
            path = self.pair / name; original = path.read_bytes()
            for content in (None, original[:-8]):
                with self.subTest(name=name, content=content is None):
                    if content is None: path.unlink()
                    else: path.write_bytes(content)
                    self.refresh()
                    code, result = self.compare()
                    self.assertEqual(code, 1)
                    self.assertEqual(result['comparison_status'], 'invalid_evidence')
                    path.write_bytes(original)

    def test_duplicate_owner_rejected(self):
        self.write_mars_parts('duplicate_owner'); self.refresh()
        self.assertEqual(self.compare()[0], 1)

    def test_reference_numbering_must_be_a_permutation(self):
        self.numbering[0] = self.numbering[1]
        self.write_reference_fields(); self.refresh()
        self.assertEqual(self.compare()[0], 1)

    def test_missing_intermediate_field_is_identified(self):
        path = self.pair / 'reference/results.e.2.0'
        with Dataset(str(path), 'a') as ds:
            ds.variables['name_nod_var'][4] = np.asarray(['PRIVATE'], dtype='S64').view('S1')
        self.refresh()
        code, result = self.compare()
        self.assertEqual(code, 1)
        self.assertEqual(result['failed_check'], 'reference_field_names')
        self.assertNotIn('PRIVATE', self.public.read_text())

    def test_nonfinite_owned_matrix_is_rejected(self):
        self.matrices['momentum'][0, 0] = np.nan
        self.write_mars_parts(); self.refresh()
        code, result = self.compare()
        self.assertEqual(code, 1)
        self.assertEqual(result['failed_check'], 'audit_nonfinite')

    def test_reference_pressure_update_must_match_supported_algebra(self):
        self.phi += .1
        self.write_reference_fields(); self.refresh()
        self.assertEqual(self.compare()[0], 1)

    def test_tampering_after_launch_is_rejected(self):
        path = self.mars / 'flow-audit-rank000000.bin'
        with path.open('ab') as stream: stream.write(b'changed')
        self.assertEqual(self.compare()[0], 1)

    def test_only_io_and_stopping_controls_change(self):
        actual = probe.load_deck((self.pair / 'reference/input.i').read_bytes())
        expected = copy.deepcopy(self.fixture.doc)
        expected['mesh']['file_path'] = str(self.fixture.mesh)
        solver = expected['simulation']['solver']
        solver['solver_control']['basic_settings']['convergence_controls'].update(min_iterations=1, max_iterations=1)
        solver['output_control'] = actual['simulation']['solver']['output_control']
        audit.enable_reference_audit(expected)
        self.assertEqual(actual, expected)
        args = probe.solver_arguments(self.pair, 'mars')
        self.assertEqual(probe.options(args)['--first-step-audit'], '1')

    def test_named_solver_lookup_preserves_settings(self):
        doc = copy.deepcopy(self.fixture.doc)
        settings = doc['simulation']['solver']['solver_control']['advanced_options']['linear_solver_settings']
        settings['coupled_navier_stokes'] = {'lookup': 'example'}
        doc['simulation']['solver']['example'] = {'family': 'PETSc', 'rtol': .001, 'diagonal_scaling': True}
        audit.enable_reference_audit(doc)
        self.assertEqual(doc['simulation']['solver']['example'],
                         {'family': 'PETSc', 'rtol': .001, 'diagonal_scaling': True, 'write_system': True})


if __name__ == '__main__':
    unittest.main()
