"""Public tetrahedra only; validate saved-field replay without private case data."""
import contextlib
import io
import json
from fractions import Fraction
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
from netCDF4 import Dataset

import simple_first_step_audit as audit
import simple_pressure_probe as pressure
from simple_snapshot_compare import EvidenceError
import test_simple_pressure_probe as probe_tests
import simple_startup_probe as startup
from simple_snapshot_compare import reference_fields


class GradientTests(unittest.TestCase):
    xyz = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    cells = np.array([[0, 1, 2, 3]])

    def apply(self, xyz, fields, cells=None):
        return audit.tet_gradient_action(xyz, [self.cells if cells is None else cells], fields)

    def test_single_tet_shifted_boundary_is_not_affine_exact(self):
        result, _ = self.apply(self.xyz, self.xyz[:, :1])
        np.testing.assert_allclose(result[:, 0], [[1., .5, .5], [2., 0., 0.], [.5, -.5, 0.], [.5, 0., -.5]])
        np.testing.assert_allclose(result[:, 0].mean(axis=0), [1., 0., 0.])

    def test_constant_translation_orientation_and_length_scale(self):
        field = np.array([[2., 7.], [3., 7.], [5., 7.], [11., 7.]])
        base, _ = self.apply(self.xyz, field)
        np.testing.assert_array_equal(base[:, 1], 0.)
        moved, _ = self.apply(8*self.xyz+1024., field, np.array([[0, 2, 1, 3]]))
        np.testing.assert_allclose(moved, base/8, rtol=1e-15, atol=0.)

    def test_pressure_difference_amplification_and_false_explanations(self):
        r = np.zeros((4, 1))
        m = self.xyz[:, :1]*1e-9
        result, size = self.apply(self.xyz*1e-6, np.column_stack([m, r, m-r]))
        gm, gr = result[:, 0].copy(), result[:, 1].copy()
        checks, _ = audit.gradient_checks(result, size, gm, gr, 1.)
        self.assertLess(np.max(abs(m-r)), 1e-5)
        self.assertEqual(checks['assessment'], 'input_difference_explains_gradient_within_replay_tolerance')
        # A common wrong gradient cancels in delta g but must fail the individual checks.
        checks, _ = audit.gradient_checks(result, size, gm+1., gr+1., 1.)
        self.assertTrue(checks['gradient_difference_matches_pressure_difference_action'])
        self.assertEqual(checks['assessment'], 'reconstruction_mismatch')
        checks, _ = audit.gradient_checks(result, size, gm*2., gr, 1.)
        self.assertFalse(checks['gradient_difference_matches_pressure_difference_action'])

    def test_loose_replay_resolution_cannot_explain_difference(self):
        values = np.column_stack([self.xyz[:, 0], self.xyz[:, 0], np.zeros(4)])
        result, size = self.apply(self.xyz, values)
        gm, gr = result[:, 0].copy(), result[:, 1].copy()
        gm[0, 0] += 2e-5
        checks, _ = audit.gradient_checks(result, size+1e8, gm, gr, 1.)
        self.assertEqual(checks['assessment'], 'insufficient_replay_resolution')
        checks, _ = audit.gradient_checks(result, size, gr, gr, 1.)
        self.assertEqual(checks['assessment'], 'gradients_within_field_tolerance')

    def test_invalid_geometry_and_missing_star_nodes_rejected(self):
        for xyz, cells in ((self.xyz, np.array([[0, 1, 2, 2]])),
                           (self.xyz, np.array([[0, 1, 2, 5]])),
                           (np.vstack([self.xyz, [2., 2., 2.]]), self.cells)):
            with self.subTest(cells=cells.tolist()), self.assertRaises(EvidenceError):
                self.apply(xyz, np.zeros((len(xyz), 1)), cells)

    def test_exodus_blocks_and_invalid_connectivity(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'public.exo'
            with Dataset(str(path), 'w') as ds:
                ds.createDimension('num_el_blk', 2); ds.createDimension('num_elem', 2)
                ds.createDimension('one', 1); ds.createDimension('four', 4)
                for number in (1, 2):
                    v = ds.createVariable('connect'+str(number), 'i4', ('one', 'four'))
                    v.elem_type = 'TET4'; v[:] = self.cells+1
            with Dataset(str(path)) as ds:
                result, _ = audit.tet_gradient_action(self.xyz, audit.tet_chunks(ds, 4), self.xyz[:, :1])
            np.testing.assert_allclose(result, self.apply(self.xyz, self.xyz[:, :1])[0])
            with Dataset(str(path), 'a') as ds: ds['connect2'][0, 0] = 0
            with Dataset(str(path)) as ds, self.assertRaises(EvidenceError):
                list(audit.tet_chunks(ds, 4))

    def test_exact_action_matches_analytic_tet_and_permutation(self):
        phi = self.xyz[:, :1]
        selected = np.arange(4)
        exact = audit.exact_gradient_action(self.xyz, [self.cells], phi, phi*0, selected)
        expected = [[1., .5, .5], [2., 0., 0.], [.5, -.5, 0.], [.5, 0., -.5]]
        for node in selected:
            self.assertEqual(exact[node][0], [Fraction(x) for x in expected[node]])
            self.assertEqual(exact[node][1], [Fraction(0)]*3)
        moved = audit.exact_gradient_action(self.xyz*8+1024, [self.cells[:, [0, 2, 1, 3]]], phi, phi*0, selected)
        for node in selected: self.assertEqual(moved[node][0], [x/8 for x in exact[node][0]])

    def test_exact_check_resolves_large_offset_without_relaxing_old_verdict(self):
        xyz = self.xyz*1e-6
        r = np.full((4, 1), 1e6)
        m = r + self.xyz[:, :1]*1e-9
        selected = np.arange(4)
        exact = audit.exact_gradient_action(xyz, [self.cells], m, r, selected)
        gm = np.array([[float(x) for x in exact[node][0]] for node in selected])
        gr = np.zeros_like(gm)
        replayed, size = self.apply(xyz, np.column_stack([m, r, m-r]))
        old, _ = audit.gradient_checks(replayed, size, gm, gr, 1.)
        self.assertEqual(old['assessment'], 'insufficient_replay_resolution')
        checks, private = audit.exact_gradient_checks(exact, selected, 4, gm, gr)
        self.assertEqual(checks['assessment'], 'input_difference_explains_checked_nodes')
        self.assertTrue(checks['all_field_tolerance_failures_checked'])
        self.assertFalse(checks['global_equivalence_verified'])
        self.assertFalse(checks['runtime_roundoff_bound_proven'])
        self.assertEqual(len(private['selected_nodes']), 4)
        # These equal errors cancel in the difference but contradict the individual reconstructions.
        common, _ = audit.exact_gradient_checks(exact, selected, 4, gm+1, gr+1)
        self.assertEqual(common['assessment'], 'individual_reconstruction_discrepancy_at_checked_nodes')
        wrong, _ = audit.exact_gradient_checks(exact, selected, 4, gm*2, gr)
        self.assertEqual(wrong['assessment'], 'captured_gradient_difference_not_explained_at_checked_nodes')

    def test_exact_conversion_precedes_subtraction(self):
        # 1 - 2^-54 rounds to 1 in binary64; preserve that small term in the rational action.
        phi = np.array([[2.**-54], [1.], [0.], [0.]])
        exact = audit.exact_gradient_action(self.xyz, [self.cells], phi, phi*0, np.array([0]))
        self.assertEqual(exact[0][0][0], Fraction(1)-Fraction(1, 2**53))
        self.assertNotEqual(exact[0][0][0], Fraction(float(1.-phi[0, 0]))-Fraction(1, 2**54))

    def test_exact_complete_star_across_chunks_and_work_limit(self):
        xyz = np.vstack([self.xyz, [1., 1., 1.]])
        cells = [self.cells, np.array([[1, 2, 3, 4]])]
        phi = xyz[:, :1]
        exact = audit.exact_gradient_action(xyz, cells, phi, phi*0, np.array([1]))
        replayed, _ = audit.tet_gradient_action(xyz, cells, phi)
        np.testing.assert_allclose([float(x) for x in exact[1][0]], replayed[1, 0])
        incomplete = audit.exact_gradient_action(xyz, cells[:1], phi, phi*0, np.array([1]))
        self.assertNotEqual(exact, incomplete)
        limited = audit.exact_gradient_action(xyz, cells, phi, phi*0, np.array([1]), max_elements=1)
        self.assertIsNone(limited)
        checks, _ = audit.exact_gradient_checks(limited, [1], 1, np.ones((5, 3)), np.zeros((5, 3)))
        self.assertEqual(checks['assessment'], 'inconclusive_work_limit')
        self.assertFalse(checks['all_field_tolerance_failures_checked'])
        with self.assertRaises(EvidenceError):
            audit.exact_gradient_checks({}, [], 0, np.zeros((5, 3)), np.zeros((5, 3)))

    def test_selection_includes_error_peaks_and_breaks_ties_by_source(self):
        gm, gr = np.ones((40, 3)), np.zeros((40, 3))
        replayed = np.zeros((40, 3, 3))
        replayed[:, 0] = gm
        replayed[:, 2] = gm-gr
        replayed[35, 0, 0] += 3
        replayed[36, 1, 0] += 4
        replayed[37, 2, 0] += 5
        selected, failing = audit.select_gradient_nodes(replayed, gm, gr, 1.)
        self.assertEqual(failing, 40)
        self.assertLessEqual(len(selected), 32)
        self.assertTrue(set(range(8)).issubset(selected))
        self.assertTrue({35, 36, 37}.issubset(selected))
        same, _ = audit.select_gradient_nodes(replayed, gm, gr, 1.)
        np.testing.assert_array_equal(selected, same)

    def test_production_geometry_kernel_parity(self):
        compiler = shutil.which('c++')
        if compiler is None: self.skipTest('C++ compiler unavailable')
        root = Path(__file__).resolve().parents[1]
        source = r'''
#include "backend/distributed/unstructured/fem/segregated/mars_segregated_geometry.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
int main() {
    int n,e; std::cin>>n>>e;
    std::vector<double> xyz(3*n), p(n), sums(3*n), vol(n);
    for (auto& x:xyz) std::cin>>x;
    for (auto& x:p) std::cin>>x;
    for (int t=0;t<e;++t) {
        int ids[4]; for (auto& i:ids) std::cin>>i;
        double x[12], v[4], out[12];
        for (int i=0;i<4;++i) { v[i]=p[ids[i]]; for(int j=0;j<3;++j) x[3*i+j]=xyz[3*ids[i]+j]; }
        mars::segregated::TetGeometry<double> g;
        if (!mars::segregated::tet_geometry(x,g)) return 2;
        mars::segregated::tet_gradient_numerator<1>(g,v,true,true,out);
        for(int i=0;i<4;++i) { vol[ids[i]]+=g.volume/4; for(int j=0;j<3;++j) sums[3*ids[i]+j]+=out[3*i+j]; }
        for(int f=0;f<4;++f) {
            double a[3], fv[3], b[9]; mars::segregated::tet_boundary_area(g,f,a);
            for(int i=0;i<3;++i) fv[i]=v[mars::segregated::tet_face_node(f,i)];
            mars::segregated::tri_gradient_numerator<1>(a,fv,true,true,b);
            for(double z:b) if(z!=0.) return 3;
        }
    }
    std::cout<<std::setprecision(17);
    for(int i=0;i<3*n;++i) std::cout<<sums[i]/vol[i/3]<<'\n';
}
'''
        with tempfile.TemporaryDirectory() as tmp:
            cpp, exe = Path(tmp)/'gate.cpp', Path(tmp)/'gate'
            cpp.write_text(source)
            subprocess.run([compiler, '-std=c++17', '-O2', '-I'+str(root), str(cpp), '-o', str(exe)], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            xyz = np.vstack([self.xyz, [1., 1., 1.]])
            cells = np.array([[0, 1, 2, 3], [1, 2, 3, 4]])
            rng = np.random.RandomState(73)
            for transform in (np.eye(3), np.array([[2., .3, 0.], [0., .1, .2], [.1, 0., 3.]])):
                coords, field = xyz @ transform + 100., rng.normal(size=(5, 1))
                data = '5 2\n'+' '.join(map(str, coords.ravel()))+'\n'+' '.join(map(str, field.ravel()))+'\n'+' '.join(map(str, cells.ravel()))
                output = subprocess.run([str(exe)], input=data, universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
                expected = np.fromstring(output.stdout, sep=' ').reshape(5, 3)
                actual, _ = self.apply(coords, field, cells)
                np.testing.assert_allclose(actual[:, 0], expected, rtol=2e-14, atol=2e-14)


class GradientEvidenceTests(unittest.TestCase):
    def geometry_capture(self):
        case = probe_tests.PressureProbeTests()
        case.setUp()
        self.addCleanup(case.doCleanups)
        capture = case.capture
        capture.fixture.add_boundary_tags()
        record_path = case.baseline / 'pair.json'
        record = startup.read_json(record_path)
        record['mesh_sha256'] = startup.digest(capture.fixture.mesh)
        record_path.write_text(json.dumps(record))
        cells = np.array([[0, 1, 2, 3], [2, 3, 4, 5], [0, 2, 4, 5]])
        result, _ = audit.tet_gradient_action(capture.fixture.xyz, [cells], capture.phi)
        capture.gradient = result[:, 0]
        capture.final = np.column_stack([capture.predictor-capture.influence*capture.gradient, .3*capture.phi])
        capture.write_reference_fields()
        capture.change_mars_final()
        case.finish()
        return case, cells

    def test_complete_capture_replay_and_unchanged_originals(self):
        case, _ = self.geometry_capture()
        before = case.hashes(case.pair)
        status = case.invoke('compare', '--pair', str(case.pair), '--output', str(case.output),
                             '--detail-dir', str(case.details), '--gradient-audit')
        result = json.loads(case.output.read_text())
        self.assertEqual(status, 0, result)
        checks = result['tightened_first_step']['pressure_gradient_checks']
        self.assertEqual(checks['assessment'], 'gradients_within_field_tolerance')
        self.assertEqual(before, case.hashes(case.pair))
        self.assertTrue((case.details/'tightened/gradient-private.json').is_file())
        self.assertNotIn('max_errors_scaled', case.output.read_text())
        self.assertNotIn(str(case.root), case.output.read_text())

    def test_exact_selected_check_preserves_failed_solve_and_private_evidence(self):
        case, cells = self.geometry_capture()
        capture = case.capture
        capture.pair = case.pair
        capture.phi[0, 0] += 1e-6
        replayed, _ = audit.tet_gradient_action(capture.fixture.xyz, [cells], capture.phi)
        capture.gradient = replayed[:, 0]
        capture.change_mars_final()
        before, baseline = case.hashes(case.pair), case.hashes(case.baseline)
        status = case.invoke('compare', '--pair', str(case.pair), '--output', str(case.output),
                             '--detail-dir', str(case.details), '--gradient-audit')
        public = json.loads(case.output.read_text())
        self.assertEqual(status, 0, public)
        self.assertEqual(public['outcome'], 'pressure_target_not_met')
        self.assertFalse(public['tighter_pressure_residuals_pass'])
        stages = public['tightened_first_step']['first_step_stage_matches']
        self.assertFalse(stages['pressure_increment_gradient'])
        checks = public['tightened_first_step']['pressure_gradient_checks']['exact_selected_node_checks']
        self.assertEqual(checks['assessment'], 'input_difference_explains_checked_nodes')
        self.assertFalse(checks['global_equivalence_verified'])
        self.assertEqual(before, case.hashes(case.pair))
        self.assertEqual(baseline, case.hashes(case.baseline))
        private = json.loads((case.details/'tightened/gradient-private.json').read_text())
        self.assertTrue(private['exact_selected_node_details']['selected_nodes'])
        for forbidden in ('source_node', 'error_fractions', str(case.root)):
            self.assertNotIn(forbidden, case.output.read_text())

    def test_reference_precision_and_duplicate_copies(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = [Path(tmp)/('results.e.2.'+str(i)) for i in (0, 1)]
            xyz = GradientTests.xyz
            for dtype, discrepancy in (('f4', 0.), ('f8', 1e-13), ('f8', 0.)):
                for rank, path in enumerate(paths):
                    with Dataset(str(path), 'w') as ds:
                        for name, count in (('num_nodes', 4), ('num_dim', 3), ('time_step', 1),
                                            ('num_nod_var', 1), ('len_name', 3)):
                            ds.createDimension(name, count)
                        ds.createVariable('node_num_map', 'i4', ('num_nodes',))[:] = np.arange(1, 5)
                        ds.createVariable('coord', 'f8', ('num_dim', 'num_nodes'))[:] = xyz.T
                        ds.createVariable('time_whole', 'f8', ('time_step',))[:] = [1.]
                        ds.createVariable('name_nod_var', 'S1', ('num_nod_var', 'len_name'))[:] = np.array([list('phi')], dtype='S1')
                        ds.createVariable('vals_nod_var', dtype, ('time_step', 'num_nod_var', 'num_nodes'))[:] = 1.+rank*discrepancy
                args = (paths, np.arange(1, 5), xyz, 1, 1e-14, np.ones(1), ('phi',))
                reference_fields(*args)
                if dtype == 'f4' or discrepancy:
                    with self.assertRaises(EvidenceError): reference_fields(*args, exact_storage=True)
                else:
                    np.testing.assert_array_equal(reference_fields(*args, exact_storage=True), 1.)

    def test_missing_mesh_fails_closed_without_private_text(self):
        case = probe_tests.PressureProbeTests()
        case.setUp()
        try:
            case.finish()
            output, detail = case.root/'gradient-public.json', case.root/'gradient-private'
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                status = pressure.main(['compare', '--pair', str(case.pair), '--output', str(output),
                                        '--detail-dir', str(detail), '--gradient-audit'])
            result = json.loads(output.read_text())
            self.assertEqual(status, 1)
            self.assertNotEqual(result['comparison_status'], 'completed')
            self.assertNotIn(str(case.root), output.read_text())
        finally:
            case.doCleanups()


if __name__ == '__main__':
    unittest.main()
