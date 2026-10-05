"""Only synthetic nodes, fields, decks and logs are used here."""

import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from netCDF4 import Dataset
import yaml

import prepare_simple_deck as bridge
import simple_snapshot_compare as compare
from test_simple_convergence_summary import fixture, serialize


ROOT = Path(__file__).resolve().parents[1]


class SnapshotTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.reference = self.root / 'private-reference'
        self.reference.mkdir()
        self.mesh = self.root / 'private-mesh.exo'
        self.prefix = self.root / 'private-flow'
        self.ids = np.array([9, 2**54 + 7, 2, 74, 16, 5], dtype=np.int64)
        self.xyz = np.array([[i * .125, i % 2, (i % 3) / 4] for i in range(6)])
        self.fields = np.array([[.1*(i + 1), .02*i, -.01*i, .03*i] for i in range(6)])
        self.write_exodus(self.mesh, list(range(6)), result=False)
        self.doc = bridge.load_deck((ROOT / 'tests/reference/openaccel/distributed_simple/preparation_fixture.i').read_bytes())
        self.doc['mesh']['file_path'] = str(self.mesh)
        settings = self.doc['simulation']['solver']['solver_control']['advanced_options']['linear_solver_settings']
        settings['pressure_correction'] = {'family': 'Hypre', 'rtol': 1e-4, 'atol': 1e-10}
        self.deck = self.reference / 'private-deck.yaml'
        self.deck.write_text(yaml.safe_dump(self.doc))
        translated, _ = bridge.translate(self.doc, 'reference')
        arguments = ['--mesh', str(self.mesh), '--mesh-format', 'exodus'] + translated + ['--reference-length', '1.0']
        self.case = self.root / 'case.json'
        self.case.write_text(json.dumps(dict(format='mars-simple-deck-v1', arguments=arguments,
                                             deck_sha256=compare.digest(self.deck), pressure_linear_policy='reference')))
        values = compare.options(arguments)
        rows, _ = fixture(count=100)
        peak = float(compare.norm3(self.fields[:, :3]).max())
        for row in rows:
            row['umax_m_s'] = peak
        log, metrics, code = serialize(rows)
        log = log.replace('high-resolution', 'upwind').replace('4 ranks', '2 ranks').replace('ranks=4', 'ranks=2')
        self.log = self.root / 'run.log'
        self.exit = self.root / 'run.exit'
        self.metrics = Path(str(self.prefix) + '-metrics.csv')
        lines = ['velocity_interpolation=' + values['--velocity-interpolation']]
        for keys in (list(compare.CONTROL_KEYS)[:5], list(compare.CONTROL_KEYS)[5:]):
            lines.append(' '.join(key + '=' + values[compare.CONTROL_KEYS[key]] for key in keys))
        lines.append('pressure_linear_rtol=0.0001 pressure_linear_atol=1e-10')
        self.log.write_text(log + '\n' + '\n'.join(lines) + '\n')
        self.metrics.write_text(metrics)
        self.exit.write_text(code)
        self.write_mars()
        self.ref_paths = [self.reference / ('results.e.2.' + str(i)) for i in range(2)]
        self.write_exodus(self.ref_paths[0], [4, 0, 2, 1])
        self.write_exodus(self.ref_paths[1], [3, 1, 5, 2], combined=True)
        self.public = self.root / 'public.json'
        self.private = self.root / 'private-report.json'

    def tearDown(self):
        self.tmp.cleanup()

    def write_exodus(self, path, source, result=True, combined=False, fields=None, ids=None, times=None):
        fields = self.fields if fields is None else fields
        times = [0., 50., 100.] if times is None else times
        with Dataset(str(path), 'w') as ds:
            ds.createDimension('num_dim', 3)
            ds.createDimension('num_nodes', len(source))
            ds.createVariable('coord', 'f8', ('num_dim', 'num_nodes'))[:] = self.xyz[source].T
            ds.createVariable('node_num_map', 'i8', ('num_nodes',))[:] = self.ids[source] if ids is None else ids
            if not result:
                return
            ds.createDimension('time_step', len(times))
            ds.createVariable('time_whole', 'f8', ('time_step',))[:] = times
            ds.createDimension('num_nod_var', 4)
            ds.createDimension('len_name', 32)
            names = np.full((4, 32), b'\0', dtype='S1')
            for i, name in enumerate(compare.FIELDS):
                names[i, :len(name)] = np.frombuffer(name.encode(), dtype='S1')
            ds.createVariable('name_nod_var', 'S1', ('num_nod_var', 'len_name'))[:] = names
            if combined:
                var = ds.createVariable('vals_nod_var', 'f8', ('time_step', 'num_nod_var', 'num_nodes'))
                for step in range(len(times)):
                    var[step] = fields[source].T if times[step] == 100 else 0.
            else:
                for i in range(4):
                    var = ds.createVariable('vals_nod_var' + str(i + 1), 'f8', ('time_step', 'num_nodes'))
                    for step in range(len(times)):
                        var[step] = fields[source, i] if times[step] == 100 else 0.

    def write_mars(self):
        for rank, nodes in enumerate(([5, 0, 2], [1, 4, 3])):
            path = Path(str(self.prefix) + '-fields-rank{:06d}.csv'.format(rank))
            with path.open('w', newline='') as stream:
                writer = csv.writer(stream)
                writer.writerow(['node', 'x', 'y', 'z', 'u', 'v', 'w', 'p'])
                for n in nodes:
                    writer.writerow([n] + list(self.xyz[n]) + list(self.fields[n]))
        Path(str(self.prefix) + '-fields.json').write_text(json.dumps(dict(
            format='mars-simple-fields-v1', nodes=6,
            parts=[self.prefix.name + '-fields-rank{:06d}.csv'.format(i) for i in range(2)])))

    def run_compare(self, expect=0, **changes):
        arguments = dict(reference_dir=self.reference, mesh=self.mesh, case=self.case,
                         mars_prefix=self.prefix, mars_log=self.log, mars_exit=self.exit,
                         iteration=100, output=self.public, private_report=self.private)
        arguments.update(changes)
        argv = [item for key, value in arguments.items() if value is not None
                for item in ('--' + key.replace('_', '-'), str(value))]
        with contextlib.redirect_stdout(io.StringIO()) as stdout, contextlib.redirect_stderr(io.StringIO()) as stderr:
            status = compare.main(argv)
        self.assertEqual(status, expect, stderr.getvalue())
        for text in (stdout.getvalue(), stderr.getvalue(), self.public.read_text()):
            self.assertNotIn(str(self.root), text)
            self.assertNotIn('private-deck', text)
            self.assertNotIn('987.654321', text)
        result = json.loads(self.public.read_text())
        self.assertTrue(all(isinstance(v, str) or type(v) is bool for v in result.values()))
        return result

    def test_distributed_permuted_ids_and_unconverged_matched_state(self):
        result = self.run_compare()
        self.assertEqual(result['comparison_status'], 'completed')
        self.assertEqual(result['mars_status'], 'iteration_limit')
        self.assertTrue(result['snapshot_fields_within_tolerance'])
        self.assertEqual(result['reference_settings_status'], 'mapped_controls_match')
        self.assertTrue(result['reference_deck_hash_matches_preparation'])
        self.assertFalse(result['full_run_provenance_verified'])
        self.assertFalse(result['identical_linear_solvers_verified'])
        detail = json.loads(self.private.read_text())
        self.assertEqual(detail['iteration'], 100)
        self.assertEqual(detail['errors']['velocity_max_scaled'], 0.)
        self.assertGreater(len(detail['sha256']), 8)

    def test_serial_reference_and_mars_output(self):
        for path in self.ref_paths:
            path.unlink()
        self.write_exodus(self.reference / 'results.e', [5, 3, 0, 2, 1, 4])
        parts = [Path(str(self.prefix) + '-fields-rank{:06d}.csv'.format(i)) for i in range(2)]
        Path(str(self.prefix) + '-fields.csv').write_text(parts[0].read_text() + '\n'.join(parts[1].read_text().splitlines()[1:]) + '\n')
        Path(str(self.prefix) + '-fields.json').unlink()
        self.assertTrue(self.run_compare()['snapshot_fields_within_tolerance'])

    def test_mesh_path_defaults_to_saved_preparation(self):
        self.assertTrue(self.run_compare(mesh=None)['snapshot_fields_within_tolerance'])

    def test_constant_pressure_offset_fails_absolute_parity(self):
        fields = self.fields.copy()
        fields[:, 3] += .02
        for path, nodes in zip(self.ref_paths, ([4, 0, 2, 1], [3, 1, 5, 2])):
            self.write_exodus(path, nodes, fields=fields)
        result = self.run_compare()
        self.assertFalse(result['snapshot_fields_within_tolerance'])
        self.assertEqual(result['pressure_max_scaled_band'], 'over_5_percent')
        self.assertEqual(result['velocity_max_scaled_band'], 'within_1e_minus_5')

    def test_matching_peak_does_not_hide_wrong_velocity_vector(self):
        fields = self.fields.copy()
        fields[0, 0] *= -1
        self.write_exodus(self.ref_paths[0], [4, 0, 2, 1], fields=fields)
        result = self.run_compare()
        self.assertEqual(result['peak_speed_relative_band'], 'within_1e_minus_5')
        self.assertFalse(result['snapshot_fields_within_tolerance'])

    def test_missing_piece_and_mixed_family_rejected(self):
        self.ref_paths[0].unlink()
        self.assertEqual(self.run_compare(1)['failed_check'], 'reference_file_family')

    def test_mixed_serial_and_partitioned_reference_rejected(self):
        self.write_exodus(self.reference / 'results.e', list(range(6)))
        self.assertEqual(self.run_compare(1)['failed_check'], 'reference_file_family')

    def test_missing_reference_node_rejected(self):
        self.write_exodus(self.ref_paths[1], [3, 1, 2])
        self.assertEqual(self.run_compare(1)['failed_check'], 'reference_node_coverage')

    def test_ghost_disagreement_rejected(self):
        fields = self.fields.copy()
        fields[1, 0] += .01
        self.write_exodus(self.ref_paths[1], [3, 1, 5, 2], fields=fields)
        self.assertEqual(self.run_compare(1)['failed_check'], 'reference_ghost_values')

    def test_wrong_reference_iteration_rejected(self):
        self.write_exodus(self.ref_paths[0], [4, 0, 2, 1], times=[0., 50., 99.])
        self.assertEqual(self.run_compare(1)['failed_check'], 'reference_saved_iteration')

    def test_wrong_requested_mars_iteration_rejected(self):
        self.assertEqual(self.run_compare(1, iteration=99)['failed_check'], 'mars_completion')

    def test_wrong_reference_global_ids_rejected(self):
        self.write_exodus(self.ref_paths[0], [4, 0, 2, 1], ids=[16, 9, 22, 2**54 + 7])
        self.assertEqual(self.run_compare(1)['failed_check'], 'reference_node_ids')

    def test_reference_coordinates_rejected(self):
        with Dataset(str(self.ref_paths[0]), 'a') as ds:
            ds['coord'][0, 0] += .1
        self.assertEqual(self.run_compare(1)['failed_check'], 'reference_coordinates')

    def test_nonfinite_field_rejected(self):
        with Dataset(str(self.ref_paths[0]), 'a') as ds:
            ds['vals_nod_var1'][-1, 0] = np.nan
        self.assertEqual(self.run_compare(1)['failed_check'], 'reference_fields')

    def test_field_name_failure_is_redacted(self):
        with Dataset(str(self.ref_paths[0]), 'a') as ds:
            ds['name_nod_var'][0, 0] = b'X'
        self.assertEqual(self.run_compare(1)['failed_check'], 'reference_field_names')

    def test_missing_mars_part_rejected(self):
        Path(str(self.prefix) + '-fields-rank000001.csv').unlink()
        self.assertEqual(self.run_compare(1)['failed_check'], 'mars_fields')

    def test_duplicate_mars_row_rejected(self):
        path = Path(str(self.prefix) + '-fields-rank000000.csv')
        path.write_text(path.read_text() + path.read_text().splitlines()[1] + '\n')
        self.assertEqual(self.run_compare(1)['failed_check'], 'mars_fields')

    def test_rank_count_must_match_manifest(self):
        self.log.write_text(self.log.read_text().replace('2 ranks', '4 ranks').replace('ranks=2', 'ranks=4'))
        self.assertEqual(self.run_compare(1)['failed_check'], 'mars_fields')

    def test_manifest_paths_cannot_escape_run_directory(self):
        path = Path(str(self.prefix) + '-fields.json')
        data = json.loads(path.read_text())
        data['parts'][0] = '../private.csv'
        path.write_text(json.dumps(data))
        self.assertEqual(self.run_compare(1)['failed_check'], 'mars_fields')

    def test_stale_mars_fields_peak_rejected(self):
        self.fields *= 2
        self.write_mars()
        self.assertEqual(self.run_compare(1)['failed_check'], 'mars_fields')

    def test_pressure_target_or_physical_controls_changed(self):
        self.log.write_text(self.log.read_text().replace('pressure_linear_rtol=0.0001', 'pressure_linear_rtol=0.01'))
        self.assertEqual(self.run_compare(1)['failed_check'], 'settings')

    def test_reformatted_deck_checks_translated_controls(self):
        self.deck.write_text(self.deck.read_text() + '# changed\n')
        result = self.run_compare()
        self.assertEqual(result['reference_settings_status'], 'mapped_controls_match')
        self.assertFalse(result['reference_deck_hash_matches_preparation'])

    def test_missing_deck_keeps_settings_unverified(self):
        self.deck.unlink()
        result = self.run_compare()
        self.assertEqual(result['reference_settings_status'], 'unavailable')
        self.assertFalse(result['snapshot_fields_within_tolerance'])
        self.assertEqual(result['velocity_max_scaled_band'], 'within_1e_minus_5')
        self.assertEqual(result['reference_deck_output_values'], 'unknown')
        self.assertFalse(result['solver_field_comparison_supported'])
        self.assertFalse(result['full_run_provenance_verified'])

    def update_deck(self):
        self.deck.write_text(yaml.safe_dump(self.doc))
        case = json.loads(self.case.read_text())
        case['deck_sha256'] = compare.digest(self.deck)
        self.case.write_text(json.dumps(case))

    def add_boundary_tags(self, elements=(1,), faces=(1,)):
        with Dataset(str(self.mesh), 'a') as ds:
            ds.createDimension('num_elem', 3)
            ds.createDimension('num_el_blk', 2)
            for number, rows in enumerate(([[1, 2, 3, 4]], [[3, 4, 5, 6], [1, 3, 5, 6]]), 1):
                suffix = str(number)
                ds.createDimension('num_el_in_blk' + suffix, len(rows))
                ds.createDimension('num_nod_per_el' + suffix, 4)
                var = ds.createVariable('connect' + suffix, 'i8', ('num_el_in_blk' + suffix, 'num_nod_per_el' + suffix))
                var.elem_type = 'TETRA4'
                var[:] = rows
            ds.createDimension('num_side_sets', 1)
            ds.createDimension('num_side_ss1', len(elements))
            ds.createVariable('elem_ss1', 'i8', ('num_side_ss1',))[:] = elements
            ds.createVariable('side_ss1', 'i8', ('num_side_ss1',))[:] = faces

    def test_corrected_reference_output_cannot_claim_solver_parity(self):
        self.doc['simulation']['solver']['output_control']['corrected_boundary_values'] = True
        self.update_deck()
        result = self.run_compare()
        self.assertEqual(result['reference_deck_output_values'], 'boundary_corrected')
        self.assertFalse(result['solver_field_comparison_supported'])
        self.assertFalse(result['snapshot_fields_within_tolerance'])
        self.assertEqual(result['velocity_max_scaled_band'], 'within_1e_minus_5')

    def test_omitted_output_correction_uses_pinned_false_default(self):
        del self.doc['simulation']['solver']['output_control']['corrected_boundary_values']
        self.update_deck()
        result = self.run_compare()
        self.assertEqual(result['reference_deck_output_values'], 'solver_values')
        self.assertTrue(result['solver_field_comparison_supported'])

    def test_malformed_output_correction_is_unknown(self):
        self.doc['simulation']['solver']['output_control']['corrected_boundary_values'] = 'private-invalid-value'
        self.update_deck()
        result = self.run_compare()
        self.assertEqual(result['reference_deck_output_values'], 'unknown')
        self.assertFalse(result['solver_field_comparison_supported'])
        self.assertNotIn('private-invalid-value', self.public.read_text())

    def test_unmatched_deck_does_not_establish_output_semantics(self):
        self.deck.write_text(self.deck.read_text() + '# unverified replacement\n')
        result = self.run_compare()
        self.assertFalse(result['solver_field_comparison_supported'])
        self.assertFalse(result['snapshot_fields_within_tolerance'])

    def test_boundary_only_difference_is_localized_without_hiding_global_error(self):
        self.add_boundary_tags()
        fields = self.fields.copy()
        fields[[0, 1, 3], :] = 0.
        for path, nodes in zip(self.ref_paths, ([4, 0, 2, 1], [3, 1, 5, 2])):
            self.write_exodus(path, nodes, fields=fields)
        result = self.run_compare()
        self.assertEqual(result['boundary_localization_status'], 'source_tags')
        self.assertEqual(result['tagged_boundary_velocity_max_scaled_band'], 'over_5_percent')
        self.assertEqual(result['other_nodes_velocity_max_scaled_band'], 'within_1e_minus_5')
        self.assertEqual(result['other_nodes_pressure_max_scaled_band'], 'within_1e_minus_5')
        self.assertFalse(result['snapshot_fields_within_tolerance'])

    def test_nonboundary_difference_is_not_attributed_to_boundary_output(self):
        self.add_boundary_tags()
        fields = self.fields.copy()
        fields[5, :] += .1
        self.write_exodus(self.ref_paths[1], [3, 1, 5, 2], fields=fields)
        result = self.run_compare()
        self.assertEqual(result['tagged_boundary_velocity_max_scaled_band'], 'within_1e_minus_5')
        self.assertEqual(result['other_nodes_velocity_max_scaled_band'], 'over_5_percent')
        self.assertFalse(result['snapshot_fields_within_tolerance'])

    def test_side_sets_use_element_rows_across_blocks_and_all_tet_ordinals(self):
        expected = [[0, 1, 3], [1, 2, 3], [0, 2, 3], [0, 1, 2]]
        self.add_boundary_tags(elements=(3,), faces=(1,))
        connectivity = np.array([0, 2, 4, 5])
        for face, nodes in enumerate(expected, 1):
            with Dataset(str(self.mesh), 'a') as ds:
                ds['side_ss1'][:] = [face]
                mask = compare.boundary_tags(ds, 6)
            np.testing.assert_array_equal(np.flatnonzero(mask), connectivity[nodes])

    def test_node_sets_join_side_sets(self):
        self.add_boundary_tags()
        with Dataset(str(self.mesh), 'a') as ds:
            ds.createDimension('num_node_sets', 1)
            ds.createDimension('num_nod_ns1', 1)
            ds.createVariable('node_ns1', 'i8', ('num_nod_ns1',))[:] = [6]
            np.testing.assert_array_equal(np.flatnonzero(compare.boundary_tags(ds, 6)), [0, 1, 3, 5])

    def test_node_set_only_and_empty_complement(self):
        with Dataset(str(self.mesh), 'a') as ds:
            ds.createDimension('num_node_sets', 1)
            ds.createDimension('num_nod_ns1', 6)
            ds.createVariable('node_ns1', 'i8', ('num_nod_ns1',))[:] = [1, 2, 3, 4, 5, 6]
        result = self.run_compare()
        self.assertEqual(result['other_nodes_status'], 'empty')
        self.assertNotIn('other_nodes_velocity_max_scaled_band', result)

    def test_repeated_unordered_sides_are_a_union(self):
        self.add_boundary_tags(elements=(3, 1, 3), faces=(2, 1, 2))
        with Dataset(str(self.mesh)) as ds:
            np.testing.assert_array_equal(np.flatnonzero(compare.boundary_tags(ds, 6)), np.arange(6))

    def test_bad_connectivity_is_rejected(self):
        self.add_boundary_tags()
        with Dataset(str(self.mesh), 'a') as ds:
            ds['connect1'][0, 3] = 0
        self.assertEqual(self.run_compare(1)['failed_check'], 'source_boundary_tags')

    def test_bad_side_element_is_rejected_without_private_values(self):
        self.add_boundary_tags(elements=(987654321,))
        self.assertEqual(self.run_compare(1)['failed_check'], 'source_boundary_tags')

    def test_bad_side_ordinal_is_rejected(self):
        self.add_boundary_tags(faces=(5,))
        self.assertEqual(self.run_compare(1)['failed_check'], 'source_boundary_tags')

    def test_unsupported_topology_keeps_regional_evidence_unavailable(self):
        self.add_boundary_tags()
        with Dataset(str(self.mesh), 'a') as ds:
            ds['connect1'].elem_type = 'PRIVATE_UNSUPPORTED_TOPOLOGY'
        result = self.run_compare()
        self.assertEqual(result['boundary_localization_status'], 'unavailable')
        self.assertNotIn('PRIVATE_UNSUPPORTED_TOPOLOGY', self.public.read_text())
        self.assertNotIn('other_nodes_velocity_max_scaled_band', result)

    def test_changed_reference_physics_are_reported(self):
        self.doc['simulation']['material_library'][0]['transport_properties']['dynamic_viscosity']['dynamic_viscosity'] = .2
        self.deck.write_text(yaml.safe_dump(self.doc))
        self.assertEqual(self.run_compare()['reference_settings_status'], 'mapped_controls_differ')

    def test_failed_run_or_mixed_completion_rejected(self):
        self.exit.write_text('137')
        self.assertEqual(self.run_compare(1)['failed_check'], 'mars_completion')

    def test_multiple_completions_rejected(self):
        self.log.write_text(self.log.read_text() + 'NOT CONVERGED: iteration limit iterations=100 ranks=2 exchange_rounds=401\n')
        self.assertEqual(self.run_compare(1)['failed_check'], 'mars_completion')

    def test_existing_public_output_never_overwritten(self):
        self.public.write_text('{"private_existing": true}')
        self.run_compare(1)
        self.assertEqual(self.public.read_text(), '{"private_existing": true}')
        self.assertFalse(self.private.exists())

    def test_existing_private_output_never_overwritten(self):
        self.private.write_text('PRIVATE EXISTING')
        self.run_compare(1)
        self.assertEqual(self.private.read_text(), 'PRIVATE EXISTING')

    def test_metrics_must_match_log(self):
        self.metrics.write_text(self.metrics.read_text().replace('0.01', '0.02'))
        self.assertEqual(self.run_compare(1)['failed_check'], 'mars_completion')

    def test_error_bands_are_fixed_and_inclusive(self):
        self.assertEqual([compare.band(x) for x in [0, 1e-5, .001, .01, .03, .05, 1]],
                         ['within_1e_minus_5']*2 + ['within_1_percent']*2 + ['within_5_percent']*2 + ['over_5_percent'])


if __name__ == '__main__':
    unittest.main()
