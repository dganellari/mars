"""Synthetic first-step fields and saved runtime records; no private data."""
import copy
import json
import os
from pathlib import Path
import shutil
import unittest
from unittest.mock import patch

import yaml
import simple_pressure_profile as profile
import simple_hypre_defaults as defaults
import simple_startup_probe as startup
import test_simple_first_step_audit as fixtures


def measured_defaults():
    values = dict.fromkeys(defaults.INT_KEYS | defaults.REAL_KEYS, 0)
    values.update(amg_coarsen_type=10, amg_interp_type=6, amg_relax_down=13, amg_relax_up=14,
        amg_relax_coarse=9, amg_sweeps_down=1, amg_sweeps_up=1, amg_sweeps_coarse=1,
        amg_p_max_elmts=4, amg_max_levels=25, amg_max_coarse_size=9, amg_num_functions=1,
        amg_cycle_type=1, amg_strong_threshold=.25, amg_jacobi_trunc_threshold=.01,
        amg_max_row_sum=.9, amg_agg_interp_type=4, amg_num_paths=1,
        gmres_restart_dimension=5, flexgmres_restart_dimension=20)
    return values


class ProfileTests(unittest.TestCase):
    def setUp(self):
        self.capture = fixtures.FirstStepTests()
        self.capture.setUp()
        self.addCleanup(self.capture.doCleanups)
        self.baseline = self.capture.pair.resolve()
        self.root = self.capture.fixture.root.resolve()
        self.capture.pair = self.baseline
        self.capture.fixture.root = self.root
        self.pair = self.root / 'profile'
        self.probe = self.root / 'defaults'
        self.probe.mkdir()
        self.exe = self.root / 'rebuilt-mars'
        self.exe.write_text('NEW SYNTHETIC EXECUTABLE')
        self.exe.chmod(0o700)
        deck = startup.load_deck((self.baseline / 'reference/input.i').read_bytes())
        config = deck['simulation']['solver']['solver_control']['advanced_options']['linear_solver_settings']['pressure_correction']
        config.update(options=dict(solver=dict(type='gmres'), precond=dict(type='boomeramg')))
        self.config = config
        (self.baseline / 'reference/input.i').write_text(yaml.safe_dump(deck))
        case = startup.read_json(self.baseline / 'case.json')
        case['deck_sha256'] = startup.digest(self.baseline / 'reference/input.i')
        (self.baseline / 'case.json').write_text(json.dumps(case))
        record = startup.read_json(self.baseline / 'pair.json')
        record.update(deck_sha256=case['deck_sha256'], case_sha256=startup.digest(self.baseline / 'case.json'))
        (self.baseline / 'pair.json').write_text(json.dumps(record))
        self.capture.refresh()
        self.libs = {'/SYNTHETIC/lib/libHYPRE.so': 'a'*64, '/SYNTHETIC/lib/libmpi.so': 'b'*64}
        for directory in ('mars', 'reference'):
            for name in ('launch-start.json', 'launch.json'):
                path = self.baseline / directory / name
                record = startup.read_json(path)
                record['libraries'] = self.libs
                path.write_text(json.dumps(record))
        raw = dict(schema=defaults.RAW_SCHEMA, version='3.1.0', build_cuda=False, memory_location=1,
            execution_policy=0, header_version_matches=True, getter_layout_checks_passed=True,
            library_path='/SYNTHETIC/lib/libHYPRE.so', mpi_library_path='/SYNTHETIC/lib/libmpi.so', defaults=measured_defaults())
        (self.probe / 'probe.log').write_text(json.dumps(raw, separators=(',', ':'))+'\n')
        public = dict(defaults.public_defaults(raw), schema=defaults.SCHEMA, comparison_status='completed', failed_check='none',
            saved_capture_identity_verified=True, captured_libraries_unchanged=True, loaded_hypre_library_matches_capture=True,
            mpi_library_identity_verified=True, mpi_enabled=True)
        self.write('public.json', public)
        self.write('private-provenance.json', dict(capture=startup.verified_launch(self.baseline, 'openaccel')))
        self.write('private-probe-libraries.json', self.libs)
        for name in ('compile.exit', 'probe.exit'): (self.probe / name).write_text('0\n')

    def write(self, name, value):
        (self.probe / name).write_text(json.dumps(value))

    def prepare(self):
        with patch.object(startup, 'runtime_libraries', return_value=self.libs):
            profile.prepare(self.baseline, self.probe, self.exe, self.pair)

    def finish(self):
        self.prepare()
        shutil.copytree(self.baseline / 'mars', self.pair / 'mars')
        for name in ('launch.json', 'launch-start.json'): (self.pair / 'mars' / name).unlink()
        log = self.pair / 'mars/run.log'
        log.write_text(log.read_text()+'\npressure_solver_profile=explicit_gpu\n')
        _, _, values, _ = profile.check_inputs(self.pair)
        actual = dict(values, hypre_release=23300, effective_levels=3, effective_relax_3=18)
        actual['effective_relax_1'] = actual.pop('relax_down')
        actual['effective_relax_2'] = actual.pop('relax_up')
        for rank in (0, 1):
            (self.pair / 'mars' / ('flow-pressure-rank-{}.settings'.format(rank))).write_text(profile.text_profile(actual))
        self.record()

    def record(self):
        old = self.capture.pair
        self.capture.pair = self.pair
        self.capture.record('mars')
        self.capture.pair = old
        for name in ('launch-start.json', 'launch.json'):
            path = self.pair / 'mars' / name
            record = startup.read_json(path)
            record.update(executable=str(self.exe), executable_sha256=startup.digest(self.exe), libraries=self.libs,
                          command=['launcher', str(self.exe)]+startup.solver_arguments(self.pair, 'mars'))
            path.write_text(json.dumps(record))

    def test_defaults_and_two_gpu_adaptations(self):
        v, changes = profile.resolve(self.config, (1e-10, 1e-4), measured_defaults())
        self.assertEqual((v['kdim'], v['miniter'], v['maxiter']), (5, 0, 20))
        self.assertEqual((v['rtol'], v['atol']), (1e-4, 1e-10))
        self.assertEqual((v['relaxtype'], v['relax_down'], v['relax_up'], v['coarserelax']), (6, 13, 14, 18))
        self.assertEqual(changes, ['hmis_to_device_pmis', 'coarse_direct_to_device_l1_jacobi'])

    def test_explicit_options_and_setter_order(self):
        config = copy.deepcopy(self.config)
        config['options']['solver'] = dict(type='FlexGMRES', KDim=31)
        config['options']['precond'].update(RelaxType=6, NumSweeps=2, CoarsenType=8)
        config.update(max_iterations=123, min_iterations=7)
        v, changes = profile.resolve(config, (1e-10, 1e-4), measured_defaults())
        self.assertEqual((v['method'], v['kdim'], v['miniter'], v['maxiter']), (1, 31, 0, 123))
        self.assertEqual((v['relax_down'], v['relax_up'], v['sweeps_1'], v['sweeps_3']), (6, 6, 2, 1))
        self.assertEqual(changes, ['coarse_direct_to_device_l1_jacobi'])

    def test_unsupported_or_malformed_controls_fail(self):
        good, _ = profile.resolve(self.config, (1e-10, 1e-4), measured_defaults())
        for key, value in [('coarsentype',10), ('coarserelax',9), ('miniter',21), ('method',3),
                           ('relax_down',8), ('rtol',float('nan')), ('kdim',True), ('maxlevels',0)]:
            with self.subTest(key=key), self.assertRaises(startup.EvidenceError):
                profile.validate(dict(good, **{key:value}))
        for fault in ({'extra':0}, {'rtol':0}):
            with self.assertRaises(startup.EvidenceError): profile.validate(dict(good, **fault))

    def test_prepare_preserves_baseline_and_original_controls(self):
        before = {str(p):startup.digest(p) for p in self.baseline.rglob('*') if p.is_file()}
        self.prepare()
        _, _, v, _ = profile.check_inputs(self.pair)
        self.assertEqual((v['rtol'],v['atol']), (1e-4,1e-10))
        self.assertEqual(before, {str(p):startup.digest(p) for p in self.baseline.rglob('*') if p.is_file()})
        self.assertFalse((self.pair / 'reference/launch.json').exists())

    def test_bad_defaults_identity_and_raw_output_fail(self):
        for name, modify in [('private-provenance.json', lambda v: v['capture'].update(executable_sha256='c'*64)),
                             ('public.json', lambda v: v.update(version='2.33.0')),
                             ('private-probe-libraries.json', lambda v: v.update({'/SYNTHETIC/lib/libHYPRE.so':'c'*64}))]:
            path = self.probe / name
            saved = path.read_text(); value=json.loads(saved); modify(value); path.write_text(json.dumps(value))
            with self.assertRaises(startup.EvidenceError): self.prepare()
            path.write_text(saved)

    def test_saved_verified_loader_alias_does_not_require_reference_mount(self):
        path=self.probe/'probe.log'
        raw=json.loads(path.read_text())
        raw['library_path']='/SYNTHETIC/real-install/lib/libHYPRE.so'
        path.write_text(json.dumps(raw,separators=(',',':'))+'\n')
        self.prepare()
        profile.check_inputs(self.pair)

    def test_tightened_baseline_rejected(self):
        path = self.baseline / 'pair.json'; v=startup.read_json(path); v['pressure_experiment']={}; path.write_text(json.dumps(v))
        self.capture.refresh()
        with self.assertRaisesRegex(startup.EvidenceError,'original_baseline_required'): self.prepare()

    def test_complete_reuses_reference_without_claiming_fresh_pair(self):
        self.finish()
        details=self.root/'details'; details.mkdir()
        public={}; profile.compare(self.pair, public, details)
        self.assertEqual(public['comparison_status'],'completed')
        self.assertTrue(public['comparison']['reference_reused'])
        self.assertFalse(public['comparison']['fresh_launch_records_verified'])
        self.assertTrue(all(public['comparison']['first_step_stage_matches'].values()))
        self.assertFalse(public['identical_linear_solvers_verified'])
        self.assertNotIn(str(self.root),json.dumps(public))

    def test_missing_or_ignored_runtime_setting_rejected(self):
        self.finish()
        path=self.pair/'mars/flow-pressure-rank-1.settings'
        path.write_text(path.read_text().replace('kdim 5\n','kdim 100\n'))
        self.record()
        details=self.root/'details'; details.mkdir()
        with self.assertRaisesRegex(startup.EvidenceError,'runtime_profile_mismatch'): profile.compare(self.pair,{},details)
        path.unlink(); self.record()
        with self.assertRaisesRegex(startup.EvidenceError,'runtime_profile_mismatch'): profile.compare(self.pair,{},details)

    def test_saved_reference_tampering_rejected(self):
        self.finish()
        (self.baseline/'reference/run.log').write_text('PRIVATE CHANGED')
        with self.assertRaises(ValueError): profile.check_inputs(self.pair)

    def test_changed_profile_and_arguments_rejected(self):
        self.prepare()
        path=self.pair/'pressure.profile'; path.write_text(path.read_text().replace('rtol 0.0001','rtol 1e-10'))
        with self.assertRaises(ValueError): profile.check_inputs(self.pair)

    def test_new_executable_allowed_but_library_change_rejected(self):
        with patch.object(startup,'runtime_libraries',return_value={'wrong':'c'*64}):
            with self.assertRaisesRegex(startup.EvidenceError,'libraries_changed'):
                profile.prepare(self.baseline,self.probe,self.exe,self.pair)

    def test_run_checks_and_restores_only_recorded_numeric_environment(self):
        self.prepare()
        with patch.dict(os.environ,{'MARS_HYPRE_FLEXGMRES':'1'},clear=True), \
             patch.object(startup,'runtime_libraries',return_value=self.libs), \
             patch.object(startup,'launch') as launch, patch.object(profile,'check_launch'):
            profile.run(self.pair,{})
            self.assertNotIn('MARS_HYPRE_FLEXGMRES',launch.call_args[0][-1])
            self.assertEqual(launch.call_args[0][1],'mars')
        with patch.dict(os.environ,{'CUDA_VISIBLE_DEVICES':'0'},clear=True), \
             patch.object(startup,'runtime_libraries',return_value=self.libs), patch.object(startup,'launch') as launch:
            with self.assertRaisesRegex(startup.EvidenceError,'solver_environment_changed'): profile.run(self.pair,{})
            launch.assert_not_called()


if __name__ == '__main__': unittest.main()
