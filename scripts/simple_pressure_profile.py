#!/usr/bin/env python3
"""Compare original-control startup fields with an explicit GPU pressure configuration."""
import copy
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import sys

import simple_hypre_defaults as defaults
import simple_pressure_probe as pressure
import simple_pressure_settings as settings
from simple_pressure_replay import reference_controls
import simple_startup_probe as startup
from simple_snapshot_compare import EvidenceError, band, digest, options, require
from simple_public_diagnostics import SafeParser

SCHEMA = 'mars-simple-pressure-profile-v1'
HISTORY = 'original_twenty_steps'
DEFAULT_KEYS = dict(coarsentype='coarsen_type', interptype='interp_type', relaxorder='relax_order',
    pmax='p_max_elmts', maxlevels='max_levels', mincoarsesize='min_coarse_size', maxcoarsesize='max_coarse_size',
    numfunctions='num_functions', coarsencutfactor='coarsen_cut_factor', cycletype='cycle_type', fcycle='fcycle',
    strongthreshold='strong_threshold', truncfactor='trunc_factor', jacobitruncthreshold='jacobi_trunc_threshold',
    maxrowsum='max_row_sum', aggnumlevels='agg_num_levels', agginterptype='agg_interp_type',
    aggtruncfactor='agg_trunc_factor', numpaths='num_paths', keeptranspose='keep_transpose', nodal='nodal',
    nodaldiag='nodal_diag', coarserelax='relax_coarse', relax_down='relax_down', relax_up='relax_up',
    sweeps_1='sweeps_down', sweeps_2='sweeps_up', sweeps_3='sweeps_coarse')
REAL_KEYS = {'rtol', 'atol', 'strongthreshold', 'truncfactor', 'jacobitruncthreshold', 'maxrowsum', 'aggtruncfactor'}
KEYS = set(DEFAULT_KEYS) | {'method', 'kdim', 'miniter', 'maxiter', 'rtol', 'atol', 'relaxtype', 'numsweeps'}
DEFAULT_FILES = ('public.json', 'probe.log', 'private-provenance.json', 'private-probe-libraries.json',
                 'compile.exit', 'probe.exit')
ERRORS = {'original_baseline_required', 'defaults_evidence', 'unsupported_gpu_profile', 'profile_identity',
          'runtime_profile_mismatch', 'executable_changed', 'libraries_changed', 'solver_environment_changed',
          'launcher_changed', 'launcher_exit', 'runtime_recovery_mismatch', 'initial_fields_mismatch',
          'reference_pressure_log', 'history_required', 'unsupported_momentum_profile'}
ERRORS |= startup.PREPARATION_ERRORS | pressure.ERRORS


def preparation_failure(error):
    # Report only fixed categories and locations in our public source, never locals or exception text.
    root = Path(__file__).resolve().parent
    sources = {str(root / name): name for name in (
        'simple_pressure_profile.py', 'simple_pressure_probe.py', 'simple_pressure_settings.py',
        'simple_startup_probe.py', 'simple_snapshot_compare.py', 'simple_hypre_defaults.py',
        'simple_pressure_replay.py', 'prepare_simple_deck.py')}
    frames = []
    trace = error.__traceback__
    while trace is not None:
        name = sources.get(str(Path(trace.tb_frame.f_code.co_filename).resolve()))
        if name is not None:
            frames.append(dict(source=name, line=trace.tb_lineno))
        trace = trace.tb_next
    kind = next((label for cls, label in (
        (EvidenceError, 'evidence_check'), (ImportError, 'import_error'),
        (FileNotFoundError, 'file_missing'), (FileExistsError, 'output_exists'),
        (PermissionError, 'permission_denied'), (OSError, 'os_error'),
        (ValueError, 'value_error'), (KeyError, 'key_error'), (TypeError, 'type_error'))
        if isinstance(error, cls)), 'other')
    label = str(error) if isinstance(error, EvidenceError) and str(error) in ERRORS else 'private_preflight_or_capture'
    return dict(failed_check=label, preparation_exception_kind=kind, preparation_source_locations=frames[-12:])


def profile_violations(v, reference_amg=False):
    if set(v) != KEYS:
        return ['profile_keys']
    invalid = []
    for k, value in v.items():
        if (type(value) not in (int, float) or not math.isfinite(value)
                or (k not in REAL_KEYS and not (value == int(value) and 0 <= value < 2**31))):
            invalid.append(k)
    if invalid:
        return sorted(invalid)
    checks = dict(method=v['method'] in (0, 1), kdim=v['kdim'] > 0, maxiter=v['maxiter'] > 0,
        miniter=v['miniter'] <= v['maxiter'], rtol=0 < v['rtol'] < 1, atol=v['atol'] >= 0,
        coarsentype=v['coarsentype'] == (10 if reference_amg else 8), interptype=v['interptype'] == 6,
        numfunctions=v['numfunctions'] == 1, cycletype=v['cycletype'] == 1,
        keeptranspose=v['keeptranspose'] in (0, 1), maxlevels=v['maxlevels'] > 0,
        maxcoarsesize=v['maxcoarsesize'] > 0, mincoarsesize=v['mincoarsesize'] <= v['maxcoarsesize'],
        numpaths=v['numpaths'] > 0, coarserelax=v['coarserelax'] == (9 if reference_amg else 18))
    checks.update((k, v[k] == 0) for k in ('relaxorder', 'aggnumlevels', 'nodal', 'nodaldiag', 'fcycle'))
    checks.update((k, v[k] in (6, 13, 14, 18)) for k in ('relaxtype', 'relax_down', 'relax_up'))
    checks.update((k, v[k] > 0) for k in ('numsweeps', 'sweeps_1', 'sweeps_2', 'sweeps_3'))
    checks.update((k, 0 <= v[k] <= 1) for k in REAL_KEYS - {'rtol', 'atol'})
    return sorted(k for k, passed in checks.items() if not passed)


def validate(v, reference_amg=False):
    require(not profile_violations(v, reference_amg), 'unsupported_gpu_profile')


def read_defaults(directory, reference):
    # Recheck saved probe evidence without requiring both uenv mounts to be visible.
    public = startup.read_json(directory / 'public.json')
    raw = defaults.raw_output(directory / 'probe.log')
    projection = defaults.public_defaults(raw)
    metadata = startup.read_json(directory / 'private-provenance.json')
    libraries = startup.read_json(directory / 'private-probe-libraries.json')
    require(public.get('schema') == defaults.SCHEMA and public.get('comparison_status') == 'completed'
        and public.get('failed_check') == 'none' and all(public.get(k) == v for k, v in projection.items())
        and all(public.get(k) is True for k in ('saved_capture_identity_verified', 'captured_libraries_unchanged',
            'loaded_hypre_library_matches_capture', 'mpi_library_identity_verified', 'mpi_enabled'))
        and all((directory / k).read_text().strip() == '0' for k in ('compile.exit', 'probe.exit')),
        'defaults_evidence')
    captured = metadata['capture']
    require(all(captured[k] == reference[k] for k in ('executable', 'executable_sha256', 'libraries')),
            'defaults_evidence')
    try:
        defaults.probe_dependencies(libraries, reference['libraries'], True)
    except ValueError:
        raise EvidenceError('defaults_evidence')
    for kind, path in (('hypre', raw['library_path']), ('mpi', raw['mpi_library_path'])):
        require(Path(path).is_absolute() and defaults.library_kind(path) == kind, 'defaults_evidence')
        # dladdr and ldd may spell the same library through different symlinks.
        # The successful saved probe checked the actual loaded path by content.
        if path in libraries:
            require(libraries[path] in {h for p, h in reference['libraries'].items()
                                       if defaults.library_kind(p) == kind}, 'defaults_evidence')
    # The omitted one-level relaxation/sweep fallback below is inspected in this version.
    require(projection['version'] == '3.1.0' and projection['build_cuda'] is False, 'defaults_evidence')
    return projection['defaults']


def resolved_values(config, target, measured, reference_amg=False):
    explicit = reference_controls(config, target)
    v = {key: measured['amg_' + name] for key, name in DEFAULT_KEYS.items()}
    method = 'flexgmres' if explicit['method'] else 'gmres'
    v.update(kdim=measured[method + '_restart_dimension'], miniter=measured[method + '_minimum_iterations'],
             relaxtype=6, numsweeps=1)
    v.update(explicit)
    # Hypre generic setters also overwrite the cycle settings, in this order.
    if 'relaxtype' in explicit:
        v.update(relax_down=v['relaxtype'], relax_up=v['relaxtype'], coarserelax=9)
    if 'numsweeps' in explicit:
        v.update(sweeps_1=v['numsweeps'], sweeps_2=v['numsweeps'], sweeps_3=1)
    changes = []
    if not reference_amg and v['coarsentype'] == 10:
        v['coarsentype'] = 8
        changes.append('hmis_to_device_pmis')
    if not reference_amg and v['coarserelax'] == 9:
        v['coarserelax'] = 18
        changes.append('coarse_direct_to_device_l1_jacobi')
    return v, changes


def resolve(config, target, measured, reference_amg=False):
    values, changes = resolved_values(config, target, measured, reference_amg)
    validate(values, reference_amg)
    return values, changes


PETSC_METHODS = {'gmres', 'fgmres', 'lgmres', 'cg', 'fcg', 'bcgs', 'ibcgs', 'bcgsl',
                 'richardson', 'chebyshev', 'preonly', 'minres', 'gcr', 'lsqr', 'lsmr'}
PETSC_PRECONDITIONERS = {'none', 'jacobi', 'bjacobi', 'pbjacobi', 'vpbjacobi', 'sor',
                         'ilu', 'icc', 'lu', 'cholesky', 'asm', 'gasm', 'gamg', 'hypre',
                         'mg', 'fieldsplit', 'redundant', 'composite', 'svd', 'ksp', 'shell'}
PETSC_ENUMS = dict(ksp_type=PETSC_METHODS, pc_type=PETSC_PRECONDITIONERS,
    sub_ksp_type=PETSC_METHODS, sub_pc_type=PETSC_PRECONDITIONERS,
    sub_sub_ksp_type=PETSC_METHODS, sub_sub_pc_type=PETSC_PRECONDITIONERS,
    ksp_norm_type={'none', 'preconditioned', 'unpreconditioned', 'natural'},
    ksp_pc_side={'left', 'right', 'symmetric'},
    pc_hypre_type={'boomeramg', 'ams', 'ads', 'euclid', 'pilut', 'ilu'},
    pc_asm_type={'basic', 'restrict', 'interpolate', 'none'},
    pc_factor_mat_solver_type={'petsc', 'mumps', 'superlu', 'superlu_dist', 'cusparse', 'mkl_pardiso'},
    sub_pc_factor_mat_solver_type={'petsc', 'mumps', 'superlu', 'superlu_dist', 'cusparse', 'mkl_pardiso'},
    mat_type={'aij', 'seqaij', 'mpiaij', 'aijcusparse', 'seqaijcusparse', 'mpiaijcusparse',
              'baij', 'seqbaij', 'mpibaij', 'aijkokkos', 'seqaijkokkos', 'mpiaijkokkos'},
    vec_type={'standard', 'seq', 'mpi', 'cuda', 'seqcuda', 'mpicuda', 'kokkos', 'seqkokkos', 'mpikokkos'})
PETSC_CONTROL_NAMES = {'ksp_rtol', 'ksp_atol', 'ksp_divtol', 'ksp_max_it', 'ksp_min_it',
    'ksp_gmres_restart', 'ksp_gmres_modifiedgramschmidt', 'ksp_gmres_classicalgramschmidt',
    'ksp_gmres_cgs_refinement_type', 'ksp_initial_guess_nonzero', 'ksp_converged_use_initial_residual_norm',
    'ksp_converged_use_min_initial_residual_norm', 'ksp_converged_maxits', 'ksp_error_if_not_converged',
    'pc_factor_levels', 'pc_factor_fill', 'pc_factor_mat_ordering_type', 'pc_factor_shift_type',
    'pc_factor_shift_amount', 'pc_asm_overlap', 'pc_bjacobi_blocks', 'sub_pc_factor_levels',
    'sub_pc_factor_fill', 'sub_pc_factor_mat_ordering_type', 'sub_pc_factor_shift_type',
    'sub_pc_factor_shift_amount', 'sub_ksp_rtol', 'sub_ksp_atol', 'sub_ksp_max_it',
    'ksp_monitor', 'ksp_monitor_true_residual', 'ksp_view', 'ksp_converged_reason', 'options_left'}


def inspect_petsc_options(config):
    # ContextPETSc lowercases flat option keys/values before PetscOptionsSetValue.
    # Omission is unresolved here: do not substitute Hypre or PETSc defaults.
    report = dict(solver='unspecified', preconditioner='unspecified',
        option_scope='declared_petsc_options_not_runtime_settings',
        runtime_settings_verified=False, tolerance_source='context_after_KSPSetFromOptions',
        minimum_iterations_forwarded_by_context=False, options_structure_valid=False)
    raw = config.get('options', {})
    if not isinstance(raw, dict):
        return report
    options = {}
    for key, value in raw.items():
        if not isinstance(key, str) or key.lower() in options:
            return report
        options[key.lower()] = value
    report.update(options_structure_valid=all(type(v) in (str, int, float, bool) for v in options.values()),
        unrecognized_option_keys_present=bool(set(options) - set(PETSC_ENUMS) - PETSC_CONTROL_NAMES),
        declared_control_names=sorted(set(options) & PETSC_CONTROL_NAMES),
        declared_petsc_options={})
    for key, allowed in PETSC_ENUMS.items():
        if key in options:
            value = options[key].lower() if isinstance(options[key], str) else ''
            report['declared_petsc_options'][key] = value if value in allowed else 'other_or_invalid'
    report['solver'] = report['declared_petsc_options'].get('ksp_type', 'unspecified')
    report['preconditioner'] = report['declared_petsc_options'].get('pc_type', 'unspecified')
    return report


def petsc_momentum_values(config):
    from prepare_simple_deck import reference_linear_target
    arguments = reference_linear_target({'momentum': config}, {}, 'momentum', 'momentum')
    target = options(arguments)
    maximum = config.get('max_iterations', 20)
    require(type(maximum) in (int, float) and math.isfinite(maximum) and maximum == int(maximum)
            and 1 <= maximum <= 2147483647, 'unsupported_momentum_profile')
    return dict(petsc_fgmres_bjacobi=1, rtol=float(target['--momentum-linear-rtol']),
                atol=float(target['--momentum-linear-atol']), maxiter=int(maximum)), arguments


def is_petsc_profile(values):
    return values.get('petsc_fgmres_bjacobi') == 1


def public_library_name(path):
    name = Path(path).name
    stem = re.sub(r'\.so(?:\.[0-9]+)*$', '', name)
    known = {'libpetsc', 'libHYPRE', 'libmpi', 'libmpich', 'libmpifort', 'libmpi_gtl_cuda',
        'libmpi_gtl_hsa', 'libc', 'libm', 'libdl', 'librt', 'libpthread', 'libstdc++', 'libgcc_s',
        'libgomp', 'libomp', 'libcuda', 'libcudart', 'libcusparse', 'libcublas', 'libcublasLt',
        'libcusolver', 'libcurand', 'libcufft', 'libnvrtc', 'libnvrtc-builtins', 'libnvJitLink',
        'libnvToolsExt', 'libnvtx3interop', 'libnetcdf', 'libhdf5', 'libhdf5_hl', 'libexodus',
        'libz', 'libzstd', 'libsz', 'libaec', 'libcurl', 'libxml2', 'libblas', 'liblapack',
        'libopenblas', 'libf2cblas', 'libf2clapack', 'libgfortran', 'libquadmath',
        'libhwloc', 'libfabric', 'libpmix', 'libucp', 'libuct', 'libucs', 'libucm',
        'libmars', 'libmars_unstructured', 'libkokkoscore', 'libkokkoscontainers', 'libkokkossimd',
        'libtrilinos_kokkoscore', 'libtrilinos_kokkoscontainers', 'libtrilinos_kokkossimd',
        'libmetis', 'libparmetis', 'libfmt', 'libyaml-cpp'}
    if stem != name and (stem in known or re.fullmatch(
            r'lib(?:mpi|mpifort)_(?:gnu(?:_[0-9]+)?|cray)|libHYPRE-[0-9.]+|ld-linux-aarch64', stem)):
        return name
    return defaults.library_kind(path) or 'other_runtime_library'


def profile_library_checks(libraries, original, petsc):
    public = dict(policy='preserve_captured_dependencies_add_petsc' if petsc else 'exact_capture',
                  matched=False, record_valid=False)
    if not all(isinstance(v, dict) and v and all(isinstance(p, str) and isinstance(h, str)
            and re.fullmatch(r'[0-9a-f]{64}', h) for p, h in v.items()) for v in (libraries, original)):
        return public
    public['record_valid'] = True
    changed = {path for path, value in original.items() if path in libraries and libraries[path] != value}
    missing = set(original) - set(libraries)
    added = set(libraries) - set(original)
    relocated = set()
    # Match one old object to one byte-identical object, never hide a changed object at its old path.
    for path in sorted(missing):
        candidates = [p for p in added if Path(p).name == Path(path).name
                      and libraries[p] == original[path]]
        if len(candidates) == 1:
            relocated.add(path)
            added.remove(candidates[0])
    missing -= relocated
    replaced = {p for p in missing if any(Path(q).name == Path(p).name and libraries[q] != original[p]
                                         for q in libraries)}
    protected = {p for p in added if defaults.library_kind(p) in ('mpi', 'hypre')}
    petsc_count = sum(Path(path).name.startswith('libpetsc.') for path in libraries)
    preserved = not changed and not missing
    public.update(captured_dependencies_preserved=preserved,
        changed_captured_libraries=sorted({public_library_name(p) for p in changed | replaced}),
        missing_captured_libraries=sorted({public_library_name(p) for p in missing - replaced}),
        relocated_identical_libraries=sorted({public_library_name(p) for p in relocated}),
        additional_libraries=sorted({public_library_name(p) for p in added}),
        additional_solver_libraries=sorted({public_library_name(p) for p in protected}),
        petsc_library_present=petsc_count > 0, petsc_library_unique=petsc_count == 1,
        matched=(preserved and not protected and petsc_count == 1) if petsc else libraries == original)
    return public


def profile_libraries_match(libraries, original, petsc):
    return profile_library_checks(libraries, original, petsc)['matched']


def inspect_configuration(config, measured, stage):
    from prepare_simple_deck import reference_linear_target

    def label(value, allowed):
        value = value.lower() if isinstance(value, str) else ''
        return value if value in allowed else 'other_or_invalid'

    report = dict(mapping_supported=False, failed_component='definition')
    if not isinstance(config, dict):
        return report
    outer = config.get('options', {})
    report.update(family=label(config.get('family'), {'hypre', 'petsc', 'trilinos', 'amgsolver', 'gmres'}),
        normalize_matrix_requested=config.get('normalize_matrix', False) is not False,
        diagonal_scaling_requested=config.get('diagonal_scaling', False) is not False)
    if report['family'] != 'hypre':
        report['failed_component'] = 'family'
        if report['family'] == 'petsc':
            report.update(inspect_petsc_options(config))
            if stage == 'momentum':
                report['failed_component'] = 'petsc_profile'
                try:
                    petsc_momentum_values(config)
                except (ValueError, KeyError, TypeError):
                    pass
                else:
                    report.update(mapping_supported=True, failed_component='none',
                        required_build='cuda_petsc', gpu_adaptations=['aij_to_device_aijcusparse', 'local_ilu_to_device_cusparse'])
        else:
            report.update(solver='unresolved', preconditioner='unresolved')
        return report
    report['failed_component'] = 'options_structure'
    if not isinstance(outer, dict):
        return report
    try:
        solver = settings.lower_options(outer.get('solver', {}))
        precond = settings.lower_options(outer.get('precond', {}))
    except ValueError:
        return report
    supported_amg = {key for key, _, _ in settings.AMG.values() if key} | set(settings.AMG_LIBRARY)
    report.update(solver=label(solver.get('type', 'gmres'), {'gmres', 'flexgmres', 'boomeramg', 'mgr'}),
        preconditioner=label(precond.get('type', 'none'), {'none', 'boomeramg', 'mgr'}),
        unrecognized_option_keys_present=bool(set(outer) - {'solver', 'precond'}
            or set(solver) - {'type', 'kdim', 'printlevel', 'logging'}
            or set(precond) - supported_amg - {'type', 'printlevel', 'logging', 'maxiter', 'tol'}))
    for component, passed in (
            ('family', report['family'] == 'hypre'),
            ('normalize_matrix', not report['normalize_matrix_requested']),
            ('diagonal_scaling', not report['diagonal_scaling_requested']),
            ('solver', report['solver'] in ('gmres', 'flexgmres')),
            ('preconditioner', report['preconditioner'] == 'boomeramg'),
            ('option_keys', not report['unrecognized_option_keys_present'])):
        if not passed:
            report['failed_component'] = component
            return report
    report['failed_component'] = 'residual_target_or_solver_definition'
    try:
        arguments = options(reference_linear_target({'equation': config}, {}, 'equation', stage))
        target = tuple(float(arguments['--' + stage + '-linear-' + key]) for key in ('atol', 'rtol'))
    except ValueError:
        return report
    report['failed_component'] = 'reference_options'
    try:
        values, changes = resolved_values(config, target, measured)
    except ValueError:
        return report
    report['unsupported_gpu_controls'] = profile_violations(values)
    report['gpu_adaptations'] = changes
    report['mapping_supported'] = not report['unsupported_gpu_controls']
    report['failed_component'] = 'none' if report['mapping_supported'] else 'gpu_controls'
    return report


def inspect_profiles(baseline, probe, public):
    before = pressure.baseline_records(baseline)
    measured = read_defaults(probe, before[2]['openaccel'])
    deck = startup.load_deck((baseline / 'reference/input.i').read_bytes())
    reports = {}
    for stage, equation in (('momentum', 'coupled_navier_stokes'), ('pressure', 'pressure_correction')):
        try:
            config = pressure.resolved_solver(deck, equation)
        except (ValueError, KeyError, TypeError):
            reports[stage] = dict(mapping_supported=False, failed_component='definition')
        else:
            reports[stage] = inspect_configuration(config, measured, stage)
    require(pressure.baseline_records(baseline) == before
            and read_defaults(probe, before[2]['openaccel']) == measured, 'profile_identity')
    public.update(comparison_status='inspection_complete', failed_check='none', profiles=reports,
        both_profiles_supported=all(r['mapping_supported'] for r in reports.values()),
        solver_launched=False, scope='saved_configuration_mapping_not_runtime_or_field_parity')


def text_profile(values):
    return ''.join('{} {:.17g}\n'.format(k, values[k]) for k in sorted(values))


def baseline_data(baseline, probe, reference_amg=False):
    record, case, launches = pressure.baseline_records(baseline)
    require('pressure_experiment' not in record and 'pressure_profile' not in record
            and launches['mars']['ranks'] == launches['openaccel']['ranks'], 'original_baseline_required')
    require(launches['mars']['environment'].get('MARS_HYPRE_SPMV_VENDOR', '0') == '0', 'unsupported_gpu_profile')
    config, arguments = settings.saved_configuration(baseline)
    target = (float(arguments['--pressure-linear-atol']), float(arguments['--pressure-linear-rtol']))
    values, adaptations = resolve(config, target, read_defaults(probe, launches['openaccel']), reference_amg)
    return record, case, launches, values, adaptations


def momentum_data(baseline, probe, reference):
    try:
        deck = startup.load_deck((baseline / 'reference/input.i').read_bytes())
        config = pressure.resolved_solver(deck, 'coupled_navier_stokes')
        if str(config.get('family', '')).lower() == 'petsc':
            values, arguments = petsc_momentum_values(config)
            return values, ['aij_to_device_aijcusparse', 'local_ilu_to_device_cusparse'], arguments + ['--momentum-backend', 'petsc']
        from prepare_simple_deck import reference_momentum_target
        solver = deck['simulation']['solver']
        arguments = reference_momentum_target(solver['solver_control']['advanced_options']['linear_solver_settings'], solver)
        target = options(arguments)
        values, adaptations = resolve(config, (float(target['--momentum-linear-atol']),
                                              float(target['--momentum-linear-rtol'])), read_defaults(probe, reference))
        return values, adaptations, arguments
    except Exception:
        # Never print a private solver name or approximate a different backend.
        raise EvidenceError('unsupported_momentum_profile')


def pair_content(baseline, pair, record, case, identity):
    case = copy.deepcopy(case)
    record = dict(record)
    if identity.get('kind') == HISTORY:
        original = startup.load_deck((baseline / 'reference/input.i').read_bytes())
        require(startup.load_deck((pair / 'reference/input.i').read_bytes()) == pressure.history_deck(original),
                'profile_identity')
        case['deck_sha256'] = digest(pair / 'reference/input.i')
        record.update(steps=startup.STEPS, first_step_audit=False, deck_sha256=case['deck_sha256'],
                      changed_controls=['iteration_limits', 'output_control', 'linear_system_output', 'gpu_pressure_profile'])
    case['arguments'] += ['--pressure-solver-profile', str(pair / 'pressure.profile')]
    case['pressure_solver_profile_sha256'] = identity['sha256']
    if identity.get('pressure_amg_mode') == 'reference':
        case['arguments'] += ['--pressure-amg-mode', 'reference']
        record['changed_controls'] = record['changed_controls'] + ['reference_pressure_amg']
    if 'momentum' in identity:
        _, _, arguments = momentum_data(baseline, pair / 'defaults', startup.verified_launch(baseline, 'openaccel'))
        case['arguments'] += arguments + ['--momentum-solver-profile', str(pair / 'momentum.profile')]
        case['momentum_linear_policy'] = 'reference'
        if '--momentum-backend' in arguments:
            case['momentum_backend'] = 'petsc'
        case['momentum_solver_profile_sha256'] = identity['momentum']['sha256']
        record['changed_controls'] = record['changed_controls'] + ['gpu_momentum_profile']
        case['notes'] = [note for note in case.get('notes', []) if not note.startswith((
            'MARS retains its momentum solver settings', 'Reference pressure residual target copied;'))] + [
            'Both linear targets are copied; explicit GPU profiles record backend adaptations.',
            'Nonlinear MARS convergence norms are unchanged; this compares saved iteration states.']
        if identity.get('pressure_amg_mode') == 'reference':
            case['notes'][-2] = ('Both linear targets are copied; reference pressure AMG permits '
                                 'Hypre host staging for HMIS and the direct coarse solve.')
    return record, case


def prepare(baseline, probe, executable, output, history=False, both_linear_solvers=False, public=None, reference_amg=False):
    require(not reference_amg or (not history and both_linear_solvers), 'unsupported_gpu_profile')
    record, case, launches, values, adaptations = baseline_data(baseline, probe, reference_amg)
    momentum = momentum_data(baseline, probe, launches['openaccel']) if both_linear_solvers else None
    require(executable.is_file() and os.access(str(executable), os.X_OK), 'executable_changed')
    libraries = startup.runtime_libraries(executable)
    checks = profile_library_checks(libraries, launches['mars']['libraries'],
                                    momentum is not None and is_petsc_profile(momentum[0]))
    if public is not None:
        public['library_checks'] = checks
    require(checks['matched'], 'libraries_changed')
    output.mkdir(mode=0o700)
    (output / 'reference').mkdir(mode=0o700)
    if history:
        import yaml
        original = startup.load_deck((baseline / 'reference/input.i').read_bytes())
        (output / 'reference/input.i').write_text(yaml.safe_dump(pressure.history_deck(original), default_flow_style=False))
    else:
        shutil.copyfile(baseline / 'reference/input.i', output / 'reference/input.i')
    (output / 'defaults').mkdir(mode=0o700)
    for name in DEFAULT_FILES:
        shutil.copyfile(probe / name, output / 'defaults' / name)
    (output / 'pressure.profile').write_text(text_profile(values))
    identity = dict(schema=SCHEMA, reference_pair=str(baseline), reference_pair_sha256=digest(baseline / 'pair.json'),
        sha256=digest(output / 'pressure.profile'), adaptations=adaptations,
        defaults_hashes={name: digest(output / 'defaults' / name) for name in DEFAULT_FILES},
        baseline_launch_hashes={s: digest(baseline / ('reference' if s == 'openaccel' else 'mars') / 'launch.json') for s in launches},
        executable=str(executable), executable_sha256=digest(executable), libraries=libraries)
    if reference_amg:
        identity['pressure_amg_mode'] = 'reference'
        if public is not None:
            public.update(pressure_amg_mode='reference', hypre_host_staging_permitted=True)
    if momentum is not None:
        (output / 'momentum.profile').write_text(text_profile(momentum[0]))
        identity['momentum'] = dict(sha256=digest(output / 'momentum.profile'), adaptations=momentum[1])
    if history:
        identity['kind'] = HISTORY
    record, case = pair_content(baseline, output, record, case, identity)
    startup.write_json(output / 'case.json', case)
    startup.write_json(output / 'pair.json', dict(record, case_sha256=digest(output / 'case.json'), pressure_profile=identity))
    check_inputs(output)


def check_inputs(pair):
    record = startup.pair_inputs(pair)
    identity = record['pressure_profile']
    require(identity['schema'] == SCHEMA and identity.get('kind') in (None, HISTORY), 'profile_identity')
    require('pressure_amg_mode' not in identity or (identity['pressure_amg_mode'] == 'reference'
            and 'kind' not in identity and 'momentum' in identity), 'profile_identity')
    reference_amg = identity.get('pressure_amg_mode') == 'reference'
    baseline = Path(identity['reference_pair'])
    require(digest(baseline / 'pair.json') == identity['reference_pair_sha256'], 'profile_identity')
    require(identity['defaults_hashes'] == {name: digest(pair / 'defaults' / name) for name in DEFAULT_FILES}, 'defaults_evidence')
    old, case, launches, values, adaptations = baseline_data(baseline, pair / 'defaults', reference_amg)
    require(identity['baseline_launch_hashes'] == {s: digest(baseline / ('reference' if s == 'openaccel' else 'mars') / 'launch.json')
                                                  for s in launches}, 'profile_identity')
    require(text_profile(values) == (pair / 'pressure.profile').read_text()
            and identity['adaptations'] == adaptations, 'profile_identity')
    petsc = False
    if 'momentum' in identity:
        momentum, changes, _ = momentum_data(baseline, pair / 'defaults', launches['openaccel'])
        petsc = is_petsc_profile(momentum)
        require(text_profile(momentum) == (pair / 'momentum.profile').read_text()
                and identity['momentum'] == dict(sha256=digest(pair / 'momentum.profile'), adaptations=changes), 'profile_identity')
    require(profile_libraries_match(identity['libraries'], launches['mars']['libraries'], petsc), 'profile_identity')
    old, expected_case = pair_content(baseline, pair, old, case, identity)
    require(startup.read_json(pair / 'case.json') == expected_case
            and record == dict(old, case_sha256=digest(pair / 'case.json'), pressure_profile=identity), 'profile_identity')
    return baseline, launches, values, identity


def check_launch(pair, baseline, old, identity, solver='mars'):
    record = startup.verified_launch(pair, solver)
    expected = identity if solver == 'mars' else old
    require(all(record[k] == expected[k] for k in ('executable', 'executable_sha256', 'libraries')), 'profile_identity')
    require(record['ranks'] == old['ranks'] and pressure.launcher(record, pair, solver) ==
            pressure.launcher(old, baseline, solver), 'launcher_changed')
    require(pressure.solver_environment(record['environment']) == pressure.solver_environment(old['environment']),
            'solver_environment_changed')
    return record


def run(pair, public, solver='mars'):
    baseline, launches, _, identity = check_inputs(pair)
    require(solver == 'mars' or (solver == 'openaccel' and identity.get('kind') == HISTORY), 'profile_identity')
    old = launches[solver]
    expected = identity if solver == 'mars' else old
    exe = Path(expected['executable'])
    require(exe.is_file() and os.access(str(exe), os.X_OK) and digest(exe) == expected['executable_sha256'], 'executable_changed')
    public['library_checks'] = profile_library_checks(startup.runtime_libraries(exe), expected['libraries'], False)
    require(public['library_checks']['matched'], 'libraries_changed')
    environment = dict(os.environ)
    for key in ('MARS_OPENACCEL_EXPORT_DIR', 'MARS_OPENACCEL_PUBLIC_FIXTURE'):
        environment.pop(key, None)
    if solver == 'openaccel':
        environment.update(OMP_NUM_THREADS='1', OMP_PROC_BIND='close', OMP_PLACES='cores')
    changes = pressure.environment_changes(environment, old['environment'])
    public['environment_check'] = pressure.environment_report(changes)
    if solver == 'mars':
        environment = pressure.restore_environment(environment, old['environment'])
    require(not pressure.environment_changes(environment, old['environment']), 'solver_environment_changed')
    launcher = pressure.launcher(old, baseline, solver)
    # Keep the complete direct command private; the wrapper also supplies the checked environment.
    with (pair / (solver + '-launch-command.txt')).open('x') as stream:
        stream.write(' '.join(shlex.quote(word) for word in launcher + [str(exe)] +
                              startup.solver_arguments(pair, solver)) + '\n')
    startup.launch(pair, solver, exe, old['ranks'], launcher, environment)
    check_launch(pair, baseline, old, identity, solver)


def check_runtime(pair, ranks, values, history=False, stage='pressure'):
    expected = dict(values)
    petsc = is_petsc_profile(values)
    if not petsc:
        expected['effective_relax_1'] = expected.pop('relax_down')
        expected['effective_relax_2'] = expected.pop('relax_up')
        expected['effective_relax_3'] = expected['coarserelax']
    prefix = 'flow-' + stage + '-'
    paths = sorted((pair / 'mars').glob(prefix + 'rank-*.settings'))
    names = [(prefix + 'rank-{}.settings').format(rank) for rank in range(ranks)]
    if history:
        paths += sorted((pair / 'mars').glob(prefix + 'step-*-rank-*.settings'))
        names = [(prefix + 'step-{}-rank-{}.settings').format(step, rank)
                 for step in range(1, startup.STEPS + 1) for rank in range(ranks)]
    require({p.name for p in paths} == set(names), 'runtime_profile_mismatch')
    for name in names:
        lines = (pair / 'mars' / name).read_text().splitlines()
        parsed = [line.split() for line in lines]
        require(all(len(row) == 2 for row in parsed), 'runtime_profile_mismatch')
        actual = {k: float(v) for k, v in parsed}
        if petsc:
            extra = {'petsc_release', 'restart', 'divergence_tolerance', 'local_blocks', 'sub_preonly',
                     'sub_ilu_levels', 'right_preconditioning', 'unpreconditioned_norm', 'device'}
            require(len(actual) == len(parsed) and set(actual) == set(expected) | extra
                    and all(math.isfinite(v) for v in actual.values()) and all(actual[k] == v for k, v in expected.items())
                    and actual['local_blocks'] == 1 and actual['sub_preonly'] == 1 and actual['sub_ilu_levels'] == 0
                    and actual['right_preconditioning'] == 1 and actual['unpreconditioned_norm'] == 1 and actual['device'] == 1
                    and actual['restart'] >= 1 and actual['divergence_tolerance'] > 1
                    and actual['petsc_release'] >= 32400, 'runtime_profile_mismatch')
            continue
        require(len(actual) == len(parsed) and set(actual) == set(expected) | {'hypre_release', 'effective_levels'}
                and all(math.isfinite(v) for v in actual.values()) and all(actual[k] == v for k, v in expected.items())
                and actual['effective_levels'] >= 1, 'runtime_profile_mismatch')


def compare(pair, public, details):
    baseline, launches, values, identity = check_inputs(pair)
    record = check_launch(pair, baseline, launches['mars'], identity)
    history = identity.get('kind') == HISTORY
    if history:
        check_launch(pair, baseline, launches['openaccel'], identity, 'openaccel')
        lines = (pair / 'mars/run.log').read_text(errors='replace').splitlines()
        require([line for line in lines if line.startswith('pressure_refinement=')] ==
                ['pressure_refinement=0 pressure_expansion=0'] and
                not any(line.startswith(('[simple-pressure-expansion]', '[simple-pressure-refinement]')) for line in lines),
                'runtime_recovery_mismatch')
    check_runtime(pair, record['ranks'], values, history)
    if 'momentum' in identity:
        momentum, _, _ = momentum_data(baseline, pair / 'defaults', launches['openaccel'])
        check_runtime(pair, record['ranks'], momentum, history, 'momentum')
    result = {}
    startup.compare(pair, result, details, reference_pair=None if history else baseline)
    check_inputs(pair)
    if history:
        require(result['initial_field_parity_verified'], 'initial_fields_mismatch')
        saved = startup.read_json(details / 'comparison-private.json')
        public.update(experiment=HISTORY, recovery_enabled=False,
            pressure_residuals_recomputed_at_every_step=False,
            outcome='twenty_step_fields_match' if result['all_twenty_snapshots_match'] else 'twenty_step_fields_differ',
            snapshot_agreement=[dict(iteration=row['iteration'], **{
                field + '_band': band(row['errors'][field]) for field in
                ('velocity_max_scaled', 'velocity_rms_scaled', 'pressure_max_scaled', 'pressure_rms_scaled')})
                for row in saved['snapshots']])
    public.update(comparison_status='completed', failed_check='none', comparison=result,
        original_pressure_targets_preserved=True, recorded_pressure_settings_verified=True,
        recorded_gpu_pressure_settings_verified=identity.get('pressure_amg_mode') != 'reference',
        pressure_amg_mode=identity.get('pressure_amg_mode', 'device'),
        hypre_host_staging_permitted=identity.get('pressure_amg_mode') == 'reference',
        gpu_adaptations=identity['adaptations'], reference_solver_launched=history,
        identical_linear_solvers_verified=False, nonlinear_convergence_verified=False)
    if 'momentum' in identity:
        public.update(original_momentum_targets_preserved=True, recorded_gpu_momentum_settings_verified=True,
                      momentum_gpu_adaptations=identity['momentum']['adaptations'])


def inspect_linear(pair, public):
    from simple_pressure_log import inspect_pressure_log
    baseline, launches, values, identity = check_inputs(pair)
    require(identity.get('kind') == HISTORY, 'history_required')
    records = {solver: check_launch(pair, baseline, launches[solver], identity, solver)
               for solver in ('mars', 'openaccel')}
    check_runtime(pair, records['mars']['ranks'], values, history=True)
    checks = inspect_pressure_log((pair / 'reference/run.log').read_text(errors='replace'),
                                  startup.STEPS, values['rtol'], values['atol'], int(values['maxiter']))
    # Recheck hashed evidence after reading the log; this command never launches a solver.
    require(check_inputs(pair) == (baseline, launches, values, identity), 'profile_identity')
    for solver in records:
        require(check_launch(pair, baseline, launches[solver], identity, solver) == records[solver], 'profile_identity')
    public.update(comparison_status='completed', failed_check='none', experiment=HISTORY,
                  solver_launched=False, reference_pressure_log_checks=checks,
                  original_pressure_targets_preserved=True, nonlinear_convergence_verified=False,
                  openaccel_field_parity_verified=False)


def main(argv=None):
    parser = SafeParser(description=__doc__)
    sub = parser.add_subparsers(dest='action')
    prep = sub.add_parser('prepare')
    for name in ('baseline-pair', 'defaults-probe', 'executable', 'output-dir'):
        prep.add_argument('--' + name, type=Path, required=True)
    prep.add_argument('--history', action='store_true')
    prep.add_argument('--both-linear-solvers', action='store_true',
                      help='Also resolve momentum from the original deck; unsupported backends/settings fail before launch')
    prep.add_argument('--reference-amg', action='store_true',
                      help='Keep reference HMIS and direct coarse solve, permitting Hypre CPU work; requires both linear solvers and a first-step audit')
    prep.add_argument('--output', type=Path,
                      help='Write a public preparation report, including a safe failure location; no solver is launched')
    inspect = sub.add_parser('inspect-profiles')
    for name in ('baseline-pair', 'defaults-probe', 'output'):
        inspect.add_argument('--' + name, type=Path, required=True)
    for action in ('run', 'compare', 'inspect-linear'):
        command = sub.add_parser(action)
        command.add_argument('--pair', type=Path, required=True)
        command.add_argument('--output', type=Path, required=True)
        if action == 'compare':
            command.add_argument('--detail-dir', type=Path, required=True)
        elif action == 'run':
            command.add_argument('--solver', choices=('mars', 'openaccel'), default='mars')
    args = parser.parse_args(argv)
    os.umask(0o077)
    public = dict(schema=SCHEMA, comparison_status='invalid_evidence', failed_check='private_preflight_or_capture')
    try:
        if args.action == 'prepare':
            inputs = (args.baseline_pair.resolve(), args.defaults_probe.resolve(), args.executable.resolve(),
                      args.output_dir.resolve(), args.history, args.both_linear_solvers)
            if args.output is not None:
                # Reserve the report first so an old result cannot be mistaken for this attempt.
                with args.output.open('x') as stream:
                    public['solver_launched'] = False
                    try:
                        prepare(*inputs, public=public, reference_amg=args.reference_amg)
                        public.update(comparison_status='preparation_complete', failed_check='none')
                    except Exception as error:
                        public.update(preparation_failure(error))
                    json.dump(public, stream, indent=2, sort_keys=True, allow_nan=False)
                    stream.write('\n')
                print('Pressure preparation status written. No solver launched; share only the public JSON.')
                return 0 if public['comparison_status'] == 'preparation_complete' else 1
            prepare(*inputs, reference_amg=args.reference_amg)
            print('Original linear profiles prepared. No solver launched; baseline preserved.')
            return 0
        require(args.action in ('run', 'compare', 'inspect-linear', 'inspect-profiles'))
        with args.output.open('x') as stream:
            try:
                if args.action == 'inspect-profiles':
                    inspect_profiles(args.baseline_pair.resolve(), args.defaults_probe.resolve(), public)
                elif args.action == 'run':
                    pair = args.pair.resolve()
                    run(pair, public, args.solver)
                    public.update(comparison_status='capture_complete', failed_check='none')
                elif args.action == 'inspect-linear':
                    inspect_linear(args.pair.resolve(), public)
                else:
                    args.detail_dir.mkdir(mode=0o700)
                    compare(args.pair.resolve(), public, args.detail_dir)
            except Exception as error:
                public['failed_check'] = str(error) if isinstance(error, EvidenceError) and str(error) in ERRORS else 'private_preflight_or_capture'
                if args.action == 'run':
                    try:
                        public['launch_diagnostics'] = pressure.failure_summary(pair, args.solver)
                    except Exception:
                        public['launch_diagnostics'] = None
            json.dump(public, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')
        print('Pressure profile status written. Share only the public JSON.')
        return 0 if public['comparison_status'] in ('completed', 'capture_complete', 'inspection_complete') else 1
    except Exception as error:
        label = str(error) if isinstance(error, EvidenceError) and str(error) in ERRORS else 'private_preflight_or_capture'
        print('ERROR: pressure profile failed (' + label + ').', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
