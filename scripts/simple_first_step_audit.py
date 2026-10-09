"""User-local first-step algebra audit. Numerical data never enters the public report."""
import math
from fractions import Fraction

import numpy as np

from prepare_simple_deck import load_deck
from simple_snapshot_compare import mars_fields, options, reference_fields, require

AUDIT_ERRORS = frozenset('audit_parts audit_format audit_rank_identity audit_dimensions audit_source_ids '
    'audit_ownership audit_offsets audit_columns audit_trailing_data audit_truncated audit_nonfinite audit_zero_row '
    'reference_matrix_count reference_matrix_index_width reference_matrix_dimensions reference_matrix_offsets '
    'reference_matrix_columns audit_residual_arithmetic audit_duplicate_owner audit_field_coverage '
    'reference_solver_ids audit_initial_fields_differ reference_pressure_update mars_predictor_update '
    'mars_velocity_update mars_pressure_update audit_field_arithmetic reference_field_names '
    'reference_ghost_values reference_coordinates reference_saved_iteration reference_node_ids reference_node_coverage '
    'audit_pressure_controls gradient_storage_precision gradient_mesh gradient_arithmetic'.split())


def enable_reference_audit(deck):
    solver = deck['simulation']['solver']
    settings = solver['solver_control']['advanced_options']['linear_solver_settings']
    for equation in ('coupled_navier_stokes', 'pressure_correction'):
        config = next((settings[key] for key in (equation, 'segregated_flow', 'default') if key in settings), None)
        require(isinstance(config, dict), 'reference_solver_configuration')
        if 'lookup' in config:
            config = solver[config['lookup']]
        config['write_system'] = True
    solver['output_control']['output_fields'] = [
        'velocity', 'pressure', 'aux', 'du', 'pressure_correction', 'pressure_correction_gradient']


class MarsPart:
    def __init__(self, path, rank, ranks, nodes):
        require(path.is_file(), 'audit_parts')
        self.data = np.memmap(str(path), dtype='u1', mode='r')
        self.position = 0
        header = self.take('<u8', 7)
        require(header[0] == 0x4d53415544495431 and header[1] == 1, 'audit_format')
        require(tuple(header[2:4]) == (rank, ranks), 'audit_rank_identity')
        n, owned, blocks = [int(x) for x in header[4:]]
        require(0 <= owned <= n <= nodes and 0 <= blocks <= n*n, 'audit_dimensions')
        self.source = self.take('<i4', n)
        self.owned = self.take('<i4', owned)
        self.offsets = self.take('<i4', n+1)
        self.columns = self.take('<i4', blocks)
        require(np.all((self.source >= 0) & (self.source < nodes)) and len(np.unique(self.source)) == n, 'audit_source_ids')
        require(np.all((self.owned >= 0) & (self.owned < n)) and len(np.unique(self.owned)) == owned, 'audit_ownership')
        require(self.offsets[0] == 0 and self.offsets[-1] == blocks and np.all(np.diff(self.offsets) >= 0), 'audit_offsets')
        require(np.all((self.columns >= 0) & (self.columns < n)), 'audit_columns')
        self.momentum = self.take('<f8', 9*blocks).reshape(blocks, 3, 3)
        self.momentum_rhs = self.take('<f8', 3*n).reshape(n, 3)
        self.increment = self.take('<f8', 3*n).reshape(n, 3)
        self.predictor = self.take('<f8', 3*n).reshape(n, 3)
        self.influence = self.take('<f8', 3*n).reshape(n, 3)
        self.pressure = self.take('<f8', blocks).reshape(blocks, 1, 1)
        self.pressure_rhs = self.take('<f8', n).reshape(n, 1)
        self.phi = self.take('<f8', n).reshape(n, 1)
        self.gradient = self.take('<f8', 3*n).reshape(n, 3)
        require(self.position == len(self.data), 'audit_trailing_data')

    def take(self, dtype, count):
        size = np.dtype(dtype).itemsize*count
        require(self.position + size <= len(self.data), 'audit_truncated')
        result = np.ndarray((count,), dtype=dtype, buffer=self.data, offset=self.position)
        self.position += size
        return result

    def row(self, stage, local, component, source_order=True):
        components = 3 if stage == 'momentum' else 1
        start, end = self.offsets[local:local+2]
        columns = self.columns[start:end]
        if source_order:
            columns = self.source[columns]
        columns = (columns.astype(np.int64)[:, None]*components + np.arange(components)).ravel()
        return canonical_row(columns, getattr(self, stage)[start:end, component, :].ravel(),
                             getattr(self, stage + '_rhs')[local, component])


def canonical_row(columns, values, rhs):
    require(len(columns) == len(values) and np.all(np.isfinite(values)) and math.isfinite(float(rhs)), 'audit_nonfinite')
    keys, inverse = np.unique(columns, return_inverse=True)
    summed = np.zeros(len(keys))
    np.add.at(summed, inverse, values)
    require(np.all(np.isfinite(summed)), 'audit_nonfinite')
    keep = summed != 0
    require(keep.any(), 'audit_zero_row')
    return keys[keep], summed[keep], float(rhs)


class ReferenceMatrix:
    def __init__(self, directory, system, components, solver_to_source):
        paths = sorted(directory.glob(system + '_????_rows.bin'))
        require(len(paths) == 1, 'reference_matrix_count')
        prefix = str(paths[0])[:-len('_rows.bin')]
        count = components*len(solver_to_source)
        width, extra = divmod(paths[0].stat().st_size, count + 1)
        require(not extra and width in (4, 8), 'reference_matrix_index_width')
        self.offsets = np.memmap(str(paths[0]), dtype='<i' + str(width), mode='r')
        self.columns = np.memmap(prefix + '_cols.bin', dtype='<i' + str(width), mode='r')
        self.values = np.memmap(prefix + '_vals.bin', dtype='<f8', mode='r')
        self.rhs = np.memmap(prefix + '_b.bin', dtype='<f8', mode='r')
        require(len(self.rhs) == count and len(self.columns) == len(self.values), 'reference_matrix_dimensions')
        require(self.offsets[0] == 0 and self.offsets[-1] == len(self.values) and np.all(np.diff(self.offsets) >= 0), 'reference_matrix_offsets')
        require(np.all((self.columns >= 0) & (self.columns < count)), 'reference_matrix_columns')
        self.components = components
        self.mapping = solver_to_source
        self.inverse = np.argsort(solver_to_source)

    def row(self, source, component):
        row = self.components*self.inverse[source] + component
        start, end = self.offsets[row:row+2]
        columns = self.columns[start:end]
        mapped = self.mapping[columns//self.components]*self.components + columns % self.components
        return canonical_row(mapped, self.values[start:end], self.rhs[row])


def row_difference(a, b, scale):
    ac, av, ar = a
    bc, bv, br = b
    an, bn = np.max(np.abs(av)), np.max(np.abs(bv))
    keys = np.union1d(ac, bc)
    difference = np.zeros(len(keys))
    difference[np.searchsorted(keys, ac)] = av/an
    difference[np.searchsorted(keys, bc)] -= bv/bn
    # Positive row scaling is allowed; the unknown and its units stay unchanged.
    return float(np.max(np.abs(difference))), float(abs(ar/an-br/bn)/(scale+max(abs(ar/an), abs(br/bn))))


def residual(row, solution, scale):
    columns, values, rhs = row
    applied = values*solution[columns]
    error = float(rhs - math.fsum(applied))
    denominator = abs(rhs) + math.fsum(np.abs(applied)) + float(np.max(np.abs(values)))*scale
    require(math.isfinite(error) and math.isfinite(denominator) and denominator > 0, 'audit_residual_arithmetic')
    return error, rhs, float(abs(error)/denominator)


def compare_system(parts, reference, stage, mars_solution, reference_solution, scale):
    components = 3 if stage == 'momentum' else 1
    matrix_error = rhs_error = 0.
    norms = {name: [0., 0., 0.] for name in ('mars', 'reference', 'mars_local')}
    ghosts_equal = True
    for part in parts:
        local_solution = getattr(part, 'increment' if components == 3 else 'phi').ravel()
        for local in part.owned:
            start, end = part.offsets[local:local+2]
            referenced = part.columns[start:end]
            saved = local_solution.reshape(-1, components)[referenced]
            owners = mars_solution[part.source[referenced]]
            ghosts_equal = ghosts_equal and np.array_equal(saved, owners)
            source = part.source[local]
            for component in range(components):
                a, b = part.row(stage, local, component), reference.row(source, component)
                da, db = row_difference(a, b, scale)
                matrix_error, rhs_error = max(matrix_error, da), max(rhs_error, db)
                local_row = part.row(stage, local, component, source_order=False)
                for name, row, solution in (('mars', a, mars_solution), ('reference', b, reference_solution),
                                            ('mars_local', local_row, local_solution)):
                    r, rhs, backward = residual(row, solution.ravel(), scale)
                    values = norms[name]
                    values[0] = math.hypot(values[0], r)
                    values[1] = math.hypot(values[1], rhs)
                    values[2] = max(values[2], backward)
    result = dict(matrix_max_row_scaled=matrix_error, rhs_max_row_scaled=rhs_error,
                  mars_referenced_copies_equal_owners=ghosts_equal)
    for name, (r, b, backward) in norms.items():
        require(all(math.isfinite(x) for x in (r, b, backward)), 'audit_residual_arithmetic')
        result[name] = dict(absolute_residual=r, rhs_norm=b, relative_residual=r/b if b else None,
                            max_row_backward_error=backward)
    return result


def pressure_accuracy(pair, controls, system):
    from simple_startup_probe import read_json
    explicit = '--pressure-linear-rtol' in controls
    environment = read_json(pair / 'mars/launch.json')['environment']
    mars_rtol = float(controls.get('--pressure-linear-rtol', 1e-12))
    mars_atol = float(controls['--pressure-linear-atol'] if explicit else environment.get('MARS_HYPRE_ABSTOL', 0))
    runtime_rtol, runtime_atol = (mars_rtol, mars_atol) if explicit else (1e-10, 1e-13)
    solver = load_deck((pair / 'reference/input.i').read_bytes())['simulation']['solver']
    settings = solver['solver_control']['advanced_options']['linear_solver_settings']
    config = next((settings[key] for key in ('pressure_correction', 'segregated_flow', 'default') if key in settings), None)
    require(isinstance(config, dict), 'audit_pressure_controls')
    if 'lookup' in config:
        config = solver[config['lookup']]
    reference_rtol, reference_atol = float(config.get('rtol', 1e-6)), float(config.get('atol', 1e-16))
    require(all(math.isfinite(x) and x >= 0 for x in
                (mars_rtol, mars_atol, reference_rtol, reference_atol)), 'audit_pressure_controls')
    family = str(config.get('family', '')).lower()
    family = family if family in ('hypre', 'petsc', 'trilinos') else 'other'
    limits = dict(mars=max(mars_atol, mars_rtol*system['mars']['rhs_norm']),
                  reference=max(reference_atol, reference_rtol*system['reference']['rhs_norm']))
    scaled = runtime_rtol*system['mars']['rhs_norm']
    limits['runtime'] = max(runtime_atol, scaled) if explicit else runtime_atol + scaled
    require(all(math.isfinite(x) for x in limits.values()), 'audit_pressure_controls')
    checks = dict(mars_runtime_target_source='explicit_pressure_options' if explicit else 'default_true_residual_check',
                  same_declared_tolerances=mars_rtol == reference_rtol and mars_atol == reference_atol,
                  mars_referenced_copies_equal_owners=system['mars_referenced_copies_equal_owners'],
                  mars_owner_residual_meets_runtime_limit=system['mars']['absolute_residual'] <= limits['runtime'],
                  mars_local_residual_meets_runtime_limit=system['mars_local']['absolute_residual'] <= limits['runtime'],
                  mars_runtime_limit_looser_than_common_relative_check=limits['runtime'] > 1e-8*system['mars']['rhs_norm'],
                  reference_family=family, reference_backend_convergence_verified=False)
    for name in ('mars', 'reference'):
        checks[name + '_residual_below_declared_unpreconditioned_limit'] = system[name]['absolute_residual'] <= limits[name]
        checks[name + '_declared_limit_looser_than_common_relative_check'] = limits[name] > 1e-8*system[name]['rhs_norm']
    # PETSc can stop on a preconditioned norm. A deck tolerance alone cannot certify its verdict.
    for key in ('normalize_matrix', 'diagonal_scaling'):
        require(type(config.get(key, False)) is bool, 'audit_pressure_controls')
        checks['reference_' + key] = config.get(key, False)
    details = dict(mars_rtol=mars_rtol, mars_atol=mars_atol, reference_rtol=reference_rtol,
                   reference_atol=reference_atol, limits=limits)
    return checks, details


def collect(parts, name, nodes, components):
    result = np.empty((nodes, components))
    seen = np.zeros(nodes, dtype=bool)
    for part in parts:
        source = part.source[part.owned]
        require(not seen[source].any(), 'audit_duplicate_owner')
        result[source] = getattr(part, name)[part.owned]
        seen[source] = True
    require(seen.all() and np.isfinite(result).all(), 'audit_field_coverage')
    return result


def compare_first_step(pair, ids, xyz, tolerance, scales, ranks, reference_paths, public, detail_dir=None, gradient_audit=False, reference_pair=None):
    reference_pair = pair if reference_pair is None else reference_pair
    from simple_startup_probe import read_json, write_json
    nodes = len(ids)
    expected = [pair / ('mars/flow-audit-rank{:06d}.bin'.format(rank)) for rank in range(ranks)]
    require(sorted((pair / 'mars').glob('*.bin')) == expected, 'audit_parts')
    parts = [MarsPart(path, rank, ranks, nodes) for rank, path in enumerate(expected)]
    public['failed_check'] = 'reference_intermediate_fields'
    fields = ('aux', 'du_x', 'du_y', 'du_z', 'pressure_correction',
              'pressure_correction_gradient_x', 'pressure_correction_gradient_y', 'pressure_correction_gradient_z')
    values = reference_fields(reference_paths, ids, xyz, 1, tolerance, np.ones(len(fields)), fields, exact_storage=gradient_audit)
    numbering = values[:, 0]
    require(np.array_equal(np.sort(numbering), np.arange(nodes)), 'reference_solver_ids')
    solver_to_source = np.argsort(numbering)
    reference_final = reference_fields(reference_paths, ids, xyz, 1, tolerance, scales)
    initial = reference_fields(reference_paths, ids, xyz, 0, tolerance, scales)
    mars_initial, _ = mars_fields(pair / 'mars/flow-step-0', xyz, tolerance, ranks)
    mars_final, _ = mars_fields(pair / 'mars/flow-step-1', xyz, tolerance, ranks)
    require(public['initial_field_parity_verified'], 'audit_initial_fields_differ')
    influence, phi, gradient = values[:, 1:4], values[:, 4:5], values[:, 5:8]
    # In the supported segregated path, the last velocity change is u = u* - d grad(p').
    predictor = reference_final[:, :3] + influence*gradient
    reference_increment = predictor - initial[:, :3]
    mars = {name: collect(parts, name, nodes, components) for name, components in (
        ('increment', 3), ('predictor', 3), ('influence', 3), ('phi', 1), ('gradient', 3))}
    controls = options(read_json(pair / 'case.json')['arguments'])
    length, alpha = float(controls['--reference-length']), float(controls['--relax-p'])
    # Check that the reference fields actually obey the update used to reconstruct u*.
    require(np.max(np.abs(reference_final[:, 3:4] - initial[:, 3:4] - alpha*phi))/scales[3] < 1e-8,
            'reference_pressure_update')
    require(np.max(np.abs(mars['predictor'] - mars_initial[:, :3] - mars['increment']))/scales[0] < 1e-10,
            'mars_predictor_update')
    require(np.max(np.abs(mars_final[:, :3] - mars['predictor'] + mars['influence']*mars['gradient']))/scales[0] < 1e-10,
            'mars_velocity_update')
    require(np.max(np.abs(mars_final[:, 3:4] - mars_initial[:, 3:4] - alpha*mars['phi']))/scales[3] < 1e-10,
            'mars_pressure_update')
    public['failed_check'] = 'first_step_matrices'
    systems = {}
    for stage, c, system, mx, rx, scale in (
            ('momentum', 3, 'coupled_navier_stokes', mars['increment'], reference_increment, scales[0]),
            ('pressure', 1, 'pressure_correction', mars['phi'], phi, scales[3])):
        public['failed_check'] = 'first_step_' + stage + '_matrix'
        matrix = ReferenceMatrix(reference_pair / 'reference', system, c, solver_to_source)
        systems[stage] = compare_system(parts, matrix, stage, mx, rx, scale)
    errors = dict(momentum_predictor=float(np.max(np.abs(mars['predictor']-predictor))/scales[0]),
                  momentum_influence=float(np.max(np.abs(mars['influence']-influence))/(length*scales[0]/scales[3])),
                  pressure_increment=float(np.max(np.abs(mars['phi']-phi))/scales[3]),
                  pressure_increment_gradient=float(np.max(np.abs(mars['gradient']-gradient))/(scales[3]/length)),
                  corrected_velocity=float(np.max(np.abs(mars_final[:, :3]-reference_final[:, :3]))/scales[0]),
                  corrected_pressure=float(np.max(np.abs(mars_final[:, 3]-reference_final[:, 3]))/scales[3]))
    require(all(math.isfinite(x) for x in errors.values()), 'audit_field_arithmetic')
    matches = {name: value <= 1e-5 for name, value in errors.items()}
    accuracy = {}
    for stage, result in systems.items():
        matches[stage + '_matrix'] = result['matrix_max_row_scaled'] <= 1e-10
        matches[stage + '_rhs'] = result['rhs_max_row_scaled'] <= 1e-10
        for name in ('mars', 'reference', 'mars_local'):
            value = result[name]['relative_residual']
            accuracy[stage + '_' + name + '_relative_residual_below_1e_8'] = value <= 1e-8 if value is not None else None
            accuracy[stage + '_' + name + '_row_backward_error_below_1e_8'] = result[name]['max_row_backward_error'] <= 1e-8
    order = ('momentum_matrix', 'momentum_rhs', 'momentum_predictor', 'momentum_influence',
             'pressure_matrix', 'pressure_rhs', 'pressure_increment', 'pressure_increment_gradient',
             'corrected_velocity', 'corrected_pressure')
    public['failed_check'] = 'first_step_pressure_controls'
    checks, target_details = pressure_accuracy(pair, controls, systems['pressure'])
    if gradient_audit:
        public['failed_check'] = 'pressure_gradient_reconstruction'
        checks_gradient, details_gradient = gradient_reconstruction(
            read_json(pair / 'pair.json')['mesh'], xyz, mars['phi'], phi, mars['gradient'], gradient, scales[3]/length)
        public['pressure_gradient_checks'] = checks_gradient
        write_json((detail_dir or pair) / 'gradient-private.json', details_gradient)
    write_json((detail_dir or pair) / 'first-step-private.json', dict(systems=systems, field_errors=errors,
        pressure_targets=target_details,
        matrix_tolerance=1e-10, field_tolerance=1e-5, common_residual_threshold=1e-8,
        reference_predictor='reconstructed_from_final_velocity_and_correction', pressure_shift_applied=False))
    public.update(first_step_stage_matches=matches, first_step_linear_accuracy=accuracy,
                  pressure_solve_checks=checks,
                  pressure_target_scope='saved_controls_and_recomputed_residuals_not_backend_convergence_status',
                  first_differing_stage=next((name for name in order if not matches[name]), 'none'),
                  reference_predictor_reconstructed=True, matrix_comparison='positive_row_scaling_equivalence',
                  linear_accuracy_scope='common_diagnostic_threshold_not_backend_stopping_test')


def tet_gradient_action(xyz, chunks, fields):
    """Replay the shifted incremental scalar operator, including complete nodal stars."""
    require(xyz.ndim == 2 and xyz.shape[1] == 3 and fields.shape[0] == len(xyz)
            and fields.ndim == 2 and np.isfinite(xyz).all() and np.isfinite(fields).all(), 'gradient_mesh')
    volume = np.zeros(len(xyz))
    numerator = np.zeros((len(xyz), fields.shape[1], 3))
    magnitude = np.zeros_like(numerator)
    with np.errstate(over='raise', invalid='raise', divide='raise'):
        for nodes in chunks:
            require(nodes.ndim == 2 and nodes.shape[1] == 4 and nodes.dtype.kind in 'iu'
                    and np.all((nodes >= 0) & (nodes < len(xyz))), 'gradient_mesh')
            vertices = xyz[nodes]
            a, b, c = (vertices[:, i]-vertices[:, 0] for i in (1, 2, 3))
            cross = np.stack((np.cross(b, c), np.cross(c, a), np.cross(a, b)), axis=1)
            det = np.sum(a*cross[:, 0], axis=1)
            require(np.all(np.isfinite(det)) and np.all(det != 0), 'gradient_mesh')
            gradients = np.empty((len(nodes), 4, 3))
            gradients[:, 1:] = cross / det[:, None, None]
            gradients[:, 0] = -np.sum(gradients[:, 1:], axis=1)
            quarter = np.abs(det)/24
            values = fields[nodes]
            for local in range(4):
                np.add.at(volume, nodes[:, local], quarter)
            for left, right in ((0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)):
                area = quarter[:, None]*(gradients[:, right]-gradients[:, left])
                term = .5*(values[:, right]-values[:, left])[:, :, None]*area[:, None, :]
                # The captures interpolate before subtracting; account for cancellation in that form.
                size = (np.abs(values[:, right])+np.abs(values[:, left]))[:, :, None]*np.abs(area[:, None, :])
                for local in (left, right):
                    np.add.at(numerator, nodes[:, local], term)
                    np.add.at(magnitude, nodes[:, local], size)
        require(np.all(volume > 0) and np.isfinite(volume).all(), 'gradient_mesh')
        result, magnitude = numerator/volume[:, None, None], magnitude/volume[:, None, None]
    require(np.isfinite(result).all() and np.isfinite(magnitude).all(), 'gradient_arithmetic')
    return result, magnitude


def tet_chunks(ds, nodes):
    require('num_el_blk' in ds.dimensions and 'num_elem' in ds.dimensions, 'gradient_mesh')
    total = 0
    for number in range(1, len(ds.dimensions['num_el_blk'])+1):
        name = 'connect' + str(number)
        if name not in ds.variables:
            require('num_el_in_blk' + str(number) not in ds.dimensions, 'gradient_mesh')
            continue
        variable = ds.variables[name]
        require(variable.ndim == 2 and variable.shape[1] == 4 and variable.dtype.kind in 'iu'
                and str(variable.getncattr('elem_type')).strip().upper() in ('TETRA', 'TETRA4', 'TET4'), 'gradient_mesh')
        total += variable.shape[0]
        for start in range(0, variable.shape[0], 65536):
            block = variable[start:start+65536]
            require(not np.any(np.ma.getmaskarray(block)) and np.all((block > 0) & (block <= nodes)), 'gradient_mesh')
            yield np.asarray(block, dtype=np.int64)-1
    require(total == len(ds.dimensions['num_elem']) and total > 0, 'gradient_mesh')


def gradient_checks(replayed, magnitude, mars_gradient, reference_gradient, scale):
    require(math.isfinite(scale) and scale > 0 and np.isfinite(mars_gradient).all()
            and np.isfinite(reference_gradient).all(), 'gradient_arithmetic')
    # This is a strict consistency tolerance, not a certified floating-point error bound.
    limits = 1e-10*(scale + magnitude)
    difference = mars_gradient-reference_gradient
    errors = np.stack((mars_gradient-replayed[:, 0], reference_gradient-replayed[:, 1],
                       difference-replayed[:, 2]), axis=1)
    limits[:, 2] += limits[:, 0]+limits[:, 1]
    require(np.isfinite(errors).all() and np.isfinite(limits).all(), 'gradient_arithmetic')
    passed = np.all(np.abs(errors) <= limits, axis=(0, 2))
    peak = float(np.max(np.abs(difference)))
    resolved = float(np.max(limits)) < .01*peak
    within = bool(peak/scale <= 1e-5)
    assessment = ('reconstruction_mismatch' if not passed.all() else
                  'gradients_within_field_tolerance' if within else
                  'input_difference_explains_gradient_within_replay_tolerance' if resolved else
                  'insufficient_replay_resolution')
    public = dict(mars_saved_gradient_matches_reconstruction=bool(passed[0]),
                  reference_saved_gradient_matches_reconstruction=bool(passed[1]),
                  gradient_difference_matches_pressure_difference_action=bool(passed[2]),
                  replay_tolerance_resolves_observed_difference=resolved,
                  gradient_difference_within_field_tolerance=within, assessment=assessment,
                  reference_float64_storage_and_exact_copies_verified=True,
                  scope='common_shifted_tet4_operator_on_saved_fields', roundoff_bound_proven=False)
    private = dict(max_errors_scaled=(np.max(np.abs(errors), axis=(0, 2))/scale).tolist(),
                   max_limits_scaled=(np.max(limits, axis=(0, 2))/scale).tolist(),
                   observed_difference_scaled=peak/scale, consistency_tolerance=1e-10,
                   field_tolerance=1e-5, required_resolution_fraction=.01)
    return public, private


def gradient_reconstruction(mesh, xyz, mars_phi, reference_phi, mars_gradient, reference_gradient, scale):
    from netCDF4 import Dataset
    # Apply G directly to delta phi; subtracting two reconstructed gradients loses accuracy.
    fields = np.column_stack((mars_phi, reference_phi, mars_phi-reference_phi))
    with Dataset(mesh) as ds:
        replayed, magnitude = tet_gradient_action(xyz, tet_chunks(ds, len(xyz)), fields)
    public, private = gradient_checks(replayed, magnitude, mars_gradient, reference_gradient, scale)
    if not public['gradient_difference_within_field_tolerance']:
        selected, failing = select_gradient_nodes(replayed, mars_gradient, reference_gradient, scale)
        with Dataset(mesh) as ds:
            exact = exact_gradient_action(xyz, tet_chunks(ds, len(xyz)), mars_phi, reference_phi, selected)
        checks, details = exact_gradient_checks(exact, selected, failing, mars_gradient, reference_gradient)
        public['exact_selected_node_checks'] = checks
        private['exact_selected_node_details'] = details
    return public, private


def select_gradient_nodes(replayed, mars_gradient, reference_gradient, scale, per_score=8):
    difference = mars_gradient-reference_gradient
    delta = np.max(np.abs(difference), axis=1)
    candidates = np.flatnonzero(delta/scale > 1e-5)
    scores = (delta, np.max(np.abs(mars_gradient-replayed[:, 0]), axis=1),
              np.max(np.abs(reference_gradient-replayed[:, 1]), axis=1),
              np.max(np.abs(difference-replayed[:, 2]), axis=1))
    selected = set()
    for score in scores:
        # Source-row order breaks ties, independent of partition or previous selection.
        order = np.lexsort((candidates, -score[candidates]))[:per_score]
        selected.update(int(node) for node in candidates[order])
    return np.array(sorted(selected), dtype=np.int64), len(candidates)


def exact_gradient_action(xyz, chunks, mars_phi, reference_phi, selected, max_elements=20000):
    """Evaluate complete selected stars exactly on the stored binary64 inputs."""
    cells = []
    count = 0
    for chunk in chunks:
        relevant = chunk[np.any(np.isin(chunk, selected), axis=1)]
        count += len(relevant)
        if count > max_elements:
            return None  # Never truncate a star and call its reconstruction complete.
        cells.extend(relevant.tolist())
    points, fields = {}, {}
    for node in {node for cell in cells for node in cell}:
        points[node] = [Fraction(float(x)) for x in xyz[node]]
        fields[node] = [Fraction(float(mars_phi[node, 0])), Fraction(float(reference_phi[node, 0]))]
    zero = Fraction(0)
    volumes = {int(node): zero for node in selected}
    sums = {int(node): [[zero]*3 for _ in range(2)] for node in selected}
    def cross(a, b):
        return [a[(j+1)%3]*b[(j+2)%3]-a[(j+2)%3]*b[(j+1)%3] for j in range(3)]
    for cell in cells:
        x = [points[node] for node in cell]
        a, b, c = ([x[k][j]-x[0][j] for j in range(3)] for k in (1, 2, 3))
        cofactors = [None, cross(b, c), cross(c, a), cross(a, b)]
        cofactors[0] = [-sum(cofactors[k][j] for k in (1, 2, 3)) for j in range(3)]
        determinant = sum(a[j]*cofactors[1][j] for j in range(3))
        require(determinant != 0, 'gradient_mesh')
        sign = 1 if determinant > 0 else -1
        for local, node in enumerate(cell):
            if node not in sums:
                continue
            volumes[node] += abs(determinant)/24
            for other, neighbor in enumerate(cell):
                if local == other:
                    continue
                for field in range(2):
                    # Convert operands before subtraction: the saved doubles are exact rational inputs.
                    delta = fields[neighbor][field]-fields[node][field]
                    for j in range(3):
                        sums[node][field][j] += sign*delta*(cofactors[other][j]-cofactors[local][j])/48
    require(bool(volumes) and all(v > 0 for v in volumes.values()), 'gradient_mesh')
    return {node: [[value/volumes[node] for value in row] for row in sums[node]] for node in sums}


def exact_gradient_checks(exact, selected, failing, mars_gradient, reference_gradient):
    public = dict(scope='selected_field_tolerance_failures_complete_stars',
                  selection='up_to_eight_per_observed_difference_and_three_replay_errors',
                  all_field_tolerance_failures_checked=False,
                  global_equivalence_verified=False, runtime_roundoff_bound_proven=False,
                  assessment='inconclusive_work_limit', complete_selected_stars_verified=False)
    if exact is None:
        return public, dict(work_limit_reached=True)
    require(len(selected) > 0 and set(int(x) for x in selected) == set(exact), 'gradient_mesh')
    ratios = [Fraction(0)]*3
    rows = []
    for node in selected:
        gm = [Fraction(float(x)) for x in mars_gradient[node]]
        gr = [Fraction(float(x)) for x in reference_gradient[node]]
        delta = [gm[j]-gr[j] for j in range(3)]
        observed = max(abs(x) for x in delta)
        require(observed > 0, 'gradient_arithmetic')
        em = [gm[j]-exact[node][0][j] for j in range(3)]
        er = [gr[j]-exact[node][1][j] for j in range(3)]
        closure = [delta[j]-(exact[node][0][j]-exact[node][1][j]) for j in range(3)]
        local = [max(abs(x) for x in values)/observed for values in (em, er, closure)]
        ratios = [max(old, new) for old, new in zip(ratios, local)]
        rows.append(dict(source_node=int(node), error_fractions_of_local_difference=[str(x) for x in local]))
    passed = [x <= Fraction(1, 100) for x in ratios]
    public.update(complete_selected_stars_verified=True,
        all_field_tolerance_failures_checked=len(selected) == failing,
        arithmetic='exact_rational_on_stored_binary64_inputs',
        mars_reconstruction_within_one_percent_of_local_difference=passed[0],
        reference_reconstruction_within_one_percent_of_local_difference=passed[1],
        closure_within_one_percent_of_local_difference=passed[2],
        assessment=('input_difference_explains_checked_nodes' if all(passed) else
                    'captured_gradient_difference_not_explained_at_checked_nodes' if not passed[2] else
                    'individual_reconstruction_discrepancy_at_checked_nodes'))
    return public, dict(work_limit_reached=False, selected_nodes=rows,
                        max_error_fractions_of_local_difference=[str(x) for x in ratios],
                        required_fraction='1/100')
