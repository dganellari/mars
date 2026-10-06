"""User-local first-step algebra audit. Numerical data never enters the public report."""
import math

import numpy as np

from simple_snapshot_compare import mars_fields, options, reference_fields, require

AUDIT_ERRORS = frozenset('audit_parts audit_format audit_rank_identity audit_dimensions audit_source_ids '
    'audit_ownership audit_offsets audit_columns audit_trailing_data audit_truncated audit_nonfinite audit_zero_row '
    'reference_matrix_count reference_matrix_index_width reference_matrix_dimensions reference_matrix_offsets '
    'reference_matrix_columns audit_residual_arithmetic audit_duplicate_owner audit_field_coverage '
    'reference_solver_ids audit_initial_fields_differ reference_pressure_update mars_predictor_update '
    'mars_velocity_update mars_pressure_update audit_field_arithmetic reference_field_names '
    'reference_ghost_values reference_coordinates reference_saved_iteration reference_node_ids reference_node_coverage'.split())


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

    def row(self, stage, local, component):
        components = 3 if stage == 'momentum' else 1
        start, end = self.offsets[local:local+2]
        columns = (self.source[self.columns[start:end]].astype(np.int64)[:, None]*components + np.arange(components)).ravel()
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
    norms = {name: [0., 0., 0.] for name in ('mars', 'reference')}
    for part in parts:
        for local in part.owned:
            source = part.source[local]
            for component in range(components):
                a, b = part.row(stage, local, component), reference.row(source, component)
                da, db = row_difference(a, b, scale)
                matrix_error, rhs_error = max(matrix_error, da), max(rhs_error, db)
                for name, row, solution in (('mars', a, mars_solution), ('reference', b, reference_solution)):
                    r, rhs, backward = residual(row, solution.ravel(), scale)
                    values = norms[name]
                    values[0] = math.hypot(values[0], r)
                    values[1] = math.hypot(values[1], rhs)
                    values[2] = max(values[2], backward)
    result = dict(matrix_max_row_scaled=matrix_error, rhs_max_row_scaled=rhs_error)
    for name, (r, b, backward) in norms.items():
        require(all(math.isfinite(x) for x in (r, b, backward)), 'audit_residual_arithmetic')
        result[name] = dict(absolute_residual=r, rhs_norm=b, relative_residual=r/b if b else None,
                            max_row_backward_error=backward)
    return result


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


def compare_first_step(pair, ids, xyz, tolerance, scales, ranks, reference_paths, public):
    from simple_startup_probe import read_json, write_json
    nodes = len(ids)
    expected = [pair / ('mars/flow-audit-rank{:06d}.bin'.format(rank)) for rank in range(ranks)]
    require(sorted((pair / 'mars').glob('*.bin')) == expected, 'audit_parts')
    parts = [MarsPart(path, rank, ranks, nodes) for rank, path in enumerate(expected)]
    public['failed_check'] = 'reference_intermediate_fields'
    fields = ('aux', 'du_x', 'du_y', 'du_z', 'pressure_correction',
              'pressure_correction_gradient_x', 'pressure_correction_gradient_y', 'pressure_correction_gradient_z')
    values = reference_fields(reference_paths, ids, xyz, 1, tolerance, np.ones(len(fields)), fields)
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
        matrix = ReferenceMatrix(pair / 'reference', system, c, solver_to_source)
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
        for name in ('mars', 'reference'):
            value = result[name]['relative_residual']
            accuracy[stage + '_' + name + '_relative_residual_below_1e_8'] = value <= 1e-8 if value is not None else None
            accuracy[stage + '_' + name + '_row_backward_error_below_1e_8'] = result[name]['max_row_backward_error'] <= 1e-8
    order = ('momentum_matrix', 'momentum_rhs', 'momentum_predictor', 'momentum_influence',
             'pressure_matrix', 'pressure_rhs', 'pressure_increment', 'pressure_increment_gradient',
             'corrected_velocity', 'corrected_pressure')
    write_json(pair / 'first-step-private.json', dict(systems=systems, field_errors=errors,
        matrix_tolerance=1e-10, field_tolerance=1e-5, common_residual_threshold=1e-8,
        reference_predictor='reconstructed_from_final_velocity_and_correction', pressure_shift_applied=False))
    public.update(first_step_stage_matches=matches, first_step_linear_accuracy=accuracy,
                  first_differing_stage=next((name for name in order if not matches[name]), 'none'),
                  reference_predictor_reconstructed=True, matrix_comparison='positive_row_scaling_equivalence',
                  linear_accuracy_scope='common_diagnostic_threshold_not_backend_stopping_test')
