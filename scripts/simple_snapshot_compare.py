#!/usr/bin/env python3
"""User-side comparison of existing steady SIMPLE snapshots; private inputs stay local."""

import csv
import glob
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import sys

from simple_convergence_summary import HEADER, COMPLETION, summarize
from simple_public_diagnostics import SafeParser


SCHEMA = 'mars-simple-snapshot-comparison-v1'
FIELDS = ('velocity_x', 'velocity_y', 'velocity_z', 'pressure')
CONTROL_KEYS = dict(rho='--rho', mu='--mu', inlet_speed='--inlet-velocity',
                    outlet_pressure='--outlet-pressure', reference_length='--reference-length',
                    alpha_u='--relax-u', alpha_p='--relax-p', alpha_mass='--relax-mass',
                    beta='--outlet-beta', pseudo_dt='--pseudo-dt')


class EvidenceError(ValueError):
    pass


def require(ok, reason=None):
    if not ok:
        if reason is not None:
            raise EvidenceError(reason)
        raise ValueError('comparison evidence check failed')


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def options(arguments):
    require(isinstance(arguments, list) and len(arguments) % 2 == 0)
    require(all(isinstance(x, str) for x in arguments))
    pairs = list(zip(arguments[::2], arguments[1::2]))
    require(all(k.startswith('--') for k, _ in pairs))
    result = dict(pairs)
    require(len(result) == len(pairs))
    return result


def controls(case_path, reference_dir, mesh, log):
    from prepare_simple_deck import load_deck, translate
    refinement = re.findall(r'(?:^|\s)pressure_refinement=([^\s]*)', '\n'.join(log))
    require(refinement in ([], ['0']) and not any('[simple-pressure-refinement]' in line for line in log),
            'pressure_refinement_not_comparable')
    case = json.loads(case_path.read_text())
    require(case['format'] == 'mars-simple-deck-v1')
    prepared = options(case['arguments'])
    require(mesh is None or Path(prepared['--mesh']).resolve() == mesh.resolve())
    require(prepared['--mesh-format'] == 'exodus')
    # The saved preparation hash selects the reference deck without printing its name.
    candidates = sorted(p for p in reference_dir.iterdir()
                        if p.is_file() and p.suffix.lower() in ('.yaml', '.yml', '.i'))
    exact = [p for p in candidates if digest(p) == case['deck_sha256']]
    deck, reference_status = None, 'unavailable'
    if exact or len(candidates) == 1:
        deck = (exact or candidates)[0]
        try:
            translated, _ = translate(load_deck(deck.read_bytes()), case.get('pressure_linear_policy', 'mars'))
            reference = options(translated)
            same = set(prepared) == set(reference) | {'--mesh', '--mesh-format', '--reference-length'}
            same = same and all(prepared[key] == value for key, value in reference.items())
            reference_status = 'mapped_controls_match' if same else 'mapped_controls_differ'
        except Exception:
            reference_status = 'unsupported_deck'
    elif candidates:
        reference_status = 'ambiguous_decks'
    headers = [HEADER.fullmatch(line.strip()) for line in log if line.startswith('SIMPLE Tet4,')]
    require(len(headers) == 1 and headers[0] is not None)
    require(headers[0].group(2) == prepared['--advection'])
    interpolation = [line.strip().split('=', 1)[1] for line in log
                     if line.startswith('velocity_interpolation=')]
    require(interpolation == [prepared['--velocity-interpolation']])
    for prefix, keys in (('rho=', list(CONTROL_KEYS)[:5]), ('alpha_u=', list(CONTROL_KEYS)[5:])):
        lines = [line for line in log if line.startswith(prefix)]
        require(len(lines) == 1)
        for key in keys:
            values = re.findall(r'(?:^|\s)' + key + r'=([^\s]+)', lines[0])
            require(len(values) == 1)
            actual, expected = float(values[0]), float(prepared[CONTROL_KEYS[key]])
            require(math.isfinite(actual) and math.isfinite(expected))
            require(math.isclose(actual, expected, rel_tol=1e-13, abs_tol=1e-15))
    target = [line for line in log if line.startswith('pressure_linear_rtol=')]
    if '--pressure-linear-rtol' in prepared:
        require(len(target) == 1)
        for label, flag in (('pressure_linear_rtol', '--pressure-linear-rtol'), ('pressure_linear_atol', '--pressure-linear-atol')):
            found = re.findall(label + r'=([^\s]+)', target[0])
            require(len(found) == 1 and float(found[0]) == float(prepared[flag]))
    else:
        require(not target)
    return prepared, deck, reference_status, bool(exact)


def result_files(reference_dir):
    paths = sorted(Path(p) for p in glob.glob(str(reference_dir / 'results.e*')))
    require(bool(paths) and all(p.is_file() for p in paths), 'reference_file_family')
    pieces = [re.fullmatch(r'(.+)\.(\d+)\.(\d+)', p.name) for p in paths]
    if not any(pieces):
        require(len(paths) == 1, 'reference_file_family')
        return paths
    require(all(pieces), 'reference_file_family')
    require(len({(m[1], int(m[2])) for m in pieces}) == 1, 'reference_file_family')
    count = int(pieces[0][2])
    require(count > 0 and len(paths) == count and {int(m[3]) for m in pieces} == set(range(count)), 'reference_file_family')
    return paths


def finite(data):
    import numpy as np
    require(not np.any(np.ma.getmaskarray(data)))
    result = np.asarray(data, dtype=np.float64)
    require(np.all(np.isfinite(result)))
    return result


def node_ids(ds, count, required=False):
    import numpy as np
    if 'node_num_map' not in ds.variables:
        require(not required)
        return np.arange(1, count + 1, dtype=np.int64)
    data = ds.variables['node_num_map'][:]
    require(not np.any(np.ma.getmaskarray(data)) and data.dtype.kind in 'iu')
    require(data.shape == (count,) and np.all(data > 0) and np.all(data <= np.iinfo(np.int64).max))
    result = np.asarray(data, dtype=np.int64)
    require(len(np.unique(result)) == count)
    return result


def coordinates(ds, start=0, stop=None):
    import numpy as np
    if 'coord' in ds.variables:
        variable = ds.variables['coord']
        require(variable.dimensions == ('num_dim', 'num_nodes') and variable.shape[0] == 3)
        return finite(variable[:, start:stop]).T
    require(all(ds.variables['coord' + axis].dimensions == ('num_nodes',) for axis in 'xyz'))
    return np.column_stack([finite(ds.variables['coord' + axis][start:stop]) for axis in 'xyz'])


def mars_fields(prefix, xyz, coordinate_tolerance, ranks):
    import numpy as np
    manifest, serial = Path(str(prefix) + '-fields.json'), Path(str(prefix) + '-fields.csv')
    require(manifest.is_file() != serial.is_file())
    paths = [serial]
    if manifest.is_file():
        record = json.loads(manifest.read_text())
        require(record['format'] == 'mars-simple-fields-v1' and type(record['nodes']) is int)
        require(record['nodes'] == len(xyz))
        parts = record['parts']
        require(isinstance(parts, list) and len(parts) == ranks)
        require(all(isinstance(p, str) and Path(p).name == p and p.endswith('.csv') for p in parts))
        require(len(parts) == len(set(parts)))
        paths = [manifest.parent / p for p in parts]
    values = np.empty((len(xyz), 4))
    seen = np.zeros(len(xyz), dtype=bool)
    for path in paths:
        with path.open(newline='') as stream:
            rows = csv.reader(stream)
            require(next(rows) == ['node', 'x', 'y', 'z', 'u', 'v', 'w', 'p'])
            for row in rows:
                require(len(row) == 8)
                node = int(row[0])
                require(0 <= node < len(xyz) and not seen[node])
                data = [float(v) for v in row[1:]]
                require(all(math.isfinite(v) for v in data))
                require(all(abs(data[k] - xyz[node, k]) <= coordinate_tolerance for k in range(3)))
                values[node] = data[3:]
                seen[node] = True
    require(seen.all())
    return values, ([manifest] if manifest.is_file() else []) + paths


def reference_fields(paths, ids, xyz, iteration, coordinate_tolerance, scales, fields=FIELDS, exact_storage=False):
    import numpy as np
    from netCDF4 import Dataset, chartostring
    order = np.argsort(ids)
    sorted_ids = ids[order]
    values = np.empty((len(ids), len(fields)))
    seen = np.zeros(len(ids), dtype=bool)
    expected_times = None
    for path in paths:
        with Dataset(str(path)) as ds:
            count = len(ds.dimensions['num_nodes'])
            mapped = node_ids(ds, count, required=len(paths) > 1)
            positions = np.searchsorted(sorted_ids, mapped)
            require(np.all(positions < len(ids)), 'reference_node_ids')
            require(np.array_equal(sorted_ids[positions], mapped), 'reference_node_ids')
            source = order[positions]
            times = finite(ds.variables['time_whole'][:])
            require(times.ndim == 1 and len(times) > 0 and np.all(np.diff(times) > 0))
            if expected_times is None:
                expected_times = times
            else:
                require(np.array_equal(times, expected_times), 'reference_saved_iteration')
            selected = np.flatnonzero(times == iteration)
            require(len(selected) == 1, 'reference_saved_iteration')
            index = int(selected[0])
            names = chartostring(np.ma.filled(ds.variables['name_nod_var'][:], b'\0'))
            names = [(s.decode() if isinstance(s, bytes) else str(s)).strip('\x00 ').lower() for s in names]
            require(all(names.count(f) == 1 for f in fields), 'reference_field_names')
            variables = []
            for name in fields:
                j = names.index(name)
                if 'vals_nod_var' in ds.variables:
                    var = ds.variables['vals_nod_var']
                    require(var.dimensions == ('time_step', 'num_nod_var', 'num_nodes'))
                    require(var.shape == (len(times), len(names), count))
                    variables.append((var, j))
                else:
                    var = ds.variables['vals_nod_var' + str(j + 1)]
                    require(var.dimensions == ('time_step', 'num_nodes') and var.shape == (len(times), count))
                    variables.append((var, None))
            if exact_storage:
                require(all(v.dtype.kind == 'f' and v.dtype.itemsize == 8 for v, _ in variables),
                        'gradient_storage_precision')
            for start in range(0, count, 65536):
                stop = min(start + 65536, count)
                nodes = source[start:stop]
                require(np.all(np.abs(coordinates(ds, start, stop) - xyz[nodes]) <= coordinate_tolerance), 'reference_coordinates')
                block = np.column_stack([finite(v[index, start:stop] if j is None else v[index, j, start:stop])
                                         for v, j in variables])
                duplicate = seen[nodes]
                if duplicate.any():
                    a, b = block[duplicate], values[nodes[duplicate]]
                    limit = 1e-12 * (scales + np.maximum(np.abs(a), np.abs(b)))
                    require(np.array_equal(a, b) if exact_storage else np.all(np.abs(a - b) <= limit),
                            'reference_ghost_values')
                values[nodes] = block
                seen[nodes] = True
    require(seen.all(), 'reference_node_coverage')
    return values


def norm3(data):
    import numpy as np
    return np.hypot(np.hypot(data[:, 0], data[:, 1]), data[:, 2])


def reference_output_values(deck):
    from prepare_simple_deck import load_deck
    if deck is None:
        return 'unknown'
    try:
        output = load_deck(deck.read_bytes())['simulation']['solver'].get('output_control', {})
        corrected = output.get('corrected_boundary_values', False)
        if type(corrected) is bool:
            return 'boundary_corrected' if corrected else 'solver_values'
    except Exception:
        pass
    return 'unknown'


def boundary_tags(ds, count):
    """Only classify stored tags; their completeness is not inferred from coordinates."""
    import numpy as np
    mask = np.zeros(count, dtype=bool)

    def indices(variable, limit):
        values = variable[:]
        require(not np.any(np.ma.getmaskarray(values)) and values.dtype.kind in 'iu')
        require(values.ndim == 1 and np.all(values > 0) and np.all(values <= limit))
        return np.asarray(values, dtype=np.int64) - 1

    sets = len(ds.dimensions['num_side_sets']) if 'num_side_sets' in ds.dimensions else 0
    nodesets = len(ds.dimensions['num_node_sets']) if 'num_node_sets' in ds.dimensions else 0
    if sets == 0 and nodesets == 0:
        return None
    for number in range(1, nodesets + 1):
        name = 'node_ns' + str(number)
        if name not in ds.variables and 'num_nod_ns' + str(number) not in ds.dimensions:
            continue  # Exodus permits empty sets.
        mask[indices(ds.variables[name], count)] = True
    if sets:
        total = len(ds.dimensions['num_elem'])
        selected, sides = [], []
        for number in range(1, sets + 1):
            suffix = str(number)
            if 'elem_ss' + suffix not in ds.variables and 'num_side_ss' + suffix not in ds.dimensions:
                continue
            elements = indices(ds.variables['elem_ss' + suffix], total)
            faces = indices(ds.variables['side_ss' + suffix], 4)
            require(elements.shape == faces.shape)
            selected.extend(elements)
            sides.extend(faces)
        selected, sides = np.asarray(selected, dtype=np.int64), np.asarray(sides, dtype=np.int64)
        # Exodus Tet4 side numbering, independent of global element IDs.
        face_nodes = np.array([[0, 1, 3], [1, 2, 3], [0, 3, 2], [0, 2, 1]])
        start = 0
        for number in range(1, len(ds.dimensions['num_el_blk']) + 1):
            name = 'connect' + str(number)
            if name not in ds.variables:
                require('num_el_in_blk' + str(number) not in ds.dimensions)
                continue
            connectivity = ds.variables[name]
            end = start + connectivity.shape[0]
            relevant = np.flatnonzero((selected >= start) & (selected < end))
            if len(relevant):
                topology = str(connectivity.getncattr('elem_type')).strip().upper()
                if topology not in ('TETRA', 'TETRA4', 'TET4') or connectivity.shape[1:] != (4,):
                    return None
                for first in range(0, len(relevant), 65536):
                    group = relevant[first:first + 65536]
                    values = connectivity[selected[group] - start, :]
                    require(not np.any(np.ma.getmaskarray(values)) and values.dtype.kind in 'iu')
                    require(np.all(values > 0) and np.all(values <= count))
                    nodes = np.asarray(values, dtype=np.int64)[np.arange(len(group))[:, None], face_nodes[sides[group]]]
                    mask[nodes.ravel() - 1] = True
            start = end
        require(start == total)
    return mask


def field_errors(difference):
    import numpy as np
    velocity = norm3(difference[:, :3])
    pressure = difference[:, 3]
    return dict(velocity_max_scaled=float(np.max(velocity)), velocity_rms_scaled=rms(velocity),
                pressure_max_scaled=float(np.max(np.abs(pressure))), pressure_rms_scaled=rms(pressure))


def rms(data):
    import numpy as np
    peak = float(np.max(np.abs(data)))
    return peak * math.sqrt(float(np.mean((data / peak)**2))) if peak else 0.


def band(error):
    for limit, label in ((1e-5, 'within_1e_minus_5'), (.01, 'within_1_percent'), (.05, 'within_5_percent')):
        if error <= limit:
            return label
    return 'over_5_percent'


def compare(args, public):
    import numpy as np
    from netCDF4 import Dataset
    public['failed_check'] = 'mars_completion'
    log = args.mars_log.read_text(errors='replace').splitlines()
    metrics_path = Path(str(args.mars_prefix) + '-metrics.csv')
    metrics_text = metrics_path.read_text()
    convergence = summarize(log, io.StringIO(metrics_text), args.mars_exit.read_text(),
                            args.residual_tol, args.mass_tol, args.change_tol)
    require(convergence['completion_matches_supplied_targets'])
    ending = [COMPLETION.fullmatch(line.strip()) for line in log if line.startswith(('CONVERGED', 'NOT CONVERGED'))]
    require(len(ending) == 1 and ending[0] is not None and int(ending[0].group(2)) == args.iteration)
    public['mars_status'] = convergence['run_status']
    public['failed_check'] = 'settings'
    prepared, deck, reference_status, deck_hash_matches = controls(args.case, args.reference_dir, args.mesh, log)
    if args.mesh is None:
        args.mesh = Path(prepared['--mesh'])
    public['logged_mars_controls_match_preparation'] = True
    public['reference_settings_status'] = reference_status
    public['reference_deck_hash_matches_preparation'] = deck_hash_matches
    output_values = reference_output_values(deck)
    public['reference_deck_output_values'] = output_values
    public['solver_field_comparison_supported'] = deck_hash_matches and output_values == 'solver_values'
    u, rho = float(prepared['--inlet-velocity']), float(prepared['--rho'])
    pressure_scale = rho * u * u
    require(all(math.isfinite(v) and v > 0 for v in (u, rho, pressure_scale)))
    scales = np.array([u, u, u, pressure_scale])
    public['failed_check'] = 'source_nodes'
    with Dataset(str(args.mesh)) as ds:
        count = len(ds.dimensions['num_nodes'])
        require(count > 0)
        ids, xyz = node_ids(ds, count), coordinates(ds)
        public['failed_check'] = 'source_boundary_tags'
        boundary = boundary_tags(ds, count)
    require(xyz.shape == (count, 3))
    coordinate_tolerance = 64 * np.finfo(float).eps * max(1., float(np.max(np.abs(xyz))))
    public['failed_check'] = 'mars_fields'
    mars, mars_paths = mars_fields(args.mars_prefix, xyz, coordinate_tolerance, int(ending[0].group(3)))
    metrics = list(csv.DictReader(io.StringIO(metrics_text)))
    peak = float(np.max(norm3(mars[:, :3])))
    require(math.isclose(peak, float(metrics[-1]['umax_m_s']), rel_tol=1e-12, abs_tol=1e-12*u))
    public['failed_check'] = 'reference_fields'
    reference_paths = result_files(args.reference_dir)
    reference = reference_fields(reference_paths, ids, xyz, args.iteration, coordinate_tolerance, scales)
    public['saved_iteration_and_node_mapping_match'] = True
    public['failed_check'] = 'field_arithmetic'
    difference = finite((mars - reference) / scales)
    reference_peak = float(np.max(norm3(reference[:, :3])))
    errors = field_errors(difference)
    regional_errors = {}
    public['boundary_localization_status'] = 'unavailable'
    if boundary is not None:
        public['boundary_localization_status'] = 'source_tags'
        for label, selected in (('tagged_boundary', boundary), ('other_nodes', ~boundary)):
            if np.any(selected):
                regional_errors[label] = field_errors(difference[selected])
                public.update((label + '_' + key + '_band', band(value))
                              for key, value in regional_errors[label].items())
            else:
                public[label + '_status'] = 'empty'
    require(all(math.isfinite(v) for v in errors.values()) and math.isfinite(reference_peak))
    # A peak alone does not establish agreement of the vector field.
    peak_relative = abs(peak - reference_peak) / reference_peak if reference_peak else None
    if peak_relative is not None:
        require(math.isfinite(peak_relative))
    public['failed_check'] = 'private_report'
    paths = [args.mesh, args.case, args.mars_log, args.mars_exit, metrics_path] + mars_paths + reference_paths
    if deck is not None:
        paths.append(deck)
    private = dict(schema=SCHEMA, iteration=args.iteration, density=rho, velocity_scale=u,
                   pressure_scale=pressure_scale, errors=errors, regional_errors=regional_errors,
                   reference_deck_output_values=output_values, mars_peak_speed=peak,
                   reference_peak_speed=reference_peak, peak_relative_difference=peak_relative,
                   pressure_mean_difference=float(np.mean(mars[:, 3] - reference[:, 3])),
                   convergence=convergence, reference_settings_status=reference_status,
                   sha256={str(p.resolve()): digest(p) for p in paths},
                   limits=['User-selected existing files, not a runtime-attested run bundle.',
                           'Checks nodal identity and coordinates, not element connectivity or side sets.',
                           'Printed MARS controls checked against preparation; reference deck availability is reported separately.',
                           'Mapped controls are only the supported deck subset; actual reference launch not attested.',
                           'Linear solvers and tolerances need not match; iteration counts are not physical time.',
                           'The saved deck declares output corrections; the actual reference launch is not attested.',
                           'Regional errors use all source side/node sets; their completeness and physical roles are not verified.',
                           'No pressure shift removed. Nodal RMS is unweighted, not a volume norm.',
                           'A field mismatch at a finite iteration does not establish a discretization defect.'])
    with args.private_report.open('x') as stream:
        json.dump(private, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')
    public.update((key + '_band', band(value)) for key, value in errors.items())
    public['peak_speed_relative_band'] = band(peak_relative) if peak_relative is not None else 'zero_reference_peak'
    public['snapshot_fields_within_tolerance'] = (public['solver_field_comparison_supported']
                                               and errors['velocity_max_scaled'] <= 1e-5
                                               and errors['pressure_max_scaled'] <= 1e-5)
    public.update(comparison_status='completed', failed_check='none')


def main(argv=None):
    parser = SafeParser(description=__doc__)
    for name in ('reference-dir', 'case', 'mars-prefix', 'mars-log', 'mars-exit', 'private-report', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--mesh', type=Path, help='Defaults to the mesh path recorded in case.json')
    parser.add_argument('--iteration', type=int, required=True)
    for name in ('residual-tol', 'mass-tol', 'change-tol'):
        parser.add_argument('--' + name, type=float, default=1e-6)
    args = parser.parse_args(argv)
    os.umask(0o077)
    result = dict(schema=SCHEMA, comparison_status='invalid_evidence', failed_check='dependencies',
                  logged_mars_controls_match_preparation=False, reference_settings_status='unavailable',
                  reference_deck_hash_matches_preparation=False, saved_iteration_and_node_mapping_match=False,
                  reference_deck_output_values='unknown', solver_field_comparison_supported=False,
                  boundary_localization_status='unavailable',
                  full_run_provenance_verified=False, identical_linear_solvers_verified=False,
                  nonlinear_convergence_required=False, mars_status='unknown')
    try:
        # Claim the public output first; never overwrite an existing report or an input.
        with args.output.open('x') as stream:
            try:
                require(args.iteration > 0 and not args.private_report.exists())
                compare(args, result)
            except EvidenceError as error:
                result['failed_check'] = str(error)  # Only fixed labels originate from require().
            except Exception:
                pass  # Exception messages may contain private paths or field values.
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')
    except Exception:
        print('ERROR: cannot create a fresh public summary; inspect paths privately.', file=sys.stderr)
        return 1
    print('Summary written. Share only the public JSON; the detailed report stays private.')
    return 0 if result['comparison_status'] == 'completed' else 1


if __name__ == '__main__':
    sys.exit(main())
