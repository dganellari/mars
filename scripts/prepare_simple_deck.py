#!/usr/bin/env python3
"""Prepare native SIMPLE arguments from a restricted OpenAccel deck; never read the mesh."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys


class Unsupported(ValueError):
    pass


def require(ok, path, reason):
    if not ok:
        raise Unsupported(path + ': ' + reason)


def mapping(value, allowed, path):
    require(isinstance(value, dict), path, 'expected a mapping')
    require(set(value) <= set(allowed.split()), path, 'contains unhandled keys; no arguments emitted')
    return value


def number(value, path, positive=False, unit=False):
    require(not isinstance(value, bool) and isinstance(value, (int, float, str)), path, 'expected a scalar number')
    try:
        result = float(value)
    except (ValueError, OverflowError):
        raise Unsupported(path + ': expected a scalar number')
    require(math.isfinite(result), path, 'expected a finite number')
    require(not positive or result > 0, path, 'must be positive')
    require(not unit or 0 < result <= 1, path, 'must lie in (0, 1]')
    return result


def equal(value, expected, path):
    require(type(value) is type(expected) and value == expected, path, 'unsupported setting')


def singleton(value, path):
    require(isinstance(value, list) and len(value) == 1, path, 'requires exactly one entry')
    return value[0]


def names(value, path):
    require(isinstance(value, list) and value, path, 'expected side-set names')
    require(all(isinstance(v, str) and v and ',' not in v and not any(ord(c) < 32 for c in v)
                for v in value), path, 'invalid side-set name')
    require(len(set(value)) == len(value), path, 'duplicate side-set names')
    return value


def zero_field(value, field, vector=False):
    path = 'initialization.' + field
    value = mapping(value, 'option ' + field, path)
    equal(value.get('option'), 'value', path + '.option')
    items = value.get(field)
    if vector:
        require(isinstance(items, list) and len(items) == 3, path, 'requires three zero components')
    else:
        items = [items]
    require(all(number(v, path) == 0 for v in items), path, 'only zero initial fields are supported')


def translate(doc):
    doc = mapping(doc, 'mesh simulation', 'document')
    mesh = mapping(doc.get('mesh'), 'file_path automatic_decomposition_type', 'mesh')
    # This selects OpenAccel's partitioner; MARS partitions through Cornerstone.
    if 'automatic_decomposition_type' in mesh:
        require(isinstance(mesh['automatic_decomposition_type'], str),
                'mesh.automatic_decomposition_type', 'expected a string')
    mesh_path = mesh.get('file_path')
    require(isinstance(mesh_path, str) and mesh_path, 'mesh.file_path', 'expected a path')
    sim = mapping(doc.get('simulation'), 'verbose physical_analysis solver material_library', 'simulation')
    physical = mapping(sim.get('physical_analysis'), 'analysis_type domains initialization', 'physical_analysis')
    analysis = mapping(physical.get('analysis_type'), 'option', 'analysis_type')
    equal(analysis.get('option'), 'steady_state', 'analysis_type.option')
    domain = mapping(singleton(physical.get('domains'), 'domains'),
                     'name location materials type domain_models fluid_models boundaries initialization', 'domain')
    equal(domain.get('type'), 'fluid', 'domain.type')
    singleton(names(domain.get('location'), 'domain.location'), 'domain.location')
    model = mapping(domain.get('domain_models', {}), 'reference_pressure', 'domain_models')
    require(number(model.get('reference_pressure', 0), 'reference_pressure') == 0,
            'reference_pressure', 'only zero reference pressure is supported')
    fluid = mapping(domain.get('fluid_models'), 'turbulence', 'fluid_models')
    turbulence = mapping(fluid.get('turbulence'), 'option', 'turbulence')
    equal(turbulence.get('option'), 'laminar', 'turbulence.option')
    init = dict(mapping(physical.get('initialization', {}), 'velocity pressure', 'global initialization'))
    init.update(mapping(domain.get('initialization', {}), 'velocity pressure', 'domain initialization'))
    zero_field(init.get('velocity'), 'velocity', True)
    zero_field(init.get('pressure'), 'pressure')

    material_name = singleton(domain.get('materials'), 'domain.materials')
    material = mapping(singleton(sim.get('material_library'), 'material_library'),
                       'name thermodynamic_properties transport_properties', 'material')
    require(material.get('name') == material_name, 'material', 'domain material does not match library')
    thermo = mapping(material.get('thermodynamic_properties'), 'equation_of_state', 'thermodynamic_properties')
    eos = mapping(thermo.get('equation_of_state'), 'option density', 'equation_of_state')
    equal(eos.get('option'), 'value', 'equation_of_state.option')
    transport = mapping(material.get('transport_properties'), 'dynamic_viscosity', 'transport_properties')
    viscosity = mapping(transport.get('dynamic_viscosity'), 'option dynamic_viscosity', 'dynamic_viscosity')
    equal(viscosity.get('option'), 'value', 'dynamic_viscosity.option')
    args = ['--rho', str(number(eos.get('density'), 'density', positive=True)),
            '--mu', str(number(viscosity.get('dynamic_viscosity'), 'dynamic_viscosity', positive=True))]

    boundaries = domain.get('boundaries')
    require(isinstance(boundaries, list) and boundaries, 'boundaries', 'missing boundaries')
    selected, wall_sets, kinds = set(), [], []
    for boundary in boundaries:
        b = mapping(boundary, 'name type location boundary_details', 'boundary')
        locations = names(b.get('location'), 'boundary.location')
        require(not selected.intersection(locations), 'boundaries', 'duplicate side-set assignment')
        selected.update(locations)
        details = mapping(b.get('boundary_details', {}), 'mass_and_momentum', 'boundary_details')
        kind = b.get('type')
        if kind == 'wall':
            mm = mapping(details.get('mass_and_momentum', {}), 'option', 'wall.mass_and_momentum')
            equal(mm.get('option', 'no_slip_wall'), 'no_slip_wall', 'wall.option')
            wall_sets.extend(locations)
            continue
        require(kind in ('inlet', 'outlet') and kind not in kinds, 'boundaries', 'requires one inlet and one outlet')
        kinds.append(kind)
        args += ['--' + kind + '-ss', singleton(locations, kind + '.location')]
        mm = mapping(details.get('mass_and_momentum'),
                     'option normal_speed' if kind == 'inlet' else 'option relative_pressure pressure_profile_blend',
                     kind + '.mass_and_momentum')
        if kind == 'inlet':
            equal(mm.get('option'), 'normal_speed', 'inlet.option')
            args += ['--inlet-velocity', str(number(mm.get('normal_speed'), 'normal_speed', positive=True))]
        else:
            option = mm.get('option')
            require(option in ('static_pressure', 'average_static_pressure'), 'outlet.option', 'unsupported setting')
            # beta=1 fixes the trace to the prescribed pressure. OpenAccel only
            # reads pressure_profile_blend for average_static_pressure.
            beta = (1.0 if option == 'static_pressure' else
                    number(mm.get('pressure_profile_blend'), 'pressure_profile_blend', unit=True))
            args += ['--outlet-pressure', str(number(mm.get('relative_pressure'), 'relative_pressure')),
                     '--outlet-beta', str(beta)]
    require(set(kinds) == {'inlet', 'outlet'} and wall_sets, 'boundaries', 'requires inlet, outlet and walls')
    args += ['--wall-ss', ','.join(wall_sets)]

    solver = mapping(sim.get('solver'), 'solver_control output_control', 'solver')
    control = mapping(solver.get('solver_control'), 'basic_settings advanced_options expert_parameters', 'solver_control')
    basic = mapping(control.get('basic_settings'),
                    'advection_scheme interpolation_scheme convergence_controls convergence_criteria', 'basic_settings')
    scheme = basic.get('advection_scheme', 'upwind')
    require(scheme in ('upwind', 'high_resolution'), 'advection_scheme', 'unsupported scheme')
    args += ['--advection', scheme.replace('_', '-')]
    interp = mapping(basic.get('interpolation_scheme', {}),
                     'velocity_interpolation_type pressure_interpolation_type velocity_gradient_interpolation_type pressure_gradient_interpolation_type', 'interpolation_scheme')
    for key in ('velocity_interpolation_type', 'pressure_interpolation_type',
                'velocity_gradient_interpolation_type', 'pressure_gradient_interpolation_type'):
        expected = 'trilinear' if key == 'velocity_interpolation_type' else 'linear_linear'
        equal(interp.get(key, expected), expected, 'interpolation_scheme.' + key)
    conv = mapping(basic.get('convergence_controls'),
                   'min_iterations max_iterations physical_timescale relaxation_parameters', 'convergence_controls')
    args += ['--pseudo-dt', str(number(conv.get('physical_timescale'), 'physical_timescale', positive=True))]
    relax = mapping(conv.get('relaxation_parameters', {}),
                    'velocity_relaxation_factor pressure_relaxation_factor relax_mass', 'relaxation_parameters')
    for key, flag in [('velocity_relaxation_factor', '--relax-u'),
                      ('pressure_relaxation_factor', '--relax-p'), ('relax_mass', '--relax-mass')]:
        args += [flag, str(number(relax.get(key, 1), key, unit=True))]
    mapping(basic.get('convergence_criteria', {}), 'residual_type residual_target', 'convergence_criteria')
    advanced = mapping(control.get('advanced_options', {}), 'equation_controls linear_solver_settings', 'advanced_options')
    eq = mapping(advanced.get('equation_controls', {}), 'sub_iterations', 'equation_controls')
    sub = mapping(eq.get('sub_iterations', {}), 'pressure_correction segregated_flow', 'sub_iterations')
    for key in ('pressure_correction', 'segregated_flow'):
        require(number(sub.get(key, 1), key) == 1, 'sub_iterations.' + key, 'requires one subiteration')
    defaults = dict(consistent=False, fractional_step_method=False, incremental_gradient_change=True,
                    limit_gradients=False, relax_gradients=True, correct_gradients=False,
                    false_mass_accumulation=True, nonlinear_stabilisation=False,
                    disable_momentum_predictor=False, high_speed_blend_damping=False)
    expert = mapping(control.get('expert_parameters', {}), ' '.join(defaults) + ' blend_factor_max', 'expert_parameters')
    for key, default in defaults.items():
        expected = False if key == 'relax_gradients' else default
        equal(expert.get(key, default), expected, 'expert_parameters.' + key)
    cap = number(expert.get('blend_factor_max', 1), 'blend_factor_max')
    require(cap == 1 if scheme == 'high_resolution' else cap in (0, 1),
            'blend_factor_max', 'unsupported limiter cap')
    return args, mesh_path


def load_deck(data):
    try:
        import yaml
    except ImportError:
        raise Unsupported('dependency: PyYAML is required in the preparation Python environment')

    class UniqueLoader(yaml.SafeLoader):
        pass

    def unique(loader, node, deep=False):
        # Reject merges too: silently overriding a physical control is unsafe here.
        result = {}
        for key_node, value_node in node.value:
            key = loader.construct_object(key_node, deep=deep)
            require(isinstance(key, str) and key not in result, 'YAML mapping', 'non-string or duplicate key')
            result[key] = loader.construct_object(value_node, deep=deep)
        return result

    UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, unique)
    try:
        return yaml.load(data, Loader=UniqueLoader)
    except yaml.YAMLError:
        raise Unsupported('YAML: invalid or unsupported syntax')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--deck', type=Path, required=True)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference-length', type=float, default=1,
                        help='MARS residual scale in metres; does not change physics')
    o = parser.parse_args()
    require(not o.output.exists(), 'output', 'choose a fresh directory')
    length = number(o.reference_length, 'reference-length', positive=True)
    data = o.deck.read_bytes()
    args, declared_mesh = translate(load_deck(data))
    expanded = os.path.expandvars(os.path.expanduser(declared_mesh))
    require('$' not in expanded, 'mesh.file_path', 'unresolved environment variable')
    reference_mesh = Path(expanded)
    if not reference_mesh.is_absolute():
        reference_mesh = o.deck.resolve().parent / reference_mesh
    require(reference_mesh.resolve() == o.mesh.resolve(), 'mesh.file_path', 'supplied mesh differs from the deck path')
    require(o.mesh.is_file(), 'mesh', 'file does not exist')
    args = ['--mesh', str(o.mesh.resolve()), '--mesh-format', 'exodus'] + args + ['--reference-length', str(length)]
    o.output.mkdir(parents=True, mode=0o700)
    record = dict(format='mars-simple-deck-v1', deck_sha256=hashlib.sha256(data).hexdigest(),
                  arguments=args, status='supported_deck_mesh_not_validated',
                  notes=['Native C++ setup must validate single-block Tet4 topology and boundary coverage.',
                         'OpenAccel automatic_decomposition_type is not translated; MARS uses Cornerstone.',
                         'Constant static_pressure uses outlet beta=1; pressure_profile_blend applies only to average_static_pressure.',
                         'Coordinates must be in metres; mesh bytes were not read.',
                         'MARS uses its own Hypre settings, zero initial fields and convergence norms.',
                         'Iteration counts and reference residual criteria are not translated.'])
    (o.output / 'case.json').write_text(json.dumps(record, indent=2) + '\n')
    (o.output / 'args.nul').write_bytes(b'\0'.join(v.encode('utf-8') for v in args) + b'\0')
    print('PASS: supported deck; native mesh setup still required. Arguments saved locally.')


if __name__ == '__main__':
    try:
        main()
    except (Unsupported, OSError, UnicodeError, RecursionError) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        sys.exit(1)
