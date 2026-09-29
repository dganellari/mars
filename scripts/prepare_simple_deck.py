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


class Checks:
    def __init__(self):
        self.errors = []

    def require(self, ok, path, reason):
        if not ok:
            message = path + ': ' + reason
            if message not in self.errors:
                self.errors.append(message)
        return ok

    def checked(self, function, *args, fallback=None, **kwargs):
        try:
            return function(*args, **kwargs)
        except Unsupported as error:
            if str(error) not in self.errors:
                self.errors.append(str(error))
            return fallback

    def mapping(self, value, allowed, path):
        if not self.require(isinstance(value, dict), path, 'expected a mapping'):
            return {}
        allowed = set(allowed.split() if isinstance(allowed, str) else allowed)
        self.require(set(value) <= allowed, path, 'unrecognized keys (names withheld); no arguments emitted')
        return value

    def number(self, value, path, **kwargs):
        return self.checked(number, value, path, fallback=0, **kwargs)

    def equal(self, value, expected, path):
        self.checked(equal, value, expected, path)

    def singleton(self, value, path):
        return self.checked(singleton, value, path)

    def names(self, value, path):
        return self.checked(names, value, path, fallback=[])

    def finish(self):
        # Fallbacks let independent checks continue, never produce runnable input.
        if self.errors:
            raise Unsupported('\n'.join(self.errors))


LINEAR_FAMILIES = {'petsc', 'hypre', 'trilinos', 'amgsolver', 'gmres'}
SOLVER_CONTROLS = {'solver_control', 'output_control', 'restart_control'}


def is_linear_solver(value):
    return (isinstance(value, dict) and isinstance(value.get('family'), str)
            and value['family'].lower() in LINEAR_FAMILIES)


def solver_mapping(value, check):
    if not isinstance(value, dict):
        return check.mapping(value, SOLVER_CONTROLS, 'solver')
    # OpenAccel linearSystem::setupSolver resolves lookup names in this mapping.
    named = {key for key, config in value.items()
             if key not in SOLVER_CONTROLS and is_linear_solver(config)}
    result = check.mapping(value, SOLVER_CONTROLS | named, 'solver')
    check.require('restart_control' not in result, 'solver.restart_control',
                  'loads saved fields; native SIMPLE currently requires a fresh start')
    return result


def check_linear_settings(value, solver, check):
    if not check.require(isinstance(value, dict), 'linear_solver_settings', 'expected a mapping'):
        return
    for config in value.values():
        path = 'linear_solver_settings.entry'
        if not check.require(isinstance(config, dict), path, 'expected a mapping'):
            continue
        if 'lookup' in config:
            name = config['lookup']
            if not check.require(isinstance(name, str) and bool(name), path + '.lookup', 'expected a name'):
                continue
            check.require(name not in SOLVER_CONTROLS and name in solver and is_linear_solver(solver[name]),
                          path + '.lookup', 'must reference a named linear-solver block with a supported family')
        else:
            check.require(is_linear_solver(config), path + '.family', 'unsupported or missing linear-solver family')


def zero_field(value, field, check, vector=False):
    path = 'initialization.' + field
    value = check.mapping(value, 'option ' + field, path)
    check.equal(value.get('option'), 'value', path + '.option')
    items = value.get(field)
    if vector:
        if not check.require(isinstance(items, list) and len(items) == 3, path, 'requires three zero components'):
            return
    else:
        items = [items]
    check.require(all(check.number(v, path) == 0 for v in items), path, 'only zero initial fields are supported')


def translate(doc):
    check = Checks()
    doc = check.mapping(doc, 'mesh simulation', 'document')
    mesh = check.mapping(doc.get('mesh'), 'file_path automatic_decomposition_type', 'mesh')
    # This selects OpenAccel's partitioner; MARS partitions through Cornerstone.
    if 'automatic_decomposition_type' in mesh:
        check.require(isinstance(mesh['automatic_decomposition_type'], str),
                'mesh.automatic_decomposition_type', 'expected a string')
    mesh_path = mesh.get('file_path')
    check.require(isinstance(mesh_path, str) and mesh_path, 'mesh.file_path', 'expected a path')
    sim = check.mapping(doc.get('simulation'), 'verbose physical_analysis solver material_library', 'simulation')
    physical = check.mapping(sim.get('physical_analysis'), 'analysis_type domains initialization', 'physical_analysis')
    analysis = check.mapping(physical.get('analysis_type'), 'option', 'analysis_type')
    check.equal(analysis.get('option'), 'steady_state', 'analysis_type.option')
    domain = check.mapping(check.singleton(physical.get('domains'), 'domains'),
                     'name location materials type domain_models fluid_models boundaries initialization', 'domain')
    check.equal(domain.get('type'), 'fluid', 'domain.type')
    check.singleton(check.names(domain.get('location'), 'domain.location'), 'domain.location')
    model = check.mapping(domain.get('domain_models', {}), 'reference_pressure', 'domain_models')
    check.require(check.number(model.get('reference_pressure', 0), 'reference_pressure') == 0,
            'reference_pressure', 'only zero reference pressure is supported')
    fluid = check.mapping(domain.get('fluid_models'), 'turbulence', 'fluid_models')
    turbulence = check.mapping(fluid.get('turbulence'), 'option', 'turbulence')
    check.equal(turbulence.get('option'), 'laminar', 'turbulence.option')
    init = dict(check.mapping(physical.get('initialization', {}), 'velocity pressure', 'global initialization'))
    init.update(check.mapping(domain.get('initialization', {}), 'velocity pressure', 'domain initialization'))
    zero_field(init.get('velocity'), 'velocity', check, True)
    zero_field(init.get('pressure'), 'pressure', check)

    material_name = check.singleton(domain.get('materials'), 'domain.materials')
    material = check.mapping(check.singleton(sim.get('material_library'), 'material_library'),
                       'name thermodynamic_properties transport_properties', 'material')
    check.require(material.get('name') == material_name, 'material', 'domain material does not match library')
    thermo = check.mapping(material.get('thermodynamic_properties'), 'equation_of_state', 'thermodynamic_properties')
    eos = check.mapping(thermo.get('equation_of_state'), 'option density', 'equation_of_state')
    check.equal(eos.get('option'), 'value', 'equation_of_state.option')
    transport = check.mapping(material.get('transport_properties'), 'dynamic_viscosity', 'transport_properties')
    viscosity = check.mapping(transport.get('dynamic_viscosity'), 'option dynamic_viscosity', 'dynamic_viscosity')
    check.equal(viscosity.get('option'), 'value', 'dynamic_viscosity.option')
    args = ['--rho', str(check.number(eos.get('density'), 'density', positive=True)),
            '--mu', str(check.number(viscosity.get('dynamic_viscosity'), 'dynamic_viscosity', positive=True))]

    boundaries = domain.get('boundaries')
    if not check.require(isinstance(boundaries, list) and boundaries, 'boundaries', 'missing boundaries'):
        boundaries = []
    selected, wall_sets, kinds = set(), [], []
    for boundary in boundaries:
        b = check.mapping(boundary, 'name type location boundary_details', 'boundary')
        locations = check.names(b.get('location'), 'boundary.location')
        check.require(not selected.intersection(locations), 'boundaries', 'duplicate side-set assignment')
        selected.update(locations)
        kind = b.get('type')
        details = check.mapping(b.get('boundary_details', {}),
                                'mass_and_momentum option' if kind in ('inlet', 'outlet') else 'mass_and_momentum',
                                'boundary_details')
        if kind == 'wall':
            mm = check.mapping(details.get('mass_and_momentum', {}), 'option', 'wall.mass_and_momentum')
            check.equal(mm.get('option', 'no_slip_wall'), 'no_slip_wall', 'wall.option')
            wall_sets.extend(locations)
            continue
        if not check.require(kind in ('inlet', 'outlet') and kind not in kinds,
                             'boundaries', 'requires one inlet and one outlet'):
            continue
        check.equal(details.get('option', 'subsonic'), 'subsonic', kind + '.boundary_details.option')
        kinds.append(kind)
        args += ['--' + kind + '-ss', check.singleton(locations, kind + '.location')]
        mm = check.mapping(details.get('mass_and_momentum'),
                     'option normal_speed' if kind == 'inlet' else 'option relative_pressure pressure_profile_blend',
                     kind + '.mass_and_momentum')
        if kind == 'inlet':
            check.equal(mm.get('option'), 'normal_speed', 'inlet.option')
            args += ['--inlet-velocity', str(check.number(mm.get('normal_speed'), 'normal_speed', positive=True))]
        else:
            option = mm.get('option')
            check.require(option in ('static_pressure', 'average_static_pressure'), 'outlet.option', 'unsupported setting')
            # beta=1 fixes the trace to the prescribed pressure. OpenAccel only
            # reads pressure_profile_blend for average_static_pressure.
            beta = (1.0 if option == 'static_pressure' else
                    check.number(mm.get('pressure_profile_blend'), 'pressure_profile_blend', unit=True))
            args += ['--outlet-pressure', str(check.number(mm.get('relative_pressure'), 'relative_pressure')),
                     '--outlet-beta', str(beta)]
    check.require(set(kinds) == {'inlet', 'outlet'} and wall_sets, 'boundaries', 'requires inlet, outlet and walls')
    args += ['--wall-ss', ','.join(wall_sets)]

    solver = solver_mapping(sim.get('solver'), check)
    control = check.mapping(solver.get('solver_control'), 'basic_settings advanced_options expert_parameters', 'solver_control')
    basic = check.mapping(control.get('basic_settings'),
                    'advection_scheme interpolation_scheme convergence_controls convergence_criteria', 'basic_settings')
    scheme = basic.get('advection_scheme', 'upwind')
    if not check.require(scheme in ('upwind', 'high_resolution'), 'advection_scheme', 'unsupported scheme'):
        scheme = 'upwind'
    args += ['--advection', scheme.replace('_', '-')]
    interp = check.mapping(basic.get('interpolation_scheme', {}),
                     'velocity_interpolation_type pressure_interpolation_type velocity_gradient_interpolation_type pressure_gradient_interpolation_type', 'interpolation_scheme')
    velocity_interp = interp.get('velocity_interpolation_type', 'trilinear')
    if check.require(velocity_interp in ('trilinear', 'linear_linear'),
                     'interpolation_scheme.velocity_interpolation_type', 'unsupported setting'):
        args += ['--velocity-interpolation', velocity_interp.replace('_', '-')]
    for key in ('pressure_interpolation_type',
                'velocity_gradient_interpolation_type', 'pressure_gradient_interpolation_type'):
        check.equal(interp.get(key, 'linear_linear'), 'linear_linear', 'interpolation_scheme.' + key)
    conv = check.mapping(basic.get('convergence_controls'),
                   'min_iterations max_iterations physical_timescale relaxation_parameters', 'convergence_controls')
    args += ['--pseudo-dt', str(check.number(conv.get('physical_timescale'), 'physical_timescale', positive=True))]
    relax = check.mapping(conv.get('relaxation_parameters', {}),
                    'velocity_relaxation_factor pressure_relaxation_factor relax_mass', 'relaxation_parameters')
    for key, flag in [('velocity_relaxation_factor', '--relax-u'),
                      ('pressure_relaxation_factor', '--relax-p'), ('relax_mass', '--relax-mass')]:
        args += [flag, str(check.number(relax.get(key, 1), key, unit=True))]
    check.mapping(basic.get('convergence_criteria', {}), 'residual_type residual_target', 'convergence_criteria')
    advanced = check.mapping(control.get('advanced_options', {}), 'equation_controls linear_solver_settings', 'advanced_options')
    check_linear_settings(advanced.get('linear_solver_settings', {}), solver, check)
    eq = check.mapping(advanced.get('equation_controls', {}), 'sub_iterations', 'equation_controls')
    sub = check.mapping(eq.get('sub_iterations', {}), 'pressure_correction segregated_flow', 'sub_iterations')
    for key in ('pressure_correction', 'segregated_flow'):
        check.require(check.number(sub.get(key, 1), key) == 1, 'sub_iterations.' + key, 'requires one subiteration')
    defaults = dict(consistent=False, fractional_step_method=False, incremental_gradient_change=True,
                    limit_gradients=False, relax_gradients=True, correct_gradients=False,
                    false_mass_accumulation=True, nonlinear_stabilisation=False,
                    disable_momentum_predictor=False, high_speed_blend_damping=False)
    expert = check.mapping(control.get('expert_parameters', {}), ' '.join(defaults) + ' blend_factor_max', 'expert_parameters')
    for key, default in defaults.items():
        expected = False if key == 'relax_gradients' else default
        check.equal(expert.get(key, default), expected, 'expert_parameters.' + key)
    cap = check.number(expert.get('blend_factor_max', 1), 'blend_factor_max')
    check.require(cap == 1 if scheme == 'high_resolution' else cap in (0, 1),
            'blend_factor_max', 'unsupported limiter cap')
    check.finish()
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
                         'OpenAccel named and inline linear-solver settings are not translated; MARS retains its Hypre solvers.',
                         'Coordinates must be in metres; mesh bytes were not read.',
                         'MARS uses its own Hypre settings, zero initial fields and convergence norms.',
                         'Iteration counts and reference residual criteria are not translated.'])
    (o.output / 'case.json').write_text(json.dumps(record, indent=2) + '\n')
    (o.output / 'args.nul').write_bytes(b'\0'.join(v.encode('utf-8') for v in args) + b'\0')
    print('PASS: supported deck; native mesh setup still required. Arguments saved locally.')


if __name__ == '__main__':
    try:
        main()
    except Unsupported as error:
        for issue in str(error).splitlines():
            print('ERROR: ' + issue, file=sys.stderr)
        sys.exit(1)
    except (OSError, UnicodeError, RecursionError):
        # OS and decoder messages can contain private paths or input fragments.
        print('ERROR: preparation: file access, encoding or YAML nesting failure; check paths and permissions locally',
              file=sys.stderr)
        sys.exit(1)
