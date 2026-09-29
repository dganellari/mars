# Synthetic public channel for deck preparation tests.
mesh:
  file_path: channel.exo
simulation:
  verbose: 1
  physical_analysis:
    analysis_type:
      option: steady_state
    domains:
      - name: channel
        location: [fluid]
        materials: [reference_fluid]
        type: fluid
        domain_models:
          reference_pressure: 0
        fluid_models:
          turbulence:
            option: laminar
        boundaries:
          - name: walls
            type: wall
            location: [walls]
          - name: inlet
            type: inlet
            location: [inlet]
            boundary_details:
              mass_and_momentum:
                option: normal_speed
                normal_speed: 0.1
          - name: outlet
            type: outlet
            location: [outlet]
            boundary_details:
              mass_and_momentum:
                option: average_static_pressure
                relative_pressure: 0
                pressure_profile_blend: 0.05
        initialization:
          velocity:
            option: value
            velocity: [0, 0, 0]
          pressure:
            option: value
            pressure: 0
  solver:
    solver_control:
      basic_settings:
        advection_scheme: upwind
        interpolation_scheme:
          velocity_interpolation_type: trilinear
          pressure_interpolation_type: linear_linear
          velocity_gradient_interpolation_type: linear_linear
          pressure_gradient_interpolation_type: linear_linear
        convergence_controls:
          min_iterations: 2
          max_iterations: 2
          physical_timescale: 0.01
          relaxation_parameters:
            velocity_relaxation_factor: 0.3
            pressure_relaxation_factor: 0.3
            relax_mass: 0.75
        convergence_criteria:
          residual_type: RMS
          residual_target: 1.0e-6
      advanced_options:
        equation_controls:
          sub_iterations:
            pressure_correction: 1
            segregated_flow: 1
        linear_solver_settings:
          default:
            family: Trilinos
            max_iterations: 200
            rtol: 1.0e-8
            atol: 1.0e-12
            options:
              belos_solver: gmres
              preconditioner: riluk
              preconditioner_parameters:
                "fact: iluk level-of-fill": 2
                "fact: drop tolerance": 1.0e-3
                "fact: absolute threshold": 1.0e-6
                "fact: relative threshold": 1.0
      expert_parameters:
        consistent: false
        fractional_step_method: false
        incremental_gradient_change: true
        limit_gradients: false
        relax_gradients: false
        blend_factor_max: 0
        nonlinear_stabilisation: false
    output_control:
      file_path: results.e
      output_frequency: 1
      output_fields: [velocity, pressure]
      corrected_boundary_values: false
  material_library:
    - name: reference_fluid
      thermodynamic_properties:
        equation_of_state:
          option: value
          density: 1
      transport_properties:
        dynamic_viscosity:
          option: value
          dynamic_viscosity: 0.1
