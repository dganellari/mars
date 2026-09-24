"""GPT/Codex, 2026-09-12: passive mean-pressure outlet host candidate.

This is a different boundary law, not a repair or reinterpretation of beta.
It reuses the public spatial replica and solves the linear block system exactly.
No production kernels, nonlinear advection, or MPI execution are tested here.
"""

import argparse
import json

import numpy as np

from outlet_coupled_time_check import CoupledStep, relative_error


class PassiveTraceStep:
    def __init__(self, dt_eff, nu, rho, resistance):
        base = CoupledStep(dt_eff, nu, 1)
        self.op = op = base.op
        self.nu = nu
        self.zeta = dt_eff * resistance / rho
        self.outlet = np.flatnonzero(op.outlet)
        area = np.zeros(op.n)
        for _, nodes, _, _, area_vector in op.faces:
            area[nodes] += np.linalg.norm(area_vector) / 3
        self.area = area[self.outlet]
        self.projector = np.eye(len(self.outlet)) - np.outer(
            np.ones(len(self.outlet)), self.area / self.area.sum())
        self.gt = op.gt[:, self.outlet]
        # N*u is raw integrated outward flux, the conjugate of trace pressure work.
        self.normal_flux = self.gt.T * np.repeat(op.mass, 3)[None, :]
        self.trace_velocity = self.zeta * self.projector @ (
            self.normal_flux / self.area[:, None])
        self.g = op.gv[op.free3]
        self.l = op.ct[:, self.outlet] + op.smooth @ self.gt
        self.s = op.cp + op.smooth @ op.gv
        feedback_free = self.trace_velocity[:, op.free3]
        feedback_fixed = self.trace_velocity[:, op.fixed3]
        self.a = base.a + self.gt[op.free3] @ feedback_free
        self.b = op.bfree + self.l @ feedback_free
        self.b_fixed = base.b_fixed + self.l @ feedback_fixed
        self.boundary_load = base.boundary_load - self.gt[op.free3] @ feedback_fixed
        self.h = np.linalg.solve(self.a, np.eye(len(op.free3)))
        self.response = self.h @ self.g
        self.schur = self.s - self.b @ self.response
        self.schur_inverse = np.linalg.solve(self.schur, np.eye(op.n))
        self.velocity_map = self.h + self.response @ self.schur_inverse @ self.b @ self.h

    def trace(self, velocity, reference=0.0):
        return reference + self.trace_velocity @ velocity.ravel()

    def residual(self, velocity, z, reference=0.0):
        op = self.op
        trace = self.trace(velocity, reference)
        gradient = op.gv @ z + self.gt @ trace
        return op.b @ velocity.ravel() + op.cp @ z + op.ct[:, self.outlet] @ trace + op.smooth @ gradient

    def step(self, velocity, previous, target, bdf):
        op = self.op
        history = velocity if bdf == 1 else (4 * velocity - previous) / 3
        fixed = target.ravel()[op.fixed3]
        rhs = history.ravel()[op.free3] + self.boundary_load @ fixed
        unconstrained = self.h @ rhs
        z = self.schur_inverse @ (-self.b @ unconstrained - self.b_fixed @ fixed)
        result = target.copy()
        result[op.free] = (unconstrained - self.response @ z).reshape(-1, 3)
        return result, z, rhs

    def check(self):
        op = self.op
        rng = np.random.default_rng(20260912)
        u, previous, target = (rng.normal(size=(op.n, 3)) for _ in range(3))
        result, z, rhs = self.step(u, previous, target, 2)
        fixed = target.ravel()[op.fixed3]
        state = np.r_[result.ravel()[op.free3], z]
        block = np.block([[self.a, self.g], [self.b, self.s]])
        block_rhs = np.r_[rhs, -self.b_fixed @ fixed]
        direct = np.linalg.solve(block, block_rhs)
        checks = {"schur_vs_block": relative_error(state, direct)}
        checks['block_backward_error'] = float(np.linalg.norm(block @ state - block_rhs) / (
            np.linalg.norm(block) * np.linalg.norm(state) + np.linalg.norm(block_rhs)))
        residual = self.residual(result, z)
        checks['continuity_relative_residual'] = float(np.linalg.norm(residual) / max(
            np.linalg.norm(op.b @ result.ravel()) + np.linalg.norm(self.s @ z), 1e-30))

        # This evaluates the physical opening law separately from the matrix map.
        normal_velocity = self.normal_flux @ u.ravel() / self.area
        fluctuation = normal_velocity - self.area @ normal_velocity / self.area.sum()
        trace = self.trace(u)
        checks['trace_law'] = relative_error(trace, self.zeta * fluctuation) if self.zeta else float(np.linalg.norm(trace))
        checks['mean_trace'] = float(abs(self.area @ trace) / max(
            np.linalg.norm(self.area) * np.linalg.norm(trace), 1e-30))
        boundary_work = float(trace @ self.normal_flux @ u.ravel())
        dissipation = float(self.zeta * (self.area @ fluctuation**2))
        checks['boundary_work_identity'] = abs(boundary_work - dissipation) / max(abs(dissipation), 1e-30)
        if boundary_work < -1e-12:
            raise RuntimeError('passive trace produced negative boundary pressure work')

        mass_free = np.repeat(op.mass[op.free], 3)
        defect = mass_free[:, None] * self.g + op.bfree.T
        checks['volume_adjoint'] = float(np.linalg.norm(defect) / np.linalg.norm(op.bfree))
        stiffness = mass_free[:, None] * (
            self.gt[op.free3] @ self.trace_velocity[:, op.free3])
        checks['impedance_symmetry'] = float(np.linalg.norm(stiffness - stiffness.T) / max(
            np.linalg.norm(stiffness), 1e-30)) if self.zeta else float(np.linalg.norm(stiffness))

        # Shift all physical pressure and the reference together: neither flux nor gradient changes.
        constant_trace = np.ones(len(self.outlet))
        constant_gradient = op.gv @ np.ones(op.n) + self.gt @ constant_trace
        checks['constant_pressure_gradient'] = float(np.max(np.abs(constant_gradient[op.free3])))
        checks['pressure_reference_shift'] = relative_error(
            self.residual(result, z + 1, 1), residual) if np.linalg.norm(residual) > 1 else float(
                np.linalg.norm(self.residual(result, z + 1, 1) - residual))

        dz = rng.normal(size=op.n)
        epsilon = 1e-5
        direction = np.zeros_like(u)
        direction[op.free] = (-self.response @ dz).reshape(-1, 3)
        plus = self.residual(result + epsilon * direction, z + epsilon * dz)
        minus = self.residual(result - epsilon * direction, z - epsilon * dz)
        checks['schur_finite_difference'] = relative_error((plus - minus) / (2 * epsilon), self.schur @ dz)

        u[op.fixed_nodes] = 0
        previous[op.fixed_nodes] = 0
        un, _, _ = self.step(u, previous, np.zeros_like(u), 2)
        mapped = self.velocity_map @ ((4 * u - previous) / 3).ravel()[op.free3]
        checks['time_map_action'] = relative_error(mapped, un.ravel()[op.free3])
        eigenvalues, vectors = np.linalg.eig(self.velocity_map)
        discriminant = np.sqrt(4 * eigenvalues**2 - 3 * eigenvalues + 0j)
        roots = np.r_[(2 * eigenvalues + discriminant) / 3, (2 * eigenvalues - discriminant) / 3]
        index = int(np.argmax(np.abs(roots)))
        vector = vectors[:, index % len(eigenvalues)]
        root = roots[index]
        state = np.r_[root * vector, vector]
        action = np.r_[self.velocity_map @ ((4 * root - 1) * vector / 3), root * vector]
        checks['bdf2_eigenpair_residual'] = relative_error(action, root * state)
        if self.zeta == 0:
            fixed_trace = CoupledStep(op.dt_eff, self.nu, 1)
            checks['zero_resistance_fixed_trace'] = relative_error(self.velocity_map, fixed_trace.velocity_map)
        for name, value in checks.items():
            limit = 1e-7 if name == 'schur_finite_difference' else 1e-10
            if not np.isfinite(value) or value > limit:
                raise RuntimeError(f'{name}={value} exceeds {limit}')
        return dict(checks=checks, spectral_radius=float(abs(root)),
                    unstable_modes=int(np.count_nonzero(abs(roots) > 1 + 1e-8)),
                    scaled_resistance=self.zeta, boundary_work=boundary_work)


def pressure_feedback_work(dt_eff, nu):
    model = CoupledStep(dt_eff, nu, .05)
    op = model.op
    values, vectors = np.linalg.eig(model.velocity_map)
    index = int(np.argmax(np.abs(values)))
    history = vectors[:, index]
    velocity = model.velocity_map @ history
    z = -model.schur_bh @ history
    mass = np.repeat(op.mass[op.free], 3)
    gradient = model.g @ z
    volume_work = float(-np.vdot(z, model.b @ velocity).real)
    trace_work = float(np.vdot(velocity, mass * (op.gt @ op.trace @ z)[op.free3]).real)
    total_work = float(np.vdot(velocity, mass * gradient).real)
    if abs(total_work - volume_work - trace_work) > 1e-12 * max(
            abs(volume_work) + abs(trace_work), 1e-30):
        raise RuntimeError('pressure work decomposition failed')
    affine = op.xyz @ np.array([2., -3., 5.])
    exact_trace_gradient = (op.gv @ affine + op.gt @ (op.outlet * affine)).reshape(-1, 3)
    gradient_error = np.max(np.abs(exact_trace_gradient - [2, -3, 5]), axis=1)
    return dict(velocity_eigenvalue=[float(values[index].real), float(values[index].imag)],
                volume_pressure_work=volume_work, trace_pressure_work=trace_work,
                total_pressure_work=total_work,
                affine_gradient_error_free_interior=float(np.max(gradient_error[~op.fixed & ~op.outlet])),
                affine_gradient_error_outlet=float(np.max(gradient_error[op.outlet])),
                normalization='unit Euclidean norm of free velocity history eigenvector')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dt', type=float, default=2e-6)
    parser.add_argument('--nu', type=float, default=1e-4)
    parser.add_argument('--rho', type=float, default=1000)
    parser.add_argument('--resistance', type=float, required=True, help='Pa s/m; a physical parameter, not beta')
    parser.add_argument('--steps', type=int, default=200)
    parser.add_argument('--audit-pressure-feedback', action='store_true')
    args = parser.parse_args()
    if (not all(np.isfinite(v) for v in (args.dt, args.nu, args.rho, args.resistance))
            or args.dt <= 0 or args.nu < 0 or args.rho <= 0 or args.resistance < 0 or args.steps < 0):
        parser.error('require finite dt>0, nu>=0, rho>0, resistance>=0, steps>=0')
    later = PassiveTraceStep(2 * args.dt / 3, args.nu, args.rho, args.resistance)
    evidence = later.check()
    evidence['model'] = 'passive velocity-feedback trace; fully implicit linear public host replica'
    evidence['parameters'] = vars(args)
    if args.audit_pressure_feedback:
        evidence['old_beta_005_pressure_work'] = pressure_feedback_work(2 * args.dt / 3, args.nu)
    if args.steps:
        first = PassiveTraceStep(args.dt, args.nu, args.rho, args.resistance)
        u = np.zeros((later.op.n, 3))
        previous = u.copy()
        trajectory = []
        for step in range(1, args.steps + 1):
            model = first if step == 1 else later
            target = np.zeros_like(u)
            target[model.op.inlet, 0] = .5 * min(1, step / 100)
            un, z, _ = model.step(u, previous, target, 1 if step == 1 else 2)
            residual = model.residual(un, z)
            if not np.all(np.isfinite(un)) or not np.all(np.isfinite(z)):
                raise RuntimeError(f'nonfinite trajectory at step {step}')
            if step in (1, 10, 50, 100, args.steps):
                trace = model.trace(un)
                physical_pressure_fluctuation = trace * args.rho / model.op.dt_eff
                trajectory.append(dict(step=step, speed_max=float(np.max(np.linalg.norm(un, axis=1))),
                    trace_fluctuation_max_pa=float(np.max(np.abs(physical_pressure_fluctuation))),
                    continuity_rms=float(np.sqrt(np.sum(residual**2 / model.op.mass) / model.op.mass.sum()))))
            previous, u = u, un
        evidence['trajectory'] = trajectory
    print(json.dumps(evidence, indent=2), flush=True)
    if evidence['unstable_modes']:
        raise RuntimeError('candidate fails the homogeneous BDF2 stability gate')


if __name__ == '__main__':
    main()
