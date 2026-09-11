"""GPT/Codex, 2026-09-11: public host model of fully implicit outlet coupling.

This counterfactual uses the spatial operators from outlet_temporal_audit.py.
It does not execute production kernels or implement a production SIMPLE solver.
"""

import argparse
import json

import numpy as np

from outlet_temporal_audit import Operators


def relative_error(actual, expected):
    return float(np.linalg.norm(actual - expected) / max(np.linalg.norm(expected), 1e-30))


class CoupledStep:
    def __init__(self, dt_eff, nu, beta, compact_only=False):
        self.op = op = Operators(dt_eff, nu, beta)
        self.h = np.kron(op.hfree, np.eye(3))
        self.g = op.gpred[op.free3]
        self.b = op.bfree
        self.b_fixed = op.b[:, op.fixed3]
        self.smooth = np.zeros_like(op.smooth) if compact_only else op.smooth
        self.s = op.cp + op.ct @ op.trace + self.smooth @ op.gpred
        self.response = self.h @ self.g
        self.schur = self.s - self.b @ self.response
        # A_hat = M^-1 * dt_eff * (M/dt_eff + nu*K), on free velocity DOFs.
        a = np.eye(len(op.free)) + dt_eff * nu * (
            op.stiffness[np.ix_(op.free, op.free)] / op.mass[op.free, None])
        self.a = np.kron(a, np.eye(3))
        self.boundary_load = np.kron(-dt_eff * nu * (
            op.stiffness[np.ix_(op.free, op.fixed_nodes)] / op.mass[op.free, None]), np.eye(3))
        self.schur_bh = np.linalg.solve(self.schur, self.b @ self.h)
        self.velocity_map = self.h + self.response @ self.schur_bh

    def pressure_flux(self, z):
        op = self.op
        trace = op.trace @ z
        gradient = op.gv @ z + op.gt @ trace
        return op.cp @ z + op.ct @ trace + self.smooth @ gradient

    def step(self, u, previous_u, target, bdf):
        op = self.op
        history = u if bdf == 1 else (4 * u - previous_u) / 3
        fixed_values = target.ravel()[op.fixed3]
        rhs = history.ravel()[op.free3] + self.boundary_load @ fixed_values
        unconstrained = self.h @ rhs
        z = np.linalg.solve(self.schur, -self.b @ unconstrained - self.b_fixed @ fixed_values)
        result = target.copy()
        result[op.free] = (unconstrained - self.response @ z).reshape(-1, 3)
        return result, z, rhs

    def check(self):
        op = self.op
        rng = np.random.default_rng(20260911)
        u, previous_u, target = (rng.normal(size=(op.n, 3)) for _ in range(3))
        result, z, rhs = self.step(u, previous_u, target, 2)
        free_result = result.ravel()[op.free3]
        fixed_values = target.ravel()[op.fixed3]

        # Solve the original block equations without Schur elimination.
        block = np.block([[self.a, self.g], [self.b, self.s]])
        block_rhs = np.r_[rhs, -self.b_fixed @ fixed_values]
        direct = np.linalg.solve(block, block_rhs)
        eliminated = np.r_[free_result, z]
        checks = {"schur_vs_block": relative_error(eliminated, direct)}
        checks["block_backward_error"] = float(np.linalg.norm(block @ eliminated - block_rhs) / (
            np.linalg.norm(block) * np.linalg.norm(eliminated) + np.linalg.norm(block_rhs)))
        checks["momentum_residual"] = relative_error(self.a @ free_result + self.g @ z, rhs)
        residual = self.b @ free_result + self.b_fixed @ fixed_values + self.pressure_flux(z)
        checks["continuity_relative_residual"] = float(np.linalg.norm(residual) / (
            np.linalg.norm(self.b @ free_result) + np.linalg.norm(self.pressure_flux(z)) + 1e-30))

        # Differentiate the separately evaluated flux after the momentum response.
        dz = rng.normal(size=op.n)
        epsilon = 1e-5
        flux_plus = self.b @ (free_result - epsilon * self.response @ dz) + self.pressure_flux(z + epsilon * dz)
        flux_minus = self.b @ (free_result + epsilon * self.response @ dz) + self.pressure_flux(z - epsilon * dz)
        checks["schur_finite_difference"] = relative_error((flux_plus - flux_minus) / (2 * epsilon), self.schur @ dz)

        zero = np.zeros_like(u)
        u[op.fixed_nodes] = 0
        previous_u[op.fixed_nodes] = 0
        un, _, _ = self.step(u, previous_u, zero, 2)
        mapped = self.velocity_map @ ((4 * u - previous_u) / 3).ravel()[op.free3]
        checks["time_map_action"] = relative_error(mapped, un.ravel()[op.free3])

        eigenvalues, vectors = np.linalg.eig(self.velocity_map)
        discriminant = np.sqrt(4 * eigenvalues**2 - 3 * eigenvalues + 0j)
        roots = np.r_[(2 * eigenvalues + discriminant) / 3,
                      (2 * eigenvalues - discriminant) / 3]
        index = int(np.argmax(np.abs(roots)))
        vector = vectors[:, index % len(eigenvalues)]
        root = roots[index]
        state = np.r_[root * vector, vector]
        action = np.r_[self.velocity_map @ ((4 * root - 1) * vector / 3), root * vector]
        checks["bdf2_eigenpair_residual"] = relative_error(action, root * state)
        for name, value in checks.items():
            limit = 1e-7 if name == "schur_finite_difference" else 1e-10
            if not np.isfinite(value) or value > limit:
                raise RuntimeError(f"{name}={value} exceeds {limit}")
        outlet = np.repeat(op.outlet[op.free], 3)
        energy = np.repeat(op.mass[op.free], 3) * np.abs(vector)**2
        return dict(checks=checks, spectral_radius=float(abs(root)),
                    unstable_modes=int(np.count_nonzero(abs(roots) > 1 + 1e-8)),
                    leading_velocity_eigenvalue=[float(eigenvalues[index % len(eigenvalues)].real),
                                                float(eigenvalues[index % len(eigenvalues)].imag)],
                    leading_mode_outlet_energy_fraction=float(energy[outlet].sum() / energy.sum()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dt', type=float, default=2e-6)
    parser.add_argument('--nu', type=float, default=1e-4)
    parser.add_argument('--beta', type=float, default=.05)
    parser.add_argument('--steps', type=int, default=200)
    parser.add_argument('--compact-only', action='store_true')
    parser.add_argument('--expect', choices=['stable', 'unstable'])
    args = parser.parse_args()
    if args.dt <= 0 or args.nu < 0 or not 0 < args.beta <= 1 or args.steps < 0:
        parser.error('require dt>0, nu>=0, 0<beta<=1 and steps>=0')
    later = CoupledStep(2 * args.dt / 3, args.nu, args.beta, args.compact_only)
    evidence = later.check()
    evidence['model'] = 'fully implicit linear host replica; advection disabled'
    evidence['parameters'] = vars(args)
    if args.steps:
        first = CoupledStep(args.dt, args.nu, args.beta, args.compact_only)
        u = np.zeros((first.op.n, 3))
        previous_u = u.copy()
        trajectory = []
        for step in range(1, args.steps + 1):
            model = first if step == 1 else later
            target = np.zeros_like(u)
            target[model.op.inlet, 0] = .5 * min(1, step / 100)
            un, z, _ = model.step(u, previous_u, target, 1 if step == 1 else 2)
            if not np.all(np.isfinite(un)) or not np.all(np.isfinite(z)):
                raise RuntimeError(f'nonfinite trajectory at step {step}')
            if step in (1, 10, 50, 100, args.steps):
                trajectory.append(dict(step=step, speed_max=float(np.max(np.linalg.norm(un, axis=1)))))
            previous_u, u = u, un
        evidence['trajectory'] = trajectory
    print(json.dumps(evidence, indent=2), flush=True)
    if args.expect and (evidence['unstable_modes'] == 0) != (args.expect == 'stable'):
        raise RuntimeError('spectrum does not match the requested classification')


if __name__ == '__main__':
    main()
