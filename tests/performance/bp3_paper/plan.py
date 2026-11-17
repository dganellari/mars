"""Explicit workloads for matched-preconditioner and best-configuration comparisons."""
import math

VARIANTS = ('matched', 'best_h0', 'best_h', 'mfem')


def cases(phase):
    if phase == 'standard':
        weak = [('weak', ranks, p, ne) for ranks, factor in ((1, 1), (8, 2), (64, 4))
                for p, base in ((3, 64), (7, 32)) for ne in (base * factor,)]
        strong = [('strong', ranks, p, ne) for ranks in (1, 4, 16, 64)
                  for p, ne in ((3, 128), (7, 64))]
        return weak + strong
    if phase == 'large':
        return [('strong', 256, 3, 128), ('strong', 256, 7, 64)]
    raise ValueError('Unknown phase: ' + phase)


def depth(variant, p):
    return (2 if p == 3 else 1) if variant == 'best_h' else 0


def mars_args(bundle, variant, p, ne, checked=False):
    if variant not in VARIANTS[:-1]:
        raise ValueError('Unknown MARS configuration')
    # h=0 alone still permits consolidation; the historical MFEM baseline disables both.
    matched = variant == 'matched'
    args = ['--ptx-dir', str(bundle / 'generated'), '--spec', 'bp3', '--p', str(p), '--ne', str(ne),
            '--deform', '0.05', '--stretch', '1', '--rhs', 'discrete', '--mg', '--amg',
            '--amg-cycles', '3', '--amg-consolidate', '0' if matched else '-1',
            '--h-levels', str(depth(variant, p)), '--pre', '2', '--post', '2',
            '--cheb-interval', 'mfem' if matched else 'mars', '--tol', '1e-13',
            '--true-tol', '1e-10', '--maxit', '2000', '--solves', '5', '--reps', '1', '--amg-print']
    return args + (['--h-check'] if checked else [])


def mfem_args(p, ne):
    return ['--p', str(p), '--ne', str(ne), '--deform', '0.05', '--stretch', '1', '--rhs', 'discrete',
            '--pc', 'mg', '--cheb-order', '2', '--amg-cycles', '3', '--amg-settings', 'mars',
            '--gpu-aware-mpi', '--device', 'cuda', '--tol', '1e-13', '--true-tol', '1e-10',
            '--maxit', '2000', '--solves', '5', '--reps', '1']


def dimensions(ranks):
    # Balanced Cartesian factorizations for the fixed power-of-two rank counts in this plan.
    if ranks not in (1, 4, 8, 16, 64, 256):
        raise ValueError('Unplanned rank count')
    dims = [1, 1, 1]
    for _ in range(int(math.log2(ranks))):
        dims[dims.index(min(dims))] *= 2
    return sorted(dims, reverse=True)


def validate_case(ranks, p, ne):
    dims = dimensions(ranks)
    factor = 1 << depth('best_h', p)
    if any(ne % (d * factor) or ne // (d * factor) < 2 for d in dims):
        raise ValueError('Unaligned h hierarchy')
    ne0 = ne
    while ne0 > 32 and ne0 % 2 == 0:
        ne0 //= 2
    if ne0 < max(dims):
        raise ValueError('MFEM serial partition is too small')
    return dict(rank_grid=dims, elements_per_gpu=ne**3 // ranks,
                element_local_dofs_per_gpu=ne**3 * (p + 1)**3 // ranks,
                global_unique_dofs=(ne * p + 1)**3, h_depth=depth('best_h', p))
