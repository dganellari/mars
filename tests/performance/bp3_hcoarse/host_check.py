"""Check the shared transfer code and the Galerkin hierarchy against an assembled oracle."""
import argparse
import ast
import ctypes
import json
from pathlib import Path
import subprocess
from types import SimpleNamespace

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]


def reference(deform):
    path = REPO / 'marsir-mlir/test/cellwise_bp_mg_ref.py'
    tree = ast.parse(path.read_text())
    definitions = [n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef)) or
                   isinstance(n, (ast.Import, ast.ImportFrom)) or
                   isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id in
                   ('corner_ref', 'SYM', 'INTERVAL') for t in n.targets)]
    scope = {'args': SimpleNamespace(deform=deform, cheb_interval='mars')}
    exec(compile(ast.Module(body=definitions, type_ignores=[]), str(path), 'exec'), scope)
    return scope


def compiled_transfer(out):
    source = out / 'transfer.cpp'
    source.write_text('''#include "backend/distributed/unstructured/solvers/mars_cellwise_h_transfer.hpp"
namespace bp3_h = mars::cellwise::h_coarse;
extern "C" void transfer(int nx, int ny, int nz, const double* x, double* y, int transpose) {
    const bp3_h::Grid g{nx,ny,nz};
    const long long n = (transpose ? g.coarse().elements() : g.elements())*8;
    for(long long i=0;i<n;++i) y[i]=transpose ? bp3_h::restrict_value(x,g,i) : bp3_h::prolong_value(x,g,i);
}
''')
    library = out / 'transfer.so'
    subprocess.run(['c++', '-std=c++17', '-O2', '-shared', '-fPIC', '-Wall', '-Wextra', '-Werror',
                    '-I' + str(REPO), str(source), '-o', str(library)], check=True)
    dll = ctypes.CDLL(str(library.resolve()))
    ptr = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
    dll.transfer.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ptr, ptr, ctypes.c_int]
    dll.transfer.restype = None

    def call(grid, x, transpose=False):
        n = np.prod(np.array(grid) // (2 if transpose else 1)) * 8
        y = np.empty(int(n))
        dll.transfer(*grid, np.ascontiguousarray(x), y, int(transpose))
        return y.reshape(-1, 8)
    return call


def stencil():
    # Independent coordinate construction, without the C++ bit indexing.
    nodes = np.array(list(np.ndindex(2, 2, 2)))
    points = (nodes[:, None, :] + nodes[None, :, :]) / 2
    return np.prod(np.where(nodes[None, None, :, :] == 1,
                            points[:, :, None, :], 1 - points[:, :, None, :]), axis=-1)


def children(grid):
    coarse = tuple(n // 2 for n in grid)
    parent = np.array(list(np.ndindex(*coarse)))
    child = np.array(list(np.ndindex(2, 2, 2)))
    xyz = 2 * parent[:, None, :] + child[None, :, :]
    return np.ravel_multi_index(xyz.transpose(2, 0, 1), grid)


def transfer_checks(call):
    rng = np.random.default_rng(1234)
    result = []
    weights = stencil()
    for grid in ((2, 2, 2), (4, 6, 8), (8, 4, 2), (8, 8, 8)):
        ids = children(grid)
        c = rng.normal(size=(len(ids), 8))
        f = rng.normal(size=(np.prod(grid), 8))
        p = call(grid, c)
        r = call(grid, f, True)
        expected_p = np.empty_like(f)
        expected_p[ids] = np.einsum('cij,ej->eci', weights, c)
        expected_r = np.einsum('cij,eci->ej', weights, f[ids])
        assert np.max(np.abs(expected_p - p)) < 3e-15
        assert np.max(np.abs(expected_r - r)) < 1e-14
        assert np.max(np.abs(call(grid, np.ones_like(c)) - 1)) == 0
        adj = abs(np.vdot(f, p) - np.vdot(r, c)) / (np.linalg.norm(f) * np.linalg.norm(p))
        assert adj < 1e-14
        # Every fine copy of the same lattice point must get the same interpolation.
        cgrid = tuple(n // 2 for n in grid)
        ci = np.array(list(np.ndindex(*cgrid)))
        nodes = np.array(list(np.ndindex(2, 2, 2)))
        cg = np.ravel_multi_index((ci[:, None, :] + nodes).transpose(2, 0, 1), tuple(n + 1 for n in cgrid))
        cv = rng.normal(size=np.prod(tuple(n + 1 for n in cgrid)))
        p = call(grid, cv[cg])
        fi = np.array(list(np.ndindex(*grid)))
        fg = np.ravel_multi_index((fi[:, None, :] + nodes).transpose(2, 0, 1), tuple(n + 1 for n in grid))
        average = np.bincount(fg.ravel(), weights=p.ravel()) / np.bincount(fg.ravel())
        cont = np.max(np.abs(p - average[fg]))
        assert cont < 2e-15
        # A missing child or multiplicity-weighted restriction must fail the same oracle.
        bad = np.einsum('cij,eci->ej', weights[:-1], f[ids[:, :-1]])
        assert np.linalg.norm(bad - expected_r) > 1e-2
        assert np.linalg.norm(expected_r / 8 - expected_r) > 1e-2
        result.append(dict(grid=grid, adjoint=adj, continuity=cont))
    return result


class MatrixLevel:
    def __init__(self, ne, matrix):
        self.ne, self.E, self.n3, self.NG = ne, ne**3, 8, ne + 1
        self.matrix = matrix
        element = np.array(list(np.ndindex(ne, ne, ne)))
        node = np.array(list(np.ndindex(2, 2, 2)))
        xyz = element[:, None, :] + node
        self.gidx = np.ravel_multi_index(xyz.transpose(2, 0, 1), (self.NG,) * 3)
        self.mask = np.all((xyz > 0) & (xyz < ne), axis=-1)
        self.inner = np.unique(self.gidx[self.mask])
        self.diag = self.dss(np.diagonal(matrix, axis1=1, axis2=2))
        self.bound = float(np.max((self.dss(np.sum(np.abs(matrix), axis=-1)) / self.diag)[self.mask]))

    def apply(self, x):
        return np.einsum('eij,ej->ei', self.matrix, x)

    def scatter(self, x):
        return np.bincount(self.gidx.ravel(), weights=x.ravel(), minlength=self.NG**3)

    def dss(self, x):
        return self.scatter(x)[self.gidx]

    def jacobi(self, x):
        return self.mask * self.dss(x) / self.diag

    def nodes(self, x):
        return self.scatter(x) / self.scatter(np.ones_like(x))

    def assembled(self):
        a = np.zeros((self.NG**3, self.NG**3))
        np.add.at(a, (self.gidx[:, :, None], self.gidx[:, None, :]), self.matrix)
        return a[np.ix_(self.inner, self.inner)]


class Hierarchy:
    def __init__(self, fine, ne, depth, call, ref):
        self.call, self.ref = call, ref
        k = np.stack([fine.apply(np.tile(np.eye(8)[j], (fine.E, 1))) for j in range(8)], axis=-1)
        self.levels = [MatrixLevel(ne, k)]
        self.galerkin_errors = []
        for _ in range(depth):
            parent = self.levels[-1]
            grid = (parent.ne,) * 3
            k = np.stack([call(grid, parent.apply(call(grid, np.tile(np.eye(8)[j], (parent.E // 8, 1)))), True)
                          for j in range(8)], axis=-1)
            oracle = np.einsum('cij,ecik,ckl->ejl', stencil(), parent.matrix[children(grid)], stencil())
            err = float(np.max(np.abs(k - oracle)) / np.max(np.abs(oracle)))
            assert err < 5e-14
            self.galerkin_errors.append(err)
            self.levels.append(MatrixLevel(parent.ne // 2, k))
        self.bottom_inverse = np.linalg.inv(self.levels[-1].assembled())

    def smooth(self, lv, b, x, zero):
        for i, (c1, c2) in enumerate(self.ref['chebyshev_steps'](lv.bound, 2)):
            z = lv.jacobi(b if i == 0 and zero else b - lv.apply(x))
            d = c2 * z if i == 0 else c1 * d + c2 * z
            x = d.copy() if i == 0 and zero else x + d
        return x

    def __call__(self, b, level=0):
        lv = self.levels[level]
        if level + 1 == len(self.levels):
            u = np.zeros(lv.NG**3)
            u[lv.inner] = self.bottom_inverse @ lv.scatter(b)[lv.inner]
            return u[lv.gidx]
        x = self.smooth(lv, b, np.zeros_like(b), True)
        r = (b - lv.apply(x)) * lv.mask
        c = self(self.call((lv.ne,) * 3, r, True), level + 1)
        x += self.call((lv.ne,) * 3, c)
        return self.smooth(lv, b, x, False)


def matrix_checks(h):
    lv = h.levels[0]
    m = len(lv.inner)
    mat = np.zeros((m, m))
    count = lv.dss(np.ones((lv.E, 8)))
    for j, g in enumerate(lv.inner):
        raw = (lv.gidx == g) / count
        mat[:, j] = lv.nodes(h(raw))[lv.inner]
    symmetry = float(np.linalg.norm(mat - mat.T) / np.linalg.norm(mat))
    eig = np.linalg.eigvalsh((mat + mat.T) / 2)
    assert symmetry < 1e-13 and eig[0] > 0
    a = lv.assembled()
    diagonal = np.diag(a)
    spectrum = np.linalg.eigvalsh(a / np.sqrt(diagonal[:, None] * diagonal[None, :]))
    assert spectrum[-1] <= lv.bound * (1 + 1e-13)
    return dict(symmetry=symmetry, min_eigenvalue=float(eig[0]), bound=lv.bound,
                actual_lmax=float(spectrum[-1]))


def run(out):
    out.mkdir(parents=True, exist_ok=True)
    call = compiled_transfer(out)
    result = dict(transfers=transfer_checks(call), solves=[])
    rng = np.random.default_rng(525)
    for ne, depth, p, deform in ((4, 1, 3, 0), (4, 1, 3, .05), (4, 1, 7, .05),
                                (8, 1, 3, .05), (8, 2, 3, .05), (8, 2, 7, .05),
                                (16, 3, 1, .05)):
        ref = reference(deform)
        corners, eidx = ref['mesh'](ne)
        degrees = [7, 3, 1] if p == 7 else [3, 1] if p == 3 else [1]
        levels = [ref['Level'](degree, ne, corners, eidx) for degree in degrees]
        h = Hierarchy(levels[-1], ne, depth, call, ref)
        entry = dict(ne=ne, depth=depth, p=p, deform=deform, galerkin=h.galerkin_errors)
        if ne == 4 and p == 3:
            entry['spd'] = matrix_checks(h)
        vc = ref['VCycle']

        class PMG(vc):
            def __call__(self, b, k=0):
                return h(b) if k + 1 == len(self.levels) else super().__call__(b, k)

        pmg = object.__new__(PMG)
        pmg.levels, pmg.pre, pmg.post = levels, 2, 2
        pmg.I = [ref['interpolation'](levels[i], levels[i + 1]) for i in range(len(levels) - 1)]
        pmg.lmax = [pmg.power(lv) for lv in levels]
        fine = levels[0]
        for rhs in ('smooth', 'random'):
            u = fine.mask * (np.prod(np.sin(np.pi * fine.xn), axis=-1) if rhs == 'smooth'
                             else rng.normal(size=fine.NG**3)[fine.gidx])
            b = fine.apply(u)
            x, hist = ref['fcg_history'](fine, pmg, b, 1e-13, 2000)
            true = fine.true_residual(b, x) / fine.true_residual(b, np.zeros_like(b))
            assert len(hist) < 2000 and hist[-1] <= 1e-13 and true <= 1e-10
            entry[rhs] = dict(iterations=len(hist), relative=float(hist[-1]), true_residual=float(true))
        # Changing the RHS and returning to it must not retain an old correction.
        b1 = rng.normal(size=(h.levels[0].E, 8)) * h.levels[0].mask
        first = h(b1)
        h(2 * b1)
        assert np.array_equal(first, h(b1))
        result['solves'].append(entry)
        (out / 'host-results.json').write_text(json.dumps(result, indent=2) + '\n')
        print(json.dumps(entry), flush=True)
    print('PASS shared transfer code, Galerkin identity, small full SPD matrices, and 14 FCG solves', flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    run(parser.parse_args().output)
