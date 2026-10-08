"""Check the production dot-product bound against exact rational arithmetic."""
import math
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest
from fractions import Fraction


class CompensatedDotTest(unittest.TestCase):
    def test_exact_enclosure(self):
        root = Path(__file__).resolve().parents[4]
        source = r'''
#include "mars_segregated_compensated_dot.hpp"
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
int main() {
    int n;
    while (std::cin >> n) {
        mars::segregated::CompensatedDot dot;
        for (int k=0;k<n;++k) {
            std::string a,b; std::cin >> a >> b;
            dot.product(std::strtod(a.c_str(),nullptr),std::strtod(b.c_str(),nullptr));
        }
        std::cout << std::hexfloat << dot.value() << ' ' << dot.error_bound() << '\n';
    }
}
'''
        cases = [[(2.**106, 1.), (2.**53, 1.), (1., 1.),
                  (-2.**106, 1.), (-2.**53, 1.)],
                 [(1. - 2.**-27, 1. + 2.**-27), (-1., 1.)],
                 [(2.**-1022, 2.**-100), (2.**-1074, .5)],
                 [(1., 1.), (-1., 1.)]]
        rng = random.Random(731)
        for _ in range(1000):
            terms = []
            for _ in range(rng.randrange(1, 45)):
                a = math.ldexp(rng.uniform(-1, 1), rng.randrange(-540, 450))
                b = math.ldexp(rng.uniform(-1, 1), rng.randrange(-540, 450))
                terms.extend([(a, b), (-a, b)])
            terms.extend([(1., 2.**rng.randrange(-100, 100))])
            rng.shuffle(terms)
            cases.append(terms)
        data = ''.join(str(len(case)) + '\n' + ''.join(
            a.hex() + ' ' + b.hex() + '\n' for a, b in case) for case in cases)
        # Keep generated files in the checkout, including on the user's laptop.
        scratch = root / '.local-worktrees/simple-performance'
        scratch.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='dot-bound-', dir=str(scratch)) as directory:
            path = Path(directory)
            (path / 'check.cpp').write_text(source)
            subprocess.run([os.environ.get('CXX', 'c++'), '-std=c++20', '-O2',
                            '-ffp-contract=fast', '-I' + str(root / 'backend/distributed/unstructured/fem/segregated'),
                            str(path / 'check.cpp'), '-o', str(path / 'check')], check=True)
            output = subprocess.check_output([str(path / 'check')], input=data, universal_newlines=True)
        rows = output.splitlines()
        self.assertEqual(len(rows), len(cases))
        for number, (case, row) in enumerate(zip(cases, rows)):
            value, bound = map(float.fromhex, row.split())
            exact = sum((Fraction(a) * Fraction(b) for a, b in case), Fraction())
            self.assertTrue(math.isfinite(value) and math.isfinite(bound), number)
            self.assertGreaterEqual(bound, 0, number)
            self.assertLessEqual(abs(exact - Fraction(value)), Fraction(bound), number)
        # The first fixture loses a unit even with compensation; its bound must expose it.
        self.assertEqual(float.fromhex(rows[0].split()[0]), 0.)
        self.assertGreaterEqual(float.fromhex(rows[0].split()[1]), 1.)


if __name__ == '__main__':
    unittest.main()
