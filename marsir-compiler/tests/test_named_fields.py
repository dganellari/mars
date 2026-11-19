"""Check the packed stored-field ABI independently of the GPU schedule."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from marsir import parse_spec, synthesize
from marsir.backends.mlir_ir import emit_galerkin


def form(body):
    return synthesize(parse_spec("operator test {\nelement = hex\nform = galerkin\n" + body + "\n}"))


class Fields(unittest.TestCase):
    def test_order(self):
        f = form("fields = growth, alpha\nv = (alpha + g2 * growth)*u")
        self.assertEqual(f.packed_fields, ("g0", "g1", "g2", "growth", "alpha"))
        self.assertIn("%G: tensor<5x5x5x5xf64>", emit_galerkin(f, p=3))

    def test_named_only(self):
        f = form("fields = sigma\nv = sigma*u")
        self.assertEqual(f.packed_fields, ("sigma",))
        self.assertIn("%G: tensor<1x5x5x5xf64>", emit_galerkin(f, p=3))

    def test_legacy_layout(self):
        self.assertEqual(form("v = g4*u").packed_fields, tuple("g%d" % i for i in range(5)))
        self.assertEqual(form("v = u").packed_fields, ())

    def test_errors(self):
        for declaration in ("u", "dv1", "g1", "g89", "wdetj", "a, a", "a,", "9sigma", "x-y"):
            with self.subTest(declaration=declaration), self.assertRaises(ValueError):
                form("fields = " + declaration + "\nv = u")
        for body in ("fields = a\nv = other*u", "fields = a\ngeometry = trilinear\nv = a*u"):
            with self.subTest(body=body), self.assertRaises(ValueError):
                form(body)


if __name__ == "__main__":
    unittest.main()
