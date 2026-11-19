"""Check source composition, coefficient identity and the single element schedule."""
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from marsir import parse_spec, parse_spec_file, synthesize
from marsir.backends.mlir_ir import emit_galerkin

SPECS = Path(__file__).resolve().parents[1] / "specs"


def body(ir):
    return [line for line in ir.splitlines() if not line.lstrip().startswith("//")]


class Composition(unittest.TestCase):
    def test_fused_kernel(self):
        direct = parse_spec("""operator reaction {
          element = hex
          form = galerkin
          fields = sigma_w
          v = sigma_w*u
          dv0 = g0*du0 + g1*du1 + g2*du2
          dv1 = g1*du0 + g3*du1 + g4*du2
          dv2 = g2*du0 + g4*du1 + g5*du2
        }""")
        fused = synthesize(parse_spec_file(SPECS / "reaction.op"))
        self.assertEqual(fused.packed_fields, ("g0", "g1", "g2", "g3", "g4", "g5", "sigma_w"))
        for p in (1, 3, 7):
            ir = emit_galerkin(fused, p=p)
            self.assertEqual(body(ir), body(emit_galerkin(synthesize(direct), p=p)))
            self.assertEqual(ir.count(":4 = mir.flux"), 1)

    def test_overlapping_outputs_and_shared_fields(self):
        def read(source):
            return parse_spec_file(SPECS / (source + ".op"))
        op = parse_spec("operator total = bp1(g0=a) + bp1(g0=b) + bp1(g0=a)", read)
        self.assertEqual(synthesize(op).packed_fields, ("a", "b"))
        self.assertEqual(op.outputs[0][1].render({}), "(((a * u) + (b * u)) + (a * u))")

    def test_rejects_incompatible_operands(self):
        def read(source):
            return parse_spec_file(SPECS / (source + ".op"))
        for expr, error in [
                ("bp3 + bp1", "ambiguous stored slots"),
                ("bp3 + bp1(g1=sigma)", "unused or repeated"),
                ("bp3 + bp1(g0=a,g0=b)", "unused or repeated"),
                ("bp3 + bp1(g0=du1)", "reserved"),
                ("bp3 + bp1(g0=g9)", "reserved"),
                ("bp3 + bp1()", "expected a stored slot"),
                ("bp3 + bp3_tri", "stored-geometry"),
                ("laplacian + bp1", "Galerkin"),
                ("bp1 + ../bp1", "invalid operator reference")]:
            with self.subTest(expr=expr), self.assertRaisesRegex(ValueError, error):
                parse_spec("operator bad = " + expr, read)

    def test_nested_and_broken_references(self):
        with TemporaryDirectory(dir=Path(__file__).resolve().parent) as tmp:
            root = Path(tmp)
            (root / "bp1.op").write_bytes((SPECS / "bp1.op").read_bytes())
            (root / "first.op").write_text("operator first = bp1(g0=a) + bp1(g0=b)")
            top = root / "top.op"
            top.write_text("operator top = first + bp1(g0=c)")
            self.assertEqual(synthesize(parse_spec_file(top)).packed_fields, ("a", "b", "c"))
            (root / "first.op").write_text("operator first = top + bp1")
            with self.assertRaisesRegex(ValueError, "top -> first -> top"):
                parse_spec_file(top)
            (root / "first.op").write_text("operator first = bp1 + missing")
            with self.assertRaisesRegex(ValueError, "missing.op"):
                parse_spec_file(top)
            (root / "first.op").write_bytes((SPECS / "bp1.op").read_bytes())
            with self.assertRaisesRegex(ValueError, "declares"):
                parse_spec_file(top)
            (root / "first.op").write_text("operator first {\nelement = hex\nform = galerkin\nv = c*u\n}")
            with self.assertRaisesRegex(ValueError, "unknown name"):
                parse_spec_file(top)


if __name__ == "__main__":
    unittest.main()
