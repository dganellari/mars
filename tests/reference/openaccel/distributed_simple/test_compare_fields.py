import contextlib
import io
from pathlib import Path
import tempfile
import unittest

import compare_fields


class FieldComparisonTest(unittest.TestCase):
    header = "node,x,y,z,u,v,w,p\n"
    row = "0,0,0,0,0.1,0,0,1\n"

    def check(self, candidate, expected=1, reference=None, tol="1e-6"):
        with tempfile.TemporaryDirectory() as directory:
            ref, got = Path(directory) / "ref.csv", Path(directory) / "got.csv"
            ref.write_text(self.header + self.row if reference is None else reference)
            got.write_text(candidate)
            with contextlib.redirect_stdout(io.StringIO()):
                result = compare_fields.main([str(ref), str(got), "--tol=" + tol])
            self.assertEqual(result, expected)

    def test_valid(self):
        self.check(self.header + self.row, 0)

    def test_nonfinite_fields_and_coordinates(self):
        for column in range(1, 8):
            for value in ("nan", "inf", "-inf"):
                fields = self.row.strip().split(",")
                fields[column] = value
                with self.subTest(column=column, value=value):
                    self.check(self.header + ",".join(fields) + "\n")

    def test_nonfinite_reference(self):
        self.check(self.header + self.row, reference=self.header + "0,nan,0,0,0.1,0,0,1\n")

    def test_empty(self):
        self.check(self.header, reference=self.header)

    def test_duplicate(self):
        self.check(self.header + "0,0,0,0,100,0,0,1\n" + self.row)

    def test_malformed(self):
        for data in ("", "node,u\n0,1\n", self.header + "0,1\n", self.header + self.row.strip() + ",9\n",
                     self.header + self.row.replace("0,0,0,0", "-1,0,0,0", 1)):
            with self.subTest(data=data):
                self.check(data)

    def test_invalid_tolerance(self):
        for tol in ("nan", "inf", "0", "-1"):
            with self.subTest(tol=tol):
                self.check(self.header + self.row, tol=tol)

    def test_pressure_level_is_not_removed(self):
        self.check(self.header + self.row.replace(",1\n", ",2\n"))

    def test_large_finite_difference(self):
        self.check(self.header + "0,0,0,0,1e308,0,0,1\n")

    def test_node_set_and_coordinates(self):
        self.check(self.header + "1,0,0,0,0.1,0,0,1\n")
        self.check(self.header + "0,0.01,0,0,0.1,0,0,1\n")


if __name__ == "__main__":
    unittest.main()
