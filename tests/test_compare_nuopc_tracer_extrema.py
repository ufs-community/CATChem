import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("compare_nuopc_tracer_extrema.py")
SPEC = importlib.util.spec_from_file_location("compare", SCRIPT)
COMPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(COMPARE)


class ExtremaParserTest(unittest.TestCase):
    def test_parser_and_relative_difference(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "out"
            path.write_text("  0:  so2 max =   0.1542806      min =   0.0000000E+00\n")
            self.assertEqual(COMPARE.parse_extrema(path)["so2"], (0.1542806, 0.0))
        self.assertAlmostEqual(COMPARE.relative_difference(2.0, 1.0), 1.0)
        self.assertEqual(COMPARE.relative_difference(0.0, 0.0), 0.0)
        self.assertEqual(COMPARE.relative_difference(1.0, 0.0), float("inf"))


if __name__ == "__main__":
    unittest.main()
