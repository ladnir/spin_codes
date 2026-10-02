from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import sparse_hill_prefix as prefix


class SparseHillPrefixTests(unittest.TestCase):
    def test_valid_limits(self):
        with TemporaryDirectory() as directory:
            prefix.validate([2, 3, 32], '.1', ['.00032', '.001'], 256, 46, 64,
                Path(directory)/'new.json')

    def test_scope_and_output_rejections(self):
        with TemporaryDirectory() as directory:
            output = Path(directory)/'new.json'
            args = [[2, 3], '.1', ['.001'], 256, 46, 64, output]
            for position, value in ((0, [1]), (0, [2, 2]), (0, [True]), (0, [33]),
                    (1, '.5'), (2, ['.001', '1/1000']), (2, ['0']),
                    (3, 64), (4, True), (5, -1)):
                candidate = args[:]
                candidate[position] = value
                with self.assertRaises(ValueError):
                    prefix.validate(*candidate)
            output.touch()
            with self.assertRaises(ValueError):
                prefix.validate(*args)


if __name__ == '__main__':
    unittest.main()
