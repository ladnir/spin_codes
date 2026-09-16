from fractions import Fraction as F
import unittest
import dense_ladder as dense


class GeometryTests(unittest.TestCase):
    def test_reshaped_tree_covers_new_root(self):
        root = dense.geometry.box(512, 8192, F(0), F(1), {'example': 1})
        a, b = dense.geometry.children(root, 'q')
        a0, a1 = dense.geometry.children(a, 'v')
        seed = dict(minimum=512, instance={'outer_rows': 8192},
                    leaves={'00': a0, '01': a1, '1': b}, splits={'': 'q', '0': 'v'})
        for rows in (32768, 131072):
            leaves, splits = dense.reshape(seed, 512, rows)
            cover = dense.geometry.check_partition(leaves, splits, 512, rows)
            self.assertTrue(cover['complete'])
            self.assertTrue(all('power' not in node for node in leaves.values()))
            self.assertEqual(leaves['1']['hi'], rows)

    def test_reject_missing_branch(self):
        seed = dict(minimum=512, instance={'outer_rows': 8192}, leaves={}, splits={})
        with self.assertRaises(Exception):
            dense.reshape(seed, 512, 32768)


if __name__ == '__main__':
    unittest.main()
