"""Exact map checks for the proposed small RS outer."""
import unittest
import rs_maps


class RsMapsTests(unittest.TestCase):
    def test_matrix_and_five_product_factor(self):
        self.assertEqual(rs_maps.verify()[4:], (
            (27, 28, 18, 20), (28, 27, 20, 18),
            (18, 20, 27, 28), (20, 18, 28, 27)))

    def test_byte_domain(self):
        for args in ((-1, 0), (0, 256)):
            with self.assertRaises(ValueError):
                rs_maps.multiply(*args)
        with self.assertRaises(ValueError):
            rs_maps.inverse(0)
        with self.assertRaises(ValueError):
            rs_maps.parity_factored((1, 2, 3))


if __name__ == '__main__':
    unittest.main()
