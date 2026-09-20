import unittest

import audit_imported_rm2sub as audit


class ConstituentAuditTests(unittest.TestCase):
    def test_literal_array(self):
        self.assertEqual(audit.array('Columns{0x01, 3,};', 'Columns', 2), [1, 3])
        for text in ('Columns{1,};', 'Columns{1,2}; Columns{1,2};', 'Columns{1, f()};'):
            with self.assertRaises(ValueError):
                audit.array(text, 'Columns', 2)

    def test_exact_spectrum_and_dual(self):
        rows = audit.generators([1, 2, 3], 2)
        self.assertEqual(rows, [5, 6])
        self.assertEqual(audit.spectrum(rows), {0: 1, 2: 3})
        self.assertEqual(audit.dual_spectrum({0: 1, 2: 3}, 3, 2), {0: 1, 3: 1})

    def test_invalid_generators(self):
        for columns in ([0], [4], [-1]):
            with self.assertRaises(ValueError):
                audit.generators(columns, 2)
        with self.assertRaises(ValueError):
            audit.spectrum([3, 3])

    def test_same_spectrum_does_not_identify_ordered_map(self):
        left = audit.generators([1, 2, 3], 2)
        right = audit.generators([3, 2, 1], 2)
        self.assertNotEqual(left, right)
        self.assertEqual(audit.spectrum(left), audit.spectrum(right))


if __name__ == '__main__':
    unittest.main()
