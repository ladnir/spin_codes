"""Exact bookkeeping regressions for BCH row-category family construction."""
from fractions import Fraction as Q
import unittest

from holder_verify import exact_types
from holder_pilot import INTERVALS, families
from row_counts import row_gamma_exact


class ExactTypes(unittest.TestCase):
    def setUp(self):
        self.caps = [0]*257
        for w, mass in ((0,1),(38,3),(62,5),(64,7),(128,11),(192,7),(194,5),(218,3),(256,1)):
            self.caps[w] = mass

    def test_total_type_mass(self):
        for p in (Q(1,4), Q(1,3)):
            rows = exact_types(self.caps, p)
            probabilities = (Q(0),p,Q(1,2),1-p,Q(1))
            gammas = [row_gamma_exact(self.caps,*domain,x) for domain,x in zip(INTERVALS,probabilities)]
            self.assertEqual(sum(h for _,h,_ in rows),sum(gammas)**4-gammas[0]**4)
            for labels,h,law in rows:
                self.assertEqual(sum(law),1)
                self.assertEqual(sum(b*mass for b,mass in enumerate(law)),sum(probabilities[label] for label in labels))
                self.assertTrue(h > 0)

    def test_twelve_families_partition(self):
        rows = exact_types(self.caps,Q(1,4))
        values = families(rows)
        names = ['middle']+[f'exactly_{n}_{side}' for n in range(1,5)
                            for side in ('low_strict','balanced','high')
                            if f'exactly_{n}_{side}' in values]
        self.assertEqual(len(names),12)
        labels = [row[0] for name in names for row in values[name]]
        self.assertEqual(len(labels),69)
        self.assertEqual(len(set(labels)),69)
        self.assertEqual(set(labels),set(row[0] for row in rows))

    def test_endpoint_types(self):
        rows = {labels:(h,law) for labels,h,law in exact_types(self.caps,Q(1,4))}
        self.assertNotIn((0,0,0,0),rows)
        self.assertEqual(rows[(4,4,4,4)],(Q(1),(Q(0),Q(0),Q(0),Q(0),Q(1))))
        self.assertEqual(rows[(0,0,0,4)],(Q(4),(Q(0),Q(1),Q(0),Q(0),Q(0))))


if __name__ == '__main__':
    unittest.main()
