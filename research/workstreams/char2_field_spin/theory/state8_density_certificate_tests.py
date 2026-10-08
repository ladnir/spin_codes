"""Exact-arithmetic regressions for the complete middle-weight replay."""
from decimal import Decimal
from fractions import Fraction as F
from math import comb
import unittest

from return_cover_certificate import Witness, interval, row_upper, terminal_power
from state8_density_certificate import (WITNESSES, cell_upper, certify, integer_span,
    interval_minimum, refined_cells, retarget_witness, transfer, weight_geometry)
from state8_joint_return import joint_return_bound
from subfield_certificate import Interval
from superregular_return_screen import local_matrix


def contains(enclosure, exact):
    return F(enclosure.lo) <= exact <= F(enclosure.hi)


def exact_packet(witness):
    result = [[F(0), F(0)], [F(0), F(0)]]
    for w in range(17):
        packet = local_matrix(w, witness.y, witness.tau)
        packet[1][0] = joint_return_bound(w, witness.y, witness.tau)
        coefficient = comb(16, w)*witness.x**w
        for i in range(2):
            for j in range(2):
                result[i][j] += coefficient*packet[i][j]
    return result


class DensityCertificateTests(unittest.TestCase):
    def test_packet_contains_independent_fraction_reference(self):
        for witness in (Witness(F(1, 3), F(2, 3), F(7), F(4, 5)),
                        Witness(F(2), F(1, 7), F(255), F(1)),
                        Witness(F(1), F(1), F(1), F(1))):
            enclosure, reference = transfer(witness), exact_packet(witness)
            for i in range(2):
                for j in range(2):
                    self.assertTrue(contains(enclosure[i][j], reference[i][j]))

    def test_complete_minimum_enclosure(self):
        first = Interval(Decimal('0.1'), Decimal('0.3'))
        second = Interval(Decimal('0.2'), Decimal('0.25'))
        result = interval_minimum(first, second)
        for a in (F(1, 10), F(1, 5), F(3, 10)):
            for b in (F(1, 5), F(9, 40), F(1, 4)):
                self.assertTrue(contains(result, min(a, b)))

    def test_terminal_power_exact(self):
        matrix = [[F(2, 3), F(1, 5)], [F(1, 7), F(4, 5)]]
        enclosed = [[interval(x) for x in row] for row in matrix]
        vector = [F(1), F(0)]
        for exponent in range(20):
            self.assertTrue(contains(terminal_power(enclosed, exponent), sum(vector)))
            vector = [sum(vector[i]*matrix[i][j] for i in range(2)) for j in range(2)]

    def test_row_and_weight_endpoints(self):
        # Exercise both monotonicity choices, not just the saved witnesses.
        for tilt in (F(1, 10), F(1), F(3)):
            for x in (F(1, 10), F(2)):
                witness = Witness(x, F(1, 2), F(2), tilt)
                row = row_upper(F(10, 100), F(13, 100), witness, 256, 80)
                total = F(0)
                for weight in range(10, 14):
                    alpha = F(weight, 100)
                    exact = (1-alpha)**256+F(2)**(-80)*(1-alpha+2*tilt*alpha)**256
                    self.assertLessEqual(exact, F(row.hi))
                    total += exact**2/(x*tilt)**weight
                bound = cell_upper(10, 13, 100, 2, witness, interval(1))
                self.assertLessEqual(total, F(bound.hi))

    def test_rational_odds_retargeting(self):
        for anchor, witness in WITNESSES:
            self.assertEqual(retarget_witness(anchor, witness, anchor), witness)
            for density in (F(1, 20), F(1, 7), F(1, 5)):
                selected = retarget_witness(anchor, witness, density)
                self.assertEqual(selected.tilt*density/(1-density),
                                 witness.tilt*anchor/(1-anchor))
                self.assertEqual((selected.x, selected.y, selected.tau),
                                 (witness.x, witness.y, witness.tau))
        with self.assertRaises(ValueError):
            retarget_witness(F(0), WITNESSES[0][1], F(1, 2))

    def test_gap_free_refinement(self):
        evaluate = lambda lo, hi: (0, interval(hi-lo+1))
        cells = list(refined_cells(3, 22, 9, evaluate, Decimal(3)))
        self.assertEqual([w for lo, hi, _, _ in cells for w in range(lo, hi+1)],
                         list(range(3, 23)))
        self.assertTrue(all(hi-lo+1 <= 3 for lo, hi, _, _ in cells))
        with self.assertRaises(ArithmeticError):
            list(refined_cells(1, 1, 1, evaluate, Decimal('0.5')))

    def test_exact_geometry(self):
        geometry = weight_geometry()
        self.assertEqual((geometry['rows'], geometry['n'], geometry['threshold']),
                         (6016, 1540096, 154009))
        self.assertEqual(integer_span(geometry['n']), (46203, 369623))
        self.assertGreaterEqual(601*81, integer_span(geometry['n'])[0])

    def test_complete_directed_replay(self):
        result = certify()
        self.assertEqual(result['covered_weights'], 323421)
        self.assertEqual(result['checked_cells'], 216)
        self.assertEqual(result['partition_sha256'],
                         '01c230076663d1bff17808ee26577d6cc1304f2761a0e599a2e124e072abecff')
        self.assertGreaterEqual(result['class_margin_bits_lower'], 344)
        self.assertLessEqual(F(result['failure_upper']), F(1, 2**344))
        self.assertFalse(result['whole_code'])


if __name__ == '__main__':
    unittest.main()
