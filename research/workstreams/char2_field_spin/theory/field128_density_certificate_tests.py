"""Exact-reference checks for the separate GRS128/87 directed middle replay."""
from decimal import Decimal
from fractions import Fraction as F
from math import comb, log
import unittest

from field128_density_certificate import (DIMENSION, LENGTH, PARITY, WITNESSES,
    cell_upper, certify, integer_span, prefactor, refined_cells, retarget_witness,
    transfer, weight_geometry)
from return_cover_certificate import Witness, interval, row_upper
from state8_density_screen import bound as original_bound
from state8_joint_return import joint_return_bound
from superregular_return_screen import local_matrix
from symmetric_density_screen import bound as generic_bound


def exact_packet(witness):
    result = [[F(0),F(0)],[F(0),F(0)]]
    for w in range(17):
        matrix = local_matrix(w,witness.y,witness.tau,state=8,width=16)
        matrix[1][0] = joint_return_bound(w,witness.y,witness.tau)
        for i in range(2):
            for j in range(2):
                result[i][j] += comb(16,w)*witness.x**w*matrix[i][j]
    return result


class Field128DensityTests(unittest.TestCase):
    def test_construction_geometry_and_class_joins(self):
        g = weight_geometry()
        self.assertEqual((LENGTH,DIMENSION,PARITY),(128,87,41))
        self.assertEqual((g['rows'],g['k'],g['n'],g['threshold']),
                         (12096,1052352,1548288,154828))
        lo,hi = integer_span(g['n'])
        self.assertEqual((lo,hi),(46449,371589))
        self.assertLess(1105*(PARITY+1),lo)
        self.assertGreaterEqual(1106*(PARITY+1),lo)
        self.assertEqual(hi+1,(6*g['n']+24)//25)
        self.assertLess(F(g['threshold'],g['n']),F(1,10))
        self.assertGreaterEqual(F(g['threshold']+1,g['n']),F(1,10))

    def test_frozen_marks_are_rational_and_valid(self):
        self.assertEqual(len(WITNESSES),8)
        for anchor,witness in WITNESSES:
            self.assertIsInstance(anchor,F)
            for name in ('x','y','tau','tilt'):
                self.assertIsInstance(getattr(witness,name),F)
            self.assertGreater(witness.x,0)
            self.assertTrue(0<witness.y<=1 and 1<=witness.tau<=255 and witness.tilt>0)
            self.assertEqual(retarget_witness(anchor,witness,anchor),witness)
            density = F(1,7)
            chosen = retarget_witness(anchor,witness,density)
            self.assertEqual(chosen.tilt*density/(1-density),witness.tilt*anchor/(1-anchor))

    def test_local_transfer_contains_independent_fraction_evaluation(self):
        for witness in (Witness(F(1,3),F(2,3),F(7),F(4,5)),
                        Witness(F(2),F(1,7),F(255),F(1))):
            enclosed,exact = transfer(witness),exact_packet(witness)
            for i in range(2):
                for j in range(2):
                    self.assertLessEqual(F(enclosed[i][j].lo),exact[i][j])
                    self.assertGreaterEqual(F(enclosed[i][j].hi),exact[i][j])

    def test_prefactor_length_state_normalization_and_threshold(self):
        witness = Witness(F(1,2),F(1,2),F(2),F(1))
        # Small exponent checks every prefactor without evaluating huge fractions.
        geometry = dict(rows=2,n=32,threshold=3)
        matrix = exact_packet(witness)
        terminal = sum(matrix[0][k]*matrix[k][j] for k in range(2) for j in range(2))
        expected = F(3)**128*terminal/(witness.tau*witness.y**3)
        result = prefactor(witness,geometry)
        self.assertLessEqual(F(result.lo),expected)
        self.assertGreaterEqual(F(result.hi),expected)

    def test_row_parity_and_both_cell_endpoint_choices(self):
        for tilt in (F(1,10),F(1),F(3)):
            for x in (F(1,10),F(2)):
                witness = Witness(x,F(1,2),F(2),tilt)
                envelope = row_upper(F(10,100),F(13,100),witness,128,41)
                total = F(0)
                for weight in range(10,14):
                    alpha = F(weight,100)
                    exact = (1-alpha)**128+F(2)**(-41)*(1-alpha+2*tilt*alpha)**128
                    self.assertLessEqual(exact,F(envelope.hi))
                    total += exact**2/(x*tilt)**weight
                self.assertLessEqual(total,F(cell_upper(10,13,100,2,witness,interval(1)).hi))

    def test_existing_and_generic_full_float_objectives_agree(self):
        g = weight_geometry()
        for anchor,witness in WITNESSES:
            point = tuple(log(float(getattr(witness,name))) for name in ('x','y','tau','tilt'))
            weight = round(anchor*g['n'])
            original = original_bound(point,weight,20,128,87,16,True)
            generic = generic_bound(point,weight,20,128,87,8,F(1,10),True)
            self.assertAlmostEqual(original,generic,places=8)

    def test_gap_free_refinement_and_rejected_singletons(self):
        evaluate = lambda lo,hi: (0,interval(hi-lo+1))
        cells = list(refined_cells(3,22,9,evaluate,Decimal(3)))
        self.assertEqual([w for lo,hi,_,_ in cells for w in range(lo,hi+1)],list(range(3,23)))
        with self.assertRaises(ArithmeticError):
            list(refined_cells(1,1,1,evaluate,Decimal('0.5')))
        with self.assertRaises(ValueError):
            certify(k_log=18)
        with self.assertRaises(ValueError):
            certify(witnesses=((F(3,100),Witness(.5,F(1,2),F(2),F(1))),))

    def test_complete_directed_replay(self):
        result = certify()
        self.assertEqual(result['covered_weights'],325141)
        self.assertEqual(result['checked_cells'],210)
        self.assertEqual(result['partition_sha256'],
                         'ea5af80dd068950da774d2daf8ed3cd5577588120652a17b7dc817952204fe83')
        self.assertGreaterEqual(result['class_margin_bits_lower'],408)
        self.assertLessEqual(F(result['failure_upper']),F(1,2**408))
        self.assertTrue(result['rigorous'])
        self.assertFalse(result['whole_code'])


if __name__ == '__main__':
    unittest.main()
