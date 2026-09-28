"""Compare proposal-only vectorized placement against exact small models."""
from itertools import product
import unittest

import numpy as np
from flint import fmpq, fmpq_mat

from typed_placement import typed_placement
from typed_screen import RESULT_SCOPE,screen_placement


def array(value):
    return np.array([[float(value[i,j]) for j in range(value.ncols())]
                     for i in range(value.nrows())])


def toy(slots):
    return {(a,b):fmpq_mat([[fmpq(1+a,3),fmpq(1+b,5)],
                            [fmpq(a+b,7),fmpq(2+a*b,11)]])
            for a in range(slots+1) for b in range(slots-a+1)}


def exact(operators,counts,epochs,slots,**kwargs):
    return typed_placement(operators,counts,epochs=epochs,slots=slots,
                           matrix=fmpq_mat,rounding=lambda value:value,**kwargs)


class TypedScreenTests(unittest.TestCase):
    def test_noncommuting_small_pairs(self):
        for slots in (1,2,3):
            operators=toy(slots)
            self.assertNotEqual(operators[1,0]*operators[0,1],operators[0,1]*operators[1,0])
            for epochs in (1,2,3):
                for n1 in range(epochs*slots+1):
                    for n2 in range(epochs*slots-n1+1):
                        actual=screen_placement(operators,(n1,n2),epochs=epochs,slots=slots)
                        reference=array(exact(operators,(n1,n2),epochs,slots))
                        np.testing.assert_allclose(actual,reference,rtol=2e-13,atol=1e-15)

    def test_vectorized_grid_matches_exact_grid(self):
        operators=toy(2)
        actual=screen_placement(operators,(2,3),epochs=3,slots=2,return_grid=True)
        reference=exact(operators,(2,3),3,2,return_grid=True)
        self.assertEqual(actual.shape,(3,4,2,2))
        self.assertEqual(set(reference),set(product(range(3),range(4))))
        for counts,value in reference.items():
            np.testing.assert_allclose(actual[counts],array(value),rtol=2e-13,atol=1e-15)
            np.testing.assert_allclose(actual[counts],
                screen_placement(operators,counts,epochs=3,slots=2),rtol=2e-13,atol=1e-15)

    def test_homogeneous_reduction_and_normalization(self):
        base=[fmpq_mat([[1,1],[0,1]]),fmpq_mat([[1,0],[1,1]]),
              fmpq_mat([[2,1],[0,1]])]
        operators={(a,b):base[a+b] for a in range(3) for b in range(3-a)}
        for total in range(7):
            reference=array(exact(operators,(total,0),3,2))
            for first in range(total+1):
                np.testing.assert_allclose(
                    screen_placement(operators,(first,total-first),epochs=3,slots=2),
                    reference,rtol=2e-13,atol=1e-13)
        identities={(a,b):np.eye(2) for a in range(4) for b in range(4-a)}
        np.testing.assert_allclose(screen_placement(identities,(10,12),epochs=12,slots=3),
                                   np.eye(2),rtol=0,atol=2e-13)

    def test_zero_and_forced_occupancy_edges(self):
        operator=fmpq_mat([[3,1],[0,2]])
        for counts,local in (((6,0),(2,0)),((0,6),(0,2)),((0,0),(0,0))):
            np.testing.assert_allclose(screen_placement({local:operator},counts,epochs=3,slots=2),
                                       array(operator**3),rtol=2e-13,atol=1e-13)
        np.testing.assert_array_equal(screen_placement({(0,0):operator},(0,0),epochs=0,slots=2),np.eye(2))
        zeros={key:np.zeros((2,2)) for key in toy(2)}
        np.testing.assert_array_equal(screen_placement(zeros,(2,2),epochs=2,slots=2),np.zeros((2,2)))

    def test_selected_zero_packet_thinning(self):
        # Operator construction averages zero outcomes but leaves the
        # original selected type counts available to the placement driver.
        no_packet=fmpq_mat([[1,1],[0,1]])
        active=fmpq_mat([[1,0],[1,1]])
        first=fmpq(1,3)*active+fmpq(2,3)*no_packet
        second=fmpq(3,4)*active+fmpq(1,4)*no_packet
        operators={(0,0):no_packet,(1,0):first,(0,1):second}
        for counts in ((1,1),(2,1),(1,2)):
            np.testing.assert_allclose(screen_placement(operators,counts,epochs=3,slots=1),
                                       array(exact(operators,counts,3,1)),rtol=2e-13,atol=1e-13)

    def test_inputs_unchanged_and_result_scope(self):
        operators={key:array(value) for key,value in toy(2).items()}
        before={key:value.copy() for key,value in operators.items()}
        result=screen_placement(operators,(2,2),epochs=3,slots=2)
        self.assertEqual(result.dtype,np.float64)
        self.assertIn('proposal only',RESULT_SCOPE.lower())
        self.assertIn('no outward certificate',RESULT_SCOPE.lower())
        for key in operators:
            np.testing.assert_array_equal(operators[key],before[key])

    def test_validation(self):
        operators=toy(2)
        with self.assertRaises(ValueError):
            screen_placement(operators,(5,0),epochs=2,slots=2)
        with self.assertRaises(ValueError):
            screen_placement({(0,0):np.array([[-1.]])},(0,0),epochs=1,slots=2)
        with self.assertRaises(ValueError):
            screen_placement({(0,0):np.array([[np.inf]])},(0,0),epochs=1,slots=2)
        with self.assertRaises(ValueError):
            screen_placement({(0,0):np.eye(2)},(1,1),epochs=1,slots=2)
        with self.assertRaises(ValueError):
            screen_placement({(2,0):np.eye(2)},(4,0),epochs=2,slots=2,return_grid=True)


if __name__=='__main__':
    unittest.main()
