"""Bounded exact checks of MDS-aware component/source counting."""

from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import unittest

import hd_mds_component as proof
from hd_affine_proof import mul
from hd_component_rank import component_moment,exact_kernels


def direct_forward(cuts,rows,dimension,kernels):
    """GRS words and all prefix diagonals, independently of source-law weights."""
    length,active = len(cuts),len(cuts[0])
    if dimension not in (1,2) or length > 4:
        raise ValueError("bounded GF4 GRS fixture required")
    code = []
    for coefficients in product(range(4),repeat=dimension):
        if any(coefficients):
            code.append(tuple(coefficients[0]^(mul(coefficients[1],x) if dimension==2 else 0)
                              for x in range(length)))
    @lru_cache(None)
    def post(occupancy):
        vector = [F(1)]+[F(0)]*(len(kernels[0])-1)
        for weight in occupancy:
            matrix = kernels[weight]
            vector = [sum(vector[i]*matrix[i][j] for i in range(len(vector)))
                      for j in range(len(vector))]
        return sum(vector)
    total = F(0)
    for words in product(code,repeat=active):
        for diagonals in product(range(1,4),repeat=active*length):
            carry = 0
            occupancy = []
            for column,positions in enumerate(cuts):
                increments = {j:mul(diagonals[r*length+column],words[r][column])
                              for r,j in positions}
                weight = 0
                for j in range(rows):
                    carry ^= increments.get(j,0)
                    weight += bool(carry)
                occupancy.append(weight)
            total += post(tuple(occupancy))
    return total/3**(active*length)


class MdsComponentTests(unittest.TestCase):
    def test_cut_degree_counts_every_field_value(self):
        for old,new in product((0,1),repeat=2):
            for previous in ((0,) if not old else (1,2,3)):
                counts = [0,0]
                for value in range(4):
                    if bool(value)==bool(new):
                        counts[value!=previous] += 1
                expected = [0,0]
                for degree,count in proof.cut_choices(old,new,4):
                    expected[degree] += count
                self.assertEqual(counts,expected)

    def test_counter_transfer_matches_field_assignments_and_forward_encoder(self):
        cases = (
            (4,2,[[(0,j)] for j in (0,1,2,3)]),
            (3,1,[[(0,0),(1,1)],[(1,0),(0,2)],[(0,1),(1,2)]]),
        )
        for rows,dimension,cuts in cases:
            kernels = exact_kernels(rows,4,F(1,4))
            exact = proof.source_moment(cuts,rows,4,dimension,kernels.__getitem__)
            self.assertEqual(exact,proof.direct_value_moment(cuts,rows,4,dimension,kernels))
            self.assertEqual(exact,direct_forward(cuts,rows,dimension,kernels))

    def test_saturated_counter_equals_exact_feasibility(self):
        rows,dimension = 3,1
        cuts = [[(0,0),(1,1)],[(1,0),(0,2)],[(0,1),(1,2)]]
        kernels = exact_kernels(rows,4,F(1,2))
        degrees = proof.degree_distribution(cuts,rows,4,kernels.__getitem__)
        exact = sum(value for degree,value in degrees.items() if all(d>=3 for d in degree))/3**4
        self.assertEqual(exact,proof.source_moment(cuts,rows,4,dimension,kernels.__getitem__,True))

    def test_positive_bounds_and_unfiltered_bridge(self):
        cases = (
            (4,2,[[(0,0)]]*4),
            (3,1,[[(0,0),(1,1)],[(1,0),(0,2)],[(0,1),(1,2)]]),
        )
        for rows,dimension,cuts in cases:
            for mark in (F(1,4),F(1,2),F(1)):
                kernels = exact_kernels(rows,4,mark)
                exact = proof.source_moment(cuts,rows,4,dimension,kernels.__getitem__)
                feasible = proof.source_moment(cuts,rows,4,dimension,kernels.__getitem__,True)
                positions = [[j for _,j in sorted(c,key=lambda item:item[1])] for c in cuts]
                relaxed = component_moment(positions,rows,4,kernels.__getitem__)/3**(
                    (len(cuts)-dimension)*len(cuts[0]))
                self.assertLessEqual(exact,feasible)
                self.assertLessEqual(feasible,relaxed)
                self.assertEqual(relaxed,proof.marked_upper(cuts,rows,4,dimension,kernels.__getitem__,1))
                for degree_mark in (2,4):
                    self.assertLessEqual(feasible,proof.marked_upper(
                        cuts,rows,4,dimension,kernels.__getitem__,degree_mark))

    def test_mark_one_has_exact_active_codeword_count(self):
        for dimension,cuts in ((2,[[(0,1)]]*4),
                               (1,[[(0,0),(1,2)]]*3)):
            kernels = exact_kernels(4,4,F(1))
            exact = proof.source_moment(cuts,4,4,dimension,kernels.__getitem__)
            self.assertEqual(exact,(4**dimension-1)**len(cuts[0]))
            unweighted = proof.degree_distribution(cuts,4,4,kernels.__getitem__)
            self.assertEqual(sum(unweighted.values()),4**(len(cuts)*len(cuts[0])))

    def test_directed_arb_kernels_preserve_exact_bounds(self):
        from flint import arb,ctx
        from hd_component_rank import arb_fraction,run_kernel_upper
        previous = ctx.prec
        ctx.prec = 160
        try:
            cuts = [[(0,0),(1,1)],[(1,0),(0,2)],[(0,1),(1,2)]]
            kernels = exact_kernels(3,16,F(9,10),collapsed=True)
            upper = lambda h:run_kernel_upper(3,h,16,F(9,10))
            tolerance = arb(2)**-120
            for feasible in (False,True):
                exact = proof.source_moment(cuts,3,16,1,kernels.__getitem__,feasible)
                bounded = proof.source_moment(cuts,3,16,1,upper,feasible)
                self.assertTrue(bounded+tolerance>arb_fraction(exact))
            for mark in (1,F(3,2),4):
                exact = proof.marked_upper(cuts,3,16,1,kernels.__getitem__,mark)
                bounded = proof.marked_upper(cuts,3,16,1,upper,mark)
                self.assertTrue(bounded+tolerance>arb_fraction(exact))
        finally:
            ctx.prec = previous


if __name__ == '__main__':
    unittest.main()
