from itertools import permutations, product
from math import comb
import unittest
from unittest.mock import patch

import numpy as np

from hd_rankone_mask_screen import packet_matrix, pair_coefficients
from hd_rankone_orbit_screen import (ORBITS,FROZEN_MEDIUM,lump_log_matrix,
    coefficients_from_unpowered,grouped_coefficients,fixed_symmetric_gate)


def logs(values):
    with np.errstate(divide='ignore'):
        return np.log2(values)


class RankOneOrbitTests(unittest.TestCase):
    def test_exact_matching_orbits_not_full_s4(self):
        edges={frozenset((0,1)),frozenset((2,3))}
        group=[p for p in permutations(range(4))
               if {frozenset((p[a],p[b])) for a,b in ((0,1),(2,3))}==edges]
        self.assertEqual(len(group),8)
        actual={frozenset(sum(((mask>>i)&1)<<p[i] for i in range(4))
                          for p in group) for mask in range(16)}
        self.assertEqual(actual,{frozenset(orbit) for orbit in ORBITS[:-1]})
        self.assertNotIn((0,2,1,3),group)
        self.assertEqual(ORBITS[0],(0,))
        self.assertEqual(ORBITS[-1],(16,))

    def test_max_rows_dominate_arbitrary_source_mixtures(self):
        matrix=np.arange(1,290,dtype=float).reshape(17,17)/1000
        grouped=np.exp2(lump_log_matrix(logs(matrix)))
        for i,sources in enumerate(ORBITS):
            for j,targets in enumerate(ORBITS):
                sums=matrix[np.ix_(sources,targets)].sum(axis=1)
                self.assertAlmostEqual(grouped[i,j],max(sums))
                self.assertLessEqual(float(np.mean(sums)),grouped[i,j]+1e-14)

    def test_conditional_normalization_does_not_fractionalize_multiplicity(self):
        count=comb(32,2)*comb(30,1)
        raw={(2,1):logs(count*np.eye(17))}
        result=np.exp2(coefficients_from_unpowered(raw,.125)[2,1])
        np.testing.assert_allclose(result,count*np.eye(7),rtol=1e-14,atol=0)
        self.assertGreater(result[0,0],count**.125)

    def test_actual_conditional_caps_lump_exactly_at_alpha_one(self):
        raw=pair_coefficients(.83,1.,2)
        grouped=coefficients_from_unpowered(raw,1.)
        for key,log_matrix in raw.items():
            matrix=np.exp2(log_matrix)
            reduced=np.exp2(grouped[key])
            for i,sources in enumerate(ORBITS):
                for j,targets in enumerate(ORBITS):
                    sums=matrix[np.ix_(sources,targets)].sum(axis=1)
                    np.testing.assert_allclose(sums,reduced[i,j],rtol=2e-13,atol=0)
        # These distinct types are ordered, not commuted in the comparison.
        a,b=(1,0),(0,1)
        original=np.exp2(raw[a])@np.exp2(raw[b])
        reduced=np.exp2(grouped[a])@np.exp2(grouped[b])
        self.assertAlmostEqual(original[0].sum(),reduced[0].sum(),delta=1e-12)

    def test_two_packet_conditional_jensen_then_path_bound(self):
        alpha=.37
        singleton=[packet_matrix(tuple(int(i==j) for i in range(4)),.73)
                   for j in range(4)]
        full=[packet_matrix(occ,.73) for occ in ((1,1,0,0),(0,0,1,1))]
        values=[(a@b)[0].sum() for a,b in product(singleton,full)]
        conditional=np.mean(np.array(values)**alpha)
        averaged=np.mean(values)**alpha
        h1=np.exp2(lump_log_matrix(logs(np.mean(singleton,axis=0))))
        h2=np.exp2(lump_log_matrix(logs(np.mean(full,axis=0))))
        grouped=(h1@h2)[0].sum()**alpha
        powered=((h1**alpha)@(h2**alpha))[0].sum()
        self.assertLessEqual(conditional,averaged+1e-14)
        self.assertLessEqual(averaged,grouped+1e-14)
        self.assertLessEqual(grouped,powered+1e-14)

    def test_budget_guards_and_fixed_screen_scope(self):
        for z,alpha,maximum in ((0.,.1,2),(.5,0.,2),(.5,.1,65),(.5,.1,True)):
            with self.assertRaises(ValueError):grouped_coefficients(z,alpha,maximum)
        for matrix in (np.zeros((16,16)),np.full((17,17),np.nan),np.full((17,17),np.inf)):
            with self.assertRaises(ValueError):lump_log_matrix(matrix)
        fixture={(0,0):logs(np.eye(7))}
        with patch('hd_rankone_orbit_screen.grouped_coefficients',return_value=fixture) as build:
            result=fixed_symmetric_gate(**FROZEN_MEDIUM)
            build.assert_called_once_with(FROZEN_MEDIUM['z'],FROZEN_MEDIUM['alpha'])
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(result['evaluations'],1)
        self.assertEqual(result['components'],7)
        self.assertEqual(result['marker_log2'],list(FROZEN_MEDIUM['marker_log2']))


if __name__=='__main__':unittest.main()
