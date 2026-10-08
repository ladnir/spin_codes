"""Exact conditioning checks and full separate K16/K18 certificate replays."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, prod
import unittest

from field128_size_certificate import (certify, combine, component_metadata,
    dense_interval_passes, geometry, middle_prefactor, modal_reciprocal, TILTS)
from return_cover_certificate import Witness
from field128_density_certificate_tests import exact_packet


class Field128SizeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipts = {size:certify(size) for size in (16,18)}

    def test_exact_geometry_and_class_joins(self):
        expected = {
            16:(768,66816,98304,9830,122,5161,23592,23593),
            18:(3072,267264,393216,39321,280,11797,94371,94372),
        }
        for size,values in expected.items():
            g = geometry(size)
            self.assertEqual(tuple(g[key] for key in ('rows','k','n','threshold',
                'sparse_last','middle_first','middle_last','dense_first')),values)
            self.assertGreaterEqual(42*(g['sparse_last']+1),g['middle_first'])
            self.assertEqual(g['dense_first'],g['middle_last']+1)
            self.assertEqual(F(2**size,g['n']),F(2,3))
            self.assertLess(F(g['threshold'],g['n']),F(1,10))
            self.assertGreaterEqual(F(g['threshold']+1,g['n']),F(1,10))
        with self.assertRaises(ValueError):
            geometry(20)

    def test_exact_modal_minimum(self):
        for rows in range(2,34,2):
            minimum = 1/modal_reciprocal(rows)
            for hits in range(rows+1):
                p = F(hits,rows)
                probability = comb(rows,hits)*p**hits*(1-p)**(rows-hits)
                self.assertGreaterEqual(probability,minimum)
            self.assertEqual(comb(rows,rows//2)*F(1,2)**rows,minimum)
        for rows in (0,1,3,-2):
            with self.assertRaises(ValueError):
                modal_reciprocal(rows)

    def test_modal_scope_excludes_nonmodal_counts(self):
        # Bin(4,1/2) at count0 is NOT the count-conditioned Bernoulli law.
        self.assertLess(F(1,2)**4,1/modal_reciprocal(4))

    def test_small_source_profile_conditioning_comparison(self):
        rows,length,parity,tau = 4,3,1,F(3)
        subsets = [tuple(combinations(range(rows),h)) for h in range(rows+1)]
        for profile in product(range(rows+1),repeat=length):
            value = F(0)
            for columns in product(*(subsets[h] for h in profile)):
                hits = [sum(row in column for column in columns) for row in range(rows)]
                value += prod(1 if h==0 else tau**(h-parity) if h>parity else 0 for h in hits)
            value /= prod(comb(rows,h) for h in profile)
            alpha = [F(h,rows) for h in profile]
            for tilt in (F(1,2),F(2)):
                bracket = prod(1-a for a in alpha)+tau**(-parity)*prod(1-a+tau*tilt*a for a in alpha)
                upper = modal_reciprocal(rows)**length*tilt**(-sum(profile))*bracket**rows
                self.assertLessEqual(value,upper)

    def test_middle_prefactor_exact_small_reference(self):
        witness = Witness(F(1,2),F(1,2),F(2),F(1))
        g = dict(rows=2,n=32,threshold=3)
        matrix = exact_packet(witness)
        terminal = sum(matrix[0][k]*matrix[k][j] for k in range(2) for j in range(2))
        exact = modal_reciprocal(2)**128*terminal/(witness.tau*witness.y**3)
        upper = middle_prefactor(witness,g)
        self.assertLessEqual(F(upper.lo),exact)
        self.assertGreaterEqual(F(upper.hi),exact)

    def test_all_dense_exponent_witnesses_exact(self):
        self.assertEqual(len(TILTS),152)
        for index,tilt in enumerate(TILTS,48):
            self.assertTrue(dense_interval_passes(index,tilt))
        for index,tilt in ((47,F(2)),(200,F(2)),(48,F(1,2))):
            with self.assertRaises(ValueError):
                dense_interval_passes(index,tilt)

    def test_complete_directed_replays(self):
        expected = {
            16:(63,[67,63,1669],28,191,
                '507d1ec6401404b4506e11436d316a0d308fd4acc81eefece4a0d4ed2fe431ab',
                'ca13f18f6dfd644eb71ca7b0e46274874f0498cf127f44f88a6c85137ad71f77'),
            18:(79,[99,79,10260],48,220,
                '1a69bcf89fcd694c299ded95210f9fc9a68cebfc73dcd2d49bad65b8af1e1bd8',
                '83ff5fee53441de52a1d0461fe1e5023017a6c429604ce54ddec631a2aa6c3dc'),
        }
        for size,(bits,components,sparse_cells,middle_cells,sparse_hash,middle_hash) in expected.items():
            r = self.receipts[size]
            self.assertEqual(r['whole_code_margin_bits_lower'],bits)
            self.assertEqual(r['component_margin_bits_lower'],components)
            self.assertEqual(r['sparse']['interval_count'],sparse_cells)
            self.assertEqual(r['middle']['checked_cells'],middle_cells)
            self.assertEqual(r['sparse']['partition_sha256'],sparse_hash)
            self.assertEqual(r['middle']['partition_sha256'],middle_hash)
            self.assertLessEqual(F(r['failure_upper']),F(1,2**bits))
            self.assertTrue(r['whole_code'] and r['rigorous'])
            for component in ('sparse','middle','dense'):
                for key,value in component_metadata(size).items():
                    self.assertEqual(r[component][key],value)

    def test_rejects_corrupt_size_geometry_or_construction(self):
        r = self.receipts[16]
        for field in component_metadata(16):
            sparse = deepcopy(r['sparse'])
            sparse[field] = None
            with self.assertRaises(ValueError):
                combine(16,sparse,r['middle'],r['dense'])
        with self.assertRaises(ValueError):
            combine(16,r['sparse'],self.receipts[18]['middle'],r['dense'])

    def test_rejects_missing_coverage_or_nonrigorous_component(self):
        r = self.receipts[16]
        for field in ('complete','rigorous'):
            middle = deepcopy(r['middle'])
            middle[field] = False
            with self.assertRaises(ValueError):
                combine(16,r['sparse'],middle,r['dense'])
        for part,field,index in (('sparse','active_range',1),('middle','weights',0),('dense','weights',0)):
            components = [deepcopy(r[key]) for key in ('sparse','middle','dense')]
            components[('sparse','middle','dense').index(part)][field][index] += 1
            with self.assertRaises(ValueError):
                combine(16,*components)

    def test_rejects_invalid_probability_bounds(self):
        r = self.receipts[16]
        for bad in ('0','-1','1'):
            middle = deepcopy(r['middle'])
            middle['failure_upper'] = bad
            with self.assertRaises(ArithmeticError):
                combine(16,r['sparse'],middle,r['dense'])
        sparse = deepcopy(r['sparse'])
        sparse['failure_upper_dyadic'] = ['0','1']
        with self.assertRaises(ArithmeticError):
            combine(16,sparse,r['middle'],r['dense'])
        dense = deepcopy(r['dense'])
        dense['failure_upper_dyadic_exponent'] = 0
        with self.assertRaises(ArithmeticError):
            combine(16,r['sparse'],r['middle'],dense)


if __name__=='__main__':
    unittest.main()
