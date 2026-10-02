import copy
from fractions import Fraction as Q
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from flint import arb,arb_mat,ctx
import shared_relaxed_strategy_obstruction as diagnostic


class ObstructionTests(unittest.TestCase):
    def setUp(self):
        self.precision = ctx.prec; ctx.prec = 256
        self.cell = (Q(1,4),Q(1,4))
        self.model = SimpleNamespace(tilt=Q(1),threshold=0,q_min=1,
            features=[Q(0),Q(1,4),Q(3,4)],active=[0,1,1],
            components=[(Q(1),Q(0),0),(Q(2),Q(1,4),1),(Q(3),Q(3,4),1)],
            family=lambda _ : ([Q(1),Q(2),Q(3)],None,[0.,1.,2.]),
            weights=lambda *_ : ([Q(3,4),Q(1,4)],None))
        parts = [dict(interval=['0','1/8'],dual=['0','0','0']),
                 dict(interval=['1/8','1/4'],dual=['0','0','0'])]
        self.witness = dict(tilt='1',parameters=['1','0','0'],variance_dual=['0','0'],
            regional_direct_counts=True,variance_partition=copy.deepcopy(parts),
            regional_count_parts=[dict(copy.deepcopy(p),mgf_witnesses=[dict(tilt='0',dual=['0','0','0'])])
                                  for p in parts])

    def tearDown(self):ctx.prec = self.precision

    def test_refinement_preserves_source_and_covers_exact_same_variance_range(self):
        original = copy.deepcopy(self.witness)
        with patch.object(diagnostic.variance,'outer_witness',return_value=(0,[Q(2),Q(0),Q(3)])) as outer, \
                patch.object(diagnostic.regional,'propose_mgf',return_value=[dict(tilt='1',dual=['0','0','0'])]) as mgf:
            refined = diagnostic.refine(self.model,self.cell,self.witness,[1],4)
        self.assertEqual(self.witness,original)
        parts = diagnostic.variance.validate(self.cell,refined['variance_partition'])
        self.assertEqual(len(parts),5)
        self.assertEqual(parts[0],((Q(0),Q(1,8)),(Q(0),Q(0),Q(0))))
        for index,((lo,hi),dual) in enumerate(parts[1:]):
            self.assertEqual((lo,hi),(Q(1,8)+Q(index,32),Q(1,8)+Q(index+1,32)))
            self.assertEqual(dual,(Q(2),Q(0),Q(3)))
        self.assertEqual(outer.call_count,4); self.assertEqual(mgf.call_count,4)
        for part,family in zip(refined['variance_partition'],refined['regional_count_parts']):
            self.assertEqual(part['interval'],family['interval'])
            self.assertEqual(part['dual'],family['dual'])
        self.assertEqual(refined['parameters'],original['parameters'])
        self.assertNotIn('upper',refined)

    def test_invalid_refinement_selection_rejected(self):
        for indices,splits in (([],4),([0,0],4),([2],4),([True],4),([0],True),([0],0),([0],17)):
            with self.subTest(indices=indices,splits=splits),self.assertRaises(ValueError):
                diagnostic.refine(self.model,self.cell,self.witness,indices,splits)

    def test_decomposition_sum_matches_hand_computed_matrix_moment(self):
        parts = diagnostic.variance.validate(self.cell,self.witness['variance_partition'])
        cache = SimpleNamespace(polynomial=lambda *_:[arb_mat([[j]]) for j in (1,2,3)])
        with patch.object(diagnostic.regional.sc,'G',2), \
                patch.object(diagnostic.regional.sc,'REGIONS',2), \
                patch.object(diagnostic.regional.sc,'PACKETS',4), \
                patch.object(diagnostic.regional,'prepare_witness',return_value=(parts,self.witness)), \
                patch.object(diagnostic.variance,'factor',return_value=Q(1)), \
                patch.object(diagnostic.regional,'count_mass_caps',return_value=[arb(1)/2,arb(1)/4,arb(1)/4]):
            upper,rows,checked = diagnostic.decompose(self.model,self.cell,self.witness,cache)
        # Each bin: (1+2+3)^2 * (1/2+2/4+3/4)^2 = 441/4.
        self.assertTrue(upper >= arb(441)/2)
        self.assertTrue(upper-arb(441)/2 < arb(2)**-200)
        self.assertEqual(len(rows),2)
        self.assertEqual([r['index'] for r in rows],[0,1])
        self.assertEqual(rows[0]['interval'],['0','1/8'])
        self.assertIs(checked,self.witness)
        self.assertAlmostEqual(sum(r['count'] for r in rows[0]['groups']),2.)


if __name__ == '__main__':unittest.main()
