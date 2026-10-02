import copy
import unittest
from fractions import Fraction as Q
from unittest.mock import patch
from flint import arb
import pairwise_cover as pc


class PairwiseCoverTests(unittest.TestCase):
    def record(self):
        return dict(schema=pc.SCHEMA,ensemble=pc.ENSEMBLE,
                    parameters=copy.deepcopy(pc.PARAMETERS),
                    leaves={'0':{'witness':{}}},unresolved={'1':{}})

    def test_claim_is_fixed(self):
        pc.validate_record(self.record())
        for key,value in (('threshold',209714),('updates',3),('minimum_groups',50),('mass_bits',259)):
            record=self.record();record['parameters'][key]=value
            with self.assertRaises(ValueError):pc.validate_record(record)
        record=self.record();record['ensemble']='independent4-gf16-r4'
        with self.assertRaises(ValueError):pc.validate_record(record)
        record=self.record();record['parameters']['mass_bits']=260
        pc.validate_record(record)
        record['parameters'].pop('prune_mixture')
        pc.validate_record(record)  # Earlier unpruned checkpoint schema.
        record['parameters']['prune_mixture']=1
        with self.assertRaises(ValueError):pc.validate_record(record)

    def test_union_objective_is_part_of_the_comparison_model(self):
        record=self.record();record['parameters'].pop('union_cost_tilt')
        pc.validate_record(record)  # Existing pairwise checkpoints remain replayable.
        record['parameters']['union_cost_tilt']='3/16'
        with self.assertRaises(ValueError):pc.validate_record(record)
        record['parameters']['prune_mixture']=True
        pc.validate_record(record)
        record['parameters']['union_cost_tilt']='1/4'
        with self.assertRaises(ValueError):pc.validate_record(record)

    def test_replay_never_accepts_a_hole_or_trusts_saved_scores(self):
        class Model:
            root=(Q(0),Q(1))
            @staticmethod
            def outward(cell,witness):return arb(2)**-100
        record=self.record();record['leaves']['0']['proposal']=-100000
        result=pc.replay(Model(),record)
        self.assertFalse(result['dense_complete'])
        record['leaves']['1']={'witness':{}};record['unresolved']={}
        result=pc.replay(Model(),record)
        self.assertTrue(result['dense_complete'])
        with patch.object(Model,'outward',return_value=arb(1)):
            self.assertFalse(pc.replay(Model(),record)['dense_complete'])


if __name__=='__main__':unittest.main()
