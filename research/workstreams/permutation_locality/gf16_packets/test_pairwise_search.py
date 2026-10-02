import copy
import unittest
from fractions import Fraction as Q
from unittest.mock import patch
from types import SimpleNamespace
import pairwise_search as ps


class PairwiseSearchTests(unittest.TestCase):
    def test_target_aware_retune_stops_at_sufficient_candidate(self):
        model=SimpleNamespace(proposal_stop_bits=66)
        base=(-75,dict(parameters=['1/10','2','3'],regional_count_parts=[{}]))
        with patch.object(ps.regional_count,'propose') as call:
            self.assertEqual(ps.retune(model,(0,0),base),base)
            call.assert_not_called()
        trial=(-67,base[1])
        with patch.object(ps.regional_count,'propose',return_value=trial) as call:
            self.assertEqual(ps.retune(model,(0,0),(-60,base[1])),trial)
            call.assert_called_once()

    def test_retune_actual_winner_and_preserve_old_witness(self):
        witness=dict(parameters=['1/10','2','3'],regional_count_parts=[{}],
                     regional_feedback_classes_through=32)
        original=copy.deepcopy(witness)
        def propose(model,cell,trial):
            self.assertEqual(trial['parameters'][1:],['2','3'])
            self.assertNotIn('regional_feedback_uniform_replace',trial)
            return float((Q(trial['parameters'][0])-Q(1,20))**2),trial
        with patch.object(ps.regional_count,'propose',side_effect=propose):
            score,result=ps.retune(None,(Q(1,32),)*2,(1,witness))
        self.assertEqual(score,0)
        self.assertEqual(result['parameters'][0],'1/20')
        self.assertEqual(witness,original)

    def test_keep_nonregional_and_already_good_bounds(self):
        with patch.object(ps.regional_count,'propose') as call:
            for base in ((1,{}),(-100,{'regional_count_parts':[]})):
                self.assertEqual(ps.retune(None,(0,0),base),base)
            call.assert_not_called()


if __name__=='__main__':unittest.main()
