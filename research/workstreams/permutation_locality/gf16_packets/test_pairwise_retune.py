import copy
from fractions import Fraction as Q
from unittest.mock import patch
import unittest
from types import SimpleNamespace
from flint import arb,ctx
import pairwise_retune as pr


class RetuneTests(unittest.TestCase):
    def test_replay_ignores_saved_numerical_success_and_rebuilds(self):
        saved=dict(schema=pr.SCHEMA,ensemble=pr.cover.ENSEMBLE,parameters=pr.cover.PARAMETERS,
            cell=['1/8','1/8'],witness={},upper=[1,-10000],proposal=-10000)
        model=SimpleNamespace(root=(Q(1,100),Q(1)),outward=lambda cell,witness:arb(8))
        prior=ctx.prec
        try:
            ctx.prec=384
            with patch.object(pr.cover,'build_model',return_value=model) as build,patch.object(pr,'validate_geometry'):
                result=pr.replay(saved,384)
            build.assert_called_once_with(384,saved['parameters'])
            m,e=result['upper'];self.assertEqual(Q(m)*Q(2)**e,8)
            self.assertTrue(result['partial_only']);self.assertEqual(result['precision'],384)
        finally:ctx.prec=prior

    def test_replay_rejects_wrong_claim_before_build(self):
        saved=dict(schema=pr.SCHEMA,ensemble=pr.cover.ENSEMBLE,parameters=pr.cover.PARAMETERS,
            cell=['1/8','1/8'],witness={})
        for change in (dict(ensemble='other'),dict(schema='screen'),dict(cell=['1/4','1/8']),dict(witness=None)):
            with patch.object(pr.cover,'build_model') as build,self.assertRaises(ValueError):
                pr.replay(dict(saved,**change),384)
            build.assert_not_called()

    def test_retunes_regional_when_nonregional_baseline_wins(self):
        witness=dict(parameters=['1','2','3'],variance_partition=[{}])
        original=copy.deepcopy(witness);calls=[]
        def propose(model,cell,trial):
            calls.append(copy.deepcopy(trial))
            lam=Q(trial['parameters'][0])
            score=-100 if lam==Q(1,2) and trial.get('regional_joint_return_through') else 20
            return score,dict(trial,regional_count_parts=[{}])
        with patch.object(pr.regional,'propose',side_effect=propose):
            best,trials=pr.search(None,(Q(1,8),)*2,(10,witness),keep=5,factors=[Q(1,2)])
        self.assertEqual(best[0],-100);self.assertEqual(best[1]['parameters'],['1/2','2','3'])
        self.assertEqual(len(trials),10);self.assertEqual(witness,original)
        self.assertNotIn('regional_count_parts',calls[0])
        self.assertTrue(all('regional_count_parts' in row for row in calls[1:]))

    def test_never_discards_better_baseline(self):
        base=(-100,dict(parameters=['1','2','3']))
        with patch.object(pr.regional,'propose',side_effect=lambda m,c,w:(10,w)):
            result,_=pr.search(None,(Q(1,8),)*2,base,keep=1,factors=[Q(2)])
        self.assertEqual(result,base)
        for keep,factors in ((0,[1]),(6,[1]),(True,[1]),(2,[0])):
            with self.assertRaises(ValueError):pr.search(None,(Q(1,8),)*2,base,keep,factors)


if __name__=='__main__':unittest.main()
