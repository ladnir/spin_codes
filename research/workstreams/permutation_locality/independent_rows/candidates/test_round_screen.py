"""Update count must reach every round-dependent four-bit local bound."""
from contextlib import ExitStack
import unittest
from unittest.mock import patch

import mass_density_screen as screen


class RoundScreen(unittest.TestCase):
    def setUp(self):self.precision=screen.ctx.prec

    def tearDown(self):screen.ctx.prec=self.precision

    def test_reject_invalid_count_before_building(self):
        for rounds in (0,33,True,2.5):
            with self.assertRaises(ValueError):
                screen.epoch_grid(['.056'],'.9',rounds=rounds)

    def test_propagate_count_and_preserve_default(self):
        for rounds in (None,2,3):
            expected=2 if rounds is None else rounds
            with ExitStack() as stack:
                helpers={}
                for name in ('prepare','local_data','fresh_census','tail_census',
                             'pair_census','zero_census','averages','fresh_refine',
                             'window_refine','multi_refine','collision_refine',
                             'zero_refine','lift','transform','tail_test'):
                    helpers[name]=stack.enter_context(patch.object(screen.baseline,name,return_value=[]))
                stack.enter_context(patch.object(screen.baseline,'prepare_inputs',return_value=(0,0,0,{})))
                for module,names in (
                    (screen.full_feedback_refinement,('census','check_pairs','refine')),
                    (screen.window_histogram,('census','moments','refine')),
                    (screen.cancellation_joint,('census','check_fresh','check_feedback','refine')),
                    (screen.feedback_density,('build',)),
                    (screen.universal_density,('refine',)),
                ):
                    for name in names:
                        helpers[module.__name__+'.'+name]=stack.enter_context(patch.object(module,name,return_value=[]))
                kwargs={} if rounds is None else dict(rounds=rounds)
                result,_=screen.epoch_grid(['.056'],'.9',**kwargs)
                self.assertEqual(set(result),{'.056'})
                for name,index in (
                    ('transform',3),('window_histogram.refine',4),
                    ('full_feedback_refinement.refine',5),('cancellation_joint.refine',4),
                    ('feedback_density.build',3),('tail_test',4),
                ):
                    self.assertEqual(helpers[name].call_args.args[index],expected,name)
                self.assertEqual(helpers['universal_density.refine'].call_args.kwargs['rounds'],expected)


if __name__=='__main__':unittest.main()
