"""Refined families use separate keys and never relabel baseline matrices."""
import unittest
from unittest.mock import patch

import replay


class RefinedMemo(unittest.TestCase):
    def test_refined_hit_skips_all_census_work(self):
        with patch('replay.local_family.source_digest',return_value='base-sources'),\
                patch('replay.refinement_sources',return_value='overlap-sources'),\
                patch('replay.local_family.load',return_value=('base','coefficients')) as load,\
                patch('replay.build') as build,patch('replay.feedback_build') as feedback:
            self.assertEqual(replay.refined_family('unused','.072','.9'),('base','coefficients'))
        build.assert_not_called();feedback.assert_not_called()
        key=load.call_args.args[1]
        self.assertEqual(key['stage'],'independent-row-overlap-family')
        self.assertEqual(key['overlap_sources'],'overlap-sources')
        self.assertFalse(key['density_chord'])

    def test_missing_zero_only_family_does_not_regenerate_baseline(self):
        with patch('replay.local_family.source_digest',return_value='base-sources'),\
                patch('replay.refinement_sources',return_value='overlap-sources'),\
                patch('replay.local_family.load',return_value=None),\
                patch('replay.build') as build:
            with self.assertRaisesRegex(ValueError,'not ready'):
                replay.refined_family('unused','.072','.9')
        build.assert_not_called()

    def test_new_refinement_saved_under_its_own_dependencies(self):
        with patch('replay.local_family.source_digest',return_value='base-sources'),\
                patch('replay.refinement_sources',return_value='overlap-sources'),\
                patch('replay.local_family.load',side_effect=[None,('original','coefficients')]),\
                patch('replay.build',return_value='overlap'),\
                patch('replay.feedback_build',return_value='feedback'),\
                patch('replay.zero_refine',return_value='changed') as refine,\
                patch('replay.local_family.save') as save:
            self.assertEqual(replay.refined_family('unused','.072','.9'),('changed','coefficients'))
        refine.assert_called_once_with('original','overlap','feedback','.072','.9')
        self.assertEqual(save.call_args.args[1]['stage'],'independent-row-overlap-family')
        self.assertEqual(save.call_args.args[2],('changed','coefficients'))

    def test_sources_changing_during_refinement_stop_replay(self):
        with patch('replay.local_family.source_digest',side_effect=['before','after']),\
                patch('replay.refinement_sources',return_value='overlap-sources'),\
                patch('replay.local_family.load',side_effect=[None,('original','coefficients')]),\
                patch('replay.build'),patch('replay.feedback_build'),\
                patch('replay.local_family.save') as save:
            with self.assertRaisesRegex(RuntimeError,'sources changed'):
                replay.refined_family('unused','.072','.9')
        save.assert_not_called()


if __name__=='__main__':unittest.main()
