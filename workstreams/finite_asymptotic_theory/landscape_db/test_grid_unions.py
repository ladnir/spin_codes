import math
import unittest

from summarize_grid_unions import best_cover


def bound(lo,hi,value,event=None,bits=None):
    return dict(lo=lo,hi=hi,log_bound=math.log(value),result_id=f'{lo}-{hi}-{event}',
                setup_event_id=event,setup_failure_bits=bits)


class GridUnionTest(unittest.TestCase):
    def test_gap_is_not_a_full_bound(self):
        self.assertIsNone(best_cover(8,[bound(1,3,.1),bound(5,8,.2)]))

    def test_overlapping_interval_is_charged_in_full(self):
        result=best_cover(8,[bound(1,4,.1),bound(3,8,.2)])
        self.assertAlmostEqual(math.exp(result['log_bound']),.3)

    def test_shared_setup_event_is_charged_once(self):
        result=best_cover(8,[bound(1,1,.01),bound(2,4,.02,'caps',4),bound(5,8,.03,'caps',4)])
        self.assertAlmostEqual(math.exp(result['log_bound']),.06+2**-4)

    def test_distinct_setup_events_cannot_be_spliced(self):
        rows=[bound(1,1,.01),bound(2,4,.02,'first',4),bound(5,8,.03,'second',4)]
        self.assertIsNone(best_cover(8,rows))

    def test_unconditional_cover_avoids_unnecessary_event_charge(self):
        rows=[bound(1,8,.1),bound(1,8,.01,'caps',1)]
        result=best_cover(8,rows)
        self.assertAlmostEqual(math.exp(result['log_bound']),.1)
        self.assertIsNone(result['setup_event_id'])

    def test_checked_event_implication_allows_one_stronger_charge(self):
        rows=[bound(1,1,.01),bound(2,4,.02,'weak',4),bound(5,8,.03,'strong',5)]
        result=best_cover(8,rows,{('strong','weak'):True},{'strong':5,'weak':4})
        self.assertAlmostEqual(math.exp(result['log_bound']),.06+2**-5)
        self.assertEqual(result['setup_event_id'],'strong')

    def test_spectrum_containment_checks_every_shell_and_parameters(self):
        from spectrum_events import implies
        weak=dict(block_bits=4,dimension=2,counts={1:3,2:3,3:3,4:1})
        strong=dict(block_bits=4,dimension=2,counts={2:2,4:1})
        self.assertTrue(implies(strong,weak))
        self.assertFalse(implies(weak,strong))
        self.assertFalse(implies(strong,dict(weak,dimension=3)))


if __name__=='__main__':unittest.main()
