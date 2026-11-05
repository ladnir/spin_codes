from fractions import Fraction as F
from itertools import product
import unittest

from flint import arb,ctx

from hd_shifted_cover_backend import histogram_vertices,profile_count,CachedBound,CONSTRUCTION
from hd_shifted_cover_ledger import Box,Ledger,refine,replay_cover
from hd_shifted_pairs_count import exact_profiles


def points(box):
    return set(product(*(range(lo,hi+1) for lo,hi in box.ranges)))


class BackendTests(unittest.TestCase):
    def test_integer_count_matches_every_small_box(self):
        for pairs in range(1,4):
            profiles = set(exact_profiles(pairs))-{(0,0,0)}
            intervals = lambda n,start=0:[(lo,hi) for lo in range(start,n+1) for hi in range(lo,n+1)]
            for ranges in product(intervals(2*pairs,1),intervals(pairs),intervals(pairs)):
                exact = sum(all(lo<=n<=hi for n,(lo,hi) in zip(point,ranges)) for point in profiles)
                self.assertEqual(profile_count(ranges,pairs),exact)

    def test_clipped_vertices_bound_convex_function(self):
        pairs = 5
        domains = (((1,10),(0,5),(0,5)),((3,7),(0,1),(2,5)),
                   ((3,3),(0,3),(0,3)),((5,8),(0,4),(1,3)))
        def objective(point):
            r,a,b = point
            hist = (pairs-r+a,r-2*a,a,pairs-r+b,r-2*b,b)
            return sum((i+1)*n*n for i,n in enumerate(hist))+r-2*a+3*b
        for ranges in domains:
            vertices = histogram_vertices(ranges,pairs)
            self.assertLessEqual(len(vertices),24)
            upper = max(map(objective,vertices))
            for p in points(Box(ranges)):
                r,a,b = p
                if min(pairs-r+a,r-2*a,a,pairs-r+b,r-2*b,b)>=0:
                    self.assertLessEqual(objective(p),upper)
        self.assertIn((F(3),F(3,2),F(3,2)),histogram_vertices(domains[2],pairs))

    def test_subtraction_partitions_including_outside_cuts(self):
        root = Box(((1,6),(0,3),(0,3)))
        for cut in (Box(((2,4),(1,2),(1,2))),Box(((0,9),(-1,4),(-1,4))),
                    Box(((4,9),(2,8),(1,2))),Box(((7,8),(0,0),(0,0)))):
            result = root.subtract(cut)
            all_points = set()
            for box in result:
                selected = points(box)
                self.assertFalse(all_points&selected)
                all_points |= selected
            self.assertEqual(all_points,points(root)-points(cut))

    @staticmethod
    def receipt(box,upper='-100'):
        return dict(ranges=box.ranges,log2_upper_rational=upper,directed=True,
                    covers_all_even_odd_splits=True,construction=CONSTRUCTION,
                    field_order=256,forbidden_weight=734003)

    def test_ledger_preserves_unresolved_and_rejects_overlap(self):
        ledger = Ledger(Box(((1,6),(0,3),(0,3))),pairs=3)
        accepted = Box(((4,6),(0,3),(0,3)))
        ledger.accept(self.receipt(accepted))
        summary = ledger.summary()
        self.assertFalse(summary['complete_local_cover'])
        self.assertFalse(summary['whole_code_certificate'])
        self.assertGreater(summary['unresolved_profile_triples'],0)
        self.assertTrue(summary['volume_accounting'])
        with self.assertRaises(ValueError):
            ledger.accept(self.receipt(accepted))
        with self.assertRaises(ValueError):
            ledger.accept(dict(self.receipt(ledger.unresolved[0]),construction='xor-outer'))
        for remaining in list(ledger.unresolved):
            ledger.accept(self.receipt(remaining))
        self.assertTrue(ledger.summary()['complete_local_cover'])
        # A toy domain can never be reported as this K20 code certificate.
        self.assertFalse(ledger.summary()['whole_code_certificate'])

    def test_scalar_half_integer_vertices_use_directed_gamma(self):
        # Construct only the scalar component; no numerical packet table is
        # needed to test the clipped-domain normalization.
        with ctx.workprec(256):
            backend = CachedBound.__new__(CachedBound)
            backend.constant,backend.slope = arb(7),arb(3)
            backend.log_markers = ((arb(-2),arb(-3)),(arb(-4),arb(-1)))
            backend.log_factorials = {}
            ranges = ((3,7),(0,3),(0,3))
            vertices = histogram_vertices(ranges)
            for matching in (0,1):
                bound = max(backend.scalar(p,matching).upper() for p in vertices)
                for p in points(Box(ranges)):
                    r,a,b = p
                    if min(4096-r+a,r-2*a,4096-r+b,r-2*b)>=0:
                        self.assertTrue(backend.scalar(p,matching).upper()<=bound)

    def test_disjoint_acceptance_can_cross_partition_slabs_and_is_snapshotted(self):
        ledger = Ledger(Box(((1,6),(0,3),(0,3))),pairs=3)
        first = Box(((2,3),(0,1),(0,1)))
        second = Box(((1,4),(2,3),(1,2)))
        receipt = self.receipt(first)
        ledger.accept(receipt)
        receipt['log2_upper_rational'] = '100'
        self.assertEqual(ledger.accepted[0]['log2_upper_rational'],'-100')
        ledger.accept(self.receipt(second))
        self.assertTrue(ledger.summary()['volume_accounting'])
        accepted_points = points(first)|points(second)
        remaining = set()
        for box in ledger.unresolved:
            self.assertFalse(remaining&points(box))
            remaining |= points(box)
        self.assertEqual(remaining,points(ledger.root)-accepted_points)

    def test_band_bound_keeps_full_leaf_mass(self):
        ledger = Ledger(Box(((1,6),(0,3),(0,3))),pairs=3)
        ledger.accept(self.receipt(Box(((4,6),(0,3),(0,3))),'-90'))
        band = ledger.band_summary(5,6)
        self.assertTrue(band['complete_local_cover'])
        self.assertGreaterEqual(F(band['accepted_union_log2_upper']),-90)
        self.assertLessEqual(F(band['accepted_union_log2_upper']),F(-90)+F(1,10**12))
        self.assertEqual(band['accepted_profile_triples'],profile_count(((5,6),(0,3),(0,3)),3))
        self.assertFalse(ledger.band_summary(3,6)['complete_local_cover'])
        self.assertEqual(ledger.summary()['complete_top_weight_range'],(4,6))
        with self.assertRaises(ValueError):
            ledger.band_summary(0,6)

    def test_top_range_never_overstates_exact_small_domain(self):
        pairs = 3
        feasible = set(exact_profiles(pairs))-{(0,0,0)}
        intervals = lambda n,start=0:[(lo,hi) for lo in range(start,n+1) for hi in range(lo,n+1)]
        for ranges in product(intervals(2*pairs,1),intervals(pairs),intervals(pairs)):
            box = Box(ranges)
            ledger = Ledger(box,pairs=pairs)
            top = ledger.summary()['complete_top_weight_range']
            if top is not None:
                self.assertFalse(any(top[0]<=r<=top[1] and all(lo<=v<=hi for v,(lo,hi) in zip((r,a,b),ranges))
                                     for r,a,b in feasible))

    def test_refinement_budget_retains_unfinished_domain(self):
        class Backend:
            def receipt(self,ranges):
                return BackendTests.receipt(Box(ranges),'1')
            def witness(self):
                return dict(test_only=True)
        ledger = Ledger(Box(((1,6),(0,3),(0,3))),pairs=3)
        before = ledger.summary()['total_profile_triples']
        result = refine(ledger,[Backend()],max_nodes=1,minimum_r=4)
        self.assertEqual(result['visited'],1)
        self.assertEqual(result['accepted'],0)
        self.assertEqual(result['summary']['unresolved_profile_triples'],before)
        self.assertFalse(result['summary']['complete_local_cover'])
        self.assertTrue(any(box.ranges[0]==(1,3) for box in ledger.unresolved))

    def test_import_band_restricts_events_not_probability_and_keeps_gaps(self):
        source = Ledger(Box(((1,6),(0,3),(0,3))),pairs=3)
        source.accept(self.receipt(Box(((2,5),(0,1),(0,2))),'-90'))
        target = Ledger(Box(((1,6),(0,3),(0,3))),pairs=3)
        imported = target.import_band(source,3,4,label='test')
        self.assertEqual(imported['imported_boxes'],1)
        self.assertFalse(imported['complete_local_cover'])
        self.assertEqual(F(target.accepted[0]['log2_upper_rational']),-90)
        self.assertEqual(target.accepted[0]['ranges'],((3,4),(0,1),(0,2)))
        self.assertEqual(target.accepted[0]['restricted_from_ranges'],((2,5),(0,1),(0,2)))
        self.assertEqual(source.accepted[0]['ranges'],((2,5),(0,1),(0,2)))
        before = target.summary()
        with self.assertRaises(ValueError):
            target.import_band(source,4,5)
        self.assertEqual(target.summary(),before)
        target.import_band(source,5,6)
        self.assertTrue(target.summary()['profile_accounting'])

    def test_import_rejects_wrong_root_and_invalid_receipt_atomically(self):
        target = Ledger(Box(((1,6),(0,3),(0,3))),pairs=3)
        source = Ledger(Box(((2,6),(0,3),(0,3))),pairs=3)
        source.accept(self.receipt(Box(((2,3),(0,3),(0,3)))))
        source.accept(self.receipt(Box(((4,5),(0,3),(0,3)))))
        with self.assertRaises(ValueError):
            target.import_band(source,1,3)
        source.accepted[-1]['construction'] = 'different-code'
        before = target.summary()
        with self.assertRaises(ValueError):
            target.import_band(source,2,5)
        self.assertEqual(target.summary(),before)

    def test_replay_cover_rejects_overlapping_assignments_before_work(self):
        for inputs in ([('a',{},1,10),('b',{},10,20)],[('a',{},7600,7650)]):
            with self.assertRaises(ValueError):
                replay_cover(inputs)


if __name__=='__main__':
    unittest.main()
