"""Memory-only completion controls; no coefficient or production-index work."""
from copy import deepcopy
import unittest
import hd_field128_gamma_completion_v31 as m


def fixture():
    common = dict(construction=m.CONSTRUCTION, field_order=1 << 128, source_dimension=1 << 20,
                  output_length=1 << 21, target_distance='2/5', whole_code_certificate=False)
    ranges = [[i,i] for i in range(1,260)]+[[260,520],[548,548],[576,916],[920,1145],[1150,6677],[6678,8192]]
    entries = [dict(path=str(i), sha256=str(i), weight_range=r, log2_upper_rational='-100') for i,r in enumerate(ranges)]
    entries[-1].update(path=m.TERMINAL, sha256=m.PINS[m.TERMINAL], log2_upper_rational='-41')
    parent = dict(common, entries=deepcopy(entries))
    entries += [dict(path='new'+str(i),sha256='new'+str(i),weight_range=list(r),log2_upper_rational='-100')
                for i,r in enumerate(m.ADDITIONS)]
    saved = dict(common, entries=entries, schema='field128-gamma-partial-cover-index-v31-composite-selected',
        completion_audit_performed=False, receipt_count=274, covered_weight_count=8192, pending_weight_count=0,
        covered_weight_ranges=[[1,8192]], pending_weight_ranges=[], added_event_count=9, added_weight_count=61,
        aggregation_precision_bits=384, covered_events_failure_bits=40, directed_aggregation=True,
        meets_covered_events_target=True, terminal_counted_once=True)
    saved['covered_events_union_log2_upper_rational'] = m.positive_union(entries)[0]
    terminal = dict(common, state_dimension=30, inner_width=64, forbidden_weight=838860,
        weight_range=[6678,8192], gamma_extension_degree=16, gamma_base_field_order=256,
        feedback='gamma A transpose', CA_zero=True, fixed_nodes='GF256 additive subspace0..63',state_nodes='128..157',
        refresh='D2 H_30 D1 with independent full-field nonzero diagonals',initial_state='zero',final_state='unrestricted',
        state_persists_across_regions=True,outer_independent_of_route_and_inner=True,full_field_random_coefficients=True,
        outer_groups=8192,outer_randomizers_per_group=64,pair_swaps=True,
        outer_randomizer='independent uniform invertible4x4 Toeplitz per group and column',
        method='inverse-image-terminal-range',directed=True,covers_all_weights_in_range_simultaneously=True,
        no_union_over_tail_weights=True,log2_upper_rational='-41')
    return saved,parent,terminal


class CompletionTests(unittest.TestCase):
    def test_complete_tail_once(self):
        s,p,t=fixture(); r=m.verify_objects(s,p,t,deepcopy(s))
        self.assertEqual(r['event_count'],274)
        self.assertEqual(r['actual_union_log2_upper_rational'],'-40999999999999/1000000000000')
        self.assertFalse(s['whole_code_certificate'])

    def test_failures(self):
        for change in ('gap','overlap','prefix','ensemble','threshold','terminal','whole','union','counter'):
            with self.subTest(change=change):
                s,p,t=fixture()
                if change=='gap':s['entries'][-1]['weight_range']=[1149,1150]
                if change=='overlap':s['entries'][-1]['weight_range']=[1148,1149]
                if change=='prefix':s['entries'][0]['log2_upper_rational']='-101'
                if change=='ensemble':s['field_order']=256
                if change=='threshold':t['forbidden_weight']=838861
                if change=='terminal':t['no_union_over_tail_weights']=False
                if change=='whole':s['whole_code_certificate']=True
                if change=='union':s['entries'][-1]['log2_upper_rational']='-39'
                if change=='counter':s['pending_weight_count']=1
                with self.assertRaises(ValueError):m.verify_objects(s,p,t,deepcopy(s))

    def test_saved_rebuild_and_zero(self):
        s,p,t=fixture(); other=deepcopy(s); other['extra']=1
        with self.assertRaises(ValueError):m.verify_objects(s,p,t,other)
        a=[dict(zero=True,log2_upper_rational=None),dict(log2_upper_rational='-50')]
        self.assertLessEqual(m.F(m.positive_union(a)[0]),-49)
        a[0]['log2_upper_rational']='-90'
        with self.assertRaises(ValueError):m.positive_union(a)

    def test_wrong_final_pins(self):
        with self.assertRaises(ValueError):m.check('0'*64,'1'*64)


if __name__=='__main__':unittest.main()
