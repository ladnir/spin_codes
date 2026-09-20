"""Grid coverage and reuse invariants, independent of numerical run progress."""
import unittest

import run_complete_q1_grid as grid


class CompleteGridTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.existing, cls.maps, _, cls.models, _ = grid.inventory()
        cls.rows = list(grid.tuples(cls.models))

    def test_complete_cartesian_grid_with_explicit_native_exclusions(self):
        self.assertEqual(len(self.rows), 3120)
        self.assertEqual(len({grid.key(r) for r in self.rows}), len(self.rows))
        invalid = [r for r in self.rows if not r['native']]
        self.assertEqual(len(invalid), 12)
        self.assertTrue(all((r['block_bits'],r['step_bits'],r['message_exponent']) == (1024,256,16)
                            for r in invalid))
        native = {grid.key(r) for r in self.rows if r['native']}
        self.assertEqual(len(native), 3108)
        self.assertEqual(len(self.existing), 1012)
        self.assertTrue(set(self.existing).issubset(native))

    def test_exact_endpoints_and_random_range(self):
        bch = sorted(r['block_bits'] for r in self.models if 'BCH' in r['series'])
        rm = sorted(r['block_bits'] for r in self.models if r['series'].startswith('RM('))
        random = sorted(r['block_bits'] for r in self.models if r['outer_model']=='random-ensemble-average')
        self.assertEqual(bch, [8,32,64,128])
        self.assertEqual(rm, [8,32,128,512])
        self.assertEqual(random, [8,16,32,64,128,256,512,1024])

    def test_existing_map_identity_retained(self):
        import json
        for (t,s),path in self.maps.items():
            payload=json.loads(path.read_text())
            full=json.loads((grid.HERE/'activation_pilot_v1/maps'/f't{t}_s20.json').read_text())
            self.assertEqual(payload['generator_words_hex'],full['generator_words_hex'][:s])
            tags={r['map_tag'] for r in self.existing.values() if (int(r['step_bits']),int(r['state_bits']))==(t,s)}
            self.assertEqual(tags, {'nested-'+grid.pilot.sha(path)[:16]})

    def test_completed_q1_grid_is_exactly_the_native_target(self):
        from read_grid_receipts import snapshot
        planned,observations,_=snapshot()
        self.assertEqual(set(observations),{grid.key(r) for r in planned if r['native']})

    def test_dense_partition_covers_all_integers_without_gaps(self):
        from run_dense_range_grid import intervals
        for length in (128,512,4194304):
            parts=intervals(length,17,64)
            self.assertEqual(parts[0][0],17)
            self.assertEqual(parts[-1][1],length)
            self.assertTrue(all(a<=b for a,b in parts))
            self.assertTrue(all(a[1]+1==b[0] for a,b in zip(parts,parts[1:])))
            self.assertEqual(sum(b-a+1 for a,b in parts),length-16)

    def test_coverage_uses_an_authenticated_common_setup_event(self):
        from register_complete_grid import common_event
        events = {
            'old': dict(block_bits=8,dimension=4,setup_failure_bits=60,counts={2:8,4:12}),
            'new': dict(block_bits=8,dimension=4,setup_failure_bits=60,counts={2:4,4:10}),
            'other': dict(block_bits=8,dimension=4,setup_failure_bits=60,counts={2:3,4:11}),
        }
        def row(event):
            return dict(block_bits='8',dimension='4',setup_failure_bits='60',setup_event_id=event)
        self.assertEqual(common_event([row('old'),row('new')],events),'new')
        with self.assertRaisesRegex(ValueError,'common spectrum event'):
            common_event([row('new'),row('other')],events)
        with self.assertRaisesRegex(ValueError,'disagrees'):
            common_event([dict(row('new'),dimension='5')],events)


if __name__ == '__main__':
    unittest.main()
