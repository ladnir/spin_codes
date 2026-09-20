"""Map identity, exact algebra, and retained calibration receipt regressions."""
import json
import math
import unittest

import numpy as np
import calibrate
import check_candidates as check
import refine_candidates as refinement


class CalibrationTests(unittest.TestCase):
    def test_map_audits_and_nested_identity(self):
        screen = json.loads((calibrate.HERE/'Q1_SCREEN.json').read_text())
        for row in screen['results']:
            path = calibrate.fixed.ROOT/row['map_source']
            self.assertEqual(calibrate.fixed.sha(path),row['map_sha256'])
            record = json.loads(path.read_text())
            t,s = record['step_bits'],record['state_bits']
            self.assertEqual(sum(record['a_counts']),1<<s)
            self.assertEqual(calibrate.fixed.outer.macwilliams(record['a_counts'],s),record['kernel_counts'])
            self.assertEqual(calibrate.algebra.kernel_weight_four_count(record['columns']),record['kernel_weight_four'])
            rows = [int(x,16) for x in record['generator_rows_hex']]
            self.assertEqual(calibrate.algebra.rank(rows),s)
            self.assertEqual(calibrate.algebra.enumerate_spectrum(rows,s,t),record['a_counts'])
            self.assertTrue(all((a&b).bit_count()%2==0 for a in rows for b in rows))
            if row['tag'].endswith('_nested') and s>14:
                prev = json.loads((calibrate.HERE/'maps'/f't{t}_s{s-1}_nested.json').read_text())
                self.assertEqual([c&((1<<(s-1))-1) for c in record['columns']],prev['columns'])

    def test_selected_baseline(self):
        record = json.loads((calibrate.HERE/'maps/t128_s19_nested.json').read_text())
        a,kernel,_ = calibrate.fixed.load_inner()
        self.assertEqual(a,{w:n for w,n in enumerate(record['a_counts']) if w and n})
        self.assertEqual(kernel,record['kernel_counts'])
        control = json.loads((calibrate.HERE/'t128_s19_nested_dense.json').read_text())
        baseline = json.loads((calibrate.fixed.HERE/'SMALLER_MARGIN_CLOSED.json').read_text())
        self.assertAlmostEqual(control['results'][1]['dense_union_margin_bits'],baseline['results'][1]['dense_union_margin_bits'],places=8)

    def test_prepared_epochs_multiple_sizes(self):
        for tag in ('t64_s16_selected','t128_s18_nested','t256_s18_nested'):
            record = json.loads((calibrate.HERE/'maps'/f'{tag}.json').read_text())
            prepared = refinement.Epochs(record)
            a = {w:n for w,n in enumerate(record['a_counts']) if w and n}
            for z in (-8.,-.5):
                direct = calibrate.fixed.general.epoch_logs(record['step_bits'],record['state_bits'],a,
                            record['kernel_counts'],math.exp(z),record['step_bits'])
                np.testing.assert_allclose(prepared.at(z),direct,atol=2e-11,rtol=2e-13)

    def test_receipt_hashes_and_union(self):
        paths = [calibrate.HERE/'Q1_SCREEN.json',*calibrate.HERE.glob('*_dense.json'),
                 *calibrate.HERE.glob('*_full.json'),*calibrate.HERE.glob('*_dense_refined.json')]
        for path in paths:
            payload = json.loads(path.read_text())
            for name,digest in payload['source_sha256'].items():
                self.assertEqual(calibrate.fixed.sha(calibrate.fixed.ROOT/name),digest,f'{path.name}: {name}')
            for row in payload['results']:
                if 'dense' not in row: continue
                dense = row['dense']
                check.check_coverage(dense['selected_boxes'],check.L,dense['occupation_min'])
                union = float(np.logaddexp.reduce([b['own_log_bound'] for b in dense['selected_boxes']]))
                self.assertAlmostEqual(union,dense['log_union_upper'],places=8)
                if 'occupation_margins_bits' in row:
                    self.assertEqual(row['covered_occupations'],[1,dense['occupation_min']-1])
                    self.assertEqual(len(row['occupation_margins_bits']),row['covered_occupations'][1])
                    sparse = float(np.logaddexp.reduce(-np.array(row['occupation_margins_bits'])*math.log(2)))
                    self.assertAlmostEqual(-sparse/math.log(2),row['sparse_union_margin_bits'],places=8)
                    self.assertAlmostEqual(-float(np.logaddexp(sparse,union))/math.log(2),row['combined_margin_bits'],places=8)

    def test_q1_replay_and_smaller_state_target(self):
        counts = {w:n for w,n in enumerate(calibrate.smaller_outer.spectrum()) if w and n}
        record = json.loads((calibrate.HERE/'maps/t64_s16_selected.json').read_text())
        rows = calibrate.q1_screen(record,counts)
        retained = json.loads((calibrate.HERE/'t64_s16_selected_full.json').read_text())['results']
        self.assertGreaterEqual(retained[0]['combined_margin_bits'],40.)
        for screen,full in zip(rows,retained):
            self.assertAlmostEqual(screen['q1_margin_bits'],full['occupation_margins_bits'][0],places=9)

    def test_candidate_dense_witness_replay(self):
        counts = {w:n for w,n in enumerate(calibrate.smaller_outer.spectrum()) if w and n}
        for path in calibrate.HERE.glob('*_full.json'):
            payload = json.loads(path.read_text())
            record = json.loads((calibrate.HERE/'maps'/f"{payload['tag']}.json").read_text())
            for row in payload['results']:
                replay = check.transport_dense(record,counts,row,row['dense']['occupation_min'])
                self.assertAlmostEqual(replay['log_union_upper'],row['dense']['log_union_upper'],places=8)
                for a,b in zip(replay['selected_boxes'],row['dense']['selected_boxes']):
                    self.assertAlmostEqual(a['own_log_bound'],b['own_log_bound'],places=8)


if __name__=='__main__': unittest.main()
