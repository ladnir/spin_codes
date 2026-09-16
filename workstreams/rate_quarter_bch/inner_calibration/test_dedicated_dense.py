"""Replay dedicated-search and diagnostic receipts without rerunning search."""
import json
import math
import unittest

import numpy as np
import dedicated_dense as dedicated
import refine_candidates as refinement

check = dedicated.check


class DedicatedDenseTests(unittest.TestCase):
    def test_tightening_preserves_simplex_not_penalty(self):
        lower,upper,total = [0,0,0],[10,2,2],10
        vertices = check.typed.vertices(lower,upper,total)
        lo,hi = vertices.min(axis=0),vertices.max(axis=0)
        np.testing.assert_array_equal(check.typed.vertices(lo,hi,total),vertices)
        self.assertLessEqual(check.typed.lattice_log_count(lo,hi),check.typed.lattice_log_count(lower,upper))

    def test_all_new_receipt_sources(self):
        for name in ('t256_s18_nested_dedicated_dense.json','t256_s18_nested_point_probe.json',
                     't256_s18_nested_path_diagnostic.json','ZERO_STATE_RANK_CEILING.json','ZERO_STATE_OBSTRUCTION.json'):
            dedicated.validate_receipt(check.calibration.HERE/name)

    def test_dedicated_dense_replay(self):
        folder = check.calibration.HERE
        data = dedicated.validate_receipt(folder/'t256_s18_nested_dedicated_dense.json')
        self.assertEqual({r['distance_target'] for r in data['results']},{'33/200','19/100'})
        record = json.loads((folder/'maps/t256_s18_nested.json').read_text())
        counts = {w:n for w,n in enumerate(check.calibration.smaller_outer.spectrum()) if w and n}
        for row in data['results']:
            replay = check.transport_dense(record,counts,row,129)
            self.assertAlmostEqual(replay['log_union_upper'],row['dense']['log_union_upper'],places=8)
            self.assertLessEqual(replay['log_union_upper'],row['search_log_union_upper_before_coordinate_tightening']+2e-6)
            combined = -float(np.logaddexp(replay['log_union_upper'],-row['sparse_union_margin_bits']*math.log(2)))/math.log(2)
            self.assertAlmostEqual(combined,row['combined_margin_bits'],places=8)

    def test_singleton_replay(self):
        folder = check.calibration.HERE
        data = dedicated.validate_receipt(folder/'t256_s18_nested_point_probe.json')
        original = dedicated.validate_receipt(folder/'t256_s18_nested_full.json')
        record = json.loads((folder/'maps/t256_s18_nested.json').read_text())
        epochs = refinement.Epochs(record)
        models = {r['distance_target']:refinement.Refiner(r,record,epochs) for r in original['results']}
        for row in data['results']:
            box = dict(lower=row['type_counts'],upper=row['type_counts'])
            bound = models[row['distance_target']].direct(box,row['witness'])[0]
            self.assertAlmostEqual(-bound/math.log(2),row['refined_margin_bits'],places=7)


if __name__=='__main__': unittest.main()
