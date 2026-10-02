import copy
import hashlib
import json
import unittest

import shared_relaxed_audit as a


class SharedAuditTests(unittest.TestCase):
    def inputs(self):
        dense = dict(schema=a.proof.DENSE_SCHEMA, updates=2, K=1 << 20, N=1 << 21,
            minimum_groups=2, zero_bits=64, variance_bins=16, base_tilt='3/16', cost_tilt='1/4',
            results=[dict(distance='1/20', threshold=104857, root=['0', '1'], cover=dict(
                leaves={'': dict(witness=dict(parameters=['1/100', '0', '0']))}, unresolved={}))])
        raw = json.dumps(dense).encode()
        report = dict(schema=a.proof.SCHEMA, ensemble=a.proof.ENSEMBLE, complete=True, precision=384,
            dense_sha256=hashlib.sha256(raw).hexdigest(), message_length=1 << 20, output_length=1 << 21,
            updates=2, threshold=104857, minimum_distance=104858, sparse_occupancies=[1, 1],
            dense_occupancies=[2, 2048], dense_leaves=1,
            dense_checked=[dict(path='', upper=[1, -40], output_tilt='1/100')],
            sparse_checked=[dict(occupancy=1, updates=2, upper=[1, -41])],
            dense_upper=[1, -40], sparse_upper=[1, -41], total_upper=[3, -41])
        return report, raw

    def test_complete_consistent_receipt(self):
        report, raw = self.inputs()
        result = a.audit(report, raw, '1/20', 20)
        self.assertTrue(result['verified_scope_and_sum'])
        self.assertTrue(result['numerical_replay_required'])
        self.assertEqual(result['minimum_distance'], 104858)

    def test_scope_tampering_rejected(self):
        report, raw = self.inputs()
        for key, value in (('ensemble', 'independent4-gf16-r2'), ('updates', 4),
                ('precision', 192), ('dense_sha256', 'wrong'), ('minimum_distance', 104859),
                ('sparse_occupancies', [1, 2]), ('dense_occupancies', [3, 2048]), ('dense_leaves', 2)):
            changed = dict(report); changed[key] = value
            with self.assertRaises(ValueError): a.audit(changed, raw, '1/20', 20)
        with self.assertRaises(ValueError): a.audit(report, raw, '3/50', 20)

    def test_leaf_and_sum_tampering_rejected(self):
        report, raw = self.inputs()
        for key, value in (('dense_checked', []), ('sparse_checked', []),
                ('dense_upper', [1, -41]), ('sparse_upper', [1, -42]), ('total_upper', [1, -40]),
                ('total_upper', [1, -20]), ('dense_checked', report['dense_checked']*2)):
            changed = dict(report); changed[key] = value
            with self.assertRaises(ValueError): a.audit(changed, raw, '1/20', 20)
        changed = copy.deepcopy(report); changed['dense_checked'][0]['output_tilt'] = '1/99'
        with self.assertRaises(ValueError): a.audit(changed, raw, '1/20', 20)


if __name__ == '__main__': unittest.main()
