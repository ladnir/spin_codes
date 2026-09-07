"""Authentication and scope checks for retained local optimization receipts."""
import math
import unittest
from pathlib import Path
import inner_candidate_screen as screen


class ReceiptTests(unittest.TestCase):
    def test_anchor_scope_and_authentication(self):
        files=list((screen.base.HERE/'generated').glob('inner_*_q*_v1.json'))
        self.assertTrue(files)
        for path in files:
            data=screen.base.read(path)
            if data['status']!='CANDIDATE_FIXED_WEIGHT_ANCHOR_OUTWARD': continue
            self.assertFalse(data['full_distance_proved'])
            self.assertEqual(data['rows'],1<<(data['message_exponent']-7))
            self.assertEqual(data['cutoff'],256*data['rows']//10)
            self.assertTrue(2<=data['occupation']<=data['rows'])
            upper=screen.base.decode(data['upper'])
            self.assertEqual(upper,screen.base.F(2)**(-data['retained_margin_bits']))
            for name,digest in data['source_sha256'].items():
                self.assertEqual(screen.base.sha(screen.base.ROOT/name),digest,name)

    def test_sampled_screens_do_not_claim_full_coverage(self):
        files=list((screen.base.HERE/'generated').glob('inner_*_dense_v1.json'))
        self.assertTrue(files)
        for path in files:
            data=screen.base.read(path)
            self.assertFalse(data['full_distance_proved'])
            self.assertEqual(data['status'],'BCH256_INNER_DENSE_SAMPLED_DIAGNOSTIC')
            self.assertTrue(all(math.isfinite(r['margin_bits']) for r in data['results']))
            for name,digest in data['source_sha256'].items():
                self.assertEqual(screen.base.sha(screen.base.ROOT/name),digest,name)


if __name__=='__main__': unittest.main()
