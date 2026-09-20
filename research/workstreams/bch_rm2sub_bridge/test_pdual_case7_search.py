"""Checks for normalized case-7 discovery and exact low-weight witnesses."""
import math
import unittest
import bridge as base
import pdual_low_weight as low
import pdual_weight30_sat as sat
import pdual_normalized_extension as extension
import pdual_rank_refinement as rank
import pdual_quintic_structure as quintic
from affine_wambach import gf_mul,gf_pow


class Case7SearchTest(unittest.TestCase):
    def test_quintic_structure_and_mobius_involution(self):
        self.assertEqual(quintic.build(),base.read(base.HERE/'generated/pdual_quintic_structure_v1.json'))
        coords=list(range(256))
        for word in (0,1,(1<<256)-1,0x123456789abcdef):
            self.assertEqual(quintic.anf(quintic.anf(word,coords),coords),word)

    def test_nonzero_normalization(self):
        self.assertEqual(math.gcd(7,255),1)
        self.assertEqual({gf_pow(a,7) for a in range(1,256)},set(range(1,256)))
        checks,rows,values=low.data()
        coords=[gf_pow(2,i) for i in range(255)]+[0];positions={x:i for i,x in enumerate(coords)}
        for word in rows:
            moved=sum(1<<positions[gf_mul(2,x)] for i,x in enumerate(coords) if word>>i&1)
            self.assertEqual(low.syndrome(moved,values),gf_mul(gf_pow(2,7),low.syndrome(word,values)))
            self.assertTrue(all((moved&r).bit_count()%2==0 for r in checks))

    def test_membership_rejects_invalid_words(self):
        with self.assertRaises(AssertionError):low.verify_word(1)
        with self.assertRaises(AssertionError):low.verify_word(0)
        with self.assertRaises(AssertionError):sat.verify_word((1<<30)-1)

    def test_search_receipts(self):
        for pattern in ('pdual_low_weight_v*.json','pdual_stern_search_v*.json','pdual_weight30_sat_v*.json','pdual_normalized_extension_v*.json'):
            files=list((base.HERE/'generated').glob(pattern));self.assertTrue(files)
            for p in files:
                saved=base.read(p)
                for name,digest in saved['source_sha256'].items():self.assertEqual(base.sha(base.ROOT/name),digest)
                if pattern.startswith(('pdual_low_weight','pdual_stern_search')):
                    for name in ('best','best_nonzero_F7'):
                        self.assertEqual(low.verify_word(int(saved['audit'][name]['word_hex'],16)),saved['audit'][name])
                    self.assertFalse(saved['status'].startswith('PROOF'))
                elif pattern.startswith('pdual_weight30_sat'):
                    self.assertFalse(saved['proof_of_absence'])
                    if saved['status']=='sat':self.assertEqual(sat.verify_word(int(saved['word_hex'],16)),saved['exact_witness_audit'])
                elif saved['found']:
                    self.assertEqual(extension.check(saved['witness']),31)
                else:self.assertFalse(saved['proof_of_impossibility'])


if __name__=='__main__':unittest.main()
