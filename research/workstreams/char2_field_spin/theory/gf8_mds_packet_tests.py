"""Exact scalar tests for the independent full-Cauchy8/4 candidate."""
import unittest

import gf8_mds_packet as m


class PacketTests(unittest.TestCase):
    def test_all_exact_receipts(self):
        result=m.check()
        self.assertEqual(result['joint_rank_cases'],1834)
        self.assertEqual(result['packed_basis_scalar_random'],6208)

    def test_public_interface(self):
        f,a,c,eta=m.matrices()
        self.assertEqual((len(f),len(f[0]),len(a),len(a[0]),len(c),len(c[0]),eta),(8,8,8,4,4,8,2))
        self.assertEqual(m.cost()['per16_gfni512'],34)

    def test_fullstate_extension(self):
        result=m.check_fullstate()
        self.assertEqual(result['packed_basis_scalar_random'],8256)
        self.assertEqual(result['four_packet_batches'],128)
        self.assertEqual(m.fullstate_cost()['batch_gfni512'],81)
        self.assertEqual(m.fullstate_factored_cost()['batch_gfni512'],62)

    def test_largerstate_extension(self):
        result=m.check_largestate()
        self.assertEqual(result['selected_rank_cases'],2048)
        self.assertEqual(result['packed_basis_scalar_random'],12352)
        self.assertEqual(result['four_packet_batches'],128)
        self.assertEqual(m.largestate_cost()['batch_gfni512'],116)


if __name__=='__main__':unittest.main()
