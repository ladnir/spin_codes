import unittest

import ladder_q1 as ladder


class LadderTests(unittest.TestCase):
    def test_old_instance_identity_preserved(self):
        for m in (16, 18, 20):
            self.assertEqual(ladder.identity.instance(m), ladder.core.instance('t128_s19', m))

    def test_extended_geometry_and_map(self):
        original = ladder.identity.instance(20)
        for m in (22, 24):
            value = ladder.identity.instance(m)
            self.assertEqual(value['rows'], 1 << (m-7))
            self.assertEqual(value['cutoff'], (1 << (m+1))//10)
            self.assertEqual(value['selection_sha256'], original['selection_sha256'])
            self.assertNotEqual(value['fingerprint'], original['fingerprint'])

    def test_reject_unintended_size(self):
        for m in (15, 25, 24.0):
            with self.assertRaises(ValueError):
                ladder.identity.instance(m)

    def test_matches_frozen_q1_at_supported_size(self):
        spec = ladder.identity.instance(16)
        expected = ladder.frozen.compute(dict(instance=spec, kind='q1', scaled_tilts=[295]))
        actual = ladder.compute(spec, [295], 256)
        self.assertEqual({str(w): ladder.base.encode(v) for w, v in actual.items()}, expected['coefficient_upper'])


if __name__ == '__main__':
    unittest.main()
