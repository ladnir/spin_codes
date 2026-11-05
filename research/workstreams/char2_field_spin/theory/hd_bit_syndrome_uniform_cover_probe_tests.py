import random
import unittest

from gf8_packet_maps import MUL, rank
import hd_bit_syndrome_uniform_cover_probe as p


class UniformCoverProbeTests(unittest.TestCase):
    def test_rank_one_formula_against_direct_field_rank(self):
        rng = random.Random(763)
        for _ in range(150):
            m, n = rng.randrange(1, 9), rng.randrange(1, 9)
            rows = [rng.randrange(1 << n) for _ in range(m)]
            u, v = [rng.randrange(256) for _ in range(m)], [rng.randrange(256) for _ in range(n)]
            direct = tuple(tuple(((rows[i] >> j) & 1) ^ MUL[u[i]][v[j]]
                                 for j in range(n)) for i in range(m))
            self.assertEqual(p.binary_plus_rank_one_rank(rows, u, v), rank(direct))
        self.assertEqual(p.binary_plus_rank_one_rank((1,), (1,), (1,)), 0)
        self.assertEqual(p.binary_plus_rank_one_rank((0,), (1,), (1,)), 1)

    def test_actual_restrictions_including_overlap(self):
        rng = random.Random(121)
        for construction in ('current', 'star-coupled'):
            oracle = p.RankOracle(construction)
            for _ in range(80):
                a, s = rng.randrange(65536), rng.randrange(65536)
                direct = tuple(tuple(oracle.matrix[i][j] for j in p.labels(a)) for i in p.labels(s))
                self.assertEqual(oracle.rank(s, a), rank(direct))
            self.assertEqual(oracle.rank(0, 65535), 0)
            self.assertEqual(oracle.rank(65535, 0), 0)

    def test_positive_witness_and_no_false_template_counterexample(self):
        oracle = p.RankOracle()
        cover = p.find_translated_cover(oracle, 1, (1,))
        self.assertEqual(len(cover), 16)
        self.assertTrue(p.verify_balanced_bases(oracle, 1, cover))
        self.assertFalse(p.verify_balanced_bases(oracle, 1, (1,)))
        self.assertIsNone(p.density_counterexample(oracle, 1, range(17)))

    def test_density_witness_semantics(self):
        class Identity:
            def rank(self, s, a):
                return (s & a).bit_count()
        result = p.density_counterexample(Identity(), 3, (12,))
        self.assertEqual(result['rank'], 0)
        self.assertGreater(result['deficit'], 0)

    def test_translation_is_balanced(self):
        for mask in (1, 7, 0x00ff, 0x0f0f, 0xffff):
            bases = tuple(set(p.translate(mask, s) for s in range(16)))
            degrees = [sum(b >> i & 1 for b in bases) for i in range(16)]
            self.assertEqual(len(set(degrees)), 1)


if __name__ == '__main__':
    unittest.main()
