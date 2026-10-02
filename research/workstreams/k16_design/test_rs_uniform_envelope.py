"""Exact finite checks of the independently randomized MDS uniform envelope."""

from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb
import unittest

from rs_uniform_envelope import UniformInputEnvelope


def image(matrix, value):
    return sum(((row & value).bit_count() & 1) << bit
               for bit, row in enumerate(matrix))


class RsUniformEnvelopeTests(unittest.TestCase):
    def test_both_target_variants(self):
        for envelope in (UniformInputEnvelope(),
                         UniformInputEnvelope(n=16, k=8, packets_per_symbol=4)):
            envelope.verify_shell_domination()
            self.assertEqual(envelope.message_bits, 128)
            self.assertEqual(envelope.output_bits, 256)
            self.assertEqual(envelope.regions, 64)
            self.assertEqual(envelope.activity, Fraction(15, 16))
            self.assertEqual(envelope.beta,
                             Fraction(1 << 256, (envelope.q - 1) ** (envelope.n - envelope.k)))
            self.assertEqual(envelope.shell_caps()[0], 0)
            self.assertEqual(sum(envelope.shell_caps()),
                             envelope.beta - envelope.pointwise_nonzero_cap)
            self.assertEqual(sum(envelope.shell_caps(include_zero=True)), envelope.beta)
            self.assertEqual(envelope.exact_pointwise_symbol_counts()[envelope.minimum_symbol_weight],
                             envelope.pointwise_nonzero_cap)
            metadata = envelope.metadata()
            self.assertTrue(metadata['iid_envelope_includes_artificial_zero'])
            self.assertFalse(metadata['whole_code_certificate'])
            self.assertFalse(metadata['deterministic_enumerator_bound'])

    def test_exact_bernoulli_interpretation(self):
        for envelope in (UniformInputEnvelope(3, 2, 1, 2),
                         UniformInputEnvelope(2, 1, 2, 1),
                         UniformInputEnvelope(16, 8, 4, 4)):
            p = envelope.activity
            for v, value in enumerate(envelope.shell_caps(include_zero=True)):
                self.assertEqual(value, envelope.beta * comb(envelope.regions, v) *
                                 p ** v * (1 - p) ** (envelope.regions - v))

    def test_all_fixed_gl_setups_pointwise(self):
        # GF4 MDS[3,2] is (u,v,u+v). Fix the whole setup before enumerating
        # messages, then average over the 6^3 independent GL2 setups.
        matrices = tuple((a, b) for a in range(1, 4) for b in range(1, 4) if a != b)
        counts = Counter()
        for setup in product(matrices, repeat=3):
            for u, v in product(range(4), repeat=2):
                if u or v:
                    counts[tuple(image(m, x) for m, x in zip(setup, (u, v, u ^ v)))] += 1
        envelope = UniformInputEnvelope(3, 2, 1, 2)
        exact_by_weight = envelope.exact_pointwise_symbol_counts()
        for output in product(range(4), repeat=3):
            h = sum(x != 0 for x in output)
            actual = Fraction(counts[output], len(matrices) ** 3)
            self.assertEqual(actual, exact_by_weight[h])
            self.assertLessEqual(actual, envelope.pointwise_nonzero_cap)
            if h < envelope.minimum_symbol_weight:
                self.assertEqual(actual, 0)
        self.assertEqual(sum(counts.values()), len(matrices) ** 3 * (4**2 - 1))
        self.assertEqual(max(Fraction(value, len(matrices) ** 3) for value in counts.values()),
                         envelope.pointwise_nonzero_cap)

    def test_common_gl_map_does_not_satisfy_the_premise(self):
        # A shared binary map preserves w=u+v, so this ensemble is supported
        # on the original code and its expected multiplicity can be one.
        matrices = tuple((a, b) for a in range(1, 4) for b in range(1, 4) if a != b)
        counts = Counter()
        for matrix in matrices:
            for u, v in product(range(4), repeat=2):
                if u or v:
                    counts[tuple(image(matrix, x) for x in (u, v, u ^ v))] += 1
        self.assertEqual(max(Fraction(c, len(matrices)) for c in counts.values()), 1)
        self.assertGreater(Fraction(counts[(1, 0, 1)], len(matrices)),
                           UniformInputEnvelope(3, 2, 1, 2).pointwise_nonzero_cap)

    def test_full_space_and_repetition_endpoints(self):
        for n, k, bits, packets in ((3, 3, 1, 2), (3, 1, 1, 2), (1, 1, 1, 1)):
            envelope = UniformInputEnvelope(n, k, bits, packets)
            envelope.verify_shell_domination()
            if n == k:
                self.assertEqual(envelope.pointwise_nonzero_cap, 1)
                self.assertEqual(envelope.beta, 1 << envelope.output_bits)
        repeated = UniformInputEnvelope(3, 1, 1, 2)
        self.assertEqual(repeated.exact_pointwise_symbol_counts(),
                         (0, 0, 0, Fraction(1, 9)))

    def test_bad_parameters(self):
        for invalid in (0, -1, 1.5, True):
            for name in ('n', 'k', 'packet_bits', 'packets_per_symbol'):
                with self.assertRaises(ValueError):
                    UniformInputEnvelope(**{name: invalid})
        for args in ((3, 4, 1, 2), (5, 2, 1, 2)):
            with self.assertRaises(ValueError):
                UniformInputEnvelope(*args)
        with self.assertRaises(ValueError):
            UniformInputEnvelope().shell_caps(include_zero=1)


if __name__ == '__main__':
    unittest.main()
