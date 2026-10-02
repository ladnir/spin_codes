import itertools
import unittest
from fractions import Fraction as Q

from canonical_counts import transport_cdf
from local_models import convolve, full_block
from monotone import suffix_expectation, transport_shells


def allocations(total, length):
    if length == 1:
        yield (total,)
    else:
        for head in range(total+1):
            for tail in allocations(total-head, length-1):
                yield (head,)+tail


class ShellTransportTests(unittest.TestCase):
    def test_suffix_bound_is_exact_linear_program_optimum(self):
        # Exhaust all small prefix caps and every kernel on this rational
        # grid. Scaling counts also exercises genuinely fractional caps.
        for scale in (Q(1), Q(1, 2)):
            for first, second in itertools.combinations_with_replacement(range(4), 2):
                caps = (Q(0), first*scale, second*scale, 3*scale)
                feasible = [tuple(scale*v for v in a) for a in allocations(3, 3)
                            if a[0] <= first and a[0]+a[1] <= second]
                for values in itertools.product((Q(0), Q(1, 3), Q(1)), repeat=3):
                    optimum = max(sum(a*v for a, v in zip(row, values)) for row in feasible)
                    self.assertEqual(suffix_expectation(caps, (Q(0),)+values), optimum)

    def test_shell_transport_against_all_small_profile_spectra(self):
        for local in (full_block(2, 1), full_block(2, 2), (Q(0), Q(1, 2), Q(0), Q(1, 2))):
            laws = [(Q(1),)]
            for _ in range(3):
                laws.append(convolve(laws[-1], local))
            for first, second in itertools.combinations_with_replacement(range(4), 2):
                caps = (0, Q(first, 2), Q(second, 2), Q(3, 2))
                shells, cdf = transport_shells(caps, local), transport_cdf(caps, local)
                self.assertTrue(all(s <= c for s, c in zip(shells, cdf)))
                for allocation in allocations(3, 3):
                    if allocation[0] > first or sum(allocation[:2]) > second:
                        continue
                    actual = [Q(0)]*len(shells)
                    for h, mass in enumerate(allocation, 1):
                        for w, probability in enumerate(laws[h]):
                            actual[w] += Q(mass, 2)*probability
                    self.assertTrue(all(a <= s for a, s in zip(actual, shells)))

    def test_direct_nonzero_block_enumeration(self):
        # A uniform nonzero two-bit block has three equiprobable values.
        caps = (0, 1, 3, 4)
        upper = transport_shells(caps, full_block(2, 1))
        direct = [Q(0)]*7
        for h, count in enumerate((0, 1, 1, 2)):
            if not count:
                continue
            for blocks in itertools.product(range(1, 4), repeat=h):
                direct[sum(x.bit_count() for x in blocks)] += Q(count, 3**h)
        self.assertTrue(all(a <= b for a, b in zip(direct, upper)))
        self.assertTrue(any(a < b for a, b in zip(upper, transport_cdf(caps, full_block(2, 1)))))

    def test_deterministic_blocks_do_not_turn_cdf_differences_into_shells(self):
        caps = (0, Q(1, 3), Q(2, 3), 1)
        result = transport_shells(caps, (0, 1))
        self.assertEqual(result, caps)
        self.assertNotEqual(sum(result), caps[-1])
        self.assertEqual(transport_shells((0, 0, 0), full_block(2)), (0,)*5)

    def test_invalid_interfaces(self):
        for caps in ((1, 2), (0, 2, 1), (0, -1, 2), (0, 1.0), (0, True), (0,)):
            with self.assertRaises(ValueError):
                transport_shells(caps)
        for local in ((1,), (Q(1, 2), Q(1, 2)), (0, 2), (0, 1.0), (0, True)):
            with self.assertRaises(ValueError):
                transport_shells((0, 1), local)
        for values in ((0,), (0, -1), (0, 0.5)):
            with self.assertRaises(ValueError):
                suffix_expectation((0, 1), values)


if __name__ == '__main__':
    unittest.main()
