"""Small exact route enumerations, not production probability screens."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
import unittest

from hd_shifted_common_source_partition import (
    IDENTITY, ZERO, good_probability, matrix_add, matrix_product,
    six_type_slot_terms, tiny_regional_partition)


def packet_occupancies(slots):
    potential, actual = [0]*4, [0]*4
    for slot, (p, a) in enumerate(slots):
        base = 2*(slot % 2)
        for lane in range(2):
            potential[base+lane] += (p >> lane) & 1
            actual[base+lane] += (a >> lane) & 1
    return tuple(potential), tuple(actual)


def caps_fixture(length=2):
    # Deliberately noncommuting and not stochastic; algebra must not rely
    # on sorted region labels, state reset, or probability row sums.
    return {o: (F(1+o[0]), F(1+2*o[1]), F(2+3*o[2]), F(1+o[3]))
            for o in product(range(length+1), repeat=4)}


class CommonSourcePartitionTests(unittest.TestCase):
    def test_exact_oriented_count_for_all_tiny_histograms(self):
        for length, packets, cap in ((2, 1, 1), (1, 2, 1), (2, 1, 0)):
            pairs = 2*length*packets
            totals, goods = Counter(), Counter()
            for source_masks in product(range(4), repeat=pairs):
                histogram = (sum(x in (1, 2) for x in source_masks), source_masks.count(3))
                totals[histogram] += 1
                good = True
                for packet in range(packets):
                    local = source_masks[2*length*packet:2*length*(packet+1)]
                    potential, _ = packet_occupancies([(x, 0) for x in local])
                    good &= max(potential) <= cap
                goods[histogram] += good
            for (s, f), total in totals.items():
                self.assertEqual(good_probability(s, f, block_length=length, packets=packets, cap=cap),
                                 F(goods[s, f], total))

    def test_unrestricted_and_small_source_events(self):
        for s in range(5):
            for f in range(5-s):
                self.assertEqual(good_probability(s, f, block_length=2, packets=1, cap=2), 1)
        # A full source pair spans two blocks, so it fits cap1 with certainty.
        self.assertEqual(good_probability(0, 1, block_length=2, packets=1, cap=1), 1)
        # Two singleton pairs can put two positions in one block.
        self.assertLess(good_probability(2, 0, block_length=2, packets=1, cap=1), 1)

    def test_slot_orientations_and_joint_masks(self):
        masses = Counter()
        for kind, potential, actual, mass in six_type_slot_terms():
            self.assertEqual(actual & ~potential, 0)
            masses[kind] += mass
        self.assertEqual(dict(masses), dict.fromkeys(range(6), F(1)))
        for position in range(16):
            for half in range(2):
                # Native route index2*position+half fills these record lanes.
                blocks = 2*half, 2*half+1
                self.assertNotEqual(*blocks)

    def test_noncommuting_conditional_operator_exact_route_enum(self):
        caps = caps_fixture()
        self.assertNotEqual(matrix_product(caps[(1, 0, 0, 0)], caps[(0, 1, 0, 0)]),
                            matrix_product(caps[(0, 1, 0, 0)], caps[(1, 0, 0, 0)]))
        target = (6, 0, 1, 0, 1, 0)
        observed = tiny_regional_partition(caps, target, block_length=2, packets=2, cap=1)
        sums, count_good, count = {True: ZERO, False: ZERO}, 0, 0
        for singleton in range(8):
            for full in range(8):
                if singleton == full:
                    continue
                for first, second in product((1, 2), repeat=2):
                    slots = [(0, 0)]*8
                    slots[singleton] = first, first
                    slots[full] = 3, second
                    result, good = IDENTITY, True
                    for packet in range(2):
                        potential, actual = packet_occupancies(slots[4*packet:4*(packet+1)])
                        good &= max(potential) <= 1
                        result = matrix_product(result, caps[actual])
                    sums[good] = matrix_add(sums[good], result)
                    count_good += good
                    count += 1
        self.assertEqual(observed['good_probability'], F(count_good, count))
        for event, key in ((True, 'good_weighted'), (False, 'bad_weighted')):
            self.assertEqual(observed[key], tuple(x/count for x in sums[event]))
        self.assertEqual(observed['good_conditional'], tuple(x/count_good for x in sums[True]))
        self.assertEqual(observed['bad_conditional'], tuple(x/(count-count_good) for x in sums[False]))
        self.assertFalse(observed['alpha_applied'])
        self.assertFalse(observed['distance_certificate'])

    def test_identity_caps_reproduce_common_event_probability(self):
        caps = dict.fromkeys(product(range(3), repeat=4), IDENTITY)
        for target in ((2, 1, 0, 0, 1, 0), (0, 0, 0, 0, 0, 4), (4, 0, 0, 0, 0, 0)):
            result = tiny_regional_partition(caps, target, block_length=2, packets=1, cap=1)
            p = result['good_probability']
            self.assertEqual(result['good_weighted'], tuple(p*x for x in IDENTITY))
            self.assertEqual(result['bad_weighted'], tuple((1-p)*x for x in IDENTITY))
            if not p:
                self.assertIsNone(result['good_conditional'])
            if p == 1:
                self.assertIsNone(result['bad_conditional'])

    def test_guards(self):
        for args in ((-1, 0), (0, -1), (5, 0)):
            with self.assertRaises(ValueError):
                good_probability(*args, block_length=2, packets=1, cap=1)
        with self.assertRaises(ValueError):
            good_probability(0, 0, cap=17)
        with self.assertRaises(ValueError):
            tiny_regional_partition({}, (1, 0, 0, 0, 0, 0))
        with self.assertRaises(ValueError):
            tiny_regional_partition({}, (6, 0, 0, 0, 0, 0), block_length=3, packets=1)


if __name__ == '__main__':
    unittest.main()
