"""Independent checks of the completed paired-map receipt, not a full replay.

The union arithmetic uses Python integers/Fractions rather than the replay's
dyadic helpers. The q=1 check uses 32 physical t64 steps, not its macro path.
"""
from fractions import Fraction
import hashlib
import json
import os
from math import comb
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'proof'))
import disjoint_pair_maps as maps
from flint import arb, ctx

RECEIPT = Path(os.environ.get('SPIN_PAIRED16_RECEIPT',
    str(HERE.parent / 'proof/disjoint-pairs-whole-v2-p256.json')))


def rational(pair):
    if (not isinstance(pair, list) or len(pair) != 2 or
            any(type(x) is not int for x in pair) or pair[0] <= 0):
        raise ValueError('positive integer dyadic pair required')
    m, e = pair
    return Fraction(m << max(e, 0), 1 << max(-e, 0))


def normalize(m, e):
    while not m & 1:
        m >>= 1
        e += 1
    return m, e


def outer_shells():
    """Inclusion-exclusion on shortened supports, then packet convolution."""
    q = 65536
    # Four packed GF16 rows are an MDS[16,8] code over a 2^16 alphabet.
    symbols = [comb(16, h) * sum((-1)**j * comb(h, j) *
               q**max(0, h-j-8) for j in range(h+1)) for h in range(1, 17)]
    assert sum(symbols) == (1 << 128)-1
    packet = [0] + [comb(4, k)*15**k for k in range(1, 5)]
    power, counts = [1], [Fraction(0)]*65
    for h, number in enumerate(symbols, 1):
        next_power = [0]*(len(power)+4)
        for i, a in enumerate(power):
            for j, b in enumerate(packet):
                next_power[i+j] += a*b
        power = next_power
        for u, numerator in enumerate(power):
            counts[u] += Fraction(number*numerator, (q-1)**h)
    assert sum(counts) == (1 << 128)-1
    return counts


class CompleteReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        ctx.prec = 256
        cls.record = json.loads(RECEIPT.read_bytes())
        cls.wrapper, cls.map_record = maps.prepare()

    def test_complete_exact_geometry(self):
        record = self.record
        expected = dict(K=65536, N=131072, threshold=13107,
            target_minimum_distance=13108, target_margin_bits=40,
            group_count=512, group_dimension=128, regions=64, packet_bits=4,
            physical_t=64, state_bits=16, physical_steps=2048,
            physical_steps_per_region=32, macro_steps_per_region=16,
            precision=256, zero_initial_state=True, final_flush=False,
            fresh_replay=True, whole_code_certificate=True,
            all_occupancies_covered=True, target_met=True,
            state_continuity='retained_between_every_step_and_region')
        for name, value in expected.items():
            self.assertEqual(record[name], value, name)
        self.assertEqual(set(record['occupancy_uppers']), set(map(str, range(1, 513))))
        self.assertEqual(set(record['tail_choices']), set(map(str, range(3, 513))))
        for first, last, tilts in record['tail_recipes']:
            for q in range(first, last+1):
                self.assertIn(record['tail_choices'][str(q)], tilts)

    def test_all_saved_source_hashes(self):
        pins = self.record['source_sha256']
        self.assertGreater(len(pins), 10)
        for filename, expected in pins.items():
            self.assertEqual(hashlib.sha256(Path(filename).read_bytes()).hexdigest(),
                             expected, filename)
        for basename in ('reproduce_disjoint_pairs.py', 'disjoint_pair_maps.py',
                         'kernel_t64.py', 'kernel_maps.py', 'kernel_birth_density.py',
                         'sparse_kernel.py', 'packet_q1.py', 'packet_q2.py',
                         'packet_uniform_tail.py', 'packet_rs_length_whole.py'):
            self.assertTrue(any(Path(p).name == basename for p in pins), basename)

    def test_receipt_map_is_fresh_declared_map(self):
        self.assertEqual(self.record['map_record'], self.map_record)
        self.assertEqual(self.wrapper['physical_steps'], 2)
        self.assertEqual(self.wrapper['physical_windows'], 16)
        self.assertEqual(self.wrapper['macro_windows'], 32)

    def test_exact_sum_single_upward_rounding_and_target(self):
        record = self.record
        endpoints = list(record['occupancy_uppers'].values())
        base = min(pair[1] for pair in endpoints)
        total = sum(m << (e-base) for m, e in endpoints)
        exact = normalize(total, base)
        self.assertEqual(exact, (int(record['union_exact_mantissa_hex'], 16),
                                 record['union_exact_exponent']))
        # ceil to a 256-bit significand, then canonicalize; independent of
        # packet_rs_length_whole.dyadic_sum/rounded_endpoint/dyadic_less.
        shift = max(0, total.bit_length()-record['precision'])
        ceiling = (total+(1 << shift)-1) >> shift
        self.assertEqual(normalize(ceiling, base+shift),
                         normalize(*record['union_upper']))
        exact_fraction = sum(map(rational, endpoints), Fraction(0))
        upper = rational(record['union_upper'])
        self.assertLessEqual(exact_fraction, upper)
        self.assertLess(upper, Fraction(1, 1 << 40))

    def test_fresh_q1_physical_geometry_and_shell_fold(self):
        q1 = maps.screen.q1
        kernel = q1.kernel_t64
        physical = kernel.authenticate(self.wrapper)
        shells = outer_shells()
        self.assertEqual(tuple(shells), tuple(maps.screen.tail.uniform_envelope(16, 8, 4)[1]))
        best = [arb(1)]*65
        for tilt in self.record['q1_tilts']:
            z = (-kernel.aq(Fraction(tilt))).exp()
            local = kernel.sparse_kernel.outward_at_z(physical, z)
            local = kernel.kernel_birth_density.refine_local(physical, local, z, Fraction(1, 2))
            regional = q1.placement(local, epochs=32, windows=16,
                                    maximum_groups=1, rounding=q1.rounded)
            moments = q1.support_moments(regional[0], regional[1], 64)
            factor = (kernel.aq(Fraction(tilt))*13107).exp()
            best = [min(old, kernel.up(factor*moment)) for old, moment in zip(best, moments)]
        fresh = kernel.up(512*sum((kernel.aq(c)*p for c, p in zip(shells, best)), arb(0)))
        saved = kernel.aq(rational(self.record['occupancy_uppers']['1']))
        self.assertLess(fresh, arb(2)**-40)
        self.assertLess(abs(fresh/saved-1), arb(2)**-200)
        exact_fold = 512*sum((count*rational(pair) for count, pair in
            zip(shells, self.record['q1_support_probability_uppers'])), Fraction(0))
        self.assertLessEqual(exact_fold, rational(self.record['occupancy_uppers']['1']))


if __name__ == '__main__':
    unittest.main()
