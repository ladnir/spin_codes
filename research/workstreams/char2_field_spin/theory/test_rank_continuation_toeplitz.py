"""Independent finite checks for the new outer law and certificate interface."""

from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import product
import unittest

from batched_transfer import multiply
from rank_continuation_toeplitz_certificate import certify, middle_upper, verify
from rank_continuation_toeplitz_witnesses import WITNESSES
from subfield_certificate import HIGH, LOW
import rank_continuation_toeplitz_sizes as sizes
from rank_continuation_toeplitz_size_witnesses import WITNESSES as SIZE_WITNESSES


def toeplitz_apply(diagonals, source, polynomial):
    size = len(source)
    output = []
    for row in range(size):
        value = 0
        for column, symbol in enumerate(source):
            value ^= multiply(diagonals[row-column+size-1], symbol, polynomial)
        output.append(value)
    return tuple(output)


def mask(vector):
    return sum((value != 0) << position for position, value in enumerate(vector))


class ToeplitzSourceLawTests(unittest.TestCase):
    def test_every_nonzero_gf4_source_has_uniform_output(self):
        q, size = 4, 3
        seeds = tuple(product(range(q), repeat=2*size-1))
        expected = Counter({word: q**(size-1) for word in product(range(q), repeat=size)})
        for source in product(range(q), repeat=size):
            if any(source):
                observed = Counter(toeplitz_apply(seed, source, 0x7) for seed in seeds)
                self.assertEqual(observed, expected, source)

    def test_exact_systematic_support_multiplicities(self):
        q, size = 4, 2
        observed = Counter()
        for seed in product(range(q), repeat=2*size-1):
            for source in product(range(q), repeat=size):
                if any(source):
                    output = source + toeplitz_apply(seed, source, 0x7)
                    observed[mask(output)] += 1
        for support in range(1 << (2*size)):
            expected = q**(size-1)*(q-1)**support.bit_count() if support & ((1 << size)-1) else 0
            self.assertEqual(observed[support], expected, support)

    def test_shared_toeplitz_does_not_factor_between_rows(self):
        q, size = 2, 2
        source = (1, 1)
        images = [toeplitz_apply(seed, source, 0x3)
                  for seed in product(range(q), repeat=2*size-1)]
        probability = Fraction(sum(not any(image) for image in images), len(images))
        independent = Fraction(sum(not any(a) and not any(b) for a in images for b in images), len(images)**2)
        self.assertEqual(probability, Fraction(1, 4))
        self.assertEqual(independent, Fraction(1, 16))
        self.assertNotEqual(probability, independent)

    def test_value_measure_identity_with_causal_weight_function(self):
        q, size = 4, 2
        def weight(word):
            state, output_weight = 0, 0
            for value in word:
                output_weight += (value ^ multiply(2, state, 0x7)) != 0
                state ^= value
            return Fraction(1, 2**output_weight)
        sources = tuple(source for source in product(range(q), repeat=size) if any(source))
        left = sum(weight(source + toeplitz_apply(seed, source, 0x7))
                   for seed in product(range(q), repeat=2*size-1) for source in sources)
        right = q**(size-1)*sum(weight(source + parity)
                               for source in sources for parity in product(range(q), repeat=size))
        envelope = q**(size-1)*sum(weight(word) for word in product(range(q), repeat=2*size))
        self.assertEqual(left, right)
        self.assertLessEqual(right, envelope)


class DirectedReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = certify(WITNESSES)

    def test_complete_replay(self):
        self.assertEqual(verify(self.record)["integer_bits"], 80)
        self.assertLess(Fraction(self.record["failure_upper"]), Fraction(1, 2**80))

    def test_reject_understated_endpoint(self):
        record = deepcopy(self.record)
        record["middle"][0]["endpoint_upper"][0] = "0"
        with self.assertRaisesRegex(ValueError, "understated endpoint"):
            verify(record)

    def test_reject_inconsistent_displayed_endpoint(self):
        record = deepcopy(self.record)
        record["middle"][0]["endpoint_upper"][0] = "1"
        with self.assertRaisesRegex(ValueError, "inconsistent with displayed endpoints"):
            verify(record)

    def test_reject_incomplete_coverage_and_nonpositive_tilt(self):
        record = deepcopy(WITNESSES)
        record["middle"][0]["first"] += 1
        with self.assertRaisesRegex(ValueError, "coverage"):
            certify(record)
        with self.assertRaises(ValueError):
            middle_upper(16, "-1/10", "1/2")

    def test_reject_margin_inconsistent_with_inflated_bound(self):
        for location in ("total", "sparse", "middle"):
            record = deepcopy(self.record)
            item = record if location == "total" else record[location][0]
            item["failure_upper"] = "1"
            with self.assertRaisesRegex(ValueError, "margin inconsistent with displayed upper bound"):
                verify(record)

    def test_reject_sparse_total_inconsistent_with_entries(self):
        record = deepcopy(self.record)
        record["sparse"][0]["failure_upper"] = "5e-25"
        with self.assertRaisesRegex(ValueError, "sparse total inconsistent"):
            verify(record)

    def test_restore_decimal_exponents_on_success_and_error(self):
        before = [(context.Emax, context.Emin) for context in (HIGH, LOW)]
        certify(WITNESSES)
        self.assertEqual([(context.Emax, context.Emin) for context in (HIGH, LOW)], before)
        bad = deepcopy(WITNESSES)
        bad["config"]["state"] = 2
        with self.assertRaises(ValueError):
            certify(bad)
        self.assertEqual([(context.Emax, context.Emin) for context in (HIGH, LOW)], before)

    def test_no_diagonal_profile_rejects_scalar_metadata(self):
        self.assertEqual(self.record["config"]["diagonal"], "none")
        self.assertEqual(self.record["format"], "gf16-systematic-toeplitz64-gl3-nodiagonal-v1")
        record = deepcopy(self.record)
        record["config"]["diagonal"] = "independent-uniform-nonzero"
        with self.assertRaisesRegex(ValueError, "incorrect construction"):
            verify(record)


class DirectedSizeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records = {k_log: sizes.certify(k_log, SIZE_WITNESSES[k_log]) for k_log in (16, 18)}

    def test_smaller_sizes_have_complete_verified_coverage(self):
        for k_log, bits in ((16, 84), (18, 82)):
            record = self.records[k_log]
            self.assertEqual(sizes.verify(record)["integer_bits"], bits)
            self.assertLess(Fraction(record["failure_upper"]), Fraction(1, 2**bits))
            self.assertEqual(record["middle"][-1]["last"], (1 << k_log)//32)
            self.assertEqual(record["config"]["diagonal"], "none")

    def test_size_configuration_and_thresholds(self):
        for k_log, threshold in ((16, 13107), (18, 52428), (20, 209715)):
            config = sizes.configuration(k_log)
            self.assertEqual(config["n"], 2 << k_log)
            self.assertEqual(config["bad_weight_threshold"], threshold)
        with self.assertRaises(ValueError):
            sizes.configuration(17)

    def test_wrapper_matches_fixed_k20_replay(self):
        fixed = certify(WITNESSES)
        wrapped = sizes.certify(20, WITNESSES)
        for field in ("config", "sparse", "middle", "sparse_upper", "failure_upper", "integer_bits"):
            self.assertEqual(wrapped[field], fixed[field], field)

    def test_size_verifier_rejects_changed_size_or_threshold(self):
        for field, value in (("k", 262144), ("bad_weight_threshold", 13108), ("diagonal", "independent-uniform-nonzero")):
            record = deepcopy(self.records[16])
            record["config"][field] = value
            with self.assertRaisesRegex(ValueError, "incorrect construction or size"):
                sizes.verify(record)

    def test_size_verifier_rejects_inconsistent_displays(self):
        record = deepcopy(self.records[16])
        record["failure_upper"] = "1"
        with self.assertRaisesRegex(ValueError, "inconsistent margin"):
            sizes.verify(record)
        record = deepcopy(self.records[16])
        record["middle"][0]["endpoint_upper"][0] = "0"
        with self.assertRaisesRegex(ValueError, "understated endpoint"):
            sizes.verify(record)
        record = deepcopy(self.records[16])
        record["middle"][0]["endpoint_upper"][0] = "1"
        with self.assertRaisesRegex(ValueError, "inconsistent with displayed endpoints"):
            sizes.verify(record)

    def test_size_replay_rejects_coverage_gaps(self):
        record = deepcopy(SIZE_WITNESSES[16])
        record["middle"][-1]["last"] -= 1
        with self.assertRaisesRegex(ValueError, "incomplete coverage"):
            sizes.certify(16, record)


if __name__ == "__main__":
    unittest.main()
