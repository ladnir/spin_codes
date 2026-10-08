"""Proof replay and deliberate-corruption checks; not encoder benchmarks."""
from contextlib import redirect_stdout
from copy import deepcopy
import io
import unittest

from return_certificate import COMPANION, generate, verify


@unittest.skipUnless(COMPANION.is_dir(), "adjacent field-universality checkout is required")
class ReturnCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with redirect_stdout(io.StringIO()):
            cls.witness = generate(1 << 17, 128, (1 << 17) // 10, 512, 40, COMPANION)

    def test_complete_k16_certificate(self):
        result = verify(self.witness)
        self.assertEqual(result["minimum_distance_at_least"], 13108)
        self.assertGreater(result["margin_bits"], 41.29)
        self.assertTrue(result["complete_distance_certificate"])

    def test_missing_bulk_interval_is_rejected(self):
        witness = deepcopy(self.witness)
        del witness["intervals"][3]
        with self.assertRaises(AssertionError):
            verify(witness)

    def test_sparse_claim_is_recomputed(self):
        witness = deepcopy(self.witness)
        witness["sparse_failure"] = "0"
        with self.assertRaises(AssertionError):
            verify(witness)

    def test_bulk_claim_is_recomputed(self):
        witness = deepcopy(self.witness)
        witness["intervals"][0]["interval_log2_upper"] = "-1000000"
        with self.assertRaises(AssertionError):
            verify(witness)

    def test_insufficient_margin_is_rejected(self):
        witness = deepcopy(self.witness)
        witness["required_margin"] = 60
        with self.assertRaises(AssertionError):
            verify(witness)


if __name__ == "__main__":
    unittest.main()
