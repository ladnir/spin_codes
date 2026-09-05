#!/usr/bin/env python3

import unittest

from prime_field_biregular_ec_certificate import validate_coverage
from prime_field_biregular_ec_singleton_certificate import (
    BLOCKS,
    EXACT_MARKERS,
    K,
)


class PrimeFieldBiregularECSingletonCertificateTests(unittest.TestCase):
    def test_partition_covers_every_nonzero_support(self) -> None:
        certificate = {
            "parameters": {"k": K},
            "exact_bands": [
                {"lo": lo, "hi": hi} for lo, hi in EXACT_MARKERS
            ],
            "blocks": [{"lo": lo, "hi": hi} for lo, hi in BLOCKS],
            "full_support": {"r": K},
        }
        validate_coverage(certificate)


if __name__ == "__main__":
    unittest.main()
