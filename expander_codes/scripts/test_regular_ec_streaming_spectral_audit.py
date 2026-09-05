import unittest

from regular_ec_streaming_spectral_audit import (
    AUXILIARY_PRIME,
    MASK64,
    auxiliary_prime_frequency_audit,
    complex_frequency_audit,
    has_distinct_differences,
    offset_digest,
    primitive_root_of_order,
    sample_topology,
    topology_digest,
)


class RegularEcStreamingSpectralAuditTests(unittest.TestCase):
    def test_deployed_seed_replays_cpp_offsets(self) -> None:
        topology = sample_topology(
            left_degree=26,
            right_degree=13,
            region_size=80659,
            seed_low=0,
            seed_high=MASK64,
        )
        self.assertEqual(topology.attempts, 1)
        self.assertTrue(has_distinct_differences(topology.offsets, 80659))
        self.assertEqual(offset_digest(topology), 0x5374E401C8EBC7EB)
        self.assertEqual(topology_digest(topology), 0xD3BA70A2C8EBC7EB)

    def test_auxiliary_root_has_exact_order(self) -> None:
        root = primitive_root_of_order(80659, AUXILIARY_PRIME)
        self.assertEqual(pow(root, 80659, AUXILIARY_PRIME), 1)
        self.assertNotEqual(pow(root, 80659 // 79, AUXILIARY_PRIME), 1)
        self.assertNotEqual(pow(root, 80659 // 1021, AUXILIARY_PRIME), 1)

    def test_reduced_all_frequency_audits(self) -> None:
        topology = sample_topology(
            left_degree=26,
            right_degree=13,
            region_size=79,
            seed_low=0,
            seed_high=MASK64,
        )
        auxiliary = auxiliary_prime_frequency_audit(
            topology, 79, batch_size=32
        )
        self.assertEqual(auxiliary["minimum_rank"], 13)
        self.assertEqual(auxiliary["deficient_frequencies"], 0)
        numerical = complex_frequency_audit(topology, 79, batch_size=32)
        self.assertEqual(numerical["counts_below"]["1e-8"], 0)


if __name__ == "__main__":
    unittest.main()
