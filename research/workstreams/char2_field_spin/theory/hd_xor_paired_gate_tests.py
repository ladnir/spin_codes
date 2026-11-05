"""Partition and wiring tests; synthetic caps are not distance certificates."""
from itertools import combinations
from math import comb, log2
from types import ModuleType
import sys
import unittest
from unittest.mock import patch

import numpy as np

import hd_xor_paired_gate as gate


class PairedGateTests(unittest.TestCase):
    def test_default_partitions_are_complete(self):
        for r in range(1, 18):
            self.assertTrue(gate.validate_plan(gate.default_plan(r), r))

    def test_partial_overlapping_and_invalid_partitions(self):
        self.assertFalse(gate.validate_plan([((0, 0), (1, 1))], 17))
        for plan in ([((0, 1), (1, 4)), ((1, 2), (3, 5))],
                     [((0, 9), (1, 17))], [((0, 1), (0, 17))]):
            with self.assertRaises(ValueError):
                gate.validate_plan(plan, 17)
        with self.assertRaises(ValueError):
            gate.validate_plan([], 18)

    def test_subset_pair_cap_for_both_matchings(self):
        # Independently enumerate fixed active sets and every output subset.
        n = 8
        def pairs(points, parity):
            return sum(i in points and (i+1) % n in points
                       for i in range(parity, n, 2))
        for r in range(1, 6):
            for source in combinations(range(n), r):
                F = max(pairs(source, 0), pairs(source, 1))
                for mask in range(1 << r):
                    output = {source[i] for i in range(r) if mask >> i & 1}
                    for parity in (0, 1):
                        self.assertLessEqual(pairs(output, parity), min(F, len(output)//2))

    def test_entrywise_envelope_and_monotonicity(self):
        r = 5
        table = {(m, f): np.array([[f-m, 2*f], [-m, m-f]], dtype=float)
                 for m in range(r+1) for f in range(m//2+1)}
        previous = None
        for F in range(r//2+1):
            result = gate.restricted_regions(table, r, F)
            for (m, f), matrix in table.items():
                if f <= F:
                    self.assertTrue(np.all(result[m] >= matrix))
            if previous is not None:
                self.assertTrue(np.all(result >= previous))
            previous = result

    def test_directed_envelope_dominates_overlapping_intervals(self):
        from flint import arb
        table = {(0, 0): (arb(1), arb(0), arb(0), arb('0.5 +/-0.01')),
                 (1, 0): (arb('1 +/-0.1'),)*4,
                 (2, 0): (arb('1.1 +/-0.2'),)*4,
                 (2, 1): (arb('1.0 +/-0.15'),)*4}
        output = gate.restricted_regions(table, 2, 1, directed=True)
        for row in range(4):
            self.assertGreaterEqual(output[2][row], table[2, 0][row].upper())
            self.assertGreaterEqual(output[2][row], table[2, 1][row].upper())

    def test_geometry_counts_enter_once_and_partial_scope_stays_partial(self):
        r = 3
        route = ModuleType('hd_xor_pair_route')
        table = {(m, f): np.zeros((2, 2)) for m in range(r+1) for f in range(m//2+1)}
        route.pair_regions_screen = lambda z, maximum: table
        route.arb_pair_regions = lambda z, maximum: None
        strata = ModuleType('hd_xor_geometry_strata')
        strata.MultiplicityCapsFloat = lambda *args: lambda mu: {}
        strata.MultiplicityCapsArb = lambda *args: lambda mu: {}
        source_counts = {0: comb(gate.REGION_LENGTH, r)-10, 1: 10}
        pair_counts = {1: 100, 2: 28**r-128, 3: 28}
        strata.geometry_stratum_counts = lambda groups: pair_counts
        geometry = ModuleType('hd_xor_group_geometry')
        geometry.matching_max_counts = lambda size, groups: source_counts
        modules = dict(hd_xor_pair_route=route, hd_xor_geometry_strata=strata,
                       hd_xor_group_geometry=geometry)
        with patch.dict(sys.modules, modules), patch.object(gate, 'joint_log_count', return_value=0.):
            result = gate.gate(r)
            prefix = -log2(gate.Q-1)-gate.LIMIT*log2(0.996548217)+63.5
            expected = prefix+log2(comb(gate.REGION_LENGTH, r))+r*log2(28)
            self.assertAlmostEqual(result['failure_log2_screen'], expected, places=9)
            self.assertTrue(result['complete_geometry_partition'])
            self.assertFalse(result['whole_code_certificate'])
            partial = gate.gate(r, plan=[((0, 0), (1, 1))])
            self.assertFalse(partial['complete_geometry_partition'])
            expected_partial = prefix+log2(source_counts[0])+log2(pair_counts[1])
            self.assertAlmostEqual(partial['failure_log2_screen'], expected_partial, places=9)


if __name__ == '__main__':
    unittest.main()
