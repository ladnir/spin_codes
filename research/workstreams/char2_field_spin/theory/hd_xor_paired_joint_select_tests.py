"""Mocked finite selection and fresh-replay tests; no expensive cap searches."""
from copy import deepcopy
from fractions import Fraction as F
from math import log2
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb, ctx

import hd_xor_paired_joint_gate as gate
import hd_xor_paired_joint_select as select


SOURCES = {(0, 0): 31, (8, 8): 7}
PAIRS = {mu: mu+1 for mu in range(1, 18)}


def seed():
    rows = []
    for source in SOURCES:
        for low, high in ((1, 8), (9, 17)):
            rows.append(dict(source_matching_counts=list(source), mu_range=[low, high],
                mark=str(F(996623217, 10**9) if source == (0, 0) and low == 9 else select.MARKS[0]),
                log2_potential='-127/2', failure_log2_screen=-100. if source == (0, 0) and low == 1 else 1e6,
                source_supports=-999, constituent_pairs=-999))
    return dict(construction=gate.CONSTRUCTION, active_groups=17, active_rows_per_group=2,
                directed=False, whole_code_certificate=False, complete_geometry_partition=True,
                source_sha256='stale-source-hash', entries=rows)


class JointSelectTests(unittest.TestCase):
    def test_scope_and_partition_validation_precede_numerics(self):
        valid = seed()
        bad = []
        for key, value in (('construction', 'another-construction'), ('active_groups', 16),
                           ('active_rows_per_group', 1), ('directed', True)):
            candidate = deepcopy(valid)
            candidate[key] = value
            bad.append(candidate)
        missing_source = deepcopy(valid)
        missing_source['entries'] = missing_source['entries'][:2]
        bad.append(missing_source)
        missing_mu = deepcopy(valid)
        missing_mu['entries'].pop()
        bad.append(missing_mu)
        overlap = deepcopy(valid)
        overlap['entries'].append(deepcopy(overlap['entries'][0]))
        bad.append(overlap)
        with patch.object(select, 'matching_joint_counts', return_value=SOURCES), \
             patch.object(select, 'geometry_stratum_counts', return_value=PAIRS), \
             patch.object(select, 'pair_regions_screen') as expensive:
            for candidate in bad:
                with self.subTest(candidate=candidate):
                    with self.assertRaises(ValueError):
                        select.select(candidate)
            expensive.assert_not_called()

    def test_finite_bank_keeps_original_mark_and_fixed_source_identity(self):
        original = seed()
        preserved = deepcopy(original)
        progress = []
        class Evaluator:
            def __init__(self, table, logp, groups):
                self.z, self.source = table
                self.logp = logp
            def __call__(self, mu):
                return self.z, self.source, self.logp, mu
        def cost(args, groups):
            z, source, logp, mu = args
            target = float(F(996623217, 10**9) if source == (0, 0) else select.MARKS[1 if mu == 8 else -1])
            # Cancel the common prefix so the synthetic winner is explicit.
            return gate.LIMIT*log2(z)+logp+1e6*abs(z-target)
        with patch.object(select, 'matching_joint_counts', return_value=SOURCES), \
             patch.object(select, 'geometry_stratum_counts', return_value=PAIRS), \
             patch.object(select, 'pair_regions_screen', side_effect=lambda z, groups: z) as raw, \
             patch.object(select, 'ordered_tables', side_effect=lambda z, a, b: (z, (a, b))), \
             patch.object(select, 'MultiplicityCapsFloat', side_effect=Evaluator) as evaluator, \
             patch.object(select, 'joint_log_count', side_effect=cost):
            result = select.select(original, progress=progress.append)
        self.assertEqual(original, preserved)
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertNotIn('source_sha256', result)
        self.assertEqual(result['selector_new_evaluations'], 22)
        self.assertEqual(raw.call_count, 8)  # Seven bank marks plus the old off-bank mark.
        self.assertEqual(evaluator.call_count, 15)  # Per-source caches never cross source strata.
        self.assertEqual(len(progress), 3)
        self.assertEqual(result['entries'][0], original['entries'][0])
        expected = (F(996623217, 10**9), select.MARKS[1], select.MARKS[-1])
        for row, winner in zip(result['entries'][1:], expected):
            self.assertEqual(F(row['mark']), winner)
            self.assertEqual(row['log2_potential'], '-127/2')
            self.assertEqual(len(row['finite_trials']), 8 if row['source_matching_counts'] == [0, 0] else 7)
            self.assertLess(row['failure_log2_screen'], 1e6)
        actual_union = float(np.logaddexp2.reduce([row['failure_log2_screen'] for row in result['entries']]))
        self.assertEqual(result['failure_log2_screen'], actual_union)

    def test_strong_cells_do_not_trigger_trials_but_remain_in_union(self):
        proposal = seed()
        for row in proposal['entries']:
            row['failure_log2_screen'] = -100.
        with patch.object(select, 'matching_joint_counts', return_value=SOURCES), \
             patch.object(select, 'geometry_stratum_counts', return_value=PAIRS), \
             patch.object(select, 'pair_regions_screen') as expensive:
            result = select.select(proposal)
        expensive.assert_not_called()
        self.assertEqual(result['selector_new_evaluations'], 0)
        self.assertEqual(len(result['entries']), 4)
        self.assertEqual(result['failure_log2_screen'], -98.)

    def test_custom_literal_mark_bank_and_invalid_banks(self):
        marks = ('998/1000', '999/1000')
        proposal = seed()
        with patch.object(select, 'matching_joint_counts', return_value=SOURCES), \
             patch.object(select, 'geometry_stratum_counts', return_value=PAIRS), \
             patch.object(select, 'pair_regions_screen', return_value=None) as raw, \
             patch.object(select, 'ordered_tables', return_value=None), \
             patch.object(select, 'MultiplicityCapsFloat', return_value=lambda mu: {}), \
             patch.object(select, 'joint_log_count', return_value=0.):
            result = select.select(proposal, marks=marks)
            self.assertEqual(result['selector_marks'], ['499/500', '999/1000'])
            self.assertEqual(result['selector_new_evaluations'], 9)
            self.assertEqual(raw.call_count, 4)
            for previous, row in zip(proposal['entries'][1:], result['entries'][1:]):
                actual = {F(trial['mark']) for trial in row['finite_trials']}
                self.assertEqual(actual, {F(mark) for mark in marks} | {F(previous['mark'])})
                self.assertEqual(F(row['mark']), F(999, 1000))
            raw.reset_mock()
            for bank, cutoff in (((), -60.), (('0',), -60.), (('1',), -60.),
                                 (('-1/3',), -60.), (marks, float('inf')), (marks, float('nan'))):
                with self.subTest(bank=bank, cutoff=cutoff):
                    with self.assertRaises(ValueError):
                        select.select(proposal, marks=bank, cutoff=cutoff)
            raw.assert_not_called()

    def test_directed_replay_ignores_saved_estimates_counts_and_trials(self):
        rows = seed()['entries']
        for row in rows:
            row.update(failure_log2_screen=-1e200, failure_log2_upper='-999999',
                       source_supports=0, constituent_pairs=0,
                       finite_trials=[dict(mark='0', failure_log2_screen=-1e300)])
        with patch.object(gate, 'matching_joint_counts', return_value=SOURCES), \
             patch.object(gate, 'geometry_stratum_counts', return_value=PAIRS), \
             patch.object(gate, 'arb_pair_regions', return_value=None) as raw, \
             patch.object(gate, 'ordered_tables', return_value=None), \
             patch.object(gate, 'MultiplicityCapsArb', return_value=lambda mu: {}), \
             patch.object(gate, 'joint_count_arb', side_effect=lambda caps, groups: arb(1)) as count:
            result = gate.run(precision=256, witness_rows=rows)
        self.assertTrue(result['directed'])
        self.assertTrue(result['complete_geometry_partition'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertEqual(count.call_count, 4)
        self.assertEqual(raw.call_count, 2)
        self.assertEqual({call.args[0] for call in raw.call_args_list}, {F(row['mark']) for row in rows})
        for row in result['entries']:
            source = tuple(row['source_matching_counts'])
            low, high = row['mu_range']
            self.assertEqual(row['source_supports'], SOURCES[source])
            self.assertEqual(row['constituent_pairs'], sum(PAIRS[mu] for mu in range(low, high+1)))
            self.assertNotIn('finite_trials', row)
            self.assertGreater(F(row['failure_log2_upper']), 0)
        from hd_shifted_coverage import rational_ceiling
        with ctx.workprec(256):
            exact = arb(0)
            for row in rows:
                z, p = F(row['mark']), F(row['log2_potential'])
                potential = arb(2)**(arb(p.numerator)/p.denominator)
                low, high = row['mu_range']
                exact += (SOURCES[tuple(row['source_matching_counts'])]*sum(PAIRS[mu] for mu in range(low, high+1))
                          *(arb(z.numerator)/z.denominator)**(-gate.LIMIT)/255/potential)
            self.assertEqual(result['failure_log2_upper'], str(rational_ceiling(exact.log()/arb(2).log())))


if __name__ == '__main__':
    unittest.main()
