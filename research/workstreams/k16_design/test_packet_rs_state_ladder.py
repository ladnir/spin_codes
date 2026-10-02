"""Light state-ladder orchestration and toy-moment equivalence tests only."""
from contextlib import ExitStack, redirect_stdout
from copy import deepcopy
from fractions import Fraction as Q
from io import StringIO
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from flint import arb, arb_mat, ctx, fmpq_mat
import packet_rs_state_ladder as ladder


class StateLadderTests(unittest.TestCase):
    def setUp(self):
        self.precision = ctx.prec
        ctx.prec = 256
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.sources = ladder.q1.source_snapshot()
        self.geometry = ladder.lengths.geometry(4096)
        self.beta, self.counts = ladder.lengths.exact_outer()
        self.count_hash = ladder.lengths.count_hash(self.counts)

    def tearDown(self):
        ctx.prec = self.precision

    def toy(self, bits=17):
        data = dict(bits=bits, birth_density='capped', map_sha256=f'toy-only-{bits}')
        record = dict(s=bits, return_denominator=(1 << bits)-1, map_sha256=f'toy-only-{bits}')
        return data, record

    def toy_context(self, bits=17):
        """Mock only state census/maps and moments; use actual proof arithmetic."""
        data, record = self.toy(bits)
        stack = ExitStack()
        adapter = ladder.adapter_for(bits)
        prepared = stack.enter_context(patch.object(adapter, 'prepare', return_value=(data, record)))
        stack.enter_context(patch.object(adapter, 'authenticate'))
        stack.enter_context(patch.object(ladder.q1, 'source_snapshot', return_value=self.sources))
        stack.enter_context(patch.object(ladder.low, 'validated', return_value=bits))
        stack.enter_context(patch.object(ladder.low, 'source_snapshot', return_value=self.sources))
        stack.enter_context(patch.object(ladder.q1.kernel_t64, 'local_operators', return_value=['toy']))
        # Positive triangular operators retain every terminal coordinate.
        tiny = arb(2)**-256
        matrix = arb_mat([[tiny, tiny], [0, tiny]])
        stack.enter_context(patch.object(ladder.q1, 'placement', return_value=[matrix]*33))
        stack.enter_context(patch.object(ladder.regional_power, 'placement_power', return_value=[matrix]*33))
        stack.enter_context(patch.object(ladder.q1, 'support_moments', return_value=[tiny]*65))
        stack.enter_context(patch.object(ladder.low.q2, 'pair_support_moments',
            return_value=[[arb(2)**-400]*(v+1) for v in range(65)]))
        stack.enter_context(patch.object(ladder.dense, 'physical_candidates', return_value={'toy': True}))
        stack.enter_context(patch.object(ladder.dense, 'select_physical', return_value=(matrix, (0,)*17)))
        stack.enter_context(redirect_stdout(StringIO()))
        return stack, prepared

    def test_explicit_adapter_selection_does_not_change_frozen_globals(self):
        for bits in (17, 18):
            self.assertIs(ladder.adapter_for(bits), ladder.small_inner)
        for bits in (19, 20, 21, 22):
            self.assertIs(ladder.adapter_for(bits), ladder.quadratic_inner)
        for bits in (16, 23, True, 17.0, '20', None):
            with self.subTest(bits=bits), self.assertRaises(ValueError):
                ladder.adapter_for(bits)
        self.assertIs(ladder.shared.inner, ladder.small_inner)
        self.assertEqual(ladder.small_inner.SUPPORTED_BITS, (17, 18))
        self.assertEqual(ladder.quadratic_inner.extension.STATE_BITS, 20)

    def test_state17_tail_equivalence_to_frozen_driver_for_both_methods(self):
        for method in ('exact', 'fugacity'):
            context, prepared = self.toy_context()
            with context:
                kwargs = dict(K=4096, bits=17, occupancies=[3, 5, 32], tilts=['.01', '.02'],
                    precision=256, method=method, marker_probabilities=['.5', '1'] if method == 'fugacity' else ())
                old = ladder.shared.run_tail(**kwargs)
                new = ladder.run_tail(**kwargs)
            self.assertEqual(prepared.call_count, 2)
            self.assertEqual(old['occupancy_uppers'], new['occupancy_uppers'])
            self.assertEqual(old['occupancy_choices'], new['occupancy_choices'])
            self.assertEqual(old['union_upper'], new['union_upper'])
            self.assertNotEqual(old['schema'], new['schema'])
            self.assertFalse(new['whole_code_certificate'])
            if method == 'exact':
                self.assertEqual(new['placement_backend'], 'linear')
            else:
                self.assertNotIn('placement_backend', new)

    def test_binary_dispatch_preserves_geometry_and_toy_endpoints(self):
        context, _ = self.toy_context()
        with context:
            linear, binary = ladder.q1.placement, ladder.regional_power.placement_power
            kwargs = dict(K=4096, bits=17, occupancies=[3, 5], tilts=['.01'], precision=256)
            new = ladder.run_tail(**kwargs, placement_backend='binary')
            linear.assert_not_called()
            binary.assert_called_once_with(['toy'], epochs=self.geometry.macros_per_region,
                windows=32, rounding=ladder.q1.rounded, maximum_groups=5)
            old = ladder.run_tail(**kwargs)
            linear.assert_called_once()
            self.assertEqual(new['placement_backend'], 'binary')
            self.assertEqual(old['placement_backend'], 'linear')
            self.assertEqual(new['occupancy_uppers'], old['occupancy_uppers'])
            self.assertEqual(new['union_upper'], old['union_upper'])

    def test_exact_rational_backends_agree_for_noncommuting_toy_operators(self):
        operators = [fmpq_mat([[1, 1], [0, 1]]), fmpq_mat([[1, 0], [1, 1]]),
                     fmpq_mat([[2, 1], [0, 1]])]
        self.assertNotEqual(operators[0]*operators[1], operators[1]*operators[0])
        for epochs in (1, 2, 3, 5):
            for degree in (0, min(3, 2*epochs), 2*epochs):
                options = dict(epochs=epochs, windows=2, matrix=fmpq_mat,
                               rounding=lambda value: value, maximum_groups=degree)
                linear = ladder.q1.placement(operators, **options)
                binary = ladder.regional_power.placement_power(operators, **options)
                self.assertEqual(linear, binary, (epochs, degree))

    def test_binary_recipe_propagates_to_exact_component_and_read_validation(self):
        context, prepared = self.toy_context()
        with context, patch.object(ladder, 'run_tail', wraps=ladder.run_tail) as tail:
            result = ladder.replay(self.folder/'binary.json', K=4096, bits=17,
                placement_backend='binary', precision=256, exact_occupancies=[3, 5], dense_min=3,
                q1_tilts=['.0024'], q2_tilts=['.0024'], exact_tilts=['.01'],
                dense_tilts=['.01'], marker_probabilities=['.5', '1'])
        prepared.assert_called_once()
        self.assertEqual(result['recipe']['placement_backend'], 'binary')
        self.assertEqual(tail.call_args_list[0].kwargs['placement_backend'], 'binary')
        self.assertNotIn('placement_backend', tail.call_args_list[1].kwargs)
        exact = json.loads((self.folder/'binary-exact.json').read_text())
        dense = json.loads((self.folder/'binary-dense.json').read_text())
        self.assertEqual(exact['placement_backend'], 'binary')
        self.assertNotIn('placement_backend', dense)
        _, record = self.toy()
        args = (self.geometry, 17, record, self.sources, 256, self.count_hash, self.beta)
        self.assertEqual(set(ladder._endpoints(exact, *args, placement_backend='binary')), {3, 5})
        with self.assertRaises(ValueError):
            ladder._endpoints(exact, *args, placement_backend='linear')
        missing = dict(exact)
        del missing['placement_backend']
        with self.assertRaises(ValueError):
            ladder._endpoints(missing, *args)
        with self.assertRaises(ValueError):
            ladder._endpoints(dict(dense, placement_backend='binary'), *args)
        self.assertTrue(result['all_occupancies_covered'])
        self.assertTrue(result['fresh_replay'])
        self.assertTrue(result['whole_code_certificate'])

    def test_state17_full_replay_equivalence_with_real_low_producer_and_toy_moments(self):
        context, prepared = self.toy_context()
        with context:
            kwargs = dict(K=4096, bits=17, precision=256, exact_occupancies=[3, 5], dense_min=3,
                q1_tilts=['.0024'], q2_tilts=['.0024'], exact_tilts=['.01'],
                dense_tilts=['.01'], marker_probabilities=['.5', '1'])
            old = ladder.shared.replay(self.folder/'old.json', **kwargs)
            new = ladder.replay(self.folder/'new.json', **kwargs)
        self.assertEqual(prepared.call_count, 2)
        for key in ('occupancy_uppers', 'selected_components', 'union_upper', 'tail_upper',
                    'q1_upper', 'q2_upper', 'target_met', 'whole_code_certificate'):
            self.assertEqual(old[key], new[key], key)
        self.assertEqual(new['occupancy_covered'], [1, 32])
        self.assertTrue(new['fresh_replay'])
        self.assertTrue(new['whole_code_certificate'])
        self.assertNotEqual(new['schema'], old['schema'])
        self.assertEqual(new['recipe']['placement_backend'], 'linear')

    def test_general_state_is_actually_passed_to_selected_adapter(self):
        for bits in (18, 19, 20, 21, 22):
            context, prepared = self.toy_context(bits)
            with context:
                result = ladder.run_tail(K=4096, bits=bits, occupancies=[3], tilts=['.01'], precision=256)
            prepared.assert_called_once_with(bits, birth_density='capped')
            self.assertEqual(result['state_bits'], bits)
            self.assertEqual(result['return_denominator'], (1 << bits)-1)
            self.assertEqual(result['inner_group'], f'GL({bits},2)')
            self.assertEqual(result['outer_symbol_group'], 'GL(16,2)')

    def test_invalid_options_and_incomplete_recipe_fail_before_preparation(self):
        base = dict(K=4096, bits=19, occupancies=[3], tilts=['.1'], precision=256,
                    output=None, method='exact', markers=())
        for key, bad in [('K', 1), ('bits', 23), ('occupancies', []), ('occupancies', [2]),
            ('occupancies', [5, 3]), ('occupancies', [3, 3]), ('tilts', []),
            ('tilts', ['0']), ('precision', True), ('method', 'old-metadata'), ('markers', ['0']),
            ('placement_backend', 'auto'), ('placement_backend', None)]:
            with self.subTest(key=key), self.assertRaises(ValueError):
                ladder._options(**dict(base, **{key: bad}))
        for kwargs in (dict(precision=192), dict(exact_occupancies=[3], dense_min=5)):
            with self.assertRaises(ValueError):
                ladder.recipe(4096, 20, **kwargs)
        with patch.object(ladder.quadratic_inner, 'prepare') as prepare, self.assertRaises(ValueError):
            ladder.run_tail(K=4096, bits=20, occupancies=[3], tilts=['.1'], method='fugacity', marker_probabilities=['1'])
        prepare.assert_not_called()
        for backend in ('auto', 'Binary', None, True, 1, {}):
            with self.subTest(backend=backend), patch.object(ladder.quadratic_inner, 'prepare') as prepare:
                with self.assertRaises(ValueError):
                    ladder.run_tail(K=4096, bits=20, occupancies=[3], tilts=['.1'], placement_backend=backend)
                with self.assertRaises(ValueError):
                    ladder.recipe(4096, 20, exact_occupancies=[], placement_backend=backend)
                prepare.assert_not_called()

    def test_preparation_requires_actual_width_denominator_and_unchanged_old_pins(self):
        data, record = self.toy(20)
        adapter = ladder.quadratic_inner
        with self.assertRaises(ValueError):
            ladder._prepare(ladder.small_inner, 20, 256, data, record)
        for bad in (dict(record, s=19), dict(record, return_denominator=65535)):
            with patch.object(adapter, 'authenticate'), self.assertRaises(ValueError):
                ladder._prepare(adapter, 20, 256, data, bad)
        for after, valid in (({'old': 'hash', 'new': 'hash2'}, True), ({'old': 'changed'}, False), ({}, False)):
            with (patch.object(adapter, 'authenticate'), patch.object(adapter, 'prepare', return_value=(data, record)),
                  patch.object(ladder.q1, 'source_snapshot', side_effect=[{'old': 'hash'}, after])):
                if valid:
                    self.assertEqual(ladder._prepare(adapter, 20, 256)[2], after)
                else:
                    with self.assertRaises(RuntimeError):
                        ladder._prepare(adapter, 20, 256)

    def test_tail_validator_rejects_other_schema_state_map_or_missing_occupancy(self):
        context, _ = self.toy_context(19)
        with context:
            result = ladder.run_tail(K=4096, bits=19, occupancies=[3, 5], tilts=['.01'], precision=256)
        data, record = self.toy(19)
        def checked(candidate):
            return ladder._endpoints(candidate, self.geometry, 19, record, self.sources,
                                     256, self.count_hash, self.beta)
        self.assertEqual(set(checked(result)), {3, 5})
        for key, bad in [('schema', ladder.shared.TAIL_SCHEMA), ('state_bits', 17),
            ('map_record', {}), ('return_denominator', 65535), ('source_sha256', {}),
            ('whole_code_certificate', True), ('every_requested_occupancy_checked', False),
            ('count_sha256', 'wrong'), ('occupancy_covered', [3, 4]),
            ('placement_backend', 'auto'), ('placement_backend', None)]:
            changed = dict(result, **{key: bad})
            with self.subTest(key=key), self.assertRaises(ValueError):
                checked(changed)

    def test_full_replay_failure_and_gaps_never_claim_a_certificate(self):
        for failure in ('bound', 'gap'):
            context, _ = self.toy_context()
            with context:
                actual = ladder.run_tail
                def altered_tail(**kwargs):
                    record = actual(**kwargs)
                    if kwargs['method'] == 'fugacity':
                        record = deepcopy(record)
                        if failure == 'gap':
                            record['occupancy_covered'].remove(32)
                            del record['occupancy_uppers']['32']
                        else:
                            record['occupancy_uppers']['32'] = [1, -39]
                    return record
                output = self.folder/f'{failure}.json'
                with patch.object(ladder, 'run_tail', side_effect=altered_tail):
                    if failure == 'gap':
                        with self.assertRaises(ValueError):
                            ladder.replay(output, K=4096, bits=17, exact_occupancies=[3],
                                dense_min=3, q1_tilts=['.0024'], q2_tilts=['.0024'], exact_tilts=['.01'],
                                dense_tilts=['.01'], marker_probabilities=['.5'])
                        self.assertFalse(output.exists())
                    else:
                        result = ladder.replay(output, K=4096, bits=17, exact_occupancies=[3],
                            dense_min=3, q1_tilts=['.0024'], q2_tilts=['.0024'], exact_tilts=['.01'],
                            dense_tilts=['.01'], marker_probabilities=['.5'])
                        self.assertFalse(result['target_met'])
                        self.assertFalse(result['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()
