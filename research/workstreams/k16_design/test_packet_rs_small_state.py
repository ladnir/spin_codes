"""Light toy algebra and fail-closed replay tests; no large map census."""
from contextlib import ExitStack, redirect_stdout
from copy import deepcopy
from fractions import Fraction as Q
from io import StringIO
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from flint import arb, arb_mat, ctx
import packet_rs_small_state as proof


def value(pair):
    return Q(pair[0])*Q(2)**pair[1]


class SmallStateProofTests(unittest.TestCase):
    def setUp(self):
        self.precision = ctx.prec
        ctx.prec = 256
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.sources = proof.q1.source_snapshot()
        self.data = dict(bits=17, birth_density='capped', map_sha256='toy-control-flow-only')
        self.map_record = dict(s=17, return_denominator=131071, full_state_census=True,
                              state_count=131072, map_sha256='toy-control-flow-only')
        self.beta, self.counts = proof.lengths.exact_outer()
        self.count_hash = proof.lengths.count_hash(self.counts)
        self.geom = proof.lengths.geometry(4096)

    def tearDown(self):
        ctx.prec = self.precision

    def component(self, q=None, *, qs=None, method='exact', endpoint=(1, -64)):
        record = proof._common(self.geom, 17, self.map_record, 256, self.sources)
        record.update(count_kind='exact_expected_shells', count_sha256=self.count_hash)
        if q is not None:
            record.update(schema=f'finite-packet-rs-state-q{q}-1', occupancy_covered=[q],
                all_supports_covered=True, all_two_group_support_pairs_covered=q == 2)
            record[f'q{q}_upper'] = list(endpoint)
        else:
            qs = list(qs)
            record.update(schema=proof.TAIL_SCHEMA, method=method,
                return_denominator=131071, birth_density='capped',
                occupancy_covered=qs, occupancy_uppers={str(q): list(endpoint) for q in qs},
                every_requested_occupancy_checked=True, every_shell_checked=True, beta=str(self.beta))
        return record

    def endpoints(self, record):
        return proof._component_endpoints(record, self.geom, 17, self.map_record,
                                          self.sources, 256, self.count_hash)

    def test_options_reject_bad_scope_before_preparation(self):
        base = dict(K=4096, bits=17, occupancies=[3, 5], tilts=['.1'], precision=256,
                    output=None, method='exact', markers=())
        for key, bad in [('K', 4097), ('bits', 16), ('bits', 20), ('bits', True),
            ('occupancies', []), ('occupancies', [5, 3]), ('occupancies', [3, 3]),
            ('occupancies', [2]), ('occupancies', [33]), ('occupancies', [True]),
            ('tilts', []), ('tilts', ['0']), ('tilts', ['.5', '1/2']), ('tilts', '.1'),
            ('precision', True), ('precision', 127), ('method', 'unknown'),
            ('markers', '.5'), ('markers', ['0']), ('markers', ['.5', '1/2'])]:
            with self.subTest(key=key, bad=bad), self.assertRaises(ValueError):
                proof._options(**dict(base, **{key: bad}))
        for markers in ([], ['1']):
            with self.assertRaises(ValueError):
                proof._options(**dict(base, method='fugacity', markers=markers))
        proof._options(**dict(base, method='fugacity', occupancies=[32], markers=['1']))
        with patch.object(proof.inner, 'prepare') as prepare, self.assertRaises(ValueError):
            proof.run_tail(K=4096, bits=16, occupancies=[3], tilts=['.1'])
        prepare.assert_not_called()

    def test_recipe_requires_complete_coverage_and_fresh_precision(self):
        plan = proof.recipe(131072, 17, exact_occupancies=list(range(3, 189)), dense_min=189)
        self.assertEqual(plan['state_bits'], 17)
        self.assertEqual(plan['dense_interval'], [189, 1024])
        self.assertIn('3/1250', plan['tilts']['q1'])  # .0024 witness retained.
        for kwargs in (dict(precision=192), dict(exact_occupancies=[3], dense_min=5)):
            with self.assertRaises(ValueError):
                proof.recipe(4096, 18, **kwargs)

    def test_prepare_allows_new_imports_not_changed_or_deleted_old_pins(self):
        before = {'old': 'digest'}
        for after, accepted in ((dict(before, new='added'), True), ({'old': 'changed'}, False), ({}, False)):
            with (patch.object(proof.q1, 'source_snapshot', side_effect=[before, after]),
                  patch.object(proof.inner, 'prepare', return_value=(self.data, self.map_record)) as prepare,
                  patch.object(proof.inner, 'authenticate') as authenticate):
                if accepted:
                    data, record, sources = proof._prepare(17, 256)
                    self.assertIs(data, self.data)
                    self.assertIs(record, self.map_record)
                    self.assertEqual(sources, after)
                    prepare.assert_called_once_with(17, birth_density='capped')
                    authenticate.assert_called_once_with(self.data, self.map_record)
                else:
                    with self.assertRaises(RuntimeError):
                        proof._prepare(17, 256)
        with self.assertRaises(ValueError):
            proof._prepare(17, 256, self.data, None)
        with patch.object(proof.inner, 'authenticate'), self.assertRaises(ValueError):
            proof._prepare(18, 256, self.data, self.map_record)

    def _tail_context(self):
        stack = ExitStack()
        stack.enter_context(patch.object(proof, '_prepare', return_value=(self.data, self.map_record, self.sources)))
        stack.enter_context(patch.object(proof.inner, 'authenticate'))
        stack.enter_context(patch.object(proof.q1, 'source_snapshot', return_value=self.sources))
        stack.enter_context(redirect_stdout(StringIO()))
        return stack

    def test_exact_tail_uses_actual_geometry_and_all_terminal_mass(self):
        regional = [arb_mat([[1, 0], [0, 1]]) for _ in range(6)]
        regional[3] = arb_mat([[1, 2], [0, 1]])
        with self._tail_context() as stack:
            local = stack.enter_context(patch.object(proof.q1.kernel_t64, 'local_operators', return_value=['toy']))
            placement = stack.enter_context(patch.object(proof.q1, 'placement', return_value=regional))
            result = proof.run_tail(K=4096, bits=17, occupancies=[3, 5], tilts=['.01'], precision=256)
        local.assert_called_once_with(self.data, Q(1, 100), activity=Q(1, 2))
        self.assertEqual(placement.call_args.kwargs['epochs'], 1)
        self.assertEqual(placement.call_args.kwargs['maximum_groups'], 5)
        self.assertEqual(result['occupancy_covered'], [3, 5])
        self.assertEqual(set(result['occupancy_uppers']), {'3', '5'})
        for q in (3, 5):
            expected = proof.lengths.occupancy_upper(regional, K=4096, occupancy=q,
                                                    beta=self.beta, tilt=Q(1, 100))
            self.assertEqual(result['occupancy_uppers'][str(q)], proof.q1.endpoint(expected))
        self.assertEqual(result['return_denominator'], 131071)
        self.assertEqual(result['inner_group'], 'GL(17,2)')
        self.assertEqual(result['outer_symbol_group'], 'GL(16,2)')
        self.assertFalse(result['whole_code_certificate'])

    def test_fugacity_tail_uses_physical_power_and_marker_endpoint(self):
        mixed = arb_mat([[1, 1], [0, 1]])
        with self._tail_context() as stack:
            candidates = stack.enter_context(patch.object(proof.dense, 'physical_candidates', return_value={'toy': 1}))
            select = stack.enter_context(patch.object(proof.dense, 'select_physical', return_value=(mixed, (0,)*17)))
            stack.enter_context(patch.object(proof.q1.kernel_t64, 'local_operators', side_effect=AssertionError('not macro path')))
            coefficient = stack.enter_context(patch.object(proof, '_fugacity_affine', wraps=proof._fugacity_affine))
            result = proof.run_tail(K=4096, bits=17, occupancies=[3, 32], tilts=['.01'],
                method='fugacity', marker_probabilities=['1/2', '1'], precision=256)
        candidates.assert_called_once_with(self.data, Q(1, 100))
        self.assertEqual([call.args[1] for call in select.call_args_list], [Q(15, 32), Q(15, 16)])
        self.assertEqual([call.args[2] for call in coefficient.call_args_list], [Q(1, 2), Q(1)])
        # [[1,1],[0,1]]^128 has total first-row mass129, not the macro power65.
        expected_moment = arb(129).log()
        for call in coefficient.call_args_list:
            self.assertTrue(call.args[3].overlaps(expected_moment))
        for q in (3, 32):
            markers = (Q(1, 2),) if q == 3 else (Q(1, 2), Q(1))
            candidates = [proof.q1.endpoint(proof.q1.kernel_t64.up(
                proof.dense.coefficient_log_bound(finite_geometry=self.geom, q=q, tilt=Q(1, 100),
                    beta=self.beta, marker_probability=p, log_moment=expected_moment).exp())) for p in markers]
            self.assertEqual(value(result['occupancy_uppers'][str(q)]), min(map(value, candidates)))
        self.assertEqual(result['count_sha256'], self.count_hash)
        self.assertFalse(result['whole_code_certificate'])

    def test_cached_fugacity_formula_matches_generic_helper(self):
        for K in (4096, 98304, 131072):
            geometry = proof.lengths.geometry(K)
            log_moment, log_beta = proof.q1.kernel_t64.aq(Q(3, 7)).log(), proof.q1.kernel_t64.aq(self.beta).log()
            for p in (Q(1, 17), Q(1, 2), Q(1)):
                base, slope = proof._fugacity_affine(geometry, Q(3, 100), p, log_moment, log_beta)
                for q in ((geometry.group_count,) if p == 1 else (3, geometry.group_count//2, geometry.group_count)):
                    generic = proof.dense.coefficient_log_bound(finite_geometry=geometry, q=q,
                        tilt=Q(3, 100), beta=self.beta, marker_probability=p, log_moment=log_moment)
                    cached = base if p == 1 else base + (1-geometry.regions)*arb(proof.comb(geometry.group_count, q)).log() + q*slope
                    self.assertTrue(cached.overlaps(generic))
                    self.assertEqual(proof.q1.endpoint(cached.exp()), proof.q1.endpoint(generic.exp()))

    def test_actual_generic_low_producer_records_pass_scope_checks_with_toy_moments(self):
        # Exercise the real producer's metadata and shell folding. Only the
        # placement moments and map authentication are mocked; no census runs.
        tiny = arb(2)**-160
        for q in (1, 2):
            with (patch.object(proof.low, 'validated', return_value=17),
                  patch.object(proof.low, 'source_snapshot', return_value=self.sources),
                  patch.object(proof.q1.kernel_t64, 'local_operators', return_value=['toy']),
                  patch.object(proof.q1, 'placement', return_value=[arb_mat([[1]])]*3),
                  patch.object(proof.q1, 'support_moments', return_value=[tiny]*65),
                  patch.object(proof.low.q2, 'pair_support_moments', return_value=[[tiny]*(v+1) for v in range(65)]),
                  redirect_stdout(StringIO())):
                record = proof.low.run_low(q, data=self.data, map_record=self.map_record,
                    geometry=self.geom, tilts=['.0024'], precision=256)
            self.assertEqual(set(self.endpoints(record)), {q})
            self.assertTrue(record['all_supports_covered'])
            self.assertEqual(record['count_sha256'], self.count_hash)

    def test_component_validation_rejects_missing_supports_and_wrong_identity(self):
        for q in (1, 2):
            self.assertEqual(self.endpoints(self.component(q)), {q: (1, -64)})
        self.assertEqual(set(self.endpoints(self.component(qs=[3, 5]))), {3, 5})
        cases = [('state_bits', 16), ('inner_group', 'GL(16,2)'), ('outer_symbol_group', 'GL(17,2)'),
            ('physical_updates_independent', False), ('zero_initial_state', 1), ('final_flush', True),
            ('threshold', 820), ('count_sha256', 'wrong'), ('all_supports_covered', False),
            ('whole_code_certificate', True), ('fresh_computation', False),
            ('source_sha256', {}), ('map_record', {}), ('all_two_group_support_pairs_covered', False)]
        for key, bad in cases:
            record = self.component(2)
            record[key] = bad
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.endpoints(record)
        for key, bad in [('return_denominator', 65535), ('every_shell_checked', False),
                         ('beta', '1'), ('occupancy_covered', [3, 4]), ('method', 'saved')]:
            record = self.component(qs=[3, 5])
            record[key] = bad
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.endpoints(record)

    def _replay(self, output, *, passing=True, gap=False, changed_source=False):
        def low_run(q, **kwargs):
            return self.component(q, endpoint=(1, -64) if passing or q == 2 else (1, -40))
        def tail_run(**kwargs):
            qs = list(kwargs['occupancies'])
            if gap and kwargs['method'] == 'fugacity':
                qs.remove(32)
            return self.component(qs=qs, method=kwargs['method'],
                                  endpoint=(1, -65) if kwargs['method'] == 'exact' else (1, -64))
        with ExitStack() as stack:
            prepared = stack.enter_context(patch.object(proof, '_prepare', return_value=(self.data, self.map_record, self.sources)))
            stack.enter_context(patch.object(proof.inner, 'authenticate'))
            stack.enter_context(patch.object(proof.q1, 'source_snapshot', return_value={} if changed_source else self.sources))
            low_run_mock = stack.enter_context(patch.object(proof.low, 'run_low', side_effect=low_run))
            tail_run_mock = stack.enter_context(patch.object(proof, 'run_tail', side_effect=tail_run))
            stack.enter_context(redirect_stdout(StringIO()))
            answer = proof.replay(output, K=4096, bits=17, exact_occupancies=[3, 5], dense_min=3)
        return answer, prepared, low_run_mock, tail_run_mock

    def test_replay_prepares_once_minimizes_overlaps_and_certifies_only_passes(self):
        for passing in (False, True):
            output = self.folder / f'whole-{passing}.json'
            result, prepared, lows, tails = self._replay(output, passing=passing)
            prepared.assert_called_once_with(17, 256)
            self.assertEqual(lows.call_count, 2)
            self.assertEqual(tails.call_count, 2)
            for call in [*lows.call_args_list, *tails.call_args_list]:
                self.assertIs(call.kwargs['data'], self.data)
                self.assertIs(call.kwargs['map_record'], self.map_record)
                self.assertEqual(call.kwargs['precision'], 256)
            self.assertEqual(result['occupancy_covered'], [1, 32])
            self.assertEqual(result['selected_components']['3'], 2)
            self.assertTrue(result['fresh_replay'])
            self.assertIs(result['target_met'], passing)
            self.assertIs(result['whole_code_certificate'], passing)
            self.assertEqual(json.loads(output.read_text()), result)
            if passing:
                self.assertEqual(value(result['union_upper']), Q(31, 2**64))

    def test_replay_rejects_holes_existing_paths_and_changed_sources(self):
        for kwargs, error in ((dict(gap=True), ValueError), (dict(changed_source=True), RuntimeError)):
            output = self.folder / 'must-not-exist.json'
            with self.subTest(kwargs=kwargs), self.assertRaises(error):
                self._replay(output, **kwargs)
            self.assertFalse(output.exists())
        output = self.folder / 'existing.json'
        for occupied in (output, self.folder/'existing-q1.json', self.folder/'existing-dense.json'):
            occupied.touch()
            with patch.object(proof, '_prepare') as prepare, self.assertRaises(ValueError):
                proof.replay(output, K=4096, bits=17)
            prepare.assert_not_called()
            occupied.unlink()


if __name__ == '__main__':
    unittest.main()
