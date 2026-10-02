import copy
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from flint import arb, ctx
import point_reuse as reuse


def witness():
    return dict(tilt='3/16', parameters=['1/20', '2', '3'], variance_dual=['0', '0'])


def point_record():
    return dict(schema=reuse.POINT_SCHEMA, updates=2, minimum_groups=33,
        maximum_groups=2048, last_lp=104, refined=True, base_tilt='3/16',
        variance_shuffle=True, variance_bins=16, regional_count=True, distance='.095',
        mixture='deliberately not consumed',
        probes=[dict(mean='1/8', witness=witness(), proposal='false', upper='false')])


def regional_point_record():
    record = point_record()
    part = dict(interval=['0', '1/8'], dual=['0', '0', '0'])
    candidate = record['probes'][0]['witness']
    candidate['variance_partition'] = [part]
    candidate['regional_count_parts'] = [dict(copy.deepcopy(part), mgf_witnesses=[])]
    return record


class FakeModel:
    q_min = 33
    tilt = Q(3, 16)
    variance_shuffle = True
    variance_bins = 16
    regional_count = True
    inner = SimpleNamespace(__name__='birth_classes')
    data = dict(windows=32, updates=2)
    root = (Q(1, 1000), Q(1))

    def __init__(self, upper=None):
        self.upper = arb(2)**-60 if upper is None else upper
        self.outward_calls = []
        self.proposal_calls = []

    def outward(self, cell, candidate):
        self.outward_calls.append((cell, copy.deepcopy(candidate), ctx.prec))
        return arb(self.upper)

    def proposal(self, cell):
        self.proposal_calls.append(cell)
        return 17., dict(fallback=True)


class PointReuseTests(unittest.TestCase):
    cell = (Q(1, 8), Q(3, 16))

    def setUp(self):
        self.precision = ctx.prec
        ctx.prec = 128

    def tearDown(self):
        ctx.prec = self.precision

    def wrap(self, underlying=None, probes=None):
        underlying = FakeModel() if underlying is None else underlying
        probes = reuse.point_proposals(point_record()) if probes is None else probes
        return reuse.PointReuseModel(underlying, probes, 40,
            regional_proposal=lambda model, cell, candidate: (-60., candidate)), underlying

    def test_bounds_and_mixture_are_not_transferred(self):
        record = point_record()
        expected = reuse.point_proposals(record)
        record.update(distance='.1', mixture=[dict(mass='bad', activity='bad')])
        record['probes'][0].update(proposal=-1e100, upper=[1, -1000000], log2_upper='lie')
        self.assertEqual(reuse.point_proposals(record), expected)
        self.assertEqual(set(expected[0]), {'mean', 'witness'})
        record['probes'][0]['witness']['parameters'][0] = '99'
        self.assertEqual(expected[0]['witness']['parameters'][0], '1/20')

    def test_wrong_source_scope_and_unknown_witness_fields_rejected(self):
        for key, value in [('schema', 'different'), ('updates', 4), ('updates', True),
                           ('minimum_groups', 1), ('regional_count', False),
                           ('base_tilt', '1/4'), ('distance', '.5')]:
            record = point_record()
            record[key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                reuse.point_proposals(record)
        record = point_record()
        record['probes'][0]['witness']['cached_upper'] = '1/100'
        with self.assertRaises(ValueError):
            reuse.point_proposals(record)

    def test_multiple_source_files_hashes_and_distances(self):
        with tempfile.TemporaryDirectory() as directory:
            paths, raws = [], []
            for index, distance in enumerate(('.095', '.1')):
                record = point_record()
                record['distance'] = distance
                raw = json.dumps(record).encode()
                path = Path(directory) / f'points{index}.json'
                path.write_bytes(raw)
                paths.append(path)
                raws.append(raw)
            probes, sources = reuse.load_points(paths)
            self.assertEqual(len(probes), 2)
            self.assertEqual([s['sha256'] for s in sources],
                             [hashlib.sha256(raw).hexdigest() for raw in raws])
            self.assertEqual([s['point_count'] for s in sources], [1, 1])

    def test_success_uses_fresh_bound_and_one_use_cache(self):
        model, underlying = self.wrap()
        score, candidate = model.proposal(self.cell)
        self.assertLess(score, -42)
        self.assertEqual(underlying.proposal_stop_bits, 42)
        self.assertEqual(len(underlying.outward_calls), 1)
        self.assertTrue(model.outward(self.cell, candidate) < arb(2)**-40)
        self.assertEqual(len(underlying.outward_calls), 1)
        model.outward(self.cell, candidate)
        self.assertEqual(len(underlying.outward_calls), 2)
        self.assertEqual(model.stats['cache_hits'], 1)
        self.assertEqual(model.root, underlying.root)

    def test_insufficient_or_borderline_bound_falls_back(self):
        for exponent in (-30, -40, -42):
            model, underlying = self.wrap(FakeModel(arb(2)**exponent))
            self.assertEqual(model.proposal(self.cell), (17., dict(fallback=True)))
            self.assertEqual(underlying.proposal_calls, [self.cell])
            self.assertEqual(model.stats['accepted'], 0)

    def test_boundary_cells_fall_back_without_rebasing(self):
        for cell in ((Q(0), Q(1, 4)), (Q(3, 4), Q(1))):
            model, underlying = self.wrap()
            self.assertEqual(model.proposal(cell), (17., dict(fallback=True)))
            self.assertEqual(underlying.outward_calls, [])

    def test_cell_witness_and_precision_cache_isolation(self):
        for changed in ('cell', 'witness', 'precision'):
            model, underlying = self.wrap()
            _, candidate = model.proposal(self.cell)
            cell = self.cell
            if changed == 'cell':
                cell = (Q(1, 8), Q(1, 4))
            if changed == 'witness':
                candidate['parameters'][1] = '100'
            if changed == 'precision':
                ctx.prec = 256
            model.outward(cell, candidate)
            self.assertEqual(len(underlying.outward_calls), 2, changed)
            self.assertEqual(model.stats['cache_hits'], 0, changed)
            self.assertEqual(model.stats['rolling_updates'], 0, changed)
            ctx.prec = 128

    def test_new_proposal_invalidates_prior_cache_and_copies_hints(self):
        probes = reuse.point_proposals(point_record())
        model, underlying = self.wrap(probes=probes)
        probes[0]['witness']['parameters'][0] = '999'
        _, first = model.proposal(self.cell)
        self.assertEqual(first['parameters'][0], '1/20')
        model.proposal((Q(1, 4), Q(3, 8)))
        model.outward(self.cell, first)
        self.assertEqual(len(underlying.outward_calls), 3)
        self.assertEqual(model.stats['cache_hits'], 0)

    def test_different_model_never_inherits_cached_bound(self):
        first, _ = self.wrap()
        _, candidate = first.proposal(self.cell)
        second, underlying = self.wrap(FakeModel(arb(2)))
        self.assertTrue(second.outward(self.cell, candidate) > 1)
        self.assertEqual(len(underlying.outward_calls), 1)

    def test_invalid_bounds_and_outward_errors_are_not_hidden(self):
        for upper in (arb(0), arb(-1), arb('nan'), arb('inf')):
            model, underlying = self.wrap(FakeModel(upper))
            with self.assertRaises(ArithmeticError):
                model.proposal(self.cell)
            self.assertEqual(underlying.proposal_calls, [])
        model, underlying = self.wrap()
        def fail(cell, candidate):
            raise ValueError('invalid witness')
        underlying.outward = fail
        model._evaluate = underlying.outward
        with self.assertRaisesRegex(ValueError, 'invalid witness'):
            model.proposal(self.cell)

    def test_outward_precision_changes_are_rejected(self):
        model, underlying = self.wrap()
        def changed_precision(cell, candidate):
            ctx.prec = 256
            return arb(2)**-60
        underlying.outward = changed_precision
        model._evaluate = underlying.outward
        with self.assertRaisesRegex(ArithmeticError, 'precision'):
            model.proposal(self.cell)

    def test_variance_and_regional_partitions_rebased_together(self):
        record = point_record()
        parts = [dict(interval=['0', '1/16'], dual=['2', '3', '4']),
                 dict(interval=['1/16', '1/8'], dual=['-2', '0', '-4'])]
        candidate = record['probes'][0]['witness']
        candidate['variance_partition'] = parts
        candidate['regional_count_parts'] = [dict(copy.deepcopy(part),
            mgf_witnesses=[dict(tilt='-1', dual=['1', '2', '-3'])]) for part in parts]
        model, underlying = self.wrap(probes=reuse.point_proposals(record))
        _, rebased = model.proposal((Q(1, 4), Q(1, 3)))
        self.assertEqual(rebased['variance_partition'][-1]['interval'][1], '1/4')
        for plain, regional in zip(rebased['variance_partition'], rebased['regional_count_parts']):
            self.assertEqual(plain['interval'], regional['interval'])
            self.assertEqual(plain['dual'], regional['dual'])
        self.assertEqual(candidate['variance_partition'][-1]['interval'][1], '1/8')
        self.assertEqual(len(underlying.outward_calls), 1)

    def test_regional_float_screen_never_accepts_without_fresh_outward(self):
        probes = reuse.point_proposals(point_record())
        # The injected screener below only tests the adapter's control flow;
        # the pure-helper tests cover actual complete-partition validation.
        probes[0]['witness']['variance_partition'] = [
            dict(interval=['0', '1/8'], dual=['0', '0', '0'])]
        probes[0]['witness']['regional_count_parts'] = [
            dict(interval=['0', '1/8'], dual=['0', '0', '0'], mgf_witnesses=[])]
        for screen, fresh, expected_calls in ((0., -60, 0), (-60., -30, 1), (-60., -60, 1)):
            underlying = FakeModel(arb(2)**fresh)
            model = reuse.PointReuseModel(underlying, probes, 40,
                regional_proposal=lambda m, c, w: (screen, w))
            score, candidate = model.proposal(self.cell)
            self.assertEqual(len(underlying.outward_calls), expected_calls)
            self.assertEqual(model.stats['accepted'], int(screen == -60. and fresh == -60))
            if model.stats['accepted']:
                self.assertLess(score, -42)
                self.assertTrue(model.outward(self.cell, candidate) < arb(2)**-40)
                self.assertEqual(len(underlying.outward_calls), 1)
            else:
                self.assertEqual(underlying.proposal_calls, [self.cell])

    def test_injected_outward_evaluator_is_used_and_checked(self):
        underlying = FakeModel(arb(2))
        calls = []
        def fresh(cell, candidate):
            calls.append((cell, copy.deepcopy(candidate), ctx.prec))
            return arb(2)**-60
        model = reuse.PointReuseModel(underlying, reuse.point_proposals(point_record()),
            40, outward=fresh)
        _, candidate = model.proposal(self.cell)
        model.outward(self.cell, candidate)
        self.assertEqual(len(calls), 1)
        self.assertEqual(underlying.outward_calls, [])

    def test_successful_fallback_supplies_next_rolling_hint(self):
        underlying = FakeModel()
        optimized = witness()
        optimized['parameters'][0] = '1/10'
        def optimize(cell):
            underlying.proposal_calls.append(cell)
            return -60., copy.deepcopy(optimized)
        def evaluate(cell, candidate):
            underlying.outward_calls.append((cell, copy.deepcopy(candidate), ctx.prec))
            return arb(2)**(-60 if candidate['parameters'][0] == '1/10' else -30)
        underlying.proposal = optimize
        model = reuse.PointReuseModel(underlying, reuse.point_proposals(point_record()),
            40, outward=evaluate)
        _, candidate = model.proposal(self.cell)
        self.assertIsNone(model._rolling)
        model.outward(self.cell, candidate)
        self.assertEqual(model.stats['rolling_updates'], 1)
        self.assertEqual(model._rolling[0], self.cell)
        self.assertEqual(len(model._rolling), 2)  # No retained numerical bound.
        candidate['parameters'][0] = '999'
        neighbor = (self.cell[1], Q(1, 4))
        _, next_candidate = model.proposal(neighbor)
        self.assertEqual(next_candidate['parameters'][0], '1/10')
        self.assertEqual(underlying.proposal_calls, [self.cell])
        self.assertEqual(model.stats['rolling_attempts'], 1)
        self.assertEqual(model.stats['rolling_accepted'], 1)
        self.assertEqual(len(underlying.outward_calls), 3)
        model.outward(neighbor, next_candidate)
        self.assertEqual(len(underlying.outward_calls), 3)
        self.assertEqual(model.stats['rolling_updates'], 2)

    def test_rolling_bound_is_fresh_and_duplicate_static_hint_skipped(self):
        model, underlying = self.wrap()
        _, candidate = model.proposal(self.cell)
        model.outward(self.cell, candidate)
        underlying.upper = arb(2)**-30
        calls = len(underlying.outward_calls)
        result = model.proposal((self.cell[1], Q(1, 4)))
        self.assertEqual(result, (17., dict(fallback=True)))
        # The rolling and nearest-point candidates have identical duals.
        self.assertEqual(len(underlying.outward_calls) - calls, 1)
        self.assertEqual(model.stats['duplicate_hints'], 1)
        self.assertEqual(model.stats['rolling_accepted'], 0)

    def test_two_distinct_hint_attempts_are_the_limit(self):
        underlying = FakeModel()
        optimized = witness()
        optimized['parameters'][0] = '1/10'
        def optimize(cell):
            underlying.proposal_calls.append(cell)
            return -60., copy.deepcopy(optimized)
        underlying.proposal = optimize
        def evaluate(cell, candidate):
            underlying.outward_calls.append((cell, copy.deepcopy(candidate), ctx.prec))
            return arb(2)**(-60 if candidate['parameters'][0] == '1/10' else -30)
        model = reuse.PointReuseModel(underlying, reuse.point_proposals(point_record()),
            40, outward=evaluate)
        _, candidate = model.proposal(self.cell)
        model.outward(self.cell, candidate)
        model._evaluate = underlying.outward
        underlying.upper = arb(2)**-30
        calls = len(underlying.outward_calls)
        attempts = model.stats['attempts']
        model.proposal((self.cell[1], Q(1, 4)))
        self.assertEqual(len(underlying.outward_calls) - calls, 2)
        self.assertEqual(model.stats['attempts'] - attempts, 2)
        self.assertEqual(len(underlying.proposal_calls), 2)

    def test_replay_outward_never_seeds_or_replaces_rolling_hint(self):
        model, underlying = self.wrap()
        other = witness()
        other['parameters'][0] = '7'
        model.outward(self.cell, other)
        self.assertIsNone(model._rolling)
        _, candidate = model.proposal(self.cell)
        model.outward(self.cell, candidate)
        retained = copy.deepcopy(model._rolling)
        model.outward((Q(1, 3), Q(1, 2)), other)
        self.assertEqual(model._rolling, retained)
        self.assertEqual(model.stats['rolling_updates'], 1)

    def test_intervening_replay_discards_issued_proposal_token(self):
        model, _ = self.wrap()
        _, candidate = model.proposal(self.cell)
        model.outward((Q(1, 3), Q(1, 2)), candidate)
        model.outward(self.cell, candidate)
        self.assertIsNone(model._rolling)
        self.assertEqual(model.stats['rolling_updates'], 0)

    def test_failed_fallback_outward_does_not_seed_rolling(self):
        underlying = FakeModel(arb(2)**-30)
        underlying.proposal = lambda cell: (-60., witness())
        model, _ = self.wrap(underlying)
        _, candidate = model.proposal(self.cell)
        self.assertTrue(model.outward(self.cell, candidate) > arb(2)**-40)
        self.assertIsNone(model._rolling)

    def test_rolling_can_be_disabled(self):
        underlying = FakeModel()
        model = reuse.PointReuseModel(underlying, reuse.point_proposals(point_record()),
            40, rolling=False)
        _, candidate = model.proposal(self.cell)
        model.outward(self.cell, candidate)
        self.assertIsNone(model._rolling)
        model.proposal((self.cell[1], Q(1, 4)))
        self.assertEqual(model.stats['rolling_attempts'], 0)
        with self.assertRaises(ValueError):
            reuse.PointReuseModel(underlying, reuse.point_proposals(point_record()), 40, rolling=1)

    def test_rolling_regional_partition_rebases_from_cell_not_point(self):
        record = point_record()
        parts = [dict(interval=['0', '1/8'], dual=['0', '0', '0'])]
        candidate = record['probes'][0]['witness']
        candidate['variance_partition'] = parts
        candidate['regional_count_parts'] = [dict(copy.deepcopy(parts[0]), mgf_witnesses=[])]
        model, _ = self.wrap(probes=reuse.point_proposals(record))
        _, first = model.proposal(self.cell)
        self.assertEqual(first['variance_partition'][0]['interval'], ['0', '3/16'])
        model.outward(self.cell, first)
        _, second = model.proposal((Q(1, 4), Q(1, 3)))
        self.assertEqual(second['variance_partition'][0]['interval'], ['0', '1/4'])
        self.assertEqual(second['regional_count_parts'][0]['interval'], ['0', '1/4'])
        self.assertEqual(model.stats['rolling_accepted'], 1)

    def test_warm_rolling_skips_only_screen_and_still_checks_fresh_bound(self):
        underlying = FakeModel()
        screens, checks = [], []
        def screen(m, cell, candidate):
            screens.append(cell)
            return -60., candidate
        def warm(candidate):
            checks.append(copy.deepcopy(candidate))
            return True
        model = reuse.PointReuseModel(underlying,
            reuse.point_proposals(regional_point_record()), 40,
            regional_proposal=screen, regional_warm=warm)
        _, first = model.proposal(self.cell)
        model.outward(self.cell, first)
        self.assertEqual(checks, [])  # Static points always keep their screen.
        neighbor = (self.cell[1], Q(1, 4))
        _, second = model.proposal(neighbor)
        self.assertEqual(screens, [self.cell])
        self.assertEqual(len(checks), 1)
        self.assertEqual(len(underlying.outward_calls), 2)
        self.assertEqual(underlying.outward_calls[-1][0], neighbor)
        self.assertEqual(model.stats['warm_exact_attempts'], 1)
        self.assertEqual(model.stats['warm_exact_accepted'], 1)
        model.outward(neighbor, second)
        self.assertEqual(len(underlying.outward_calls), 2)

    def test_cold_rolling_keeps_float_screen(self):
        underlying = FakeModel()
        screens = []
        def screen(m, cell, candidate):
            screens.append(cell)
            return -60., candidate
        model = reuse.PointReuseModel(underlying,
            reuse.point_proposals(regional_point_record()), 40,
            regional_proposal=screen, regional_warm=lambda candidate: False)
        _, first = model.proposal(self.cell)
        model.outward(self.cell, first)
        neighbor = (self.cell[1], Q(1, 4))
        model.proposal(neighbor)
        self.assertEqual(screens, [self.cell, neighbor])
        self.assertEqual(model.stats['warm_exact_attempts'], 0)
        self.assertEqual(len(underlying.outward_calls), 2)

    def test_poor_warm_exact_hint_falls_through_distinct_static_hint(self):
        underlying = FakeModel()
        probes = reuse.point_proposals(regional_point_record())
        optimized = reuse.rebase_witness((Q(1, 8), Q(1, 8)), self.cell, probes[0]['witness'])
        optimized['parameters'][0] = '1/10'
        phase, screens, checks = [0], [], []
        def optimize(cell):
            underlying.proposal_calls.append(cell)
            return -60., copy.deepcopy(optimized)
        def evaluate(cell, candidate):
            underlying.outward_calls.append((cell, copy.deepcopy(candidate), ctx.prec))
            good = candidate['parameters'][0] == ('1/10' if phase[0] == 0 else '1/20')
            return arb(2)**(-60 if good else -30)
        def screen(m, cell, candidate):
            screens.append(candidate['parameters'][0])
            return -60., candidate
        def warm(candidate):
            checks.append(candidate['parameters'][0])
            return True
        underlying.proposal = optimize
        model = reuse.PointReuseModel(underlying, probes, 40,
            outward=evaluate, regional_proposal=screen, regional_warm=warm)
        _, first = model.proposal(self.cell)
        model.outward(self.cell, first)
        phase[0] = 1
        _, second = model.proposal((self.cell[1], Q(1, 4)))
        self.assertEqual(second['parameters'][0], '1/20')
        self.assertEqual(checks, ['1/10'])
        self.assertEqual(screens, ['1/20', '1/20'])
        self.assertEqual(len(underlying.outward_calls), 4)
        self.assertEqual(underlying.proposal_calls, [self.cell])
        self.assertEqual(model.stats['warm_exact_attempts'], 1)
        self.assertEqual(model.stats['warm_exact_rejected'], 1)
        self.assertEqual(model.stats['warm_exact_accepted'], 0)

    def test_warmth_predicate_changes_or_nonboolean_results_fail(self):
        for kind in ('precision', 'witness', 'nonboolean', 'exception'):
            model, underlying = self.wrap(probes=reuse.point_proposals(regional_point_record()))
            _, first = model.proposal(self.cell)
            model.outward(self.cell, first)
            def bad(candidate):
                if kind == 'precision':
                    ctx.prec += 1
                if kind == 'witness':
                    candidate['parameters'][0] = '1/9'
                if kind == 'exception':
                    raise ValueError('invalid selector')
                return 1 if kind == 'nonboolean' else True
            model._regional_warm = bad
            with self.subTest(kind=kind), self.assertRaises((ArithmeticError, ValueError)):
                model.proposal((self.cell[1], Q(1, 4)))
            self.assertEqual(len(underlying.outward_calls), 1)
            ctx.prec = 128

    def test_nonregional_rolling_never_asks_warmth(self):
        model, _ = self.wrap()
        def unexpected(candidate):
            raise AssertionError('nonregional warmth queried')
        model._regional_warm = unexpected
        _, first = model.proposal(self.cell)
        model.outward(self.cell, first)
        model.proposal((self.cell[1], Q(1, 4)))
        self.assertEqual(model.stats['warm_exact_attempts'], 0)

    def test_failed_wide_hints_request_bisection_without_optimization(self):
        underlying = FakeModel(arb(2)**-30)
        model = reuse.PointReuseModel(underlying, reuse.point_proposals(point_record()),
            40, maximum_optimization_width=Q(1, 16))
        wide = (Q(1, 8), Q(1, 4))
        self.assertEqual(model.proposal(wide), (0., {}))
        self.assertEqual(underlying.proposal_calls, [])
        self.assertEqual(len(underlying.outward_calls), 1)
        self.assertEqual(model.stats['wide_splits'], 1)
        self.assertEqual(model.stats['fallback'], 0)
        self.assertIsNone(model._pending)
        self.assertIsNone(model._awaiting)

    def test_width_boundary_and_narrower_cells_keep_fallback(self):
        underlying = FakeModel(arb(2)**-30)
        model = reuse.PointReuseModel(underlying, reuse.point_proposals(point_record()),
            40, maximum_optimization_width=Q(1, 16))
        for cell in (self.cell, (Q(1, 8), Q(5, 32))):
            self.assertEqual(model.proposal(cell), (17., dict(fallback=True)))
        self.assertEqual(len(underlying.proposal_calls), 2)
        self.assertEqual(model.stats['fallback'], 2)
        self.assertEqual(model.stats['wide_splits'], 0)

    def test_fresh_wide_hint_can_pass_without_subdivision(self):
        underlying = FakeModel()
        model = reuse.PointReuseModel(underlying, reuse.point_proposals(point_record()),
            40, maximum_optimization_width=Q(1, 16))
        wide = (Q(1, 8), Q(1, 4))
        score, candidate = model.proposal(wide)
        self.assertLess(score, -42)
        self.assertTrue(model.outward(wide, candidate) < arb(2)**-40)
        self.assertEqual(len(underlying.outward_calls), 1)
        self.assertEqual(model.stats['accepted'], 1)
        self.assertEqual(model.stats['wide_splits'], 0)
        self.assertEqual(underlying.proposal_calls, [])

    def test_width_cap_does_not_hide_invalid_math(self):
        underlying = FakeModel(arb('nan'))
        model = reuse.PointReuseModel(underlying, reuse.point_proposals(point_record()),
            40, maximum_optimization_width=Q(1, 16))
        with self.assertRaises(ArithmeticError):
            model.proposal((Q(1, 8), Q(1, 4)))
        self.assertEqual(model.stats['wide_splits'], 0)
        self.assertEqual(underlying.proposal_calls, [])

    def test_optimization_width_validation_and_default(self):
        for cap in (0, -1, 2, True):
            with self.subTest(cap=cap), self.assertRaises(ValueError):
                reuse.PointReuseModel(FakeModel(), reuse.point_proposals(point_record()),
                    40, maximum_optimization_width=cap)
        model, underlying = self.wrap(FakeModel(arb(2)**-30))
        self.assertEqual(model.proposal((Q(1, 8), Q(3, 4))), (17., dict(fallback=True)))
        self.assertEqual(len(underlying.proposal_calls), 1)
        self.assertEqual(model.stats['wide_splits'], 0)


if __name__ == '__main__':
    unittest.main()
