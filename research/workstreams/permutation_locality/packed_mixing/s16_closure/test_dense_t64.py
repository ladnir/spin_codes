"""Dense-wrapper regressions using tiny geometry and mocked outer preparation.

These tests do not run a production count, optimizer, or proof preparation.
The two-physical-step kernel itself is checked in test_kernel_t64.py.
"""
from contextlib import ExitStack, redirect_stdout
from copy import deepcopy
from fractions import Fraction as Q
import hashlib
import io
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from flint import arb, arb_mat, ctx
import dense_t64 as driver


def identity(path):
    return dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def scope_fixture():
    return dict(schema='packed-gl32-t64-s16-dense-context-1', root=['1/8', '1'],
        physical_t=64, state_bits=16, physical_steps=32768, macro_t=128,
        macro_windows=32, macro_steps=16384, physical_steps_per_macro=2,
        regions=256, macros_per_region=64, birth_density='capped',
        maps=dict(spectrum={'0': 1}), maps_sha256='fresh map',
        outer_premises=dict(authentication='fresh'), threshold=209715)


def record_fixture(scope=None):
    scope = deepcopy(scope if scope is not None else scope_fixture())
    leaves = {path: dict(cell=list(map(str, driver.geometry.path_cell(scope['root'], path))),
        witness=dict(parameters=['1/4'], selected=path), upper=['not', 'evidence'])
        for path in ('0', '1')}
    return dict(schema=driver.SCHEMA, scope=scope,
        cover=dict(leaves=leaves, unresolved={}, visited=2), aggregate_upper='ignored')


class DenseT64Tests(unittest.TestCase):
    def setUp(self):
        self.precision = ctx.prec
        ctx.prec = 256
        # Accidental fresh outer preparation is a test error, not a long job.
        guard = patch.object(driver.prior, 'authenticated_bch_cdf',
            side_effect=AssertionError('production outer preparation forbidden in unit tests'))
        guard.start()
        self.addCleanup(guard.stop)
        self.addCleanup(setattr, ctx, 'prec', self.precision)

    def test_private_adapter_keeps_legacy_globals_untouched(self):
        self.assertIs(driver._regional_propose.__globals__, driver._bindings)
        self.assertIs(driver._regional_outward.__globals__, driver._bindings)
        self.assertIs(driver._bindings['local_operators'], driver.local_operators)
        self.assertIsNot(driver.regional.propose.__globals__, driver._bindings)
        self.assertIsNot(driver.regional.local_operators, driver.local_operators)

    def test_private_outward_calls_macro_factory_with_tiny_geometry(self):
        local = [arb_mat([[1]])]*33
        model = SimpleNamespace(data=dict(map_sha256='tiny', birth_density='capped'),
            tilt=Q(3, 16), threshold=0, q_min=1, features=[Q(1, 4)], active=[1],
            weights=lambda *args: ([Q(1), Q(1)], None),
            family=lambda *args: ([Q(1)], None, None))
        witness = dict(parameters=['1/4'], regional_direct_counts=True, variance_dual=[])
        checked = dict(regional_count_parts=[dict(mgf_witnesses=[])])
        bindings = dict(sc=SimpleNamespace(G=1, REGIONS=1, PACKETS=1),
            prepare_witness=Mock(return_value=([((Q(0), Q(1, 4)), (Q(0), Q(0), Q(0)))], checked)),
            placement=Mock(return_value=[arb_mat([[1]])]*2),
            count_mass_caps=Mock(return_value=[arb(0), arb(3)/16]))
        with patch.dict(driver._bindings, bindings), \
             patch.object(driver.kernel_t64, 'authenticate') as authenticate, \
             patch.object(driver.kernel_t64, 'local_operators', return_value=local) as factory, \
             patch.object(driver.regional, 'local_operators', side_effect=AssertionError('legacy operator used')), \
             patch.object(driver.variance_partition, 'factor', return_value=Q(1)), redirect_stdout(io.StringIO()):
            upper, result = driver._regional_outward(model, (Q(1, 8), Q(1, 4)), witness)
        authenticate.assert_called_once_with(model.data)
        factory.assert_called_once_with(model.data, Q(1, 4), Q(1, 2))
        bindings['placement'].assert_called_once_with(local)
        self.assertIs(result, checked)
        self.assertTrue(upper >= 1 and upper < 1+arb(2)**-100)

    def test_local_cache_binds_precision_map_density_tilt_and_activity(self):
        model = SimpleNamespace(data=dict(map_sha256='first', birth_density='capped'))
        witness = dict(parameters=['1/4'], t64_birth_density_activity='1/3')
        scratch = {}
        with patch.object(driver.kernel_t64, 'authenticate') as authenticate, \
             patch.object(driver.kernel_t64, 'local_operators', side_effect=lambda *args: object()) as factory:
            first = driver.local_operators(model, witness, _proposal_scratch=scratch)
            self.assertIs(first, driver.local_operators(model, witness, _proposal_scratch=scratch))
            self.assertEqual(factory.call_count, 1)
            # Authentication is still required on cache hits.
            self.assertEqual(authenticate.call_count, 2)
            for mutate in (lambda: setattr(ctx, 'prec', 192),
                           lambda: model.data.update(map_sha256='second'),
                           lambda: model.data.update(birth_density='classes'),
                           lambda: witness.update(parameters=['1/5']),
                           lambda: witness.update(t64_birth_density_activity='2/3')):
                previous = factory.call_count
                mutate()
                driver.local_operators(model, witness, _proposal_scratch=scratch)
                self.assertEqual(factory.call_count, previous+1)
            # Outward calls without proposal scratch always rebuild.
            driver.local_operators(model, witness)
            driver.local_operators(model, witness)
            self.assertEqual(factory.call_count, 8)
            self.assertEqual(authenticate.call_count, 9)

    def test_invalid_local_parameters_fail_before_factory(self):
        model = SimpleNamespace(data=dict(map_sha256='tiny'))
        with patch.object(driver.kernel_t64, 'authenticate'), \
             patch.object(driver.kernel_t64, 'local_operators') as factory:
            for lam, activity in (('0', '1/2'), ('-1', '1/2'), ('1/4', '-1'), ('1/4', '2')):
                with self.subTest(lam=lam, activity=activity), self.assertRaises(ValueError):
                    driver.local_operators(model, dict(parameters=[lam], t64_birth_density_activity=activity))
            factory.assert_not_called()

    def fresh_case(self, *, changed=None, fail_mixture=False, mutate_source=False):
        caps, comparison = [0, 1], [0, 2]
        premises = dict(canonical_cdf=[0, 1], authenticated='new count')
        fingerprint = driver.prior.cell_search.dense.fingerprint
        old = dict(K=1 << 20, N=1 << 21, block_width=8, minimum_groups=33,
            maximum_groups=2048, comparison='direct-expected-shell-majorant',
            outer_premises=premises, expected_cdf_sha256=fingerprint(caps),
            comparison_caps_sha256=fingerprint(comparison), mixture=['rational proposal'])
        if changed:
            old[changed[0]] = changed[1]
        with TemporaryDirectory() as directory, ExitStack() as stack:
            source = Path(directory)/'proposal.json'
            source.write_text(json.dumps(dict(scope=old, upper='ignored numerical endpoint')))
            events = []
            def authenticate():
                events.append('outer'); return caps, premises
            def mixture(*args):
                events.append('mixture')
                if fail_mixture:
                    raise ValueError('unverified shell majorant')
                return [(Q(2), Q(1, 3))], dict(verified='every shell')
            data, maps = dict(new_kernel=True), dict(spectrum={0: 1}, physical_t=64)
            def prepare(**kwargs):
                events.append('maps')
                if mutate_source:
                    source.write_text(source.read_text()+' ')
                return data, maps
            stack.enter_context(patch.object(driver.prior, 'authenticated_bch_cdf', side_effect=authenticate))
            stack.enter_context(patch.object(driver.prior, 'full_block', return_value='GL32'))
            transport = stack.enter_context(patch.object(driver.prior, 'transport_shells', return_value=comparison))
            exact = stack.enter_context(patch.object(driver.prior.cell_search.dense, 'exact_mixture', side_effect=mixture))
            prepare_mock = stack.enter_context(patch.object(driver.kernel_t64, 'prepare', side_effect=prepare))
            fake = SimpleNamespace(root=(Q(1, 8), Q(1)), tilt=Q(3, 16))
            constructor = stack.enter_context(patch.object(driver, 'Model', return_value=fake))
            stack.enter_context(redirect_stdout(io.StringIO()))
            if changed or fail_mixture or mutate_source:
                with self.assertRaises(ValueError):
                    driver.fresh_model(source, precision=192, variance_bins=8)
                if not mutate_source:
                    prepare_mock.assert_not_called()
                return
            model, scope, metadata = driver.fresh_model(source, precision=192, variance_bins=8)
            self.assertIs(model, fake)
            self.assertEqual(events, ['outer', 'mixture', 'maps'])
            transport.assert_called_once_with(premises['canonical_cdf'], 'GL32')
            exact.assert_called_once_with(comparison, old['mixture'])
            prepare_mock.assert_called_once_with(birth_density='capped')
            constructor.assert_called_once_with([(Q(2), Q(1, 3))], data, 209715, 8)
            self.assertEqual(scope['maps']['spectrum'], {'0': 1})
            self.assertEqual(scope['maps_sha256'], driver.prior.fingerprint(scope['maps']))
            self.assertEqual(scope['physical_t']*scope['physical_steps'], scope['N'])
            self.assertEqual(scope['macro_t']*scope['macro_steps'], scope['N'])
            self.assertEqual(scope['regions']*scope['macros_per_region'], scope['macro_steps'])
            self.assertEqual(scope['physical_steps_per_macro']*scope['macro_steps'], scope['physical_steps'])
            self.assertEqual(scope['macro_windows'], 32)
            self.assertEqual(scope['state_bits'], 16)
            self.assertTrue(scope['outer_authenticated'])
            self.assertEqual(metadata, identity(source))
            self.assertEqual(ctx.prec, 192)

    def test_fresh_count_mixture_map_binding_and_macro_geometry(self):
        self.fresh_case()

    def test_stale_count_geometry_or_unverified_mixture_is_rejected(self):
        for changed in (('K', 1 << 19), ('N', 1 << 22), ('block_width', 4),
                        ('minimum_groups', 1), ('maximum_groups', 1024),
                        ('comparison', 'legacy bound'), ('outer_premises', {}),
                        ('expected_cdf_sha256', 'stale'), ('comparison_caps_sha256', 'stale')):
            with self.subTest(changed=changed):
                self.fresh_case(changed=changed)
        self.fresh_case(fail_mixture=True)
        self.fresh_case(mutate_source=True)

    def test_invalid_preparation_options_are_rejected_without_preparation(self):
        for option in (dict(precision=127), dict(precision=True), dict(variance_bins=0),
                       dict(variance_bins=65), dict(distance=0), dict(distance=Q(1, 2)),
                       dict(birth_density='old_s19')):
            with self.subTest(option=option), self.assertRaises(ValueError):
                driver.fresh_model('not read', **option)

    def test_exact_cover_validation_rejects_gaps_overlap_and_stale_cells(self):
        good = record_fixture()
        self.assertEqual(set(driver.validate_cover_record(good, complete=True)), {'0', '1'})
        changes = []
        gap = deepcopy(good); del gap['cover']['leaves']['1']; changes.append(gap)
        overlap = deepcopy(good); overlap['cover']['leaves'][''] = dict(witness={}); changes.append(overlap)
        stale = deepcopy(good); stale['cover']['leaves']['0']['cell'][1] = '1/2'; changes.append(stale)
        missing = deepcopy(good); del missing['cover']['leaves']['0']['witness']; changes.append(missing)
        pending = deepcopy(good); pending['cover']['unresolved']['1'] = pending['cover']['leaves'].pop('1'); changes.append(pending)
        schema = deepcopy(good); schema['schema'] = 'legacy s19'; changes.append(schema)
        for record in changes:
            with self.assertRaises(ValueError):
                driver.validate_cover_record(record, complete=True)

    def replay_case(self, *, change_scope=None, bounds=(-60, -61), mutation=None, bad_value=None):
        with TemporaryDirectory() as directory, ExitStack() as stack:
            outer, source, output = [Path(directory)/name for name in ('outer.json', 'search.json', 'fresh.json')]
            outer.write_text('mocked outer premise source')
            scope = scope_fixture(); old = record_fixture(scope)
            if change_scope:
                old['scope'][change_scope] = 'different scope'
            source.write_text(json.dumps(old))
            seen = []
            def outward(cell, witness):
                seen.append((cell, witness))
                if mutation == 'input':
                    source.write_text(source.read_text()+' ')
                elif mutation == 'outer':
                    outer.write_text(outer.read_text()+' ')
                elif mutation == 'precision':
                    ctx.prec = 160
                return bad_value if bad_value is not None else arb(2)**bounds[len(seen)-1]
            model = SimpleNamespace(root=tuple(map(Q, scope['root'])), outward=outward)
            stack.enter_context(patch.object(driver, 'fresh_model', return_value=(model, scope, identity(outer))))
            stack.enter_context(patch.object(sys, 'argv', ['dense_t64.py', 'replay', '--input', str(source),
                '--source', str(outer), '--precision', '192', '--output', str(output)]))
            stack.enter_context(redirect_stdout(io.StringIO()))
            if change_scope or mutation or bad_value is not None:
                with self.assertRaises((ValueError, ArithmeticError)):
                    driver.main()
                if change_scope:
                    self.assertEqual(seen, [])
                if output.exists():
                    result = json.loads(output.read_text())
                    self.assertFalse(result['complete_dense'])
                    self.assertNotIn('aggregate_below_2_minus_40', result)
                return
            driver.main()
            result = json.loads(output.read_text())
            self.assertEqual([item[0] for item in seen], [driver.geometry.path_cell(scope['root'], p) for p in ('0', '1')])
            self.assertEqual([item[1] for item in seen], [old['cover']['leaves'][p]['witness'] for p in ('0', '1')])
            expected = sum((arb(2)**power for power in bounds), arb(0))
            self.assertEqual(result['aggregate_upper'], driver.endpoint(expected))
            self.assertTrue(result['complete_dense'])
            self.assertEqual(result['aggregate_below_2_minus_40'], bool(expected < arb(2)**-40))
            self.assertFalse(result['whole_code_certificate'])
            self.assertEqual(result['input_source'], identity(source))
            self.assertEqual(len(result['cells']), 2)

    def test_replay_rebuilds_and_sums_bounds_not_saved_endpoints(self):
        self.replay_case()
        # A complete replay is not automatically a successful bound.
        self.replay_case(bounds=(-39, -39))

    def test_replay_requires_fresh_exact_scope_before_any_bound(self):
        for field in ('physical_t', 'physical_steps', 'macro_steps', 'state_bits',
                      'physical_steps_per_macro', 'birth_density', 'maps', 'outer_premises', 'threshold'):
            with self.subTest(field=field):
                self.replay_case(change_scope=field)

    def test_replay_fails_closed_on_source_precision_or_endpoint_changes(self):
        for mutation in ('input', 'outer', 'precision'):
            with self.subTest(mutation=mutation):
                self.replay_case(mutation=mutation)
        for value in (arb(0), arb(-1), arb('nan'), arb('inf')):
            with self.subTest(value=str(value)):
                self.replay_case(bad_value=value)

    def test_search_reuses_only_partition_geometry(self):
        with TemporaryDirectory() as directory, ExitStack() as stack:
            outer, source, output = [Path(directory)/name for name in ('outer.json', 'geometry.json', 'new.json')]
            outer.write_text('mocked outer premise source')
            scope = scope_fixture(); old = record_fixture(scope)
            for row in old['cover']['leaves'].values():
                row['witness'] = 'poisoned old witness'
            source.write_text(json.dumps(old))
            model = SimpleNamespace(root=tuple(map(Q, scope['root'])))
            def search(received_model, assigned, **options):
                self.assertIs(received_model, model)
                self.assertEqual(assigned, ['0', '1'])
                self.assertEqual(options['precision'], 192)
                cover = record_fixture(scope)['cover']
                for row in cover['leaves'].values():
                    row['witness'] = dict(fresh=True)
                options['checkpoint'](cover)
            stack.enter_context(patch.object(driver, 'fresh_model', return_value=(model, scope, identity(outer))))
            stack.enter_context(patch.object(driver.prior.cell_search, 'search', side_effect=search))
            stack.enter_context(patch.object(sys, 'argv', ['dense_t64.py', 'search', '--geometry', str(source),
                '--precision', '192', '--output', str(output)]))
            stack.enter_context(redirect_stdout(io.StringIO()))
            driver.main()
            result = json.loads(output.read_text())
            self.assertTrue(result['final_replay_required'])
            self.assertFalse(result['whole_code_certificate'])
            self.assertTrue(result['complete_search_partition'])
            self.assertEqual(result['geometry_source'], identity(source))
            self.assertEqual([r['witness'] for r in result['cover']['leaves'].values()], [dict(fresh=True)]*2)


if __name__ == '__main__':
    unittest.main()
