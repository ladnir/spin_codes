import copy
from contextlib import redirect_stderr, redirect_stdout
from fractions import Fraction as Q
import io
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

from flint import arb, ctx
import hill_probe as hill
from test_cell_search import source


class HillTests(unittest.TestCase):
    def test_retarget_discards_old_bounds_and_hints(self):
        old = source()
        old['comparison'] = 'direct-expected-shell-majorant'
        old['probes'] = [{'upper': [1, -999]}]
        old['point_hint_sources'] = ['stale']
        before = copy.deepcopy(old)
        new = hill.retarget(old, '.1')
        self.assertEqual(new['threshold'], 209715)
        self.assertEqual(new['distance'], '1/10')
        self.assertEqual(new['cover']['leaves'], {})
        self.assertEqual(set(new['cover']['unresolved']), {''})
        self.assertNotIn('probes', new)
        self.assertNotIn('point_hint_sources', new)
        self.assertEqual(old, before)
        new['mixture'][0]['mass'] = '2'
        self.assertEqual(old['mixture'][0]['mass'], '1')

    def test_no_cdf_comparison_or_invalid_distance(self):
        old = source()
        with self.assertRaises(ValueError):
            hill.retarget(old, '.1')
        old['comparison'] = 'direct-expected-shell-majorant'
        for distance in ('0', '-.01', '.5'):
            with self.assertRaises(ValueError):
                hill.retarget(old, distance)

    def test_refined_factory_rebuilds_comparison_and_separates_receipt_scope(self):
        old = source()
        old['comparison'] = 'direct-expected-shell-majorant'
        previous = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', previous)
        premises = dict(canonical_cdf=[0, 2], authenticated=True)
        created = []

        def model_factory(*args, **kwargs):
            created.append((args, kwargs))
            return SimpleNamespace(root=(Q(1, 100), Q(1)))

        with patch('outer_hill_incidence.authenticated_bch_cdf', return_value=((0, 2), premises)) as auth, \
                patch('monotone.transport_shells', return_value=(0, 3)) as shells, \
                patch.object(hill.cell_search.dense, 'exact_mixture',
                    return_value=([(Q(7), Q(1, 4))], {'fresh': True})) as mixture, \
                patch.dict('sys.modules', birth_classes=SimpleNamespace(actual=lambda r: {'updates': r}),
                    shared_mixture=SimpleNamespace(as_components=lambda values: values),
                    scalar_cover=SimpleNamespace(Model=model_factory)):
            _, scope = hill.refined_model(old, '.1', 256, 64)
        auth.assert_called_once_with()
        self.assertEqual(shells.call_args.args[0], [0, 2])
        mixture.assert_called_once_with((0, 3))
        self.assertEqual(created[0][0][2], 209715)
        self.assertEqual(created[0][1]['variance_bins'], 64)
        self.assertEqual(scope['mixture'], [dict(mass='7', activity='1/4')])
        self.assertEqual(scope['variance_bins'], 64)
        self.assertEqual(scope['outer_premises'], premises)
        self.assertEqual(scope['root'], ['1/100', '1'])
        with self.assertRaises(ValueError):
            hill.cell_search.dense.validate_record(scope)

    def test_all_requested_refinements_run_even_after_passing(self):
        original = dict(tilt='3/16', parameters=['2', '0', '0'])
        model = SimpleNamespace(tilt=Q(3, 16),
            propose_with=lambda *_: (-1000., copy.deepcopy(original)))
        calls, events = [], []

        def propose(_model, _cell, witness):
            calls.append(copy.deepcopy(witness))
            return -100.-float(Q(witness['parameters'][0])), witness

        variance = SimpleNamespace(propose=lambda _m, _c, b: b)
        regional = SimpleNamespace(prepare_witness=lambda _m, _c, w: ([], w), propose=propose)
        with patch.dict('sys.modules', variance_partition=variance, regional_count=regional):
            trials, best = hill.candidates(model, (Q(1, 10),)*2, ['1/4', '1'],
                ['classified', 'uniform'], lambda row, _: events.append(row))
        self.assertEqual(len(calls), 4)
        self.assertEqual(trials, events)
        self.assertEqual([w['parameters'][0] for w in calls], ['1/2', '1/2', '2', '2'])
        self.assertNotIn('regional_feedback_uniform_classes', calls[0])
        self.assertTrue(calls[1]['regional_feedback_uniform_classes'])
        self.assertNotIn('regional_feedback_uniform_classes', calls[2])
        self.assertEqual(best[0], -1000.)
        self.assertNotIn('regional_count_parts', best[1])
        self.assertEqual(original, dict(tilt='3/16', parameters=['2', '0', '0']))

    def test_actual_updates_and_ekr_factory_are_propagated_without_r2_claims(self):
        old = source()
        old.update(comparison='direct-expected-shell-majorant',
                   whole_code_certificate=True, certified_upper=[1, -100])
        before = copy.deepcopy(old)
        previous = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', previous)
        for updates in (2, 3, 4):
            for intersection in (False, True):
                with self.subTest(updates=updates, intersection=intersection):
                    actual = Mock(side_effect=lambda r: {'updates': r, 'fresh': True})
                    factory = Mock(return_value=SimpleNamespace(root=(Q(1, 100), Q(1))))
                    incidence = dict(canonical_cdf=[0, 2], schema='synthetic-incidence')
                    ekr = dict(canonical_cdf=[0, 3], schema='synthetic-ekr')
                    with patch('outer_hill_incidence.authenticated_bch_cdf',
                               return_value=((0, 2), incidence)) as auth_incidence, \
                            patch('outer_hill_intersection.authenticated_bch_cdf',
                                  return_value=((0, 3), ekr)) as auth_ekr, \
                            patch('monotone.transport_shells', return_value=(0, 4)) as transport, \
                            patch.object(hill.cell_search.dense, 'exact_mixture',
                                return_value=([(Q(7), Q(1, 4))], {'fresh': True})), \
                            patch.dict(sys.modules,
                                birth_classes=SimpleNamespace(actual=actual),
                                shared_mixture=SimpleNamespace(as_components=lambda values: values),
                                scalar_cover=SimpleNamespace(Model=factory)):
                        _, scope = hill.refined_model(old, '.1', 256, 16,
                            updates=updates, intersection=intersection)
                    actual.assert_called_once_with(updates)
                    self.assertEqual(factory.call_args.args[1], {'updates': updates, 'fresh': True})
                    self.assertEqual(scope['updates'], updates)
                    self.assertEqual(scope['ensemble'], f'canonical-gl32-width8-shared4-r{updates}')
                    self.assertEqual(scope['outer_premises'], ekr if intersection else incidence)
                    transport.assert_called_once()
                    self.assertEqual(transport.call_args.args[0], [0, 3] if intersection else [0, 2])
                    (auth_ekr if intersection else auth_incidence).assert_called_once_with()
                    (auth_incidence if intersection else auth_ekr).assert_not_called()
                    self.assertEqual(scope['cover']['leaves'], {})
                    self.assertEqual(set(scope['cover']['unresolved']), {''})
                    self.assertNotIn('certified_upper', scope)
                    self.assertNotIn('whole_code_certificate', scope)
                    with self.assertRaises(ValueError):
                        hill.cell_search.dense.validate_record(scope)
        self.assertEqual(old, before)

    def test_invalid_actual_update_and_intersection_flags_fail_before_authentication(self):
        old = source()
        old['comparison'] = 'direct-expected-shell-majorant'
        with patch('outer_hill_incidence.authenticated_bch_cdf') as incidence, \
                patch('outer_hill_intersection.authenticated_bch_cdf') as ekr:
            for updates in (True, 1, 5, '3', 3.0):
                with self.assertRaises(ValueError):
                    hill.refined_model(old, '.1', 256, 16, updates=updates)
            for intersection in (0, 1, 'true', None):
                with self.assertRaises(ValueError):
                    hill.refined_model(old, '.1', 256, 16, intersection=intersection)
        incidence.assert_not_called()
        ekr.assert_not_called()

    def test_birth_class_toy_operators_use_requested_refresh_probability(self):
        # Exercise the real preparation/evaluation chain on an eight-state
        # toy, not the expensive production census. R changes the actual
        # lazy-state coefficient, not merely the receipt's label.
        parent = Path(hill.__file__).resolve().parents[1]
        previous = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', previous)
        ctx.prec = 192
        with patch.object(sys, 'path', [str(parent/'gf16_packets'), str(parent), *sys.path]):
            import birth_classes as birth
            import occupancy_birth_classes as occupancy
        for updates in (2, 3, 4):
            with self.subTest(updates=updates):
                data = birth.prepare(list(range(8)), [1, 2, 4, 3, 5, 7, 6, 1], 3, updates)
                self.assertEqual(data['updates'], updates)
                zero_input = occupancy.outward_at_z(data, arb(3)/4)[0]
                for index, level in enumerate(data['birth_class_levels'], 3):
                    expected = arb(2)**(-updates)*(arb(3)/4)**int(level)
                    self.assertGreaterEqual(zero_input[index, index], expected)
                    self.assertLess(float(zero_input[index, index]-expected), 1e-40)

    def test_cli_requires_fresh_factory_for_changed_construction(self):
        for flags in (['--updates', '3'], ['--updates', '4'], ['--intersection'],
                      ['--updates', '1', '--outer-refinement'], ['--variance-bins', '65']):
            with self.subTest(flags=flags), \
                    patch.object(sys, 'argv', ['hill_probe.py', 'unread-source.json',
                        '--output', 'unused-output.json', *flags]), \
                    patch.object(Path, 'exists', return_value=False), \
                    patch.object(Path, 'read_bytes') as read, \
                    redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
                hill.main()
            self.assertEqual(error.exception.code, 2)
            read.assert_not_called()

    def test_quick_mode_checks_fresh_bound_and_labels_changed_point_only_scope(self):
        old = source()
        old['comparison'] = 'direct-expected-shell-majorant'
        raw = json.dumps(old).encode()
        scope = dict(schema='packed-canonical-gl32-hill-dense-context-1', updates=3,
                     ensemble='canonical-gl32-width8-shared4-r3')
        witness = {'parameters': ['1/2', '0', '0']}
        model = SimpleNamespace(root=(Q(1, 100), Q(1)),
            proposal=Mock(return_value=(-999., witness)),
            outward=Mock(return_value=arb(2)**-70))
        writes = []
        previous = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', previous)
        with patch.object(sys, 'argv', ['hill_probe.py', 'synthetic-source.json',
                    '--output', 'synthetic-output.json', '--means', '.032', '--quick',
                    '--updates', '3', '--outer-refinement', '--intersection']), \
                patch.object(Path, 'exists', return_value=False), \
                patch.object(Path, 'read_bytes', return_value=raw), \
                patch.object(Path, 'mkdir'), patch.object(Path, 'replace'), \
                patch.object(Path, 'write_text', autospec=True,
                    side_effect=lambda _path, text: writes.append(json.loads(text))), \
                patch.object(hill, 'refined_model', return_value=(model, scope)) as fresh, \
                patch.object(hill, 'candidates', side_effect=AssertionError('quick must use adaptive proposal')), \
                redirect_stdout(io.StringIO()):
            hill.main()
        fresh.assert_called_once_with(old, '.1', 256, 16, updates=3, intersection=True)
        model.proposal.assert_called_once_with((Q('.032'), Q('.032')))
        model.outward.assert_called_once_with((Q('.032'), Q('.032')), witness)
        self.assertEqual(model.proposal_stop_bits, 54)
        record = writes[-1]
        self.assertTrue(record['changed_construction'])
        self.assertTrue(record['partial_only'])
        self.assertFalse(record['whole_code_certificate'])
        self.assertEqual(record['search'], 'adaptive')
        self.assertEqual(record['scope']['updates'], 3)
        self.assertEqual(record['count_refinement'], 'H5 exact + incidence + EKR')
        self.assertEqual(record['points'][0]['trials'], [])
        pair = record['points'][0]['checked']['upper']
        self.assertEqual(Q(pair[0])*Q(2)**pair[1], Q(2)**-70)

    def test_outward_is_fresh_and_requires_finite_positive_value(self):
        previous = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', previous)
        calls = []
        model = SimpleNamespace(outward=lambda cell, witness:
            calls.append((cell, witness)) or arb(2)**-60)
        cell, witness = (Q(1, 10),)*2, {'new': True}
        checked = hill.check_point(model, cell, (-1e6, witness), 192)
        self.assertEqual(calls, [(cell, witness)])
        self.assertEqual(Q(checked['upper'][0])*Q(2)**checked['upper'][1], Q(2)**-60)
        for bad in (arb(0), arb(-1), arb('nan')):
            model.outward = lambda *_: bad
            with self.assertRaises(ArithmeticError):
                hill.check_point(model, cell, (-100., witness), 192)


if __name__ == '__main__':
    unittest.main()
