"""Scope and fresh-replay tests; no production proof job is run here."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from flint import arb, arb_mat
import sparse_t64 as driver
import kernel_t64
import numpy as np


def fixture():
    return dict(schema='t64-s16-uniform-gl-sparse-search-1', K=1 << 20, N=1 << 21,
        group_count=2048, outer='BCH256128', block_rows=4, block_columns=8,
        inner=driver.INNER.copy(), kernel_variant=driver.KERNEL_VARIANT, **driver.SCOPE, distance='1/10', threshold=209715,
        requested=[1], tilts=['.001'], support_min=6, count_sha256='fresh',
        count_premises=dict(scope='fresh'), map_record=dict(spectrum={'0': 1}),
        results=[dict(occupancy=1, passed=True, upper=['untrusted', 'not a bound'],
            details=dict(method='all-support exact placement', support_min=6,
                support_choices=['.001']*257, support_probability_uppers='ignored'))])


class T64SparseReplayTests(unittest.TestCase):
    def test_scope_rejects_old_state_or_update_distribution(self):
        for field, value in (('K', 1 << 18), ('inner', dict(t=128, s=19, updates=4)),
                             ('kernel_variant', 'old-t128'), ('threshold', 209716)):
            record = fixture(); record[field] = value
            with self.assertRaises(ValueError):
                driver.validate_record(record)

    def test_new_schema_requires_explicit_mix_route_count_scope(self):
        record = fixture()
        old = deepcopy(record); old['schema'] = 's16-uniform-gl-sparse-search-2'
        with self.assertRaises(ValueError):
            driver.validate_record(old)
        self.assertEqual(len(driver.validate_record(record)), 1)
        for name in driver.SCOPE:
            changed = deepcopy(record); changed[name] = 'different construction'
            with self.assertRaises(ValueError):
                driver.validate_record(changed)

    def test_rational_tilt_aliases_share_fresh_operators(self):
        data = dict(bits=16, birth_density='capped', windows=32,
                    macro_step_bits=128, physical_step_bits=64, physical_steps=2)
        with patch.object(kernel_t64, 'authenticate'), \
             patch.object(driver.kernel_t64, 'local_operators', return_value=[arb_mat([[1]])]*33) as local, \
             patch.object(driver, 'placement', return_value=[arb_mat([[1]])]) as place:
            operators = driver.build_operators(data, ['.001', '1/1000'], 192, 1)
            self.assertIs(operators['.001', '1'], operators['1/1000', '1'])
            self.assertEqual(local.call_count, 1)
            self.assertEqual(place.call_args.kwargs['epochs'], 64)
            self.assertEqual(place.call_args.kwargs['windows'], 32)

    def test_physical_and_macro_geometry_are_distinct_and_required(self):
        self.assertEqual(driver.INNER['physical_t']*driver.INNER['physical_steps'], 1 << 21)
        self.assertEqual(driver.INNER['macro_bits']*driver.INNER['macro_steps'], 1 << 21)
        for name in driver.INNER:
            changed = fixture(); changed['inner'][name] = 'wrong'
            with self.assertRaises(ValueError):
                driver.validate_record(changed)
        for field in ('bits', 'birth_density', 'physical_step_bits', 'macro_step_bits', 'physical_steps', 'windows'):
            data = dict(bits=16, birth_density='capped', windows=32,
                        macro_step_bits=128, physical_step_bits=64, physical_steps=2)
            data[field] = -1
            with self.assertRaises(ValueError):
                driver.build_operators(data, ['.001'], 192, 1)

    def test_saved_bounds_are_not_selected_as_evidence(self):
        record = fixture()
        self.assertEqual(len(driver.validate_record(record)), 1)
        record['results'][0]['details']['support_choices'][4] = '.5'
        with self.assertRaises(ValueError):
            driver.validate_record(record)

    def test_fresh_replay_ignores_saved_endpoints_and_normalizes_map_keys(self):
        self.run_replay_case()

    def test_fresh_map_mismatch_is_rejected(self):
        self.run_replay_case(changed_map=True)

    def test_source_mutation_is_rejected(self):
        self.run_replay_case(mutate_source=True)

    def test_extension_reuses_fresh_operators_only_after_screen_passes(self):
        import support_cover
        for fail_screen in (False, True):
            with self.subTest(fail_screen=fail_screen), TemporaryDirectory() as directory:
                fake_sparse = SimpleNamespace(count_hash=lambda counts: 'fresh',
                    endpoint=lambda value: [int(x) for x in value.upper().man_exp()],
                    up=lambda value: arb(value.upper()),
                    save=lambda path, data: path.write_text(json.dumps(data)))
                operators = {('.001', '1'): ([arb_mat(8, 8)]*33, [])}
                seen = []
                def cover(args, *unused, **kwargs):
                    seen.append(args.groups)
                    return (None if fail_screen and args.groups == 32 else arb(2)**-60), {}
                with patch.object(kernel_t64, 'prepare', return_value=({}, {})), \
                     patch.object(driver, 'legacy_sparse', return_value=fake_sparse), \
                     patch.object(driver.old_variant, 'authenticated_counts',
                                  return_value=([Q(0)]*257, 6, {})), \
                     patch.object(driver, 'build_operators', return_value=operators) as build, \
                     patch.object(driver.old_variant, 'one_group', return_value=(arb(2)**-60, {})), \
                     patch.object(support_cover, 'cover', side_effect=cover):
                    result = driver.run([1, 2, 8, 32], ['.001'],
                        output=Path(directory)/'new.json', precision=192, extend_full=True)
                self.assertEqual(build.call_count, 1)
                self.assertEqual(build.call_args.args[-1], 32)
                if fail_screen:
                    self.assertEqual(result['requested'], [1, 2, 8, 32])
                    self.assertEqual(seen, [2, 8, 32])
                    self.assertFalse(result['complete_sparse_prefix'])
                else:
                    self.assertEqual(set(result['requested']), set(range(1, 33)))
                    self.assertEqual(len(result['results']), 32)
                    self.assertEqual(result['aggregate_upper'], [1, -55])
                    self.assertTrue(result['complete_sparse_prefix'])
                self.assertFalse(result['whole_code_certificate'])

    def test_two_receipts_preserve_aliases_and_sum_fresh_bounds(self):
        with TemporaryDirectory() as directory:
            paths = [Path(directory)/f'{i}.json' for i in range(2)]
            first, second = fixture(), fixture()
            second.update(requested=[2], tilts=['1/1000'])
            second['results'] = [dict(occupancy=2, passed=True, upper=['ignored'], details=dict(
                method='full support-box CDF cover', support_min=6, leaves=[dict(tilt='1/1000')]))]
            for path, record in zip(paths, (first, second)):
                path.write_text(json.dumps(record))
            fake_sparse = SimpleNamespace(count_hash=lambda counts: 'fresh',
                endpoint=lambda value: [int(x) for x in value.upper().man_exp()],
                up=lambda value: arb(value.upper()), save=lambda path, data: path.write_text(json.dumps(data)))
            pair = ([arb_mat(8, 8)]*3, [])
            operators = {('.001', '1'): pair, ('1/1000', '1'): pair}
            import support_cover
            with patch.object(kernel_t64, 'prepare', return_value=({}, dict(spectrum={0: 1}))), \
                 patch.object(driver, 'legacy_sparse', return_value=fake_sparse), \
                 patch.object(driver.old_variant, 'authenticated_counts',
                              return_value=([Q(0)]*257, 6, dict(scope='fresh'))), \
                 patch.object(driver, 'build_operators', return_value=operators) as build, \
                 patch.object(driver.old_variant, 'one_group', return_value=(arb(2)**-60, {})), \
                 patch.object(support_cover, 'replay', return_value=arb(2)**-61):
                result = driver.replay(paths, Path(directory)/'fresh.json', precision=192)
                self.assertEqual(set(build.call_args.args[1]), {'.001', '1/1000'})
                self.assertEqual(result['aggregate_upper'], [3, -61])
                self.assertEqual(result['requested'], [1, 2])
                self.assertTrue(result['complete_requested'])
                self.assertFalse(result['complete_sparse_prefix'])

    def test_q2_malformed_partition_never_reaches_a_bound(self):
        import support_cover
        box = [[6, 256], [6, 256]]
        good = dict(method='full support-box CDF cover', support_min=6, domain_volume=251**2,
            tree=[dict(id=1, parent=None, children=[], box=box, multiplicity=1)],
            selected_ids=[1], leaves=[dict(box=box, multiplicity=1, tilt='.001', numerators=[1, 1])])
        changes = []
        duplicate = deepcopy(good); duplicate['tree'].append(deepcopy(duplicate['tree'][0])); changes.append(duplicate)
        mult = deepcopy(good); mult['leaves'][0]['multiplicity'] = 2; changes.append(mult)
        gap = deepcopy(good); gap['tree'][0]['box'] = [[7, 256], [7, 256]]; changes.append(gap)
        boundary = deepcopy(good); boundary['leaves'][0]['numerators'][1] = 0; changes.append(boundary)
        counts = [Q(0)]*6+[Q(1)]*251
        operators = {('.001', '1'): ([arb_mat(8, 8)]*3, [])}
        for details in changes:
            with self.assertRaises(ValueError):
                support_cover.replay(details, operators, counts, np.ones(8), groups=2,
                                     cutoff=209715, precision=192, target_bits=46)

    def run_replay_case(self, changed_map=False, mutate_source=False):
        with TemporaryDirectory() as directory:
            source, output = Path(directory)/'source.json', Path(directory)/'fresh.json'
            record = fixture()
            if changed_map:
                record['map_record']['spectrum']['0'] = 2
            source.write_text(json.dumps(record))
            fake_sparse = SimpleNamespace(count_hash=lambda counts: 'fresh',
                endpoint=lambda value: [int(x) for x in value.upper().man_exp()],
                up=lambda value: arb(value.upper()),
                save=lambda path, data: path.write_text(json.dumps(data)))
            operators = {('.001', '1'): ([arb_mat(8, 8), arb_mat(8, 8)], [])}
            def fresh_bound(*args):
                if mutate_source:
                    source.write_text(source.read_text()+' ')
                return arb(2)**-60, {}
            with patch.object(kernel_t64, 'prepare', return_value=({}, dict(spectrum={0: 1}))), \
                 patch.object(driver, 'legacy_sparse', return_value=fake_sparse), \
                 patch.object(driver.old_variant, 'authenticated_counts',
                              return_value=([Q(0)]*257, 6, dict(scope='fresh'))), \
                 patch.object(driver, 'build_operators', return_value=operators), \
                 patch.object(driver.old_variant, 'one_group', side_effect=fresh_bound):
                if changed_map or mutate_source:
                    with self.assertRaises(ValueError):
                        driver.replay([source], output, precision=192)
                else:
                    result = driver.replay([source], output, precision=192)
                    self.assertEqual(result['results'][0]['upper'], [1, -60])
                    self.assertTrue(result['complete_requested'])
                    self.assertFalse(result['complete_sparse_prefix'])
                    self.assertFalse(result['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()
