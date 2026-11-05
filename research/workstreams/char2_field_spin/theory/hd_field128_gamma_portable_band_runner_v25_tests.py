"""Pure controller/resource tests: no AliasBank or coefficient calculation."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

import hd_field128_gamma_portable_band_runner_v25 as runner


def seed():
    return dict(schema='gamma-portable-exact-witness-v1', field_order=1 << 128,
        construction=runner.cert.local.construction(30), active_groups=376,
        state_dimension=30, target_distance='2/5', directed=False, whole_code_certificate=False,
        witness=dict(z='1/4', alpha='1/8', markers=[['1/2', '3/4'], ['2', '3']],
                     grid=128, state_dimension=30),
        source=dict(path='C:/old/unavailable/source.json', sha256='1'*64))


def synthetic(roots=None, *, zero=False, fail=None):
    roots = roots or (runner.bands.Rectangle(0, ((1, 2), (0, 3), (0, 2))),
                      runner.bands.Rectangle(1, ((1, 2), (0, 3), (0, 2))))
    calls = []

    def receipt(cell, bits, budget):
        calls.append(cell)
        if fail:
            raise fail('synthetic failure')
        singleton = all(lo == hi for lo, hi in cell.ranges)
        ceiling = '-200' if singleton else '0'
        count = 1
        for lo, hi in cell.ranges:
            count *= hi-lo+1
        return dict(cell=cell.descriptor(), status='complete', complete_cell=True,
            count_upper=0 if zero else count, zero=zero,
            log2_upper_rational=None if zero else ceiling,
            log2_upper_enclosure=None if zero else ceiling,
            failure_bits=bits, meets_cell_budget=zero or F(ceiling) <= -bits)

    table = SimpleNamespace(roots=roots, root_by_parity={x.parity: x for x in roots},
                            precision=128, receipt=receipt)
    return SimpleNamespace(tables=table, _active=False, check_pins=Mock(), calls=calls,
        context=lambda: dict(requested_weight_range=[2, 5], requested_weights=[2, 3, 4, 5],
                             synthetic_test_only=True))


def points(cell):
    return {(cell.parity, *point) for point in product(*(range(a, b+1) for a, b in cell.ranges))}


class ControllerTests(unittest.TestCase):
    def test_absolute_current_width_and_ties(self):
        cell = lambda widths: runner.bands.Rectangle(0, tuple((0, x) for x in widths))
        for widths, expected in (((7, 191, 191), 1), ((7, 95, 191), 2),
                                 ((7, 3, 3), 0), ((7, 7, 7), 1),
                                 ((0, 0, 5), 2), ((3, 0, 0), 0)):
            self.assertEqual(runner.split_axis(cell(widths)), expected)
        with self.assertRaises(ValueError):
            runner.split_axis(object())

    def test_complete_partition_compact_checkpoints_and_union(self):
        value, checkpoints, contexts = synthetic(), [], []
        result = runner.cover(value, max_nodes=256, checkpoint=checkpoints.append,
                              context_checkpoint=contexts.append)
        self.assertTrue(result['complete_requested_band'])
        self.assertFalse(result['directed'])
        self.assertEqual(result['numerical_backend'], 'synthetic-test-only')
        self.assertEqual(result['pending_frontier'], [])
        self.assertEqual(result['checkpoints'], len(result['nodes'])+1)
        self.assertEqual(value.check_pins.call_count, 2)
        replay = {}
        for i, record in enumerate(checkpoints):
            self.assertNotIn('pending_frontier', record)
            self.assertEqual(record['sequence'], i)
            self.assertEqual(record['context_sha256'], result['context_sha256'])
            for key in ('directed', 'complete_requested_band', 'covers_all_supports_in_band',
                        'batch_source_verified_at_exit', 'whole_code_certificate'):
                self.assertIs(record[key], False)
            self.assertIsNone(record['log2_upper_rational'])
            replay.update((row['id'], row) for row in record['changed_nodes'])
        self.assertEqual([replay[i] for i in sorted(replay)], result['nodes'])
        expected = set().union(*(points(cell) for cell in value.tables.roots))
        seen = set()
        for identity in result['accepted_leaves']:
            row = result['nodes'][identity]
            desc = row['cell']
            cell = runner.bands.Rectangle(desc['parity'], tuple(map(tuple, desc['u_deficit_ranges'])))
            self.assertFalse(seen & points(cell))
            seen |= points(cell)
        self.assertEqual(seen, expected)
        self.assertEqual(len(result['accepted_leaves']), len(expected))
        self.assertLess(F(result['log2_upper_rational']), -190)
        self.assertEqual(result['controller']['path'], runner.NAME)
        self.assertEqual(result['controller_source_sha256'][runner.NAME], result['controller']['sha256'])
        self.assertEqual(result['partition_policy'], runner.POLICY)
        self.assertEqual(result['context'], contexts[0]['context'])
        # Duplicate/omitted/misparented children are still rejected by frozen audit_tree.
        for mode in ('duplicate', 'omit', 'parent'):
            bad = deepcopy(result['nodes'])
            if mode == 'duplicate':
                bad.append(deepcopy(bad[-1]))
            elif mode == 'omit':
                bad.pop()
            else:
                bad[-1]['parent'] = None
            with self.assertRaises(ValueError):
                runner.bands.audit_tree(bad, value.tables.roots, 72, result['context_sha256'])

    def test_unchanged_receipts_union_and_all_other_checkpoint_fields(self):
        roots = (runner.bands.Rectangle(0, ((1, 1), (0, 5), (0, 0))),)
        old_records, new_records = [], []
        old = runner.bands.cover(synthetic(roots), max_nodes=64, cell_failure_bits=72,
                                 aggregate_bits=54, checkpoint=old_records.append)
        new = runner.cover(synthetic(roots), max_nodes=64, checkpoint=new_records.append)
        for record in old_records:
            record.pop('pending_frontier')
        self.assertEqual(old_records, new_records)
        for key in ('nodes', 'accepted_leaves', 'root_cells', 'context', 'context_sha256',
                    'zero', 'log2_upper_rational', 'log2_upper_enclosure', 'pending_frontier',
                    'complete_requested_band', 'meets_requested_band_target'):
            self.assertEqual(old[key], new[key], key)

    def test_partial_timeout_arithmetic_zero_and_singleton(self):
        partial = runner.cover(synthetic(), max_nodes=1)
        self.assertFalse(partial['complete_requested_band'])
        self.assertTrue(partial['pending_frontier'])
        self.assertIsNone(partial['log2_upper_rational'])
        self.assertFalse(partial['directed'])
        self.assertEqual(partial['stop_reason'], 'max_nodes')
        for error in (TimeoutError, ArithmeticError):
            value = synthetic(fail=error)
            result = runner.cover(value)
            self.assertFalse(result['complete_requested_band'])
            self.assertEqual(len(result['pending_frontier']), 2)
            self.assertTrue(all(row['status'] == 'unresolved' for row in result['nodes']))
        zero = runner.cover(synthetic(zero=True))
        self.assertTrue(zero['zero'])
        self.assertIsNone(zero['log2_upper_rational'])
        value = synthetic((runner.bands.Rectangle(0, ((0, 0),)*3),))
        original = value.tables.receipt
        value.tables.receipt = lambda *args: dict(original(*args), log2_upper_rational='0',
                                                   log2_upper_enclosure='0', meets_cell_budget=False)
        result = runner.cover(value)
        self.assertEqual(result['pending_frontier'], [0])
        self.assertEqual(result['nodes'][0]['status'], 'unresolved')

    def test_limits_source_and_exception_guards(self):
        for options in (dict(max_nodes=4097), dict(max_nodes=True), dict(max_seconds=601),
                        dict(max_seconds=float('nan')), dict(cell_failure_bits=True)):
            with self.assertRaises(ValueError):
                runner.cover(synthetic(), **options)
        value = synthetic()
        with self.assertRaises(OSError):
            runner.cover(value, checkpoint=Mock(side_effect=OSError('output budget')))
        self.assertFalse(value._active)
        self.assertEqual(value.check_pins.call_count, 2)
        value = synthetic()
        def reenter(row):
            runner.cover(value)
        with self.assertRaisesRegex(RuntimeError, 'already active'):
            runner.cover(value, checkpoint=reenter)
        self.assertFalse(value._active)
        value = synthetic()
        check = runner.check_source_pins
        changed = [False]
        def guard(pins):
            if changed[0]:
                raise RuntimeError('Pinned v25 source changed')
            check(pins)
        with patch.object(runner, 'check_source_pins', side_effect=guard):
            with self.assertRaisesRegex(RuntimeError, 'Pinned v25'):
                runner.cover(value, checkpoint=lambda row: changed.__setitem__(0, True))
        self.assertFalse(value._active)


class PortableTests(unittest.TestCase):
    def test_assignment_seed_and_authentic_metadata(self):
        for args in ((368, 383, 4096), (376, 376, 512), (368, 375, 2048)):
            runner.validate_band(*args)
        for args in ((367, 383, 512), (368, 384, 512), (376, 368, 512),
                     (368, 383, 4097), (368, 383, True)):
            with self.assertRaises(ValueError):
                runner.validate_band(*args)
        with tempfile.TemporaryDirectory() as directory, patch.object(runner, 'HERE', Path(directory)):
            raw = json.dumps(seed()).encode()
            Path(directory, 'seed.json').write_bytes(raw)
            value, parameters, source = runner.load_seed('seed.json', sha256(raw).hexdigest())
            self.assertEqual(value, seed())
            self.assertEqual(parameters[:2], (F(1, 4), F(1, 8)))
            self.assertEqual(source['path'], 'seed.json')
            with self.assertRaises(ValueError):
                runner.load_seed('../seed.json', sha256(raw).hexdigest())
            with self.assertRaises(ValueError):
                runner.load_seed('seed.json', '0'*64)
        pins = runner.source_pins()
        runner.check_source_pins(pins)
        self.assertEqual(pins[Path(runner.bands.__file__).name], runner.FROZEN_BANDS)
        record = runner.driver_record(pins, source, 368, 383, 4096, 376)
        self.assertEqual(record['partition_controller'], runner.NAME)
        self.assertEqual(record['numerical_provider'], Path(runner.bands.__file__).name)
        self.assertEqual(record['requested_weight_range'], [368, 383])
        self.assertEqual((record['cell_failure_bits'], record['aggregate_failure_bits']), (72, 54))
        self.assertEqual(record['max_total_compute_seconds'], 600)
        self.assertEqual(record['max_address_space_bytes'], 4 << 30)
        self.assertEqual(record['max_table_bytes'], 64 << 20)
        self.assertEqual(record['max_output_bytes'], 32 << 20)
        self.assertTrue(record['final_frontier_preserved'])
        self.assertFalse(record['construction_changed'])
        self.assertFalse(record['fresh_bank_rebuilt'])

    def test_output_limit_reserve_exact_bytes_and_no_completed_summary(self):
        checkpoint = dict(complete_requested_band=False, directed=False,
                          log2_upper_rational=None, validation_only='\u03b1')
        size = len((json.dumps(checkpoint, indent=2, allow_nan=False)+'\n').encode())
        usage = SimpleNamespace(free=runner.RESERVE+4096)
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(runner.shutil, 'disk_usage', return_value=usage), \
                patch.object(runner, 'MAX_OUTPUT_BYTES', 2*size):
            root = Path(directory)
            writer = runner.OutputBudget(root)
            usage.free = runner.RESERVE
            with self.assertRaises(RuntimeError):
                writer.write_new(root/'first.json', checkpoint)
            self.assertEqual(writer.used_bytes, 0)
            usage.free += size
            writer.write_new(root/'node_000000.json', checkpoint)
            writer.write_new(root/'node_000001.json', checkpoint)
            self.assertEqual(writer.used_bytes, 2*size)
            self.assertEqual(sum(path.stat().st_size for path in root.iterdir()), 2*size)
            with self.assertRaises(RuntimeError):
                writer.write_new(root/'summary.json', {'complete_requested_band': True})
            self.assertFalse((root/'summary.json').exists())
            with self.assertRaises(FileExistsError):
                writer.write_new(root/'node_000000.json', checkpoint)
            with self.assertRaises(ValueError):
                writer.write_new(root/'..'/'outside.json', checkpoint)
        usage.free = runner.RESERVE+runner.MAX_OUTPUT_BYTES-1
        with patch.object(runner.shutil, 'disk_usage', return_value=usage):
            with self.assertRaises(RuntimeError):
                runner.reserve_guard(runner.HERE, runner.MAX_OUTPUT_BYTES)
            usage.free += 1
            runner.reserve_guard(runner.HERE, runner.MAX_OUTPUT_BYTES)

    def test_memory_limit_mock_only_and_shared_deadline(self):
        limits = [(8 << 30, 16 << 30)]
        fake = SimpleNamespace(RLIMIT_AS=9, RLIM_INFINITY=-1,
            getrlimit=lambda unused: limits[0],
            setrlimit=lambda unused, value: limits.__setitem__(0, value))
        with patch.object(runner.sys, 'platform', 'linux'), patch.dict(runner.sys.modules, {'resource': fake}):
            result = runner.install_memory_limit()
            self.assertEqual(limits[0], (4 << 30, 16 << 30))
            self.assertEqual(result['applied_bytes'], 4 << 30)
            limits[0] = (2 << 30, 3 << 30)
            runner.install_memory_limit()
            self.assertEqual(limits[0], (2 << 30, 3 << 30))
        with patch.object(runner.sys, 'platform', 'win32'):
            with self.assertRaises(RuntimeError):
                runner.install_memory_limit()
        clock = SimpleNamespace(check=Mock(), max_seconds=600, started=10.)
        with patch.object(runner, 'monotonic', return_value=110.):
            self.assertEqual(runner.remaining_seconds(clock), 500)
            self.assertEqual(runner.remaining_seconds(clock, 180), 180)
        clock.check.side_effect = TimeoutError('expired')
        with self.assertRaises(TimeoutError):
            runner.remaining_seconds(clock)

    def test_seed_and_controller_byte_mutation(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(runner, 'HERE', Path(directory)):
            root = Path(directory)
            (root/'controller.py').write_bytes(b'source')
            (root/'seed.json').write_bytes(b'{}')
            pins = {'controller.py': sha256(b'source').hexdigest()}
            witness = dict(path='seed.json', sha256=sha256(b'{}').hexdigest())
            runner.check_pins(pins, witness)
            (root/'seed.json').write_bytes(b'{"changed":true}')
            with self.assertRaises(RuntimeError):
                runner.check_pins(pins, witness)
            (root/'seed.json').write_bytes(b'{}')
            (root/'controller.py').write_bytes(b'changed')
            with self.assertRaises(RuntimeError):
                runner.check_pins(pins, witness)


if __name__ == '__main__':
    unittest.main()
