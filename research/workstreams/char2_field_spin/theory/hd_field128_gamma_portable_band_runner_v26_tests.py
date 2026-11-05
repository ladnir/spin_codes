"""Pure v26 validation/provenance tests; no coefficients, banks or proof jobs."""
from hashlib import sha256
from pathlib import Path, PurePosixPath, PureWindowsPath
import unittest
from unittest.mock import patch

import hd_field128_gamma_portable_band_runner_v26 as runner


class WrapperTests(unittest.TestCase):
    def test_explicit_band_boundaries_and_terminal_rejection(self):
        for args in ((1, 64, 4096), (385, 439, 4096), (576, 639, 2048),
                     (640, 703, 512), (704, 704, 1), (8128, 8191, 4096)):
            runner.validate_band(*args)
            self.assertEqual(runner.controller.bands.checked_weights(4096, args[0], args[1]),
                             tuple(range(args[0], args[1]+1)))
        for args in ((0, 1, 4096), (1, 65, 4096), (640, 704, 4096),
                     (8192, 8192, 4096), (8129, 8192, 4096), (4, 3, 4096),
                     (True, 3, 4096), (1, 3, True), (1, 3, 0), (1, 3, 4097)):
            with self.subTest(args=args), self.assertRaisesRegex(ValueError, '8192'):
                runner.validate_band(*args)

    def test_only_orchestration_identity_and_assignment_fields_change(self):
        pins = {runner.NAME: '1'*64, runner.controller.NAME: '2'*64}
        witness = dict(path='portable-seed.json', sha256='3'*64)
        old = runner.controller.driver_record(pins, witness, 385, 439, 4096, 400)
        new = runner.driver_record(pins, witness, 385, 439, 4096, 400)
        differing = {key for key in old.keys() & new.keys() if old[key] != new[key]}
        self.assertEqual(differing, {'schema', 'driver'})
        self.assertEqual(new['partition_controller'], runner.controller.NAME)
        self.assertEqual(new['partition_controller_sha256'], pins[runner.controller.NAME])
        self.assertEqual(new['source_sha256'], pins)
        self.assertEqual(new['witness_source'], witness)
        self.assertEqual(new['requested_weight_range'], [385, 439])
        self.assertEqual(new['source_calibration_weight'], 400)
        self.assertFalse(new['terminal_weight_8192_supported'])
        for link in (new['driver'], new['witness_source']):
            self.assertFalse(PurePosixPath(link['path']).is_absolute())
            self.assertFalse(PureWindowsPath(link['path']).is_absolute())
        for key in ('directed', 'whole_code_certificate', 'construction_changed',
                    'old_numerical_bounds_reused', 'fresh_bank_rebuilt'):
            self.assertIs(new[key], False)

    def test_unchanged_limits_and_frozen_controller_pin(self):
        pins = runner.source_pins()
        self.assertEqual(pins[runner.controller.NAME], runner.FROZEN_CONTROLLER)
        self.assertEqual(pins[runner.NAME], sha256(Path(runner.__file__).read_bytes()).hexdigest())
        runner.controller.check_source_pins(pins)
        record = runner.driver_record(pins, dict(path='seed.json', sha256='0'*64), 385, 439, 4096, 400)
        self.assertEqual((record['cell_failure_bits'], record['aggregate_failure_bits']), (72, 54))
        self.assertEqual((record['precision'], record['coefficient_bits'], record['grid']), (256, 80, 128))
        self.assertEqual(record['max_total_compute_seconds'], 600)
        self.assertEqual(record['max_output_bytes'], 32 << 20)
        self.assertEqual(record['max_table_bytes'], 64 << 20)
        self.assertEqual(record['max_address_space_bytes'], 4 << 30)
        self.assertEqual(record['disk_reserve_bytes'], 20 << 30)
        self.assertEqual(record['partition_policy'], runner.controller.POLICY)
        self.assertEqual(record['checkpoint_format'], runner.controller.CHECKPOINT_FORMAT)

    def test_frozen_controller_mutation_fails_before_any_work(self):
        read = Path.read_bytes
        target = Path(runner.controller.__file__)
        def changed(path):
            value = read(path)
            return value+b'\n# simulated source mutation\n' if path == target else value
        with patch.object(Path, 'read_bytes', changed):
            with self.assertRaisesRegex(RuntimeError, 'Frozen v25 controller'):
                runner.source_pins()

    def test_invalid_run_never_loads_a_seed_or_builds_a_bank(self):
        # Validation is the first operation, before a seed path, output directory,
        # memory limit or numerical producer can be touched.
        with patch.object(runner, 'source_pins', side_effect=AssertionError('unexpected work')):
            for first, last in ((0, 1), (1, 65), (8192, 8192)):
                with self.assertRaises(ValueError):
                    runner.run('unopened.json', '0'*64, first, last, 'unused-output')


if __name__ == '__main__':
    unittest.main()
