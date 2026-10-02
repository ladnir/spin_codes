"""Provenance and measurement-scope regression tests; never run benchmarks."""
from contextlib import ExitStack
from copy import deepcopy
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import implementation_t64 as binding


FILES = ('GfniT64.h', 'FusedR4Gfni.h', 'gfni_t64_probe.cpp',
         'gfni_t64_codegen.py', 'gfni_t64_run.sh')
CPP_HASH = 'c'*64
SOURCE_HASH = 'd'*64
ROWS, COLUMNS = [1, 2], [1, 2, 3]
T64_CHECK = ('checks PASS: explicit A64 zeta basis, C=A^T, grouped emission/raw feedback, '
    'every64-step update/transpose basis, independent original updates, inner forward/reverse, '
    'scalar outer/direct forward BCH, inner+FULL adjoint, complete output/suffix; '
    'step=64; state=16; refresh=uniform-gl16; map_sha256='+CPP_HASH+
    '; physical_epochs={epochs}; hot_update_bytes=1048576; sampled_matrix_draws=113315; tile_mode=1\n')
CONTROL_CHECK = ('checks PASS: every-epoch basis, all3 inner modes/original dense+forward adjoint, '
    'original sequential/full scalar GL32+production BCH reference,suffix; '
    'gfni_row_bytes=48; packed_state_bytes=384; gfni_setup_bytes=786432; tile_mode=1\n')


class Fixture:
    def __init__(self, directory):
        self.local, self.archive = Path(directory)/'local', Path(directory)/'archive'
        self.local.mkdir(); self.archive.mkdir()
        for name in FILES:
            value = '// generated fixture\n' if name == 'GfniT64.h' else name+'\n'
            (self.local/name).write_text(value)
            (self.archive/name).write_text(value)
        (self.archive/'gfni-t64').write_bytes(b'fixture binary, not executable')
        self.build = self.archive/'measurements/t64-build-one'
        self.check = self.archive/'measurements/t64-check-one'
        self.confirm = self.archive/'measurements/t64-confirm-one'
        for path in (self.build, self.check, self.confirm):
            path.mkdir(parents=True)
        names = ('GfniT64.h', 'FusedR4Gfni.h', 'gfni_t64_probe.cpp', 'gfni-t64')
        self.hashes = self.build/'hashes.txt'
        self.hashes.write_text(''.join(binding.digest(self.archive/name)+'  /build/'+name+'\n' for name in names))
        for name, epochs in (('full-14-1.txt', 512), ('full-14-17.txt', 512), ('full-20-1.txt', 32768)):
            (self.check/name).write_text(T64_CHECK.format(epochs=epochs))
        for i, (seed, rep) in enumerate(((1, 0), (1, 1), (17, 0), (17, 1))):
            (self.confirm/f'{seed}-{rep}-control.csv').write_text(CONTROL_CHECK+
                'exponent,mode,seed,calls,median_ms,p10_ms,p90_ms,checksum,updates,tile_mode\n'+
                f'20,1,{seed},101,{5+i/10:.1f},4.9,5.5,deadbeefdeadbeef,4,1\n')
            (self.confirm/f'{seed}-{rep}-t64.csv').write_text(T64_CHECK.format(epochs=32768)+
                'exponent,step,state,refresh,seed,calls,median_ms,p10_ms,p90_ms,checksum,tile_mode\n'+
                f'20,64,16,uniform-gl16,{seed},101,{6+i/10:.1f},5.9,6.5,deadbeefdeadbeef,1\n')
        self.wrapper = dict(map_sha256='p'*64)
        self.maps = dict(expansion_rows_hex=list(map(hex, ROWS)), feedback_columns=COLUMNS,
                         source=dict(sha256=SOURCE_HASH))
        self.prepared = (ROWS, COLUMNS, None, None, None, CPP_HASH, SOURCE_HASH)

    def verify(self, *, maps=None, prepared=None):
        with ExitStack() as stack:
            stack.enter_context(patch.object(binding, 'PERMUTATION', self.local))
            stack.enter_context(patch.object(binding.kernel_t64, 'prepare',
                return_value=(self.wrapper, self.maps if maps is None else maps)))
            stack.enter_context(patch.object(binding.codegen, 'prepare',
                return_value=self.prepared if prepared is None else prepared))
            stack.enter_context(patch.object(binding.codegen, 'main',
                side_effect=lambda: print('// generated fixture')))
            return binding.verify(self.archive)


class ImplementationBindingTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.fixture = Fixture(self.temporary.name)

    def assertRejected(self):
        with self.assertRaises((ValueError, ArithmeticError, FileNotFoundError)):
            self.fixture.verify()

    def test_fresh_explicit_map_binding_and_measurement_scope(self):
        result = self.fixture.verify()
        self.assertEqual(result['map_proof_sha256'], 'p'*64)
        self.assertEqual(result['map_cpp_sha256'], CPP_HASH)
        self.assertTrue(result['map_hashes_bind_identical_explicit_matrices'])
        self.assertEqual(result['element_bits'], 128)
        self.assertEqual(result['physical_step_bits'], 64)
        self.assertEqual(result['state_bits'], 16)
        self.assertAlmostEqual(result['median_of_process_medians_ms']['control'], 5.15)
        self.assertAlmostEqual(result['median_of_process_medians_ms']['t64'], 6.15)
        self.assertAlmostEqual(result['relative_latency'], 6.15/5.15)
        self.assertFalse(result['whole_code_certificate'])
        self.assertIn('setup excluded', result['timing_scope'])
        self.assertIn('deterministic setup', result['limitation'])
        self.assertIn('ideal independent uniform setup', result['limitation'])

    def test_changed_expansion_feedback_or_map_source_rejected(self):
        for field, value in (('expansion_rows_hex', ['0x2', '0x1']),
                             ('feedback_columns', [1, 3, 2]),
                             ('source', dict(sha256='e'*64))):
            with self.subTest(field=field):
                changed = deepcopy(self.fixture.maps); changed[field] = value
                with self.assertRaises(ArithmeticError):
                    self.fixture.verify(maps=changed)

    def test_changed_generated_header_rejected_even_if_source_copy_matches(self):
        for directory in (self.fixture.local, self.fixture.archive):
            (directory/'GfniT64.h').write_text('// different circuit\n')
        self.assertRejected()

    def test_changed_archived_or_local_source_rejected(self):
        for directory in (self.fixture.local, self.fixture.archive):
            path = directory/'gfni_t64_probe.cpp'; old = path.read_text()
            path.write_text(old+'// changed\n')
            with self.subTest(directory=directory):
                self.assertRejected()
            path.write_text(old)

    def test_changed_binary_rejected(self):
        (self.fixture.archive/'gfni-t64').write_bytes(b'changed binary')
        self.assertRejected()

    def test_missing_build_hash_rejected(self):
        rows = self.fixture.hashes.read_text().splitlines()
        self.fixture.hashes.write_text('\n'.join(row for row in rows if not row.endswith('/gfni-t64'))+'\n')
        self.assertRejected()

    def test_duplicate_build_hash_rejected(self):
        self.fixture.hashes.write_text(self.fixture.hashes.read_text()+
            binding.digest(self.fixture.archive/'gfni-t64')+'  /other/gfni-t64\n')
        self.assertRejected()

    def test_missing_check_or_check_scope_rejected(self):
        path = self.fixture.check/'full-20-1.txt'; original = path.read_text()
        for token in ('checks PASS:', 'inner+FULL adjoint', 'complete output/suffix',
                      'step=64; state=16; refresh=uniform-gl16',
                      'map_sha256='+CPP_HASH, 'physical_epochs=32768'):
            with self.subTest(token=token):
                path.write_text(original.replace(token, 'omitted'))
                self.assertRejected()
        path.unlink()
        self.assertRejected()

    def test_missing_timing_or_correctness_line_rejected(self):
        path = self.fixture.confirm/'1-0-t64.csv'; original = path.read_text()
        path.write_text(original.replace('checks PASS:', 'checks missing:'))
        self.assertRejected()
        path.unlink()
        self.assertRejected()

    def test_wrong_seed_call_count_or_refresh_rejected(self):
        path = self.fixture.confirm/'1-0-t64.csv'; original = path.read_text()
        for before, after in (('uniform-gl16,1,101,', 'uniform-gl16,17,101,'),
                              ('uniform-gl16,1,101,', 'uniform-gl16,1,31,'),
                              ('20,64,16,uniform-gl16,', '20,64,16,transvections,')):
            with self.subTest(after=after):
                path.write_text(original.replace(before, after))
                self.assertRejected()

    def test_wrong_timed_map_or_epoch_scope_rejected(self):
        path = self.fixture.confirm/'1-0-t64.csv'; original = path.read_text()
        for before, after in ((CPP_HASH, 'e'*64), ('physical_epochs=32768', 'physical_epochs=16384')):
            with self.subTest(after=after):
                path.write_text(original.replace(before, after))
                self.assertRejected()

    def test_duplicate_timing_row_rejected(self):
        path = self.fixture.confirm/'1-0-t64.csv'; text = path.read_text()
        path.write_text(text+text.splitlines()[-1]+'\n')
        self.assertRejected()

    def test_extra_timing_file_rejected(self):
        (self.fixture.confirm/'1-2-t64.csv').write_text((self.fixture.confirm/'1-0-t64.csv').read_text())
        self.assertRejected()

    def test_nonfinite_decimal_timing_rejected(self):
        path = self.fixture.confirm/'1-0-t64.csv'
        path.write_text(path.read_text().replace(',101,6.0,', ',101,'+'9'*400+','))
        self.assertRejected()

    def test_wrong_csv_update_count_or_tile_mode_rejected(self):
        for kind, before, after in (('control', ',4,1\n', ',8,1\n'),
                                   ('control', ',4,1\n', ',4,2\n'),
                                   ('t64', ',deadbeefdeadbeef,1\n', ',deadbeefdeadbeef,2\n')):
            with self.subTest(kind=kind, after=after):
                path = self.fixture.confirm/f'1-0-{kind}.csv'; original = path.read_text()
                path.write_text(original.replace(before, after))
                self.assertRejected()
                path.write_text(original)

    def test_missing_or_extra_csv_suffix_fields_rejected(self):
        path = self.fixture.confirm/'1-0-t64.csv'; original = path.read_text()
        for replacement in (',deadbeefdeadbeef\n', ',deadbeefdeadbeef,1,extra\n'):
            with self.subTest(replacement=replacement):
                path.write_text(original.replace(',deadbeefdeadbeef,1\n', replacement))
                self.assertRejected()

    def test_checksum_must_match_repetitions_for_same_seed_and_candidate(self):
        for kind in ('control', 't64'):
            with self.subTest(kind=kind):
                path = self.fixture.confirm/f'1-1-{kind}.csv'; original = path.read_text()
                path.write_text(original.replace('deadbeefdeadbeef', 'cafebabecafebabe'))
                self.assertRejected()
                path.write_text(original)

    def test_nonhexadecimal_checksum_rejected(self):
        path = self.fixture.confirm/'1-0-t64.csv'
        path.write_text(path.read_text().replace('deadbeefdeadbeef', 'not-a-checksum'))
        self.assertRejected()

    def test_multiple_build_or_confirmation_batches_rejected(self):
        extra = self.fixture.archive/'measurements/t64-confirm-extra'; extra.mkdir()
        self.assertRejected()
        extra.rmdir()
        extra = self.fixture.archive/'measurements/t64-build-extra'; extra.mkdir()
        (extra/'hashes.txt').write_text(self.fixture.hashes.read_text())
        self.assertRejected()


ARCHIVE = binding.PERMUTATION.parents[2]/'tmp/packed-hill/gfni-t64-performance'


class ActualArchiveTests(unittest.TestCase):
    @unittest.skipUnless(ARCHIVE.is_dir(), 'optional local research performance archive absent')
    def test_current_generator_and_fresh_proof_maps_bind_actual_archive(self):
        result = binding.verify(ARCHIVE)
        self.assertAlmostEqual(result['median_of_process_medians_ms']['control'], 5.3241165)
        self.assertAlmostEqual(result['median_of_process_medians_ms']['t64'], 5.3445045)
        self.assertAlmostEqual(result['relative_latency'], 1.0038293677458032)
        self.assertFalse(result['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()
