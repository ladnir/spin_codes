"""Lightweight schema/provenance tests; synthetic bounds are not certificates."""
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr
from io import StringIO
import runpy
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from hd_shifted_cover_backend import CONSTRUCTION, profile_count
from hd_shifted_cover_ledger import terminal_receipt
from hd_shifted_outer_density import DEFAULT_OUTER_MEASURE
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION
import hd_shifted_verified_cover as check


def proposal(first=1, last=2):
    return dict(construction=SWAP_CONSTRUCTION, requested_weight_range=[first, last],
        fixed_witness=dict(construction=SWAP_CONSTRUCTION, pair_swaps=True,
            orientation_aggregation='independent-fair-swaps', grid=64,
            source_exact_witness=dict(z='9/10', alpha='1/10', potential='1/8', markers=['1/4 1/16'])))


def receipt(entry, source, precision):
    first, last = entry['weight_range']
    ranges = ((first, last), (0, 4096), (0, 4096))
    return dict(construction=SWAP_CONSTRUCTION, field_order=256, forbidden_weight=734003,
        directed=True, whole_code_certificate=False, covers_all_even_odd_splits=True,
        ranges=ranges, weight_range=[first, last], integer_profile_triples=profile_count(ranges),
        precision=precision, pair_swaps=True, source_dimension=1048576,
        output_length=2097152, outer_groups=8192, every_integer_weight_replayed=True,
        covers_all_supports_at_each_weight=True, support_scope='scalar-code active input groups',
        support_count_method='common-potential-run-convolution',
        coefficient_method='fixed-normalized-directed-alias-bank',
        proposal_sha256=entry['proposal']['sha256'], weights_replayed=last-first+1,
        witness=check.expected_potential_witness(source, entry, precision),
        results=[dict(active_groups=r, regional_histograms=r//2-max(0, r-4096)+1,
                      log2_upper_rational='-100') for r in range(first, last+1)],
        log2_upper_rational='-999999')


class VerifiedCoverTests(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.base = Path(self.directory.name)

    def tearDown(self):
        self.directory.cleanup()

    def save(self, name, value):
        path = self.base/name
        path.write_text(json.dumps(value, separators=(',', ':')), encoding='utf-8')
        return dict(path=name, sha256=sha256(path.read_bytes()).hexdigest())

    def make(self, first=1, last=2):
        source = proposal(first, last)
        entry = dict(id='potential', kind='potential-band', weight_range=[first, last],
                     observed_fresh_dual_replay=True, outer_measure=DEFAULT_OUTER_MEASURE,
                     proposal=self.save('proposal.json', source), receipts={})
        for precision in (256, 384):
            entry['receipts'][str(precision)] = self.save(f'{precision}.json', receipt(entry, source, precision))
        return dict(schema_version=1, construction=SWAP_CONSTRUCTION,
                    trust_basis=check.TRUST_BASIS, entries=[entry])

    def verify(self, manifest):
        ref = self.save('manifest.json', manifest)
        return check.verify_manifest(self.base/ref['path'], ref['sha256'])

    def alter(self, manifest, precision, callback):
        entry = manifest['entries'][0]
        ref = entry['receipts'][str(precision)]
        value = json.loads((self.base/ref['path']).read_text())
        callback(value)
        entry['receipts'][str(precision)] = self.save(ref['path'], value)

    def test_external_pin_and_receipt_byte_pin(self):
        manifest = self.make()
        ref = self.save('manifest.json', manifest)
        with self.assertRaisesRegex(ValueError, 'Byte hash mismatch'):
            check.verify_manifest(self.base/ref['path'], '0'*64)
        (self.base/'256.json').write_text('{}', encoding='utf-8')
        with self.assertRaisesRegex(ValueError, 'Byte hash mismatch'):
            check.verify_manifest(self.base/ref['path'], ref['sha256'])

    def test_forged_summary_is_unused_and_incomplete_cannot_certify(self):
        result = self.verify(self.make())
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['arithmetic_validity_from_hashes'])
        self.assertTrue(result['cached_arithmetic_requires_trusted_observed_replays'])
        self.assertGreaterEqual(F(result['summary']['accepted_union_log2_upper']), -99)
        self.assertLess(F(result['summary']['accepted_union_log2_upper']), F(-989, 10))
        self.assertEqual(result['missing_weights'][0], 3)

    def test_every_integer_weight_and_precision_agreement(self):
        for mutation in (lambda x: x['results'].pop(),
                         lambda x: x['results'][1].update(active_groups=1),
                         lambda x: x['results'][0].update(log2_upper_rational='-101')):
            manifest = self.make()
            self.alter(manifest, 384, mutation)
            with self.assertRaises(ValueError):
                self.verify(manifest)

    def test_exact_witness_mode_and_scalar_scope(self):
        for mutation in (lambda x: x['witness'].update(z='89/100'),
                         lambda x: x.update(outer_measure='invertible-toeplitz'),
                         lambda x: x.update(covers_all_supports_at_each_weight=False),
                         lambda x: x.update(source_dimension=524288),
                         lambda x: x.update(precision=256)):
            manifest = self.make()
            self.alter(manifest, 384, mutation)
            with self.assertRaises(ValueError):
                self.verify(manifest)

    def test_explicit_legacy_missing_proposal_hash_never_allows_wrong_hash(self):
        manifest = self.make()
        for precision in (256, 384):
            self.alter(manifest, precision, lambda x: x.pop('proposal_sha256'))
        with self.assertRaisesRegex(ValueError, 'Legacy receipt'):
            self.verify(manifest)
        manifest['entries'][0]['legacy_missing_proposal_pin'] = True
        self.assertFalse(self.verify(manifest)['whole_code_certificate'])
        self.alter(manifest, 384, lambda x: x.update(proposal_sha256='0'*64))
        with self.assertRaisesRegex(ValueError, 'different proposal'):
            self.verify(manifest)

    def test_overlap_rejected_before_loading_artifacts(self):
        manifest = self.make()
        second = deepcopy(manifest['entries'][0])
        second['id'] = 'overlap'
        manifest['entries'].append(second)
        with self.assertRaisesRegex(ValueError, 'overlap'):
            self.verify(manifest)

    def test_full_geometry_fixture_and_fresh_directed_sum(self):
        # These are deliberately synthetic trusted arithmetic fixtures. The
        # test proves accounting behavior, not any encoder distance statement.
        result = self.verify(self.make(1, 8192))
        self.assertTrue(result['whole_code_certificate'])
        self.assertEqual(result['missing_weights'], [])
        self.assertEqual(result['summary']['accepted_profile_triples'], 11465830401)
        self.assertGreaterEqual(F(result['summary']['accepted_union_log2_upper']), -87)
        self.assertLess(F(result['summary']['accepted_union_log2_upper']), F(-869, 10))

    def test_sparse_current_family_and_legacy_jsonl(self):
        from hd_shifted_swap_sparse import family_receipt
        entry = dict(id='sparse', kind='sparse-family', weight_range=[1, 2],
                     observed_fresh_dual_replay=True, receipts={})
        for precision in (256, 384):
            value = family_receipt(2, precision=precision)
            rows = value.pop('weight_receipts')
            path = self.base/f'sparse{precision}.jsonl'
            path.write_text('\n'.join(json.dumps(row) for row in rows+[value]), encoding='utf-8')
            entry['receipts'][str(precision)] = dict(path=path.name, format='sparse-jsonl',
                sha256=sha256(path.read_bytes()).hexdigest())
        manifest = dict(schema_version=1, construction=SWAP_CONSTRUCTION,
                        trust_basis=check.TRUST_BASIS, entries=[entry])
        result = self.verify(manifest)
        self.assertEqual(result['summary']['accepted_profile_triples'], 4)
        self.assertFalse(result['whole_code_certificate'])

    def test_terminal_is_recomputed_and_inherited_explicitly(self):
        entry = dict(id='terminal', kind='terminal', weight_range=[7649, 8192],
                     observed_fresh_dual_replay=True,
                     receipts={str(p): self.save(f'terminal{p}.json', terminal_receipt(p)) for p in (256, 384)})
        manifest = dict(schema_version=1, construction=SWAP_CONSTRUCTION,
                        trust_basis=check.TRUST_BASIS, entries=[entry])
        self.assertFalse(self.verify(manifest)['whole_code_certificate'])
        ref = entry['receipts']['384']
        value = json.loads((self.base/ref['path']).read_text())
        value['log2_upper_rational'] = '-999999'
        entry['receipts']['384'] = self.save(ref['path'], value)
        with self.assertRaisesRegex(ValueError, 'fresh computation'):
            self.verify(manifest)

    def plan_fixture(self, swapped=True):
        construction = SWAP_CONSTRUCTION if swapped else CONSTRUCTION
        ranges = ((1, 1), (0, 4096), (0, 4096))
        source = dict(construction=construction, requested_weight_range=[1, 1],
            profile_coordinates=['r', 'a', 'b'], all_even_odd_splits=True,
            accepted_boxes=[dict(ranges=ranges, witness_id=0)],
            witnesses=[dict(mark=.9, alpha=.1, marker_log_odds=[[0., 0.], [0., 0.]],
                            pair_swaps=True, orientation_aggregation='independent-fair-swaps')])
        entry = dict(id='plan', kind='swap-plan' if swapped else 'unsigned-plan', weight_range=[1, 1],
                     observed_fresh_dual_replay=True, proposal=self.save('plan.json', source), receipts={})
        for precision in (256, 384):
            leaf = dict(ranges=ranges, construction=construction, field_order=256, forbidden_weight=734003,
                directed=True, whole_code_certificate=False, covers_all_even_odd_splits=True,
                log2_upper_rational='-100', support_count_method='matching-profile-envelope', witness_index=0)
            value = dict(proposal_sha256=entry['proposal']['sha256'], accepted=[leaf], failed=[],
                         all_proposed_leaves_replayed=True,
                         witnesses=[check.expected_plan_witness(source, 0, precision, swapped)])
            entry['receipts'][str(precision)] = self.save(f'plan{precision}.json', value)
        return dict(schema_version=1, construction=SWAP_CONSTRUCTION,
                    trust_basis=check.TRUST_BASIS, entries=[entry])

    def test_plan_geometry_witnesses_and_unsigned_transfer(self):
        for swapped in (False, True):
            manifest = self.plan_fixture(swapped)
            self.assertFalse(self.verify(manifest)['whole_code_certificate'])
            self.alter(manifest, 384, lambda x: x['witnesses'][0].update(alpha='1/2'))
            with self.assertRaisesRegex(ValueError, 'witness identity'):
                self.verify(manifest)

    def test_swap_plan_must_have_explicit_construction(self):
        manifest = self.plan_fixture()
        entry = manifest['entries'][0]
        source = json.loads((self.base/'plan.json').read_text())
        source.pop('construction')
        entry['proposal'] = self.save('plan.json', source)
        with self.assertRaisesRegex(ValueError, 'plan identity'):
            self.verify(manifest)

    def test_summary_only_plan_and_unknown_inheritance_fail(self):
        manifest = self.plan_fixture(False)
        self.alter(manifest, 384, lambda x: x['accepted'][0].update(support_count_method='fixed-oriented'))
        with self.assertRaisesRegex(ValueError, 'Unestablished'):
            self.verify(manifest)

    def test_checkpoint_reconstructs_inputs_but_checks_every_source_leaf(self):
        manifest = self.plan_fixture()
        entry = manifest['entries'][0]
        entry.update(kind='checkpoint', checkpoint_kind='dense')
        entry['plan_bands'] = [dict(label='plan.json', kind='swap-plan', weight_range=[1, 1],
                                    proposal=entry.pop('proposal'))]
        source = Path(check.__file__).with_name('hd_shifted_cover_ledger.py')
        entry['replay_source_pin'] = dict(path=str(source), sha256=sha256(source.read_bytes()).hexdigest())
        for precision in (256, 384):
            ref = entry['receipts'][str(precision)]
            value = json.loads((self.base/ref['path']).read_text())
            value.update(kind='dense', precision=precision, construction=SWAP_CONSTRUCTION,
                         fresh_arithmetic=True, whole_code_certificate=False,
                         proposal_sha256={'plan.json': entry['plan_bands'][0]['proposal']['sha256']})
            value['accepted'][0]['source_label'] = 'plan.json'
            value.pop('witnesses')  # Explicitly exercise deterministic reconstruction.
            entry['receipts'][str(precision)] = self.save(ref['path'], value)
        result = self.verify(manifest)
        self.assertTrue(result['entries'][0]['witness_metadata_reconstructed'])
        self.alter(manifest, 384, lambda x: x['accepted'][0].update(source_label='another-plan.json'))
        with self.assertRaisesRegex(ValueError, 'source identity'):
            self.verify(manifest)
        # A valid byte pin on a same-basename different implementation does
        # not authorize the checker's hardcoded current replay conversion.
        alternate = self.base/'hd_shifted_cover_ledger.py'
        alternate.write_text('# different implementation\n', encoding='utf-8')
        entry['replay_source_pin'] = dict(path=alternate.name, sha256=sha256(alternate.read_bytes()).hexdigest())
        with self.assertRaisesRegex(ValueError, 'reviewed implementation'):
            self.verify(manifest)

    def test_checkpoint_requires_replay_source_pin(self):
        manifest = self.plan_fixture()
        entry = manifest['entries'][0]
        entry.update(kind='checkpoint', checkpoint_kind='dense', plan_bands=[])
        for precision in (256, 384):
            self.alter(manifest, precision, lambda x: x.update(kind='dense', precision=precision,
                construction=SWAP_CONSTRUCTION, fresh_arithmetic=True, whole_code_certificate=False))
        with self.assertRaises(KeyError):
            self.verify(manifest)
        manifest = self.plan_fixture()
        self.alter(manifest, 384, lambda x: x.pop('accepted'))
        with self.assertRaises(KeyError):
            self.verify(manifest)

    def test_strict_json_and_invalid_grid(self):
        for raw in ('{"x":1,"x":2}', '{"x":NaN}'):
            with self.assertRaises(ValueError):
                check.strict_json(raw)
        source = proposal()
        source['fixed_witness']['grid'] = 7
        entry = dict(weight_range=[1, 2], outer_measure=DEFAULT_OUTER_MEASURE)
        with self.assertRaisesRegex(ValueError, 'grid'):
            check.expected_potential_witness(source, entry, 256)
        manifest = self.make()
        manifest['target'] = dict(field_order=65536)
        with self.assertRaisesRegex(ValueError, 'Manifest target'):
            self.verify(manifest)

    def test_fresh_uses_existing_replay_not_cache_and_marks_incomplete(self):
        manifest = self.make()
        ref = self.save('manifest.json', manifest)
        class Empty:
            accepted, unresolved = [], []
        result = dict(summary=dict(whole_code_certificate=True))
        with patch('hd_shifted_swap_ledger.replay_swap_cover', return_value=(Empty(), result)) as replay:
            value = check.fresh_manifest(self.base/ref['path'], ref['sha256'])
        self.assertEqual(replay.call_args.kwargs['threshold'], -40)
        self.assertFalse(value['summary']['whole_code_certificate'])
        self.assertFalse(value['cached_bounds_used'])

    def test_cli_exclusive_receipt_and_preflight_overwrite_guard(self):
        manifest = self.make()
        ref = self.save('manifest.json', manifest)
        output = self.base/'verified.json'
        argv = ['hd_shifted_verified_cover.py', '--manifest', str(self.base/ref['path']),
                '--expected-manifest-sha256', ref['sha256'], '--output', str(output)]
        with patch('sys.argv', argv), redirect_stdout(StringIO()):
            runpy.run_module('hd_shifted_verified_cover', run_name='__main__')
        value = json.loads(output.read_text())
        self.assertFalse(value['whole_code_certificate'])
        original = output.read_bytes()
        argv[2] = str(self.base/'does-not-exist.json')
        with patch('sys.argv', argv), redirect_stderr(StringIO()), self.assertRaises(SystemExit) as caught:
            runpy.run_module('hd_shifted_verified_cover', run_name='__main__')
        self.assertNotEqual(caught.exception.code, 0)
        self.assertEqual(output.read_bytes(), original)


if __name__ == '__main__':
    unittest.main()
