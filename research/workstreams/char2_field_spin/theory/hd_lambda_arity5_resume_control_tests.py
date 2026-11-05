"""Exact ordering and mocked resume tests; no provider or proof control."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import hd_lambda_arity5_cover_control as base
import hd_lambda_arity5_resume_control as resume
from hd_pair_scheduled_outer import column_permutations


def pins():
    result = {name: name for name in resume.SOURCE_FILES}
    result['hd_lambda_arity5_cover_control.py'] = resume.BASE_SHA256
    return result


def row(source, prefix=0., tail=0., exact_tail=None, timeout=False):
    root = tail/18-base.RS_ROOT_CHARGE
    result = dict(original_source=list(source), original_mask=base.source_mask(source),
        marked_tail_log2=tail, marked_normalized_root_log2=root,
        normalized_fixed_A_root_log2=root, fixed_A_repeated_r_bound_log2=prefix+18*root,
        selected_method='marked-fibre', exact_fallback_status='not_needed')
    if exact_tail is not None:
        exact_root = exact_tail/18-base.RS_ROOT_CHARGE
        result.update(exact_tail_log2=exact_tail, exact_normalized_root_log2=exact_root,
            normalized_fixed_A_root_log2=min(root, exact_root), exact_fallback_status='complete',
            selected_method='exact-GF4' if exact_root < root else 'marked-fibre',
            fixed_A_repeated_r_bound_log2=prefix+18*min(root, exact_root))
    if timeout:
        result['exact_fallback_status'] = 'timeout_marked_retained'
    return result


def fixture(rows=None, prefix=0.):
    saved, measure = dict(mark='9927/10000', source_count=123), {'log2_atom00': -1.}
    rows = deepcopy(rows if rows is not None else [row(base.original_sources()[0], prefix)])
    attempts = []
    for record in rows:
        if record['exact_fallback_status'] not in ('complete', 'timeout_marked_retained'):
            continue
        attempt = dict(original_source=record['original_source'], original_mask=record['original_mask'],
            status=record['exact_fallback_status'],
            marked_fixed_A_bound_log2=prefix+18*record['marked_normalized_root_log2'])
        if attempt['status'] == 'complete':
            attempt.update({key: record[key] for key in (
                'exact_tail_log2', 'exact_normalized_root_log2', 'selected_method')})
            attempt['selected_fixed_A_bound_log2'] = record['fixed_A_repeated_r_bound_log2']
        attempts.append(attempt)
    base_pins = {name: pins()[name] for name in base.SOURCE_FILES}
    seed = dict(resume.SCOPE, mark=saved['mark'], source_count=saved['source_count'],
        common_measure=measure, common_prefix_log2=prefix, records=rows, exact_fallbacks=attempts,
        exact_fallback_limit=128, exact_fallback_cap_skipped=0,
        input_artifact=base.INPUT, input_sha256=base.INPUT_SHA,
        raw_artifact_sha256_after=base.INPUT_SHA, source_pins_unchanged=True,
        source_pins_before=base_pins, source_pins_after=base_pins,
        complete=False, target_met=False, coverage={'untrusted_cached_summary': True})
    return seed, saved, measure, prefix


class ResumeControlTests(unittest.TestCase):
    def test_full_orbits_and_scheduled_cache_identity(self):
        orbits = resume.translation_orbits()
        self.assertEqual(len(orbits), 273)
        flattened = tuple(source for orbit in orbits for source in orbit)
        self.assertEqual(len(flattened), 4368)
        self.assertEqual(set(flattened), set(base.original_sources()))
        for source in (orbits[0][0], orbits[101][0], orbits[-1][0]):
            mapped = set()
            for j in range(16):
                before, _ = column_permutations(resume.SCOPE['schedule'], j)
                mapped.add(tuple(sorted(before[i] for i in source)))
                for d in range(16):
                    shifted, _ = column_permutations(resume.SCOPE['schedule'], j ^ d)
                    self.assertEqual(sorted(before[i ^ d] for i in source),
                                     sorted(shifted[i] for i in source))
            self.assertEqual(len(mapped), 16)
        first = row(flattened[0])
        self.assertEqual(resume.remaining_sources([first]), flattened[1:])

    def test_seed_checks_scope_pins_root_algebra_and_attempt_ledger(self):
        source = base.original_sources()[0]
        seed, saved, measure, prefix = fixture([row(source, tail=50., exact_tail=10.)])
        validated, summary = resume.validate_seed(seed, pins(), saved, measure, prefix)
        self.assertEqual(len(validated['exact_fallbacks']), 1)
        self.assertEqual(summary['unique_original_A'], 1)
        self.assertFalse(summary['all_original_A_covered'])
        mutations = [
            lambda x: x.update(source_class='Fmax=1'),
            lambda x: x.update(common_prefix_log2=1.),
            lambda x: x['source_pins_before'].update({'hd_lambda_outer.py': 'wrong'}),
            lambda x: x['records'].append(deepcopy(x['records'][0])),
            lambda x: x['records'][0].update(marked_tail_log2=51.),
            lambda x: x['records'][0].update(normalized_fixed_A_root_log2=-999.),
            lambda x: x.update(exact_fallbacks=[]),
            lambda x: x['exact_fallbacks'].append(deepcopy(x['exact_fallbacks'][0])),
        ]
        for mutate in mutations:
            bad = deepcopy(seed)
            mutate(bad)
            with self.assertRaises((ValueError, KeyError)):
                resume.validate_seed(bad, pins(), saved, measure, prefix)

    def test_zero_logs_and_timed_out_attempt_remain_charged(self):
        sources = base.original_sources()
        seed, saved, measure, prefix = fixture([
            row(sources[0], exact_tail=float('-inf')), row(sources[1], timeout=True)])
        serialized = json.loads(json.dumps(base._jsonable(seed), allow_nan=False))
        validated, summary = resume.validate_seed(serialized, pins(), saved, measure, prefix)
        self.assertEqual(validated['records'][0]['normalized_fixed_A_root_log2'], float('-inf'))
        self.assertEqual(len(validated['exact_fallbacks']), 2)
        self.assertEqual(summary['unique_original_A'], 2)

    def run_mocked(self, seed, saved, measure, prefix, processor, *, pins_after=None):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'seed.json'
            data = json.dumps(base._jsonable(seed), allow_nan=False).encode()
            path.write_bytes(data)
            digest = sha256(data).hexdigest()
            with patch.object(resume, '_pins', side_effect=[pins(), pins_after or pins()]), \
                 patch.object(base, '_input_hash', return_value=base.INPUT_SHA), \
                 patch.object(base, 'load_measure', return_value=(saved, measure, prefix)), \
                 patch.object(base, 'moment_table', return_value=None), \
                 patch.object(base, '_process_sources', side_effect=processor) as process:
                result = resume.run(path, seed_sha256=digest)
            self.assertEqual(path.read_bytes(), data)
            return result, process

    def test_resume_filters_seed_and_charges_cumulative_attempts(self):
        sources = base.original_sources()
        seed, saved, measure, prefix = fixture([row(source, tail=50., exact_tail=10.)
                                                for source in sources[:47]])
        def process(remaining, table, table2, actual_prefix, **kwargs):
            self.assertEqual(actual_prefix, prefix)
            self.assertEqual(kwargs['fallback_limit'], 81)
            self.assertEqual(remaining, resume.remaining_sources(seed['records']))
            self.assertEqual(len(remaining), 4321)
            return dict(records=[row(source) for source in remaining], exact_fallbacks=[],
                exact_fallback_cap_skipped=0, pending_preparations_without_tail=0, stop_reason=None)
        result, _ = self.run_mocked(seed, saved, measure, prefix, process)
        self.assertTrue(result['complete'])
        self.assertEqual(result['coverage']['unique_original_A'], 4368)
        self.assertEqual(result['seed_exact_attempts'], 47)
        self.assertEqual(len(result['exact_fallbacks']), 47)
        self.assertFalse(result['orbit_quotient_used'])
        self.assertEqual(result['records'][:47], seed['records'])
        self.assertTrue(result['source_pins_unchanged'])

    def test_exhausted_fallback_cap_does_not_stop_remaining_marked_rows(self):
        sources = base.original_sources()
        seed, saved, measure, prefix = fixture([row(source, timeout=True) for source in sources[:128]])
        def process(remaining, *unused, **kwargs):
            self.assertEqual(kwargs['fallback_limit'], 0)
            rows = [row(source) for source in remaining]
            for record in rows:
                record['exact_fallback_status'] = 'attempt_cap_skipped'
            return dict(records=rows, exact_fallbacks=[], exact_fallback_cap_skipped=len(rows),
                pending_preparations_without_tail=0, stop_reason=None)
        result, _ = self.run_mocked(seed, saved, measure, prefix, process)
        self.assertTrue(result['complete'])
        self.assertTrue(result['fallback_budget_exhausted'])
        self.assertEqual(len(result['exact_fallbacks']), 128)
        # A new pinned continuation can itself be validated as a later seed.
        validated, summary = resume.validate_seed(result, pins(), saved, measure, prefix)
        self.assertEqual(len(validated['records']), 4368)
        self.assertTrue(summary['all_original_A_covered'])

    def test_timeout_and_pin_change_withhold_full_claim(self):
        seed, saved, measure, prefix = fixture()
        def partial(remaining, *unused, **kwargs):
            return dict(records=[row(remaining[0])], exact_fallbacks=[],
                exact_fallback_cap_skipped=0, pending_preparations_without_tail=3,
                stop_reason=dict(type='TimeoutError', message='toy'))
        result, _ = self.run_mocked(seed, saved, measure, prefix, partial)
        self.assertFalse(result['complete'])
        self.assertEqual(len(result['records']), 2)
        self.assertIsNone(result['coverage']['full_selected_arity_bound_log2'])
        changed = pins()
        changed['hd_lambda_arity5_resume_control.py'] = 'changed'
        def full(remaining, *unused, **kwargs):
            return dict(records=[row(source) for source in remaining], exact_fallbacks=[],
                exact_fallback_cap_skipped=0, pending_preparations_without_tail=0, stop_reason=None)
        result, _ = self.run_mocked(seed, saved, measure, prefix, full, pins_after=changed)
        self.assertTrue(result['coverage']['all_original_A_covered'])
        self.assertFalse(result['complete'])
        self.assertFalse(result['source_pins_unchanged'])
        self.assertIsNone(result['coverage']['full_selected_arity_bound_log2'])

    def test_external_seed_hash_and_resource_guards(self):
        for kwargs in ({'batch_size': 64}, {'max_seconds': 361}, {'max_seconds': 0}):
            with self.assertRaises(ValueError):
                resume.run('unused', seed_sha256='a'*64, **kwargs)
        with self.assertRaises(ValueError):
            resume.run('unused', seed_sha256='unobserved')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'seed.json'
            path.write_text('{}', encoding='utf-8')
            with self.assertRaises(ValueError):
                resume.run(path, seed_sha256='0'*64)


if __name__ == '__main__':
    unittest.main()
