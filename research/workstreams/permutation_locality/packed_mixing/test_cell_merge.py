"""Synthetic provenance/geometry tests; no numerical proof jobs."""
import copy
from fractions import Fraction as Q
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import cell_merge as merge


def fixture():
    root = ['1/1000', '1']
    source = dict(schema=merge.dense_cover.SCHEMA, ensemble=merge.dense_cover.ENSEMBLE,
        K=1 << 20, N=1 << 21, updates=2, block_width=8, distance='19/200', threshold=199229,
        minimum_groups=33, maximum_groups=2048, root=root,
        comparison='direct-expected-shell-majorant', expected_cdf_sha256='a'*64,
        comparison_caps_sha256='b'*64, last_lp=104, refined=True,
        base_tilt='3/16', variance_bins=16, regional_count=True,
        outer_premises=dict(schema='synthetic-test-only'), mixture=[dict(mass='100', activity='1/2')],
        precision=256, target_bits=36,
        cover=dict(leaves={'00': dict(witness={}, proposal=-50, old_marker='preserve')},
            unresolved={'01': {}, '1': {}}, visited=7),
        fresh_replay=dict(complete_dense=True), verified=True, upper=[1, -999])
    provenance = dict(path=str(Path('synthetic-base.json').resolve()), sha256='c'*64)
    fragments = []
    for index, (assigned, accepted, unresolved) in enumerate((('01', ['010'], ['011']), ('1', ['1'], []))):
        cover = {name: {path: dict(cell=list(map(str, merge.geometry.path_cell(root, path))),
                    **(dict(witness={}, proposal='NOT A BOUND', upper='NOT A BOUND') if name == 'leaves' else {}))
                for path in paths} for name, paths in (('leaves', accepted), ('unresolved', unresolved))}
        cover['visited'] = index+1
        fragment = {key: copy.deepcopy(source[key]) for key in merge.SCOPE_FIELDS}
        fragment.update(schema=merge.PARTIAL_SCHEMA, source=copy.deepcopy(provenance),
            assigned_roots=[assigned], precision=256, target_bits=36, cover=cover)
        metadata = dict(path=str(Path(f'synthetic-partial-{index}.json').resolve()), sha256=str(index)*64)
        fragments.append((fragment, metadata))
    return source, provenance, fragments


class CellMergeTests(unittest.TestCase):
    def setUp(self):
        self.source, self.provenance, self.partials = fixture()

    def merged(self, partials=None):
        return merge.merge_records(self.source, self.provenance, self.partials if partials is None else partials)

    def test_exact_frontier_preserves_old_leaves_and_ignores_worker_scores(self):
        before = copy.deepcopy((self.source, self.provenance, self.partials))
        result = self.merged()
        self.assertEqual(result['cover']['leaves']['00'], self.source['cover']['leaves']['00'])
        self.assertEqual(set(result['cover']['leaves']), {'00', '010', '1'})
        self.assertEqual(set(result['cover']['unresolved']), {'011'})
        self.assertEqual(set(result['cover']['leaves']['010']), {'cell', 'witness'})
        self.assertEqual(result['cover']['visited'], 10)
        self.assertEqual(before, (self.source, self.provenance, self.partials))
        self.assertTrue(result['final_replay_required'])
        for key in ('fresh_replay', 'verified', 'complete', 'upper'):
            self.assertNotIn(key, result)
        merge.dense_cover.validate_record(result)
        merge.geometry.partition(result['root'], result['cover']['leaves'], result['cover']['unresolved'])

    def test_partial_with_only_unresolved_children_is_valid(self):
        part = self.partials[0][0]
        part['cover']['unresolved'].update(part['cover']['leaves'])
        part['cover']['leaves'] = {}
        result = self.merged()
        self.assertEqual(set(result['cover']['unresolved']), {'010', '011'})
        self.assertEqual(set(result['cover']['unresolved']['010']), {'cell'})

    def test_full_search_coverage_still_has_no_certificate_claim(self):
        part = self.partials[0][0]
        part['cover']['leaves']['011'] = dict(part['cover']['unresolved'].pop('011'), witness={})
        result = self.merged()
        self.assertEqual(result['cover']['unresolved'], {})
        self.assertTrue(result['final_replay_required'])
        self.assertNotIn('fresh_replay', result)

    def test_missing_duplicate_and_unknown_assignments_rejected(self):
        for parts in ([], self.partials[:1], self.partials+[self.partials[0]]):
            with self.subTest(parts=len(parts)), self.assertRaises(ValueError):
                self.merged(parts)
        for roots in ([], ['01', '01'], ['00'], ['0'], ['010'], ['01', '1']):
            parts = copy.deepcopy(self.partials)
            parts[0][0]['assigned_roots'] = roots
            with self.subTest(roots=roots), self.assertRaises(ValueError):
                self.merged(parts)

    def test_scope_and_source_tampering_rejected(self):
        changes = dict(ensemble='other', distance='1/10', threshold=199230,
            root=['0', '1'], minimum_groups=32, expected_cdf_sha256='d'*64,
            comparison_caps_sha256='e'*64, outer_premises={}, mixture=[dict(mass='99', activity='1/2')],
            base_tilt='1/4', precision=384, target_bits=35, updates=True, schema='other')
        for key, value in changes.items():
            parts = copy.deepcopy(self.partials)
            parts[0][0][key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.merged(parts)
        for field, value in (('path', str(Path('other.json').resolve())), ('sha256', '0'*64)):
            parts = copy.deepcopy(self.partials)
            parts[0][0]['source'][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.merged(parts)

    def test_gaps_overlaps_cross_assignment_and_bad_coordinates_rejected(self):
        def missing(cover):
            cover['unresolved'].clear()
        def overlap(cover):
            cover['leaves']['01'] = dict(cell=list(map(str, merge.geometry.path_cell(self.source['root'], '01'))), witness={})
        def escape(cover):
            cover['leaves']['00'] = dict(cell=list(map(str, merge.geometry.path_cell(self.source['root'], '00'))), witness={})
        def stale(cover):
            cover['leaves']['010']['cell'] = ['0', '1']
        def duplicate(cover):
            cover['unresolved']['010'] = copy.deepcopy(cover['leaves']['010'])
        def missing_witness(cover):
            del cover['leaves']['010']['witness']
        for mutate in (missing, overlap, escape, stale, duplicate, missing_witness):
            parts = copy.deepcopy(self.partials)
            mutate(parts[0][0]['cover'])
            with self.subTest(mutation=mutate.__name__), self.assertRaises(ValueError):
                self.merged(parts)

    def test_cross_fragment_duplicate_roots_rejected_even_with_new_provenance(self):
        extra = copy.deepcopy(self.partials[0])
        extra[1]['path'] = str(Path('different-partial.json').resolve())
        with self.assertRaises(ValueError):
            self.merged([*self.partials, extra])

    def test_multiple_assigned_roots_in_one_fragment(self):
        first = copy.deepcopy(self.partials[0])
        second = self.partials[1][0]
        first[0]['assigned_roots'] += second['assigned_roots']
        first[0]['cover']['leaves'].update(second['cover']['leaves'])
        result = self.merged([first])
        self.assertEqual(set(result['cover']['leaves']), {'00', '010', '1'})

    def test_bad_source_and_noninteger_work_counts_fail(self):
        for change in (dict(root=['0', '0']), dict(mixture=[dict(mass='1', activity=.5)]),
                dict(precision=True), dict(expected_cdf_sha256='bad')):
            saved = copy.deepcopy(self.source)
            saved.update(change)
            with self.subTest(change=change), self.assertRaises(ValueError):
                merge.merge_records(saved, self.provenance, self.partials)
        for visits in (-1, True, 1.0):
            parts = copy.deepcopy(self.partials)
            parts[0][0]['cover']['visited'] = visits
            with self.assertRaises(ValueError):
                self.merged(parts)

    def test_filesystem_provenance_no_overwrite_and_no_output_on_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source_path, output = directory/'source.json', directory/'merged.json'
            source_path.write_text(json.dumps(self.source), encoding='utf-8')
            _, provenance = merge.read_source(source_path)
            partial_paths = []
            for index, (record, _) in enumerate(self.partials):
                record = copy.deepcopy(record)
                record['source'] = provenance
                path = directory/f'partial{index}.json'
                path.write_text(json.dumps(record), encoding='utf-8')
                partial_paths.append(path)
            with self.assertRaises(ValueError):
                merge.run(source_path, partial_paths[:1], output)
            self.assertFalse(output.exists())
            result = merge.run(source_path, partial_paths, output)
            self.assertEqual(json.loads(output.read_text()), result)
            self.assertEqual(result['source'], provenance)
            self.assertEqual(len(result['partial_sources']), 2)
            original = output.read_bytes()
            with self.assertRaises(ValueError):
                merge.run(source_path, partial_paths, output)
            self.assertEqual(output.read_bytes(), original)
            with self.assertRaises(ValueError):
                merge.run(source_path, [*partial_paths, partial_paths[0]], directory/'duplicate.json')

    def test_source_mutation_during_merge_fails_before_output(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source_path, output = directory/'source.json', directory/'merged.json'
            source_path.write_text(json.dumps(self.source), encoding='utf-8')
            _, provenance = merge.read_source(source_path)
            paths = []
            for index, (record, _) in enumerate(self.partials):
                record = copy.deepcopy(record)
                record['source'] = provenance
                path = directory/f'p{index}.json'
                path.write_text(json.dumps(record), encoding='utf-8')
                paths.append(path)
            original = merge.merge_records
            def changed(*args):
                result = original(*args)
                source_path.write_text(source_path.read_text()+'\n', encoding='utf-8')
                return result
            with patch.object(merge, 'merge_records', side_effect=changed), self.assertRaises(ValueError):
                merge.run(source_path, paths, output)
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
