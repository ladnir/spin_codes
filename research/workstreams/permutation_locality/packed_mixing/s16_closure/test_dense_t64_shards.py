"""Exact partition/provenance tests with tiny models, not proof jobs."""
from copy import deepcopy
from fractions import Fraction as Q
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from flint import arb
import dense_t64_shards as shards


class TinyModel:
    root = (Q(0), Q(1))
    def proposal(self, cell):
        return -80., dict(parameters=['1/100', '0', '0'])
    def outward(self, cell, witness):
        return arb(2)**-70


def row(path, accepted=False):
    result = dict(cell=list(map(str, shards.geometry.path_cell((0, 1), path))))
    if accepted:
        result.update(witness=dict(parameters=['1/100', '0', '0']), upper=['poisoned old endpoint'])
    return result


class T64ShardTests(unittest.TestCase):
    def setUp(self):
        self.temporary = TemporaryDirectory(); self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.proposal = self.directory/'proposal.json'; self.proposal.write_text('{}')
        _, proposal_meta = shards.read(self.proposal)
        self.record = dict(schema=shards.dense.SCHEMA, precision=192, cell_target_bits=52,
            scope=dict(schema='packed-gl32-t64-s16-dense-context-1', root=['0', '1'],
                maps=dict(map='fixed-t64'), distance='1/10', variance_bins=16, birth_density='capped'),
            source=proposal_meta, whole_code_certificate=False,
            cover=dict(leaves={'00': row('00', True)},
                unresolved={p: row(p) for p in ('01', '10', '11')}, visited=7))
        self.source = self.directory/'snapshot.json'; self.source.write_text(json.dumps(self.record))
        _, self.source_meta = shards.read(self.source)

    def partial(self, index, *, split=False):
        roots = shards.assigned_paths(self.record, index, 2)
        leaves = {p: row(p, True) for p in roots}; unresolved = {}
        if split:
            p = roots[0]; del leaves[p]
            leaves[p+'0'] = row(p+'0', True); unresolved[p+'1'] = row(p+'1')
        record = dict(schema=shards.SCHEMA, snapshot_source=self.source_meta,
            source=self.record['source'], scope=deepcopy(self.record['scope']), precision=192,
            cell_target_bits=52, shard_index=index, shard_count=2, assigned_roots=roots,
            cover=dict(leaves=leaves, unresolved=unresolved, visited=3))
        path = self.directory/f'partial-{index}.json'; path.write_text(json.dumps(record))
        return shards.read(path)

    def test_interleaved_assignments_and_invalid_limits(self):
        self.assertEqual(shards.assigned_paths(self.record, 0, 2), ['01', '11'])
        self.assertEqual(shards.assigned_paths(self.record, 1, 2), ['10'])
        for index, count in ((-1, 2), (2, 2), (0, 0), (0, 4), (True, 2)):
            with self.assertRaises(ValueError):
                shards.assigned_paths(self.record, index, count)

    def test_worker_freshly_prepares_only_its_assigned_subtrees(self):
        output = self.directory/'result.json'
        with patch.object(shards.dense, 'fresh_model',
                return_value=(TinyModel(), deepcopy(self.record['scope']), self.record['source'])) as prepare:
            result = shards.run_shard(self.source, 0, 2, output, max_cells=5, checkpoint_every=1)
        self.assertEqual(set(result['cover']['leaves']), {'01', '11'})
        self.assertFalse(result['cover']['unresolved'])
        self.assertEqual(prepare.call_args.args, (str(self.proposal.resolve()),))
        self.assertEqual(prepare.call_args.kwargs['precision'], 192)
        self.assertEqual(result['snapshot_source'], self.source_meta)
        self.assertFalse(result['whole_code_certificate'])
        self.assertTrue(result['final_replay_required'])

    def test_fresh_map_or_source_mismatch_rejected(self):
        for change in ('maps', 'source'):
            scope, source = deepcopy(self.record['scope']), deepcopy(self.record['source'])
            if change == 'maps':
                scope['maps']['map'] = 'different'
            else:
                source['sha256'] = '0'*64
            with self.subTest(change=change), patch.object(shards.dense, 'fresh_model',
                    return_value=(TinyModel(), scope, source)):
                with self.assertRaises(ValueError):
                    shards.run_shard(self.source, 0, 2, self.directory/(change+'.json'))

    def test_snapshot_changed_during_search_rejected(self):
        original = self.source.read_text()
        class MutatingModel(TinyModel):
            def proposal(inner, cell):
                self.source.write_text(original+' ')
                return super().proposal(cell)
        with patch.object(shards.dense, 'fresh_model',
                return_value=(MutatingModel(), self.record['scope'], self.record['source'])):
            with self.assertRaises(ValueError):
                shards.run_shard(self.source, 0, 2, self.directory/'changed.json', checkpoint_every=1)

    def test_changed_proposal_source_rejected_before_fresh_model(self):
        self.proposal.write_text('{"changed":true}')
        with patch.object(shards.dense, 'fresh_model') as prepare:
            with self.assertRaises(ValueError):
                shards.run_shard(self.source, 0, 2, self.directory/'changed.json')
            prepare.assert_not_called()

    def test_merge_retains_only_search_witnesses_and_exact_partition(self):
        partials = [self.partial(0, split=True), self.partial(1)]
        merged = shards.merge_records(self.record, self.source_meta, partials)
        self.assertEqual(set(merged['cover']['leaves']), {'00', '010', '10', '11'})
        self.assertEqual(set(merged['cover']['unresolved']), {'011'})
        self.assertEqual(merged['source'], self.record['source'])
        self.assertEqual(merged['scope'], self.record['scope'])
        self.assertEqual(merged['cover']['visited'], 13)
        self.assertFalse(merged['complete_search_partition'])
        self.assertTrue(merged['final_replay_required'])
        self.assertNotIn('upper', merged['cover']['leaves']['00'])
        self.assertEqual(merged['cover']['leaves']['00']['witness'],
                         self.record['cover']['leaves']['00']['witness'])
        self.assertIn('upper', self.record['cover']['leaves']['00'])
        shards.dense.validate_cover_record(merged)

    def test_complete_merge_accepted_as_normal_dense_search(self):
        merged = shards.merge_records(self.record, self.source_meta, [self.partial(0), self.partial(1)])
        self.assertTrue(merged['complete_search_partition'])
        self.assertEqual(set(shards.dense.validate_cover_record(merged, complete=True)),
                         {'00', '01', '10', '11'})
        self.assertFalse(merged['whole_code_certificate'])

    def test_missing_or_duplicate_assignment_rejected(self):
        first, second = self.partial(0), self.partial(1)
        duplicate = deepcopy(first); duplicate[1]['path'] = str(self.directory/'duplicate.json')
        for partials in ([first], [first, first, second], [first, duplicate, second]):
            with self.assertRaises(ValueError):
                shards.merge_records(self.record, self.source_meta, partials)

    def test_stale_source_map_scope_or_target_rejected(self):
        for field, value in (('snapshot_source', dict(path=str(self.source), sha256='0'*64)),
                             ('source', dict(path=str(self.proposal), sha256='0'*64)),
                             ('precision', 256), ('cell_target_bits', 51)):
            first, second = self.partial(0), self.partial(1)
            first[0][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                shards.merge_records(self.record, self.source_meta, [first, second])
        first, second = self.partial(0), self.partial(1)
        first[0]['scope']['maps']['map'] = 'changed'
        with self.assertRaises(ValueError):
            shards.merge_records(self.record, self.source_meta, [first, second])

    def test_gaps_overlap_and_wrong_global_coordinates_rejected(self):
        for change in ('gap', 'overlap', 'coordinate', 'escape', 'witness', 'assignment'):
            first, second = self.partial(0), self.partial(1)
            cover = first[0]['cover']
            if change == 'gap':
                del cover['leaves']['01']
            elif change == 'overlap':
                cover['leaves']['010'] = row('010', True)
            elif change == 'coordinate':
                cover['leaves']['01']['cell'] = ['0', '1/4']
            elif change == 'escape':
                cover['leaves']['00'] = row('00', True)
            elif change == 'witness':
                del cover['leaves']['01']['witness']
            else:
                first[0]['assigned_roots'] = ['01']
            with self.subTest(change=change), self.assertRaises(ValueError):
                shards.merge_records(self.record, self.source_meta, [first, second])

    def test_changed_file_during_merge_rejected(self):
        parts = [self.partial(0), self.partial(1)]
        original = shards.merge_records
        def mutate(*args):
            result = original(*args)
            self.source.write_text(self.source.read_text()+' ')
            return result
        with patch.object(shards, 'merge_records', side_effect=mutate):
            with self.assertRaises(ValueError):
                shards.run_merge(self.source, [p[1]['path'] for p in parts], self.directory/'merged.json')

    def test_duplicate_json_keys_rejected(self):
        path = self.directory/'bad.json'; path.write_text('{"cover":{},"cover":{}}')
        with self.assertRaises(ValueError):
            shards.read(path)


if __name__ == '__main__':
    unittest.main()
