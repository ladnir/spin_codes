import copy
from fractions import Fraction as Q
import io
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from contextlib import redirect_stderr, redirect_stdout
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from flint import arb, ctx
import cell_search as search


def source():
    return dict(schema=search.dense.SCHEMA, ensemble=search.dense.ENSEMBLE,
        K=1 << 20, N=1 << 21, updates=2, block_width=8, distance='19/200',
        threshold=199229, minimum_groups=33, maximum_groups=2048,
        root=['0', '1'], comparison='expected-cdf-shell-majorant',
        expected_cdf_sha256='cdf', comparison_caps_sha256='cdf', last_lp=104,
        refined=True, base_tilt='3/16', variance_bins=16, regional_count=True,
        outer_premises={}, mixture=[dict(mass='1', activity='1/2')],
        precision=256, target_bits=36,
        cover=dict(leaves={'00': dict(witness={'old': True}, proposal=-100)},
            unresolved={p: dict(cell=list(map(str, search.geometry.path_cell((0, 1), p))))
                for p in ('01', '100', '101', '11')}, visited=300))


class ToyModel:
    root = (Q(0), Q(1))

    def __init__(self, accept_width=Q(1), value=None):
        self.accept_width = accept_width
        self.value = value
        self.proposed, self.checked = [], []

    def proposal(self, cell):
        self.proposed.append(cell)
        return (-80 if cell[1]-cell[0] <= self.accept_width else 0), dict(parameters=['1', '0', '0'])

    def outward(self, cell, witness):
        self.checked.append(cell)
        return arb(2)**-60 if self.value is None else self.value


class CellSearchTests(unittest.TestCase):
    def setUp(self):
        previous = ctx.prec
        self.addCleanup(setattr, ctx, 'prec', previous)

    def test_cold_authentication_keeps_legacy_dense_cover_namespace_free(self):
        # Run the real cold import chain that formerly failed in
        # screen_dense's class definition. Stop at outer_caps(), before
        # expensive certificate arithmetic or the actual inner census.
        script = """
import builtins
import json
from pathlib import Path
import sys
from unittest.mock import patch
sys.path.insert(0, sys.argv[1])
before = sys.path[:]
import cell_search
assert sys.path == before
assert cell_search.dense.__name__ == 'packed_gl32_dense_interface'
assert 'dense_cover' not in sys.modules
class ReachedAuthentication(Exception):
    pass
def stop_before_arithmetic():
    legacy = sys.modules['dense_cover']
    assert callable(legacy.Model)
    assert Path(legacy.__file__).resolve() != Path(cell_search.dense.__file__).resolve()
    raise ReachedAuthentication
original = builtins.__import__
def importing(name, *args, **kwargs):
    module = original(name, *args, **kwargs)
    if name == 'verify_progress':
        module.outer_caps = stop_before_arithmetic
    return module
with patch('builtins.__import__', side_effect=importing):
    try:
        cell_search.fresh_model(json.loads(sys.argv[2]), 256, 36)
    except ReachedAuthentication:
        print('actual BCH authentication import chain reached')
    else:
        raise AssertionError('authentication entry point was not reached')
"""
        result = subprocess.run([sys.executable, '-B', '-c', script,
            str(Path(__file__).resolve().parent), json.dumps(source())],
            capture_output=True, text=True, timeout=45)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertIn('actual BCH authentication import chain reached', result.stdout)

    def test_contiguous_shards_partition_only_pending_frontier(self):
        original = source()
        shards = [search.assigned_paths(original, shard_index=i, shard_count=3) for i in range(3)]
        self.assertEqual(shards, [['01'], ['100'], ['101', '11']])
        self.assertEqual(sum(shards, []), sorted(original['cover']['unresolved']))
        for args in (dict(paths=['00']), dict(paths=['01', '01']), dict(paths=['0']),
                dict(paths=['01'], shard_index=0, shard_count=1),
                dict(shard_index=True, shard_count=2), dict(shard_index=3, shard_count=3)):
            with self.assertRaises(ValueError):
                search.assigned_paths(original, **args)

    def test_assigned_only_global_root_and_source_unchanged(self):
        original = source()
        before = copy.deepcopy(original)
        assigned = search.assigned_paths(original, paths=['101', '01'])
        model = ToyModel()
        result = search.search(model, assigned, max_cells=10)
        self.assertEqual(set(result['leaves']), {'01', '101'})
        self.assertEqual(result['visited'], 2)
        self.assertFalse(result['unresolved'])
        self.assertEqual(model.root, (0, 1))
        self.assertEqual(model.proposed, [(Q(1, 4), Q(1, 2)), (Q(5, 8), Q(3, 4))])
        self.assertEqual(original, before)
        search.validate_cover(model.root, assigned, result)

    def test_subdivision_and_budget_preserve_every_assigned_interval(self):
        model = ToyModel(accept_width=Q(1, 16))
        checkpoints = []
        result = search.search(model, ['01', '11'], max_cells=5, checkpoint_every=2,
            checkpoint=checkpoints.append)
        self.assertEqual(result['visited'], 5)
        self.assertEqual([row['visited'] for row in checkpoints], [0, 2, 4, 5])
        for state in checkpoints:
            search.validate_cover(model.root, ['01', '11'], state)
        self.assertTrue(result['unresolved'])

    def test_depth_cap_leaves_failed_paths_unresolved(self):
        result = search.search(ToyModel(accept_width=Q(0)), ['01'], max_cells=20, max_depth=3)
        self.assertEqual(set(result['unresolved']), {'010', '011'})
        self.assertEqual(result['visited'], 3)
        self.assertFalse(result['leaves'])

    def test_exact_gate_not_proposal_accepts_and_precision_checked(self):
        result = search.search(ToyModel(value=arb(1)), ['01'], max_cells=1, max_depth=2)
        self.assertFalse(result['leaves'])
        self.assertEqual(set(result['unresolved']), {'01'})
        for value in (arb(0), arb('nan'), arb('inf')):
            with self.assertRaises(ArithmeticError):
                search.search(ToyModel(value=value), ['01'], max_cells=1)
        class Mutated(ToyModel):
            def outward(self, *args):
                ctx.prec = 128
                return arb(2)**-80
        with self.assertRaises(ArithmeticError):
            search.search(Mutated(), ['01'], max_cells=1)

    def test_partial_coverage_rejects_gap_overlap_wrong_cell_and_outside(self):
        good = search.search(ToyModel(), ['01'], max_cells=1)
        for change in ('gap', 'wrong', 'outside', 'duplicate'):
            bad = copy.deepcopy(good)
            if change == 'gap':
                bad['leaves'].clear()
            elif change == 'wrong':
                bad['leaves']['01']['cell'] = ['0', '1']
            elif change == 'outside':
                bad['leaves']['11'] = bad['leaves'].pop('01')
            else:
                bad['unresolved']['01'] = dict(cell=['1/4', '1/2'])
            with self.assertRaises(ValueError):
                search.validate_cover((0, 1), ['01'], bad)

    def test_invalid_work_limits_rejected(self):
        for changed in (dict(max_cells=-1), dict(max_cells=True), dict(precision=127),
                dict(target_bits=19), dict(max_depth=1), dict(checkpoint_every=0)):
            with self.assertRaises(ValueError):
                search.search(ToyModel(), ['01'], **changed)

    def test_fresh_premise_or_comparison_mismatch_stops_before_model(self):
        import canonical_counts
        original_path = sys.path[:]
        self.addCleanup(lambda: sys.path.__setitem__(slice(None), original_path))
        with patch.object(canonical_counts, 'authenticated_bch_cdf', return_value=([0, 1], {})), \
                patch.object(search.dense, 'fingerprint', return_value='different'), \
                patch.object(search.dense, 'exact_mixture') as mixture:
            with self.assertRaises(ValueError):
                search.fresh_model(source(), 256, 36)
            mixture.assert_not_called()

    def test_cli_receipt_has_source_scope_and_partial_status(self):
        with TemporaryDirectory() as directory:
            snapshot, output = Path(directory)/'source.json', Path(directory)/'partial.json'
            snapshot.write_text(json.dumps(source()))
            argv = ['cell_search', str(snapshot), '--paths', '01', '--output', str(output), '--max-cells', '1']
            with patch.object(sys, 'argv', argv), patch.object(search, 'fresh_model', return_value=ToyModel()), \
                    redirect_stdout(io.StringIO()):
                search.main()
            receipt = json.loads(output.read_text())
            self.assertEqual(receipt['schema'], search.SCHEMA)
            self.assertEqual(receipt['assigned_roots'], ['01'])
            self.assertTrue(receipt['complete_assigned'])
            self.assertNotIn('fresh_replay', receipt)
            self.assertEqual(receipt['source']['sha256'], search.hashlib.sha256(snapshot.read_bytes()).hexdigest())
            self.assertEqual(receipt['source']['path'], str(snapshot.resolve()))
            self.assertEqual(receipt['root'], ['0', '1'])
            self.assertEqual(receipt['cover']['visited'], 1)
            with patch.object(sys, 'argv', argv), patch.object(search, 'fresh_model') as build, \
                    redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                search.main()
            build.assert_not_called()

    def test_changed_source_aborts_checkpoint(self):
        with TemporaryDirectory() as directory:
            snapshot, output = Path(directory)/'source.json', Path(directory)/'partial.json'
            snapshot.write_text(json.dumps(source()))
            class ChangesSource(ToyModel):
                def proposal(self, cell):
                    snapshot.write_text('{}')
                    return super().proposal(cell)
            argv = ['cell_search', str(snapshot), '--paths', '01', '--output', str(output)]
            with patch.object(sys, 'argv', argv), \
                    patch.object(search, 'fresh_model', return_value=ChangesSource()), \
                    redirect_stdout(io.StringIO()), self.assertRaises(ValueError):
                search.main()
            receipt = json.loads(output.read_text())
            self.assertFalse(receipt['complete_assigned'])
            self.assertEqual(receipt['cover']['visited'], 0)


if __name__ == '__main__':
    unittest.main()
