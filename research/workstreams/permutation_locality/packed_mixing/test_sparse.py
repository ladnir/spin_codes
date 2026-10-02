"""Small exact checks for the packed sparse proof boundary."""
from fractions import Fraction as Q
import importlib.util
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
spec = importlib.util.spec_from_file_location('packed_sparse', HERE/'sparse.py')
packed = importlib.util.module_from_spec(spec)
spec.loader.exec_module(packed)
import support_cover
from flint import arb, arb_mat, ctx


class SparseTests(unittest.TestCase):
    def setUp(self):
        ctx.prec = 256

    def test_fractional_mass_is_not_rounded_away(self):
        counts = [Q(0)]*5+[Q(1, 1 << 100)]*252
        weights = [arb(0)]*257
        weights[5] = arb(1)
        result = packed.fold_cdf(counts, weights)
        self.assertGreaterEqual(result, arb(2)**-100)
        self.assertLess(result, arb(2)**-99)

    def test_abel_fold_bounds_every_compatible_toy_shell(self):
        from itertools import product
        caps = [Q(0), Q(1, 3), Q(3, 4), Q(2)]
        weights = [arb(0), arb(4), arb(1), arb(3)]
        upper = packed.fold_cdf(caps, weights)
        for shell in product(range(5), repeat=3):
            values = [Q(x, 6) for x in shell]
            if all(sum(values[:u]) <= caps[u] for u in range(1, 4)):
                actual = sum((packed.aq(values[u-1])*weights[u] for u in range(1, 4)), arb(0))
                self.assertLessEqual(actual, upper)

    def test_domain_rejects_omitted_sub38_tail(self):
        counts = [Q(0)]*5+[Q(1, 1 << 100)]*252
        self.assertEqual(support_cover.domain(counts, 2, 5), ((5, 256),)*2)
        self.assertRaises(ValueError, support_cover.domain, counts, 2, 38)
        self.assertRaises(ValueError, support_cover.domain, counts, 2, 4)
        self.assertRaises(ValueError, support_cover.domain, counts, True, 5)

    def test_legacy_domain_is_unchanged_when_requested(self):
        counts = [0]*38+[1]*219
        self.assertEqual(support_cover.domain(counts, 3, 38), ((38, 256),)*3)

    def test_exact_rational_binomial_fold(self):
        counts = [Q(0)]*5+[Q(1, 1 << 100)]*252
        actual = support_cover.fold_arb(counts, 5, 5, arb(1)/2)
        from math import comb
        exact = Q(1 << 156, comb(256, 5))
        self.assertGreaterEqual(actual, packed.aq(exact))
        self.assertLess(actual, packed.aq(exact)*(1+arb(2)**-200))

    def test_partition_checks_tree_not_only_volume(self):
        old = support_cover.old
        tree = old.RetainedCover()
        root = ((5, 8),)*2
        item = lambda i, box, mult: (0., i, box, mult, ('.01', '1'), [1]*2)
        tree.add(item(1, root, 1))
        children = list(old.split(root, 1))
        selected = []
        for i, (box, mult) in enumerate(children, 2):
            value = item(i, box, mult)
            tree.add(value, 1)
            selected.append(value)
        self.assertTrue(support_cover.checked_partition(tree, selected, root))
        self.assertRaises(ValueError, support_cover.checked_partition, tree, selected[:-1], root)
        first = tree.nodes[2]['item']
        tree.nodes[2]['item'] = item(2, ((6, 7),)*2, first[3])
        self.assertRaises(ValueError, support_cover.checked_partition, tree, selected, root)

    def test_checked_counts_authenticates_total_and_first_support(self):
        counts = [Q(0)]*5+[Q(1, 1 << 100)]*251+[Q((1 << 512)-1)]
        self.assertEqual(packed.checked_counts(counts, 5)[1], 5)
        self.assertRaises(ValueError, packed.checked_counts, counts, 38)
        self.assertRaises(ValueError, packed.checked_counts, counts[:-1]+[Q(1)], 5)

    def test_saved_replay_checks_degree_and_geometry(self):
        counts = [Q(0)]*5+[Q(1, 7)]*252
        box = [[5, 256], [5, 256]]
        details = dict(support_min=5, domain_volume=252**2,
            tree=[dict(id=1, parent=None, children=[], box=box, multiplicity=1)], selected_ids=[1],
            leaves=[dict(box=box, multiplicity=1, tilt='.0001',
                         numerators=[support_cover.old.DENOMINATOR//2]*2)])
        exact = [arb_mat([[arb(1)/1024]]) for _ in range(3)]
        operators = {('.0001', '1'): (exact, None)}
        args = dict(groups=2, cutoff=0, precision=256, target_bits=32)
        upper = support_cover.replay(details, operators, counts, [1], **args)
        self.assertLess(upper, arb(2)**-32)
        self.assertRaises(ValueError, support_cover.replay, details,
            {('.0001', '1'): (exact[:2], None)}, counts, [1], **args)
        details['leaves'][0]['box'] = [[6, 256], [6, 256]]
        self.assertRaises(ValueError, support_cover.replay, details, operators, counts, [1], **args)

    def test_receipt_scope_rejects_other_ensembles(self):
        import copy
        record = dict(schema='packed-gl32-sparse-1', K=1 << 20, N=1 << 21,
            group_count=2048, outer='BCH256128', block_rows=4, block_columns=8,
            mixing='independent uniform GL32 per group/block; no additional GF16 stage',
            routing='independent shared column shuffle per group and independent regional shuffles',
            inner=dict(t=128, s=19, updates=2), support_min=5,
            distance='12/125', threshold=201326, results=[dict(occupancy=1, passed=True, upper=[1, -40])])
        self.assertEqual(len(packed.validate_record(record)), 1)
        for key, value in (('support_min', 38), ('block_columns', 4), ('threshold', 201325),
                           ('mixing', 'one matrix shared across groups'), ('inner', dict(t=128, s=19, updates=1))):
            changed = copy.deepcopy(record)
            changed[key] = value
            self.assertRaises(ValueError, packed.validate_record, changed)

    def test_existing_receipt_is_never_overwritten(self):
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'receipt.json'
            path.write_text('preserve')
            self.assertRaises(ValueError, packed.run, [1], '.095', ['.001'], output=path)
            self.assertEqual(path.read_text(), 'preserve')


if __name__ == '__main__':
    unittest.main()
