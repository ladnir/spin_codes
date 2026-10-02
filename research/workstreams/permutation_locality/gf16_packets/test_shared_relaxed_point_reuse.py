import copy
from fractions import Fraction as Q
import unittest

import shared_relaxed_point_reuse as reuse


class PointReuseTests(unittest.TestCase):
    def witness(self):
        parts = [dict(interval=['0', '1/16'], dual=['2', '3', '4']),
                 dict(interval=['1/16', '1/8'], dual=['-2', '0', '-4'])]
        return dict(tilt='3/16', parameters=['1/20', '2', '3'], variance_dual=['0', '0'],
            variance_partition=parts, regional_count_parts=[dict(part,
                mgf_witnesses=[dict(tilt='-1', slopes=['1', '2', '-3'])]) for part in parts])

    def test_rescale_both_directions_without_changing_source_or_duals(self):
        source = (Q(1,8), Q(1,8))
        witness = self.witness()
        original = copy.deepcopy(witness)
        for cell, upper in (((Q(1,4), Q(1,3)), Q(1,4)),
                            ((Q(1,64), Q(1,32)), Q(1,32))):
            result = reuse.rebase_witness(source, cell, witness)
            self.assertEqual(witness, original)
            parts = reuse.variance.validate(cell, result['variance_partition'])
            self.assertEqual(parts[-1][0][1], upper)
            self.assertEqual(parts[0][0][1], upper/2)
            self.assertEqual(result['parameters'], witness['parameters'])
            for old, new, part in zip(witness['regional_count_parts'],
                    result['regional_count_parts'], result['variance_partition']):
                self.assertEqual(old['dual'], new['dual'])
                self.assertEqual(old['mgf_witnesses'], new['mgf_witnesses'])
                self.assertEqual(new['interval'], part['interval'])
                self.assertEqual(new['dual'], part['dual'])
            self.assertNotIn('upper', result)

    def test_incomplete_or_mismatched_partitions_rejected(self):
        point = (Q(1,8), Q(1,8))
        for mutation in ('gap', 'missing', 'interval', 'dual', 'regional-only'):
            witness = self.witness()
            if mutation == 'gap': witness['variance_partition'][0]['interval'][0] = '1/128'
            if mutation == 'missing': witness['regional_count_parts'].pop()
            if mutation == 'interval': witness['regional_count_parts'][0]['interval'] = ['0', '1/32']
            if mutation == 'dual': witness['regional_count_parts'][0]['dual'] = ['2', '4', '4']
            if mutation == 'regional-only': witness.pop('variance_partition')
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                reuse.rebase_witness(point, (Q(1,4), Q(1,3)), witness)

    def test_invalid_endpoints_rejected(self):
        for cell in ((0,Q(1,2)), (Q(1,2),1), (Q(1,3),Q(1,4))):
            with self.assertRaises(ValueError):
                reuse.rebase_witness((Q(1,8), Q(1,8)), cell, self.witness())

    def test_nearest_point_does_not_consult_saved_bounds(self):
        probes = [dict(mean='1/8', witness=self.witness(), upper='false', proposal='false'),
                  dict(mean='1/2', witness=dict(tilt='3/16', parameters=['1','0','0']), upper=None)]
        index, mean, witness = reuse.nearest_witness(probes, (Q(1,5), Q(1,4)))
        self.assertEqual(index, 0)
        self.assertEqual(mean, Q(1,8))
        reuse.variance.validate((Q(1,5), Q(1,4)), witness['variance_partition'])
        self.assertNotIn('upper', witness)
        probes[0]['upper'] = [1,100000]
        self.assertEqual(reuse.nearest_witness(probes, (Q(1,5), Q(1,4)))[2], witness)

    def test_plain_witness_copied_not_aliased(self):
        witness = dict(tilt='3/16', parameters=['1','0','0'])
        result = reuse.rebase_witness((Q(1,8), Q(1,8)), (Q(1,4), Q(1,3)), witness)
        self.assertEqual(result, witness)
        result['parameters'][0] = '2'
        self.assertEqual(witness['parameters'][0], '1')

    def test_point_descendants_have_exact_tree_geometry(self):
        parent = (Q(1,8), Q(1,4)); point = Q(1,5)
        for depth in range(8):
            suffix, cell = reuse.point_descendant(parent, point, depth)
            self.assertEqual(len(suffix), depth)
            self.assertLessEqual(cell[0], point)
            self.assertGreaterEqual(cell[1], point)
            self.assertEqual(cell[1]-cell[0], (parent[1]-parent[0])/2**depth)
            replay = parent
            for digit in suffix:
                mid = sum(replay)/2
                replay = (replay[0], mid) if digit == '0' else (mid, replay[1])
            self.assertEqual(cell, replay)
        for point, depth in ((Q(1,2), 1), (Q(1,5), True), (Q(1,5), -1), (Q(1,5), 17)):
            with self.assertRaises(ValueError): reuse.point_descendant(parent, point, depth)


if __name__ == '__main__':
    unittest.main()
