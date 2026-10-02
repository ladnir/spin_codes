import unittest
from fractions import Fraction as Q
from types import SimpleNamespace
from shared_relaxed_dense import cutoff, reuse_exact_census, starting_mixture


class SharedRelaxedDenseTests(unittest.TestCase):
    def test_cutoff_is_exact_floor(self):
        for distance in (Q(1, 20), Q(7, 100), Q(2, 25), Q(9, 100)):
            value = cutoff(distance)
            self.assertLessEqual(value, distance * (1 << 21))
            self.assertGreater(value + 1, distance * (1 << 21))
        self.assertEqual(cutoff('0.05'), 104857)

    def test_invalid_distance(self):
        for distance in ('0', '-.01', '1/2', '1'):
            with self.assertRaises(ValueError):
                cutoff(distance)

    def test_starting_mixture_extracts_only_exact_coefficients(self):
        rows = [dict(mass='5', activity='1/2')]
        for record in (dict(schema='shared-gf16-relaxed-dense-1', mixture=rows, upper=[1,-999]),
                dict(schema='shared-relaxed-alternative-mixtures-1', rows=[dict(mixture=rows)])):
            result = starting_mixture(record)
            self.assertEqual(result, rows)
            self.assertIsNot(result[0], rows[0])
        for record in (dict(schema='unknown', mixture=rows),
                dict(schema='shared-relaxed-alternative-mixtures-1', rows=[]),
                dict(schema='shared-relaxed-alternative-mixtures-1', rows=[dict(mixture=rows)]*2)):
            with self.assertRaises(ValueError): starting_mixture(record)
        for rows in ([], [dict(mass='-1', activity='1/2')], [dict(mass=True, activity='1/2')],
                [dict(mass='1', activity=.5)], [dict(mass='1', activity='0')]):
            with self.assertRaises(ValueError):
                starting_mixture(dict(schema='shared-gf16-relaxed-dense-1', mixture=rows))

    def test_only_exact_census_from_identical_inner_is_shared(self):
        data = {}
        old = SimpleNamespace(data=data, joint_return_cache={3: [1, 2]},
            fiber_density_data={'all_fiber_trimmed': [4]}, cache={'weighted': 5},
            regional_proposal_scratch={'old': 6}, variance_cache={'dual': 7})
        new = SimpleNamespace(data=data)
        reuse_exact_census(old, new)
        self.assertIs(new.joint_return_cache, old.joint_return_cache)
        self.assertIs(new.fiber_density_data, old.fiber_density_data)
        for name in ('cache', 'regional_proposal_scratch', 'variance_cache'):
            self.assertFalse(hasattr(new, name))
        with self.assertRaises(ValueError):
            reuse_exact_census(old, SimpleNamespace(data={}))


if __name__ == '__main__':
    unittest.main()
