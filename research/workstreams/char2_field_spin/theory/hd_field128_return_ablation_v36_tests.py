import unittest
import numpy as np
import hd_field128_return_ablation_v36 as subject


class Tests(unittest.TestCase):
    def test_only_return_changes(self):
        logs = np.arange(20, dtype=float).reshape(5, 2, 2)-20
        copy = logs.copy()
        altered = subject.ablate_returns(logs)
        np.testing.assert_array_equal(logs, copy)
        np.testing.assert_array_equal(altered[:, 0], logs[:, 0])
        np.testing.assert_array_equal(altered[:, 1, 1], logs[:, 1, 1])
        self.assertTrue(np.isneginf(altered[:, 1, 0]).all())

    def test_tiny_explicit_polynomial(self):
        # One pair-slot perpacket,2packets perregion,2regions,alpha1.
        matrices = np.array([[[1., 0.], [0., .25]],
                             [[.2, .3], [.1, .4]], [[.1, .5], [.2, .5]]])
        with np.errstate(divide='ignore'):
            logs = np.log2(matrices)
        for alter in (False, True):
            current = subject.ablate_returns(logs) if alter else logs
            m = matrices.copy()
            if alter:
                m[:, 1, 0] = 0
            packet = m.sum(axis=0)  # odds=(0,0),all3histograms count1
            regional = packet@packet/2  # 2slots,weight1,full0
            total = regional@regional
            result = subject.evaluate(current, 1., [['1', '1']], 1, 0, 0,
                                      pairs=2, blocks=1, repetitions=1)
            self.assertAlmostEqual(result['transfer'], np.log2(total[0].sum()), places=12)

    def test_bad_shape_and_nonfinite(self):
        for bad in (np.zeros((2, 2)), np.full((3, 2, 2), np.nan), np.full((3, 2, 2), np.inf)):
            with self.assertRaises(ValueError):
                subject.ablate_returns(bad)

    def test_json_structural_zero(self):
        self.assertEqual(subject.finite_json([1., -np.inf]), [1., '-inf'])
        with self.assertRaises(ValueError):
            subject.finite_json(np.nan)


if __name__ == '__main__':
    unittest.main()
