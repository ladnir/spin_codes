import copy
import unittest

import quarter_binding as binding


class QuarterBindingTests(unittest.TestCase):
    def test_rejects_ambiguous_source_replacement(self):
        for original in ('nothing', 'xx'):
            with self.assertRaises(AssertionError):
                binding.replace_once(original, 'x', 'y')
        self.assertEqual(binding.replace_once('axb', 'x', 'y'), 'ayb')

    def test_timing_policy_and_identity(self):
        runs = [binding.prior.read(binding.DATA / f'asymmetric_on-m20-r{i}.jsonl') for i in (1, 2, 3)]
        self.assertEqual(binding.timing_cell(runs, 20)['median_ms'], 15.254337)
        for key, bad in (('inplace', False), ('m', 18), ('outer_dimension', 64),
                         ('trials', 99), ('configuration', 'old-rm2sub'),
                         ('output_hash', 'different-map'), ('workspace_bytes', 0)):
            changed = copy.deepcopy(runs)
            changed[0][key] = bad
            with self.assertRaises(AssertionError, msg=key):
                binding.timing_cell(changed, 20)


if __name__ == '__main__':
    unittest.main()
