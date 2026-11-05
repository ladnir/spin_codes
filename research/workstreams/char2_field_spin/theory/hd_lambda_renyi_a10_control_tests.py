"""Inventory validation and same-formula regression tests; no regional run."""
from copy import deepcopy
from fractions import Fraction as F
import json
import unittest

from flint import ctx
import hd_lambda_renyi_a10_control as control


class ArityTenControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((control.DIRECTORY/'hd_lambda_renyi_a10_inventory234.json').read_text())

    def test_complete_labelled_inventory(self):
        moments = control.checked_a10(self.data)
        self.assertEqual(len(moments), 8008)
        self.assertTrue(all(len(key) == 6 for key in moments))

    def test_reject_missing_row(self):
        data = deepcopy(self.data)
        data['rows'].pop()
        with self.assertRaises(ValueError):
            control.checked_a10(data)

    def test_reject_changed_coordinate_witness(self):
        data = deepcopy(self.data)
        data['rows'][0]['representative_to_complement_permutation'].reverse()
        with self.assertRaises(ValueError):
            control.checked_a10(data)

    def test_reject_wrong_scope(self):
        data = dict(self.data, arities=[11])
        with self.assertRaises(ValueError):
            control.checked_a10(data)

    def test_reject_impossible_moment(self):
        data = deepcopy(self.data)
        data['rows'][0]['moments']['4'] = '1'
        with self.assertRaises(ValueError):
            control.checked_a10(data)

    def test_same_existing_scheduled_weights(self):
        moments = control.inputs()
        self.assertEqual(len(moments), 14893)
        with ctx.workprec(128):
            for order in (2, 3, 4):
                weights, diagnostics = control.scheduled_weights(moments, order)
                prior, prior_diagnostics = control.activity.scheduled_weights(moments, order)
                self.assertEqual([str(x) for x in weights[1:]], [str(x) for x in prior])
                self.assertEqual(diagnostics[1:], prior_diagnostics)
                self.assertEqual(diagnostics[0]['labelled_rowsets'], 8008)

    def test_reject_unsupported_order(self):
        with self.assertRaises(ValueError):
            control.scheduled_weights({}, 5)


if __name__ == '__main__':
    unittest.main()
