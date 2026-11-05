"""Small controller/scope tests; regional operators are mocked here."""
from copy import deepcopy
from fractions import Fraction as F
from unittest.mock import patch
import unittest
from flint import arb
import hd_lambda_high_arity_certificate as target


def bank(mark,alpha,maximum,bits):
    assert (mark,alpha,maximum,bits)==(F(691,700),F(1,5),18,60)
    return {(18,f):tuple(arb(x)/arb(2)**128 for x in (2,1,3,2)) for f in range(10)}


class LambdaHighTests(unittest.TestCase):
    def test_binding_and_fresh_call(self):
        target.checked_dependencies()
        with patch.object(target.arithmetic.uniform,'regional_directed',side_effect=bank) as mock:
            data=target.receipt(256)
        mock.assert_called_once_with(target.MARK,target.ALPHA,18,60)
        rows,ceiling=target.checked_receipt(data,256)
        self.assertEqual(len(rows),99)
        self.assertNotEqual(data['construction'],target.arithmetic.CONSTRUCTION)
        self.assertFalse(data['crossed_outer_layer'])
        self.assertFalse(data['whole_code_certificate'])

    def test_scope_tampering_rejected(self):
        with patch.object(target.arithmetic.uniform,'regional_directed',side_effect=bank):
            data=target.receipt(256)
        mutations=(lambda x:x.update(lambda_aes_byte=2),lambda x:x.update(total_arity_minimum=246),
            lambda x:x.update(input_schedule='static'),lambda x:x.update(whole_code_certificate=True),
            lambda x:x['profiles'].pop(),lambda x:x['profiles'].append(deepcopy(x['profiles'][0])))
        for mutate in mutations:
            changed=deepcopy(data);mutate(changed)
            with self.assertRaises(ValueError):target.checked_receipt(changed,256)

    def test_replay_guards(self):
        for kwargs in ({'precision':True},{'precision':512},{'max_seconds':float('nan')},
                       {'max_seconds':True},{'max_seconds':61}):
            with self.assertRaises(ValueError):target.receipt(**kwargs)


if __name__=='__main__':unittest.main()
