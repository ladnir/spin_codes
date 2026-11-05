"""Pure AST-isolated controls; no numerical module import or search.

Root runs these remotely. The tested functions come from the held source,
but only decode_point and the first two run guards are executed.
"""
import ast
from fractions import Fraction
from hashlib import sha256
from math import isfinite, log2
from pathlib import Path
import unittest

SOURCE_SHA = '88e2450b79e6293a39bc89b552ceded61cd63474979f60c922dbf760b0d65e96'
REMOTE = Path('/tmp/spin-field128-proof-52GBmlZD/bundle/theory')


def isolated():
    path = Path(__file__).with_name('hd_field128_wide_state_retune_v63.py')
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != SOURCE_SHA:
        raise AssertionError('Held V63 source changed')
    tree = ast.parse(raw)
    functions = {node.name: node for node in tree.body
                 if isinstance(node, ast.FunctionDef)}
    guard = functions['run']
    if not all(isinstance(node, ast.If) for node in guard.body[:2]):
        raise AssertionError('Run guard layout changed')
    guard.body = guard.body[:2] + [ast.Return(value=ast.Constant(True))]
    module = ast.fix_missing_locations(ast.Module(
        body=[functions['decode_point'], guard], type_ignores=[]))
    namespace = dict(isfinite=isfinite, log2=log2, Path=Path, HERE=REMOTE)
    exec(compile(module, str(path), 'exec'), namespace)
    return namespace


class Tests(unittest.TestCase):
    def setUp(self):
        self.scope = isolated()
        self.decode = self.scope['decode_point']

    def test_origin_and_separate_marker_coordinates(self):
        z, alpha, left, right = self.decode(
            [0, 0, 0, 0, 0, 0], 2.**-16, 2.**-8, (1., 2.), (3., 4.))
        self.assertEqual((z, alpha, left, right),
                         (2.**-16, 2.**-8, (1., 2.), (3., 4.)))
        value = self.decode([1, 1, .5, -.5, -.25, .25],
                            z, alpha, left, right)
        self.assertEqual(value, (2.**-32, 2.**-7, (1.5, 1.5), (2.75, 4.25)))
        self.assertEqual(float(Fraction(value[0])), value[0])
        self.assertEqual(float(Fraction(value[1])), value[1])

    def test_reject_nonfinite_or_outside_chart(self):
        for point in ([0]*5, [0]*7, [float('nan')]+[0]*5,
                      [float('inf')]+[0]*5, [4.01]+[0]*5):
            with self.subTest(point=point), self.assertRaises(ValueError):
                self.decode(point, .25, .25, (0., 0.), (0., 0.))
        for point, z, alpha, left in (
                ([4, 0, 0, 0, 0, 0], 2.**-100, .25, (0., 0.)),
                ([0, 4, 0, 0, 0, 0], .25, .25, (0., 0.)),
                ([0, 0, 1, 0, 0, 0], .25, .25, (32., 0.))):
            with self.assertRaises(ValueError):
                self.decode(point, z, alpha, left, (0., 0.))

    def test_output_round_and_remote_guards_only(self):
        guard = self.scope['run']
        self.assertTrue(guard('result.json', rounds=1))
        self.assertTrue(guard('result.json', rounds=2))
        for name, rounds in (('../result.json', 1), ('result.txt', 1),
                             ('result.json', 0), ('result.json', 3)):
            with self.assertRaises(ValueError):
                guard(name, rounds=rounds)
        self.scope['HERE'] = Path('/not-the-approved-remote-directory')
        with self.assertRaises(RuntimeError):
            guard('result.json')


if __name__ == '__main__':
    unittest.main()
