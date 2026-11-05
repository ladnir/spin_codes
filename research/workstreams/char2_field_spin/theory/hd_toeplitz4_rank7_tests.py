"""Small exact algebra checks; no native execution or performance timing."""
import ast
import inspect
from itertools import product
import unittest
from unittest.mock import patch

import hd_toeplitz4_rank7 as circuit


def basis(length):
    for bit in range(8*length):
        yield tuple((1 << (bit % 8)) if i == bit//8 else 0 for i in range(length))


class RankSevenTests(unittest.TestCase):
    def test_actual_subfields_and_explicit_constants(self):
        w, g = circuit.OMEGA, circuit.GAMMA
        self.assertEqual(circuit.mul(w, w) ^ w, 1)
        self.assertEqual(circuit.mul(g, g) ^ g, w)
        self.assertEqual(circuit.power(g, 4) ^ g, 1)
        self.assertEqual(circuit.power(g, 16), g)
        self.assertEqual(sum(circuit.power(x, 8) == x for x in range(256)), 2)
        self.assertEqual(sum(circuit.power(x, 4) == x for x in range(256)), 4)
        self.assertEqual(sum(circuit.power(x, 16) == x for x in range(256)), 16)
        self.assertEqual(len(set(circuit.POINTS)), 7)

    def test_explicit_evaluations_and_algebraic_transposes(self):
        evaluation = tuple(circuit.evaluation_row(point, 4) for point in circuit.POINTS)
        c = circuit.congruence_matrix()
        for x in basis(4):
            self.assertEqual(circuit.eval7_novel(circuit.monomial_to_novel(x)),
                             circuit.apply(evaluation, x))
            self.assertEqual(circuit.eval7_congruent(x),
                             circuit.apply(evaluation, circuit.apply(c, x)))
        for r in basis(7):
            self.assertEqual(circuit.novel_transpose_to_monomial(circuit.eval7_novel_transpose(r)),
                             circuit.apply(circuit.transpose(evaluation), r))
            self.assertEqual(circuit.eval7_congruent_transpose(r),
                             circuit.apply(circuit.transpose(c), circuit.apply(circuit.transpose(evaluation), r)))
        for x in basis(4):
            for r in basis(7):
                self.assertEqual(circuit.dot(circuit.quadratic_eval(x), r),
                                 circuit.dot(x, circuit.quadratic_eval_transpose(r)))

    def test_key_transform_reconstructs_all_diagonals(self):
        v_transpose = circuit.transpose(tuple(circuit.evaluation_row(point, 7)
                                             for point in circuit.POINTS))
        for t in basis(7):
            self.assertEqual(circuit.apply(v_transpose, circuit.prepare7(t)), t)

    def test_every_binary_bilinear_basis_pair(self):
        # 56 key bits x32 input bits per circuit proves each GF2-bilinear map.
        for t in basis(7):
            keys7, keys8 = circuit.prepare7(t), circuit.prepare8(t)
            for x in basis(4):
                self.assertEqual(circuit.original_from_keys(keys7, x), circuit.dense_original(t, x))
                expected = circuit.dense_congruent(t, x)
                self.assertEqual(circuit.congruent_from_keys(keys7, x), expected)
                self.assertEqual(circuit.quadratic_from_keys(keys8, x), expected)

    def test_operation_counts_are_from_explicit_schedules(self):
        functions = (
            (circuit.eval7_novel, 4, 9), (circuit.eval7_novel_transpose, 4, 12),
            (circuit.monomial_to_novel, 0, 2), (circuit.novel_transpose_to_monomial, 0, 2),
            (circuit.eval7_congruent, 3, 8), (circuit.eval7_congruent_transpose, 3, 11),
            (circuit.quadratic_eval, 2, 6), (circuit.quadratic_eval_transpose, 2, 9))
        for function, fixed, xors in functions:
            tree = ast.parse(inspect.getsource(function))
            self.assertEqual(sum(isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                                 and node.func.id == 'mul' for node in ast.walk(tree)), fixed)
            self.assertEqual(sum(isinstance(node, ast.BinOp) and isinstance(node.op, ast.BitXor)
                                 for node in ast.walk(tree)), xors)
        keys7 = circuit.prepare7((1, 2, 3, 4, 5, 6, 7))
        keys8 = circuit.prepare8((1, 2, 3, 4, 5, 6, 7))
        for function, keys, name in (
            (circuit.original_from_keys, keys7, 'original_toeplitz_rank7'),
            (circuit.congruent_from_keys, keys7, 'congruent_rank7'),
            (circuit.quadratic_from_keys, keys8, 'congruent_quadratic_rank8')):
            with patch.object(circuit, 'mul', wraps=circuit.mul) as multiply:
                function(keys, (8, 9, 10, 11))
            count = circuit.COUNTS[name]
            self.assertEqual(multiply.call_count, count['variable_products']+count['fixed_maps'])

    def test_only_five_GF4_rank_one_Hankel_directions(self):
        subfield = (0, 1, circuit.OMEGA, circuit.OMEGA2)
        directions = set()
        for v in product(subfield, repeat=4):
            if not any(v):
                continue
            matrix = tuple(tuple(circuit.mul(a, b) for b in v) for a in v)
            diagonals = []
            for degree in range(7):
                entries = {matrix[i][degree-i] for i in range(4) if 0 <= degree-i < 4}
                if len(entries) != 1:
                    break
                diagonals.append(entries.pop())
            if len(diagonals) == 7:
                first = next(x for x in diagonals if x)
                scale = circuit.power(first, 254)
                directions.add(tuple(circuit.mul(scale, x) for x in diagonals))
        expected = {circuit.evaluation_row(x, 7) for x in (*subfield, None)}
        self.assertEqual(directions, expected)
        self.assertEqual(len(directions), 5)

    def test_congruence_is_not_the_old_coordinate_map(self):
        t, x = (1, 2, 3, 4, 5, 6, 7), (1, 0, 0, 0)
        self.assertNotEqual(circuit.dense_original(t, x), circuit.dense_congruent(t, x))
        self.assertEqual(circuit.original_from_keys(circuit.prepare7((0,)*7), x), (0,)*4)
        self.assertEqual(circuit.congruent_from_keys(circuit.prepare7((0,)*7), x), (0,)*4)


if __name__ == '__main__':
    unittest.main()
