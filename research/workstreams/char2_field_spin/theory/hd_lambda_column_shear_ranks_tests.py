from collections import Counter
from itertools import product
from unittest.mock import patch
import unittest

import hd_lambda_column_shear_ranks as shear
import hd_lambda_scheduled_relations as scheduled

MUL = ((0,0,0,0), (0,1,2,3), (0,2,3,1), (0,3,1,2))


def generic_rank(matrix):
    matrix = [list(row) for row in matrix]
    rank = 0
    for column in range(len(matrix[0]) if matrix else 0):
        pivot = next((i for i in range(rank, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = (0,1,3,2)[matrix[rank][column]]
        matrix[rank] = [MUL[inverse][v] for v in matrix[rank]]
        for i in range(len(matrix)):
            if i != rank and matrix[i][column]:
                factor = matrix[i][column]
                matrix[i] = [a ^ MUL[factor][b] for a,b in zip(matrix[i],matrix[rank])]
        rank += 1
    return rank


def apply(rows, values):
    answer = []
    for row in rows:
        value = 0
        for i, x in enumerate(values):
            value ^= MUL[(row >> (2*i)) & 3][x]
        answer.append(value)
    return tuple(answer)


class ColumnShearRankTests(unittest.TestCase):
    def test_packed_rank_against_independent_elimination(self):
        for height, width in ((2,3),(5,7),(16,24),(32,32)):
            matrix = [[((i+1)*(j+3) ^ (i >> 1) ^ (j >> 2)) & 3 for j in range(width)]
                      for i in range(height)]
            self.assertEqual(shear.rank_gf4(tuple(map(shear.pack,matrix)),width),generic_rank(matrix))
        self.assertEqual(shear.rank_gf4((),0),0)

    def test_column_label_order_matches_existing_exact_adapter(self):
        for source, j in product(((0,), (3,12), (0,5,9,14)), (0,1,6,15)):
            self.assertEqual(shear.column_rows(source,j), scheduled.row_forms(source,shear.SCHEDULE,j))

    def test_exact_toy_zero_counts_and_joint_injectivity(self):
        for source, j, n in product(((0,), (0,1), (0,8)), (0,3), (1,7)):
            left,right = shear.paired_rows(source,j,n)
            data = shear.pair_ranks(source,j,n,verify_joint=True)
            counts = Counter()
            a = len(source)
            for values in product(range(4),repeat=2*a):
                u,v = apply(left,values),apply(right,values)
                counts[not any(u),not any(v)] += 1
            self.assertEqual(counts[True,True],1)
            self.assertEqual(sum(n for (z,_),n in counts.items() if z),4**(2*a-data['first_zero_rank']))
            self.assertEqual(sum(n for (_,z),n in counts.items() if z),4**(2*a-data['second_zero_rank']))

    def test_singleton_spreading_and_support_bounds(self):
        for source, j, normal in product(((0,), (7,), (15,)), (0,5), (1,6,15)):
            data = shear.pair_ranks(source,j,normal)
            self.assertEqual((data['first_zero_rank'],data['second_zero_rank']),(2,2))
            left,right = shear.paired_rows(source,j,normal)
            for x,y in product(range(4),repeat=2):
                a,b = apply(left,(x,y)),apply(right,(x,y))
                weight = sum(v != 0 for v in a+b)
                if bool(x) != bool(y):
                    self.assertEqual(weight,24)
                elif x and y:
                    self.assertGreaterEqual(weight,16)
                    self.assertTrue(any(a) and any(b))

    def test_pair_order_and_full_arity(self):
        for source in ((0,1,3,7,9), tuple(range(16))):
            even = shear.pair_ranks(source,4,9,verify_joint=True)
            odd = shear.pair_ranks(source,5,9,verify_joint=True)
            self.assertEqual(even['first_zero_rank'],odd['second_zero_rank'])
            self.assertEqual(even['second_zero_rank'],odd['first_zero_rank'])
            if len(source) == 16:
                self.assertEqual((even['first_zero_rank'],even['second_zero_rank']),(16,16))

    def test_plan_and_bounded_mock_screen(self):
        plan = shear.source_plan()
        counts = Counter(row['arity'] for row in plan)
        self.assertEqual(counts[1],16)
        self.assertEqual(counts[2],120)
        self.assertTrue(all(1 <= counts[a] <= 4 for a in range(3,17)))
        fake = [dict(arity=1,source=(0,),all_labelled_arity_covered=False)]
        with patch.object(shear,'source_plan',return_value=fake):
            data = shear.screen(normals=(1,),columns=(0,1),max_pairs=2,max_seconds=1)
        self.assertEqual(data['completed_pairs'],2)
        self.assertTrue(data['complete'])
        self.assertFalse(data['RS_word_probabilities_claimed'])
        self.assertFalse(data['whole_code_certificate'])

    def test_guards_and_pins(self):
        self.assertEqual(shear.dependencies()['hd_pair_scheduled_outer.py'],shear.PINS['hd_pair_scheduled_outer.py'])
        with self.assertRaises(ValueError):
            shear.pair_ranks((0,0),0)
        with self.assertRaises(ValueError):
            shear.pair_ranks((0,),0,0)
        with self.assertRaises(ValueError):
            shear.rank_gf4((1,),33)
        with self.assertRaises(ValueError):
            shear.screen(max_pairs=1)


if __name__ == '__main__':
    unittest.main()
