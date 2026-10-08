"""Independent GF4 enumeration of the proposed local rank lemmas.

This tests local conditional inequalities, not an active-input operator
or a whole-code rank certificate.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from arithmetic import multiply
from test_return_audit import rank

Q, MODULUS = 4, 0x7


def mul(x, y):
    return multiply(x, y, MODULUS)


def add_vectors(vectors, dimension=2):
    result = [0] * dimension
    for vector in vectors:
        for j, value in enumerate(vector):
            result[j] ^= value
    return tuple(result)


def scale(value, vector):
    return tuple(mul(value, x) for x in vector)


def matmul(a, b):
    return [add_vectors((scale(value, row) for value, row in zip(left, b)))
            for left in a]


def return_gate(active, zeros, width=4):
    canceled = tuple(j for j in active if j in zeros)
    uncanceled = tuple(j for j in active if j not in zeros)
    inactive_zeros = tuple(j for j in zeros if j not in active)
    diagonal = [[1, 0], [0, 1]]
    for j in canceled:
        diagonal[0][0] ^= 1
        diagonal[0][1] ^= j
        diagonal[1][0] ^= j
        diagonal[1][1] ^= mul(j, j)
    if not uncanceled:
        annihilator = ((1, 0), (0, 1))
    elif len(uncanceled) == 1:
        annihilator = ((uncanceled[0], 1),)
    else:
        annihilator = ()
    gate = [(1, j) for j in inactive_zeros] + matmul(annihilator, diagonal)
    return gate, len(canceled) + min(2, len(uncanceled))


def evaluate(form, message):
    value = 0
    for coefficient, entry in zip(form, message):
        value ^= mul(coefficient, entry)
    return value


def alive(equations, required_nonzero=(), state=()):
    for message in product(range(Q), repeat=2):
        if not any(message):
            continue
        if any(evaluate(equation, message) for equation in equations):
            continue
        if any(evaluate(form, message) == 0 for form in required_nonzero):
            continue
        if state and not any(evaluate(form, message) for form in state):
            continue
        return True
    return False


SL2 = [((a, b), (c, d)) for a, b, c, d in product(range(Q), repeat=4)
       if mul(a, d) ^ mul(b, c) == 1]


class RankLemmaAudit(unittest.TestCase):
    def test_diagonal_pivots_with_dependent_required_forms(self):
        forms = [form for form in product(range(Q), repeat=2) if any(form)]
        coefficient_matrices = (
            ((1, 0), (0, 1)), ((1, 1), (0, 0)),
            ((1, 2), (2, 3)), ((0, 1), (0, 0)),
        )
        cases = [(matrix, selected) for matrix in coefficient_matrices
                 for selected in product(forms, repeat=2)]
        # Nonpivot exposure and pivot forms killed by prior equations.
        cases += [(matrix, selected)
                  for matrix in (((1, 0, 1), (0, 1, 1)),
                                 ((1, 1, 0), (0, 1, 1)))
                  for selected in (((1, 0), (1, 0), (1, 0)),
                                   ((1, 0), (0, 1), (1, 1)),
                                   ((1, 0), (1, 0), (0, 1)))]
        for matrix, selected in cases:
            histogram = Counter()
            for diagonal in product(range(1, Q), repeat=len(selected)):
                weighted = [scale(d, form) for d, form in zip(diagonal, selected)]
                equations = matmul(matrix, weighted)
                if alive(equations, selected):
                    histogram[rank(equations, MODULUS)] += 1
            for tau in (F(1), F(2), F(4)):
                mean = sum(count * tau ** -r for r, count in histogram.items())
                mean /= (Q - 1) ** len(selected)
                g = F(1, Q - 1) + F(Q - 2, Q - 1) / tau
                self.assertLessEqual(mean, g ** rank(matrix, MODULUS))

    def test_required_nonzero_condition_is_essential(self):
        # Requiring f != 0 makes f=0 a dead branch. Merely requiring
        # the message nonzero would leave witnesses (0,1), incorrectly
        # permitting the second dependent pivot to consume another g.
        tau = F(2)
        matrix = ((1, 0), (0, 1))
        selected = ((1, 0), (1, 0))
        mean_without_requirement = F(0)
        for diagonal in product(range(1, Q), repeat=2):
            equations = matmul(matrix, [scale(d, form) for d, form in zip(diagonal, selected)])
            self.assertTrue(alive(equations))
            self.assertFalse(alive(equations, selected))
            mean_without_requirement += tau ** -rank(equations, MODULUS)
        mean_without_requirement /= (Q - 1) ** 2
        g = F(1, Q - 1) + F(Q - 2, Q - 1) / tau
        self.assertGreater(mean_without_requirement, g ** 2)

    def test_sl2_return_rank_factors_for_every_fixed_target(self):
        states = (((1, 0), (0, 0)), ((1, 0), (0, 1)))
        for state_rank, state in enumerate(states, 1):
            for entries in product(range(Q), repeat=4):
                target = (entries[:2], entries[2:])
                histogram = Counter()
                for matrix in SL2:
                    refreshed = matmul(matrix, state)
                    equations = [add_vectors((row, expected))
                                 for row, expected in zip(refreshed, target)]
                    histogram[rank(equations, MODULUS)] += 1
                self.assertLessEqual(histogram[0],
                                     len(SL2) // (Q * Q - 1) if state_rank == 1 else 1)
                if state_rank == 2:
                    self.assertLessEqual(F(histogram[0] + histogram[1], len(SL2)),
                                         F(1, Q - 1))
                for tau in (F(1), F(2), F(4)):
                    mean = sum(count * tau ** -r for r, count in histogram.items()) / len(SL2)
                    if state_rank == 1:
                        upper = 1 / tau + (1 - 1 / tau) / (Q * Q - 1)
                    else:
                        upper = (tau ** -2 + (tau ** -1 - tau ** -2) / (Q - 1)
                                 + (1 - tau ** -1) / (Q * (Q * Q - 1)))
                    self.assertLessEqual(mean, upper)

    def test_empty_input_rank_operator(self):
        width, z, tau = 4, F(3, 4), F(2)
        states = (((1, 0), (0, 0)), ((1, 0), (0, 1)))
        g0 = F(Q + 1 - width, Q + 1) * z ** width
        g0 += F(width, Q + 1) * z ** (width - 1)
        for state_rank, state in enumerate(states, 1):
            totals = [F(0), F(0), F(0)]
            for matrix in SL2:
                refreshed = matmul(matrix, state)
                outputs = matmul([(1, j) for j in range(width)], refreshed)
                for count in range(width + 1):
                    for zeros in combinations(range(width), count):
                        equations = [outputs[j] for j in zeros]
                        nonzeros = [outputs[j] for j in range(width) if j not in zeros]
                        if not alive(equations, nonzeros, state):
                            continue
                        lost = rank(equations, MODULUS)
                        remaining = rank(list(equations) + refreshed, MODULUS) - lost
                        totals[remaining] += tau ** -lost * z ** (width - count) / len(SL2)
            expected = [F(0), g0, F(0)] if state_rank == 1 else [
                F(0), width / tau * z ** (width - 1), z ** width]
            self.assertEqual(totals, expected)

    def test_return_gate_is_necessary_and_low_p_never_free(self):
        width = 4
        for inputs in product(range(Q), repeat=width):
            active = tuple(j for j, x in enumerate(inputs) if x)
            refreshed = add_vectors(scale(x, (1, j)) for j, x in enumerate(inputs))
            if not any(refreshed):
                continue
            outputs = tuple(x ^ evaluate((1, j), refreshed) for j, x in enumerate(inputs))
            zeros = tuple(j for j, y in enumerate(outputs) if not y)
            gate, _ = return_gate(active, zeros)
            self.assertFalse(any(evaluate(equation, refreshed) for equation in gate))
        for mask in range(16):
            active = tuple(j for j in range(width) if (mask >> j) & 1)
            for zero_mask in range(16):
                zeros = tuple(j for j in range(width) if (zero_mask >> j) & 1)
                uncanceled = tuple(j for j in active if j not in zeros)
                if len(uncanceled) <= 1:
                    gate, _ = return_gate(active, zeros)
                    self.assertGreaterEqual(rank(gate, MODULUS), 1)

    def test_sharper_return_bound_with_fresh_matrix_then_diagonal(self):
        width, tau = 4, F(2)
        g = F(1, Q - 1) + F(Q - 2, Q - 1) / tau
        cases = []
        for j in range(width):
            cases += [((j,), ((1, 0),), zeros)
                      for zeros in ((), (0,), (1,), (2,), (3,), (0, 1))]
        cases += [((0, 1), forms, zeros)
                  for forms in (((1, 0), (1, 0)), ((1, 0), (0, 1)))
                  for zeros in ((), (0,), (1,), (0, 1))]
        cases += [((1, 2, 3), ((1, 0), (0, 1), (1, 1)), zeros)
                  for zeros in ((), (1,), (1, 2, 3))]
        for state_rank, state in enumerate((
                ((1, 0), (0, 0)), ((1, 0), (0, 1))), 1):
            for active, forms, zeros in cases:
                gate, pivots = return_gate(active, zeros)
                gate_rank = rank(gate, MODULUS)
                if gate_rank == 2:
                    upper = F(0)
                elif gate_rank == 1:
                    upper = (F(1, Q + 1) if state_rank == 1 else 1 / tau) * g ** pivots
                else:
                    upper = g ** pivots
                total = F(0)
                for matrix in SL2:
                    refreshed = matmul(matrix, state)
                    for diagonal in product(range(1, Q), repeat=len(active)):
                        inputs = {j: scale(d, form) for j, d, form in zip(active, diagonal, forms)}
                        outputs = [add_vectors((inputs.get(j, (0, 0)),
                                                refreshed[0], scale(j, refreshed[1])))
                                   for j in range(width)]
                        next_state = [add_vectors((refreshed[0],) + tuple(inputs.values())),
                                      add_vectors((refreshed[1],) + tuple(
                                          scale(j, form) for j, form in inputs.items()))]
                        equations = [outputs[j] for j in zeros] + next_state
                        required = list(forms) + [outputs[j] for j in range(width) if j not in zeros]
                        if alive(equations, required, state):
                            total += tau ** -rank(equations, MODULUS)
                total /= len(SL2) * (Q - 1) ** len(active)
                self.assertLessEqual(total, upper)


if __name__ == '__main__':
    unittest.main()
