from fractions import Fraction as F
from itertools import product
from unittest import TestCase, main

import hd_coupled_state_return as gate


class ThreePacketReturnTests(TestCase):
    def test_gf4_and_actual_refresh_law(self):
        self.assertEqual(gate.multiply(2, 2), 3)
        self.assertEqual(gate.multiply(2, 3), 1)
        actions = gate.refresh_actions()
        self.assertEqual(len(actions), 81)
        for action in actions:
            self.assertEqual(action[0], 0)
            self.assertEqual(set(action), set(range(16)))
        kernel = gate.kernel_counts()
        self.assertEqual(max(F(n, 81) for row in kernel[1:] for n in row[1:]), F(1, 9))
        for i in range(16):
            self.assertEqual(sum(kernel[i]), 81)
            self.assertEqual(sum(row[i] for row in kernel), 81)

    def test_all_conditional_triples_before_refresh_identity(self):
        counts, kernel = gate.direct_three_return_counts('before'), gate.kernel_counts()
        self.assertEqual(len(counts), 15**3)
        for (c1, c2, c3), count in counts.items():
            expected = kernel[c1][c3]*81**2 if c2 == c3 else 0
            self.assertEqual(count, expected)
            self.assertLessEqual(F(count, 81**3), F(1, 9))

    def test_all_conditional_triples_after_refresh_identity(self):
        counts, kernel = gate.direct_three_return_counts('after'), gate.kernel_counts()
        for (c1, c2, c3), count in counts.items():
            self.assertEqual(count, kernel[c1][c3]*kernel[c2][c3]*81)
            self.assertLessEqual(F(count, 81**3), F(1, 81))

    def test_every_first_two_key_pair_obeys_symbolic_return_condition(self):
        # One fixed third key suffices here because the other test enumerates
        # all third keys. Each of these input triples covers all81**2 first keys.
        actions = gate.refresh_actions()
        for c1, c2, c3 in ((1, 5, 5), (1, 1, 5), (7, 9, 14)):
            for first, second in product(actions, repeat=2):
                for order in ('before', 'after'):
                    state = (0, 0)
                    for c, action in zip((c1, c2, c3), (first, second, actions[0])):
                        state = gate.transition(state, c, action, order)
                    expected = ((c2 == c3 and first[c1] == c3) if order == 'before'
                                else first[c1] == second[c2] == c3)
                    self.assertEqual(state == (0, 0), expected)

    def test_worst_case_gain_is_not_a_uniform_input_average_gain(self):
        result = gate.summary()
        self.assertEqual(result['before_maximum'], '1/9')
        self.assertEqual(result['after_maximum'], '1/81')
        self.assertEqual(result['before_uniform_feedback_average'], '1/225')
        self.assertEqual(result['after_uniform_feedback_average'], '1/225')


if __name__ == '__main__':
    main()
