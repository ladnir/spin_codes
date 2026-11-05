"""Tiny exact noncommuting controls for the pure bit-syndrome contraction."""
from fractions import Fraction as F
from itertools import product
from math import comb, log2
from unittest import TestCase, main
from unittest.mock import patch

import numpy as np

import hd_bit_syndrome_gate as gate
from hd_bit_syndrome_kernel import information_cover, information_sets


def elementary(values):
    answer = [F(1)]
    for value in values:
        answer.append(F(0))
        for j in range(len(answer) - 1, 0, -1):
            answer[j] += value * answer[j - 1]
    return answer


def atom_words(size, constants):
    p, b00, b01, b10, a01, a10, a11 = constants
    edges = ((0, 0, 1-p, F(0)), (0, 0, p, b00), (0, 1, a01, b01),
             (1, 0, a10, b10), (1, 1, a11, F(1)))
    def visit(initial, state, markers, mass):
        if len(markers) == size:
            yield initial, state, markers, mass
        else:
            for incoming, outgoing, weight, marker in edges:
                if state == incoming and weight:
                    yield from visit(initial, outgoing, markers + (marker,), mass * weight)
    for initial in range(2):
        yield from visit(initial, initial, (), F(1))


def exact_nonquiet_column(source, info, counts, groups, constants, *, size, q):
    answer = [[F(0), F(0)], [F(0), F(0)]]
    outside = tuple(i for i in range(size) if i not in source)
    for start, end, markers, mass in atom_words(size, constants):
        e = elementary(tuple(markers[i] for i in outside))
        value = sum(c * x for c, x in zip(counts, e)) ** groups
        for i in source:
            value *= ((1 + (q-1)*markers[i]) if i in info else max(1, markers[i])) ** groups
        answer[start][end] += mass * value
    return answer


def exact_cover_column(source, cover, counts, groups, constants, *, size, q):
    """Toy exponents have integral 1/c and groups*c, so this is rational."""
    answer = [[F(0), F(0)], [F(0), F(0)]]
    outside = tuple(i for i in range(size) if i not in source)
    for start, end, markers, mass in atom_words(size, constants):
        e = elementary(tuple(markers[i] for i in outside))
        value = sum(c*x for c, x in zip(counts, e)) ** groups
        for i, c in zip(source, cover):
            if not c:
                value *= max(1, markers[i]) ** groups
            else:
                if (1/c).denominator != 1 or (groups*c).denominator != 1:
                    raise AssertionError('Exact toy exponent is not integral')
                value *= (1+(q-1)*markers[i]**int(1/c))**int(groups*c)
        answer[start][end] += mass * value
    return answer


class BitSyndromeGateTests(TestCase):
    def setUp(self):
        # 01 atom deliberately exceeds1: ED must not be replaced by E0.
        self.constants = tuple(map(F, ('2/5', '1/3', '2', '3/5', '3/7', '1/11', '2/13')))
        self.measure = dict(zip(('log2_probability', 'log2_atom00', 'log2_atom01',
            'log2_atom10', 'log2_scale01', 'log2_scale10', 'log2_scale11'), map(log2, self.constants)))

    def test_lexicographic_information_and_labelled_gap_codes(self):
        rows, infos, ranks, gaps = gate.row_geometry()
        self.assertEqual(len(rows), comb(16, 8))
        self.assertEqual(sum(ranks == 4), 30)
        self.assertEqual(sum(ranks == 5), 12840)
        for j in (0, 1, 11, 511, 6000, len(rows)-1):
            self.assertEqual(infos[j], min(information_sets(rows[j])))
            outside = [i for i in range(16) if i not in rows[j]]
            decoded, position = [], 0
            for k, code in enumerate(gaps[j]):
                labels = bin(int(code))[3:]
                for label in labels:
                    self.assertIn(position, rows[j])
                    self.assertEqual(label == '1', position in infos[j])
                    decoded.append(position)
                    position += 1
                if k < len(outside):
                    self.assertEqual(position, outside[k])
                    position += 1
            self.assertEqual(tuple(decoded), rows[j])
            self.assertEqual(position, 16)
        self.assertFalse(gaps.flags.writeable)

    def test_weighted_gaps_and_noncommuting_words(self):
        bank = gate.weighted_gap_matrices(self.measure, 2, q=4)
        expected_i = [[F(0), F(0)], [F(0), F(0)]]
        expected_d = [[F(0), F(0)], [F(0), F(0)]]
        for i, j, markers, mass in atom_words(1, self.constants):
            expected_i[i][j] += mass * (1 + 3*markers[0]) ** 2
            expected_d[i][j] += mass * max(1, markers[0]) ** 2
        np.testing.assert_allclose(2**bank['information'], np.asarray(expected_i, dtype=float), rtol=2e-14)
        np.testing.assert_allclose(2**bank['constrained'], np.asarray(expected_d, dtype=float), rtol=2e-14)
        self.assertGreater(bank['constrained'][0, 1], bank['empty'][0, 1])
        products = gate.ordered_gap_products(bank['information'], bank['constrained'], 3)
        i, d = 2**bank['information'], 2**bank['constrained']
        np.testing.assert_allclose(2**products[0b101], d @ i, rtol=2e-14)
        np.testing.assert_allclose(2**products[0b110], i @ d, rtol=2e-14)
        self.assertFalse(np.allclose(i @ d, d @ i))

    def test_exact_cover_transport_preserves_labelled_coordinates(self):
        for size, arity in ((2, 1), (4, 2), (8, 4), (16, 8)):
            rows = gate.row_geometry(size, arity)[0]
            covers = gate.information_covers(size, arity)
            self.assertEqual(len(covers), len(rows))
            indices = range(len(rows)) if size <= 8 else (0, 1, 11, 31, 511, 6000, len(rows)-1)
            for index in indices:
                direct = information_cover(rows[index], bits=size.bit_length()-1)
                self.assertEqual(covers[index], tuple(direct[i] for i in rows[index]))
        covers, values, parents, labels, codes = gate.cover_gap_geometry()
        rows = gate.row_geometry()[0]
        words = [()]
        for index in range(1, len(parents)):
            words.append(words[parents[index]] + (values[labels[index]],))
        for index in (0, 1, 31, 511, 6000, len(rows)-1):
            outside = [i for i in range(16) if i not in rows[index]]
            restored, previous = {}, -1
            for k, following in enumerate((*outside, 16)):
                word = words[int(codes[index, k])]
                self.assertEqual(len(word), following-previous-1)
                restored.update(zip(range(previous+1, following), word))
                previous = following
            self.assertEqual(tuple(restored[i] for i in rows[index]), covers[index])
        self.assertFalse(codes.flags.writeable)

    def test_cover_dp_matches_exact_labelled_atoms(self):
        counts = (0, 2, 1, 3, 4)
        logs = np.array([-np.inf] + [log2(x) for x in counts[1:]])
        result = gate.nonquiet_column_matrices({3: logs, 4: logs}, 4, self.measure,
            size=8, arity=4, q=4, batch_size=7, information_cover=True)
        self.assertTrue(result['information_cover'])
        for index in (0, 1, 8, len(result['rowsets'])-1):
            expected = exact_cover_column(result['rowsets'][index],
                result['information_cover_exponents'][index], counts, 4, self.constants, size=8, q=4)
            np.testing.assert_allclose(2**result['column_log2_matrices'][index],
                                      np.asarray(expected, dtype=float), rtol=9e-14)
        bank = result['gap_log2_matrices']
        zero = bank['cover_values'].index(F(0))
        np.testing.assert_array_equal(bank['cover_matrices'][zero], bank['constrained'])
        self.assertGreater(bank['cover_matrices'][zero, 0, 1], bank['empty'][0, 1])

    def test_cover_matrices_and_noncommuting_prefixes(self):
        values = (F(0), F(1, 4), F(1))
        bank = gate.weighted_cover_matrices(self.measure, 4, values, q=4)
        for c, matrix in zip(values, bank['cover_matrices']):
            exact = [[F(0), F(0)], [F(0), F(0)]]
            for i, j, markers, mass in atom_words(1, self.constants):
                t = markers[0]
                value = max(1, t)**4 if not c else (1+3*t**int(1/c))**int(4*c)
                exact[i][j] += mass*value
            np.testing.assert_allclose(2**matrix, np.asarray(exact, dtype=float), rtol=3e-14)
        parents, labels = (-1, 0, 1, 0, 3), (-1, 0, 1, 1, 0)
        products = gate.ordered_cover_products(bank['cover_matrices'], parents, labels)
        d, c = 2**bank['cover_matrices'][:2]
        np.testing.assert_allclose(2**products[2], d@c, rtol=3e-14)
        np.testing.assert_allclose(2**products[4], c@d, rtol=3e-14)
        self.assertFalse(np.allclose(c@d, d@c))

    def test_nonquiet_dp_matches_exact_labelled_atoms(self):
        counts = (0, 2, 1, 3, 4)
        logs = np.array([-np.inf] + [log2(x) for x in counts[1:]])
        result = gate.nonquiet_column_matrices({3: logs, 4: logs}, 2, self.measure,
            size=8, arity=4, q=4, batch_size=7)
        for index in (0, 1, 8, len(result['rowsets'])-1):
            expected = exact_nonquiet_column(result['rowsets'][index], result['information_sets'][index],
                counts, 2, self.constants, size=8, q=4)
            np.testing.assert_allclose(2**result['column_log2_matrices'][index],
                                      np.asarray(expected, dtype=float), rtol=7e-14)

    def test_rank_specific_coefficients_and_batching(self):
        logs3 = np.array([-np.inf, 0., 1., log2(3), 2.])
        logs4 = logs3 + 1
        first = gate.nonquiet_column_matrices({3: logs3, 4: logs4}, 2, self.measure,
            size=8, arity=4, q=4, batch_size=1)
        second = gate.nonquiet_column_matrices({3: logs3, 4: logs4}, 2, self.measure,
            size=8, arity=4, q=4, batch_size=13)
        np.testing.assert_array_equal(first['column_log2_matrices'], second['column_log2_matrices'])
        base = gate.nonquiet_column_matrices({3: logs3, 4: logs3}, 2, self.measure,
            size=8, arity=4, q=4, batch_size=13)
        for rank in (3, 4):
            chosen = first['outside_ranks'] == rank
            np.testing.assert_allclose(first['column_log2_matrices'][chosen] - base['column_log2_matrices'][chosen],
                                       0 if rank == 3 else 2., atol=3e-14)

    def test_wholeword_fixed_rowset_minkowski_and_young(self):
        activity, lambdas = F(3, 2), (F(3, 4), F(1, 4))
        result = gate.contract([0., 0.], {1: np.array([-np.inf, 1.])}, self.measure,
            2, log2(activity), size=2, columns=2, q=4, lambdas=lambdas, batch_size=1)
        # Direct common positive functional across both columns, with the
        # fixed source label inside its complete-word function.
        actual = F(0)
        for initial, end, markers, mass in atom_words(4, self.constants):
            if initial:
                continue
            fixed_sum = F(0)
            for source in result['rowsets']:
                i = source[0]
                value = F(1)
                for column in range(2):
                    t = markers[2*column:2*column+2]
                    quiet = 1 + activity*t[i]
                    nonquiet = activity*2*t[1-i]*max(1, t[i])
                    value *= quiet + nonquiet
                fixed_sum += value
            actual += mass * fixed_sum**2
        self.assertGreaterEqual(result['combined_word_log2']+1e-12, log2(actual))
        quiet, nonquiet = (result['branch_column_log2_matrices'][key] for key in ('quiet', 'nonquiet'))
        expected = np.logaddexp2(quiet-log2(lambdas[0]), nonquiet-log2(lambdas[1])+2*log2(activity))
        np.testing.assert_allclose(expected, result['column_log2_matrices'], rtol=2e-14, atol=2e-14)
        words = gate.log_matrix_power(expected, 2)
        endpoints = np.logaddexp2.reduce(words[:, 0, :], axis=1)
        np.testing.assert_allclose(endpoints, result['per_rowset_word_log2'], rtol=2e-14, atol=2e-14)
        self.assertFalse(result['extra_rowset_factor'])

    def test_production_coefficients_and_wrapper_scope_without_run(self):
        quiet, outside = gate.coefficient_tables()
        self.assertEqual(quiet, (1, 0, 0, 0, 255, 0, 64770, 16322040, 4162573845))
        self.assertEqual(tuple(rank for rank, _ in outside), (4, 5))
        with patch.object(gate, 'contract', return_value={}) as contraction:
            result = gate.full_row_bounds(self.measure, 18, 1.)
        self.assertEqual(contraction.call_count, 1)
        self.assertEqual(contraction.call_args.kwargs['columns'], 16)
        self.assertEqual(result['metadata']['intended_inner'], 'original-one-state')
        self.assertFalse(result['metadata']['global_prefix_included'])
        self.assertFalse(result['metadata']['exact_sparse_restoration'])
        self.assertFalse(result['metadata']['whole_code_certificate'])
        self.assertFalse(contraction.call_args.kwargs['information_cover'])
        with patch.object(gate, 'contract', return_value={}) as contraction:
            result = gate.full_row_bounds(self.measure, 18, 1., information_cover=True)
        self.assertTrue(contraction.call_args.kwargs['information_cover'])
        self.assertTrue(result['metadata']['information_cover'])
        self.assertFalse(result['metadata']['lexicographic_information_sets'])

    def test_guards_and_budget(self):
        with self.assertRaises(ValueError):
            gate.checked_lambdas((F(1), F(1)))
        with self.assertRaises(ValueError):
            gate.row_geometry(15, 8)
        with self.assertRaises(ValueError):
            gate.nonquiet_column_matrices({4: [0.]*9}, 18, self.measure)
        with self.assertRaises(ValueError):
            gate.contract([0., 0.], {1: [-np.inf, 0.]}, self.measure, 2, -1., size=2)
        with self.assertRaises(ValueError):
            gate.weighted_cover_matrices(self.measure, 2, (F(-1),))
        def stop():
            raise TimeoutError('toy budget')
        with self.assertRaises(TimeoutError):
            gate.full_row_bounds(self.measure, 18, 0., check_time=stop)


if __name__ == '__main__':
    main()
