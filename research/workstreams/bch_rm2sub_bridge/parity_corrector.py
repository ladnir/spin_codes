"""A linear right inverse using three actual BCH rows of the unchanged encoder.

This selects messages after row permutations are fixed. It is not a new
encoding step or a modification of the code or setup distribution.
"""
import sys
import bridge as base
sys.path.insert(0, str(base.BCH/'code'))
from bch_quotient import generator_polynomial
from affine_wambach import gf_pow, ALPHA


def fixed_c_basis():
    generator = generator_polynomial(37)
    assert 255-(generator.bit_length()-1) == 131
    coordinates = [gf_pow(ALPHA, j) for j in range(255)]+[0]
    powers = [gf_pow(x, 37) for x in coordinates]
    pivots, basis = {}, []
    for j in range(131):
        word = generator << j
        word |= (word.bit_count() % 2) << 255
        syndrome = 0
        for i in range(256):
            if word >> i & 1:
                syndrome ^= powers[i]
        label = syndrome >> 5
        while label:
            pivot = label.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = label, word
                break
            old, old_word = pivots[pivot]
            label ^= old
            word ^= old_word
        if label == 0:
            basis.append(word)
    assert len(pivots) == 3 and len(basis) == 128
    assert all(word.bit_count() % 2 == 0 for word in basis)
    return basis


def permute(word, permutation):
    return sum(((word >> i) & 1) << j for i, j in enumerate(permutation))


def right_inverse(basis, permutations):
    assert len(basis) == 128 and len(permutations) == 3
    for pi in permutations:
        assert sorted(pi) == list(range(256))
    columns = [permute(word, pi) for pi in permutations for word in basis]
    pivots = {}
    for j, column in enumerate(columns):
        value, combination = column, 1 << j
        while value:
            pivot = value.bit_length()-1
            if pivot not in pivots:
                pivots[pivot] = value, combination
                break
            old, old_combination = pivots[pivot]
            value ^= old
            combination ^= old_combination
    if len(pivots) != 255:
        return len(pivots), None
    result = []
    for i in range(255):
        value, combination = (1 << i) | (1 << 255), 0
        while value:
            pivot = value.bit_length()-1
            old, old_combination = pivots[pivot]
            value ^= old
            combination ^= old_combination
        check = 0
        for j, column in enumerate(columns):
            if combination >> j & 1:
                check ^= column
        assert check == (1 << i) | (1 << 255)
        result.append(combination)
    return 255, result
