"""Exact three-packet return control over GF(4); no regional proof bank.

Feedbacks c1,c2,c3 are fixed nonzero vectors, independent of fresh sampled
refreshes. Before-refresh coupling is (V,M(U+V+c)); after-refresh coupling
is (V,M(U+c)+V). Both start at (0,0). This module enumerates the actual
D2 H D1 key law, not a uniform-GL surrogate. Emission weights are absent.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product


def multiply(a, b):
    """GF(4) in the polynomial basis with x**2+x+1=0."""
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


def apply(matrix, vector):
    x, y = vector & 3, vector >> 2
    first = multiply(matrix[0], x)^multiply(matrix[1], y)
    second = multiply(matrix[2], x)^multiply(matrix[3], y)
    return first+4*second


@lru_cache(None)
def refresh_actions():
    """All 3**4 equally likely diagonal keys, retaining duplicate matrices.

    H=[[1,1],[1,2]] is MDS and invertible over GF(4). Vectors are encoded
    as x+4*y; coordinatewise addition is bitwise XOR.
    """
    actions = []
    for d0, d1, e0, e1 in product(range(1, 4), repeat=4):
        matrix = (multiply(e0, d0), multiply(e0, d1),
                  multiply(e1, d0), multiply(e1, multiply(2, d1)))
        actions.append(tuple(apply(matrix, x) for x in range(16)))
    return tuple(actions)


def transition(state, feedback, action, order):
    u, v = state
    if order == 'before':
        return v, action[u^v^feedback]
    if order == 'after':
        return v, action[u^feedback]^v
    raise ValueError('Explicit before/after-refresh coupling order required')


def advance(distribution, feedback, actions, order):
    result = Counter()
    for state, count in distribution.items():
        for action in actions:
            result[transition(state, feedback, action, order)] += count
    return result


@lru_cache(None)
def kernel_counts():
    actions = refresh_actions()
    return tuple(tuple(sum(action[x] == y for action in actions) for y in range(16))
                 for x in range(16))


@lru_cache(None)
def direct_three_return_counts(order):
    """Exhaust all nonzero feedback triples and refresh keys via exact counting.

    Two steps retain the full joint state distribution. At step three,
    states with V!=0 cannot return, because the next first coordinate is V.
    For the remaining states every third refresh key is still enumerated.
    The denominator is 81**3 for EACH fixed feedback triple.
    """
    if order not in ('before', 'after'):
        raise ValueError('Explicit coupling order required')
    actions = refresh_actions()
    result = {}
    for c1 in range(1, 16):
        first = advance({(0, 0): 1}, c1, actions, order)
        for c2 in range(1, 16):
            second = advance(first, c2, actions, order)
            candidates = [(state, count) for state, count in second.items() if state[1] == 0]
            for c3 in range(1, 16):
                result[c1, c2, c3] = sum(count
                    for state, count in candidates for action in actions
                    if transition(state, c3, action, order) == (0, 0))
    return result


def summary():
    denominator = len(refresh_actions())**3
    before, after = direct_three_return_counts('before'), direct_three_return_counts('after')
    maximum_before, maximum_after = max(before.values()), max(after.values())
    return dict(field_order=4, state_dimension=2, nonzero_states=15,
        diagonal_key_tuples=len(refresh_actions()), fixed_nonzero_feedback_triples=15**3,
        per_triple_key_denominator=denominator, beta=str(F(15, 9)),
        maximum_single_atom=str(max(F(n, 81) for row in kernel_counts()[1:] for n in row[1:])),
        before_maximum=str(F(maximum_before, denominator)),
        after_maximum=str(F(maximum_after, denominator)),
        before_maximizer=next(key for key, value in before.items() if value == maximum_before),
        after_maximizer=next(key for key, value in after.items() if value == maximum_after),
        before_uniform_feedback_average=str(F(sum(before.values()), denominator*15**3)),
        after_uniform_feedback_average=str(F(sum(after.values()), denominator*15**3)),
        fixed_feedback_conditioning=True, emission_weights_included=False,
        regional_certificate=False, whole_code_certificate=False)


if __name__ == '__main__':
    import json
    print(json.dumps(summary(), indent=2))
