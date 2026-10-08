"""Positive zero-set expansion of the common-rank existence potential.

This is NOT a fixed-word first moment.  A witness contributes the sum of
weights over all subsets of its actual zero outputs; that sum is y**wt(Y).
The source multiplier remains tau**(D-1), not the number of field words.
The construction is GRS[32,22] with Cauchy emission, eta=2, independent
nonzero input diagonals, and four fresh state diagonals in every packet.
"""

from fractions import Fraction
from functools import lru_cache
from math import comb

from diagonal_refresh_rank import rank_potential_bound


@lru_cache(maxsize=None)
def root_factors(tau, q=256, state=4):
    """Exact for rational tau; uniform over every conditional state image."""
    return tuple(max(rank_potential_bound(r, z, tau, q, state)
                     for r in range(1, state+1)) for z in range(state))


def local_rank_transfer(weight, y, tau, q=256, width=16, state=4,
                        joint_one=True):
    """Positive 2x2 matrix for one FIXED input support of given weight.

    Scalar arguments can be Fraction for exact replay, or float for screens.
    Rows/columns mean actual zero/nonzero witness state. The current state
    diagonal is consumed only in the root factor, except for the separate
    joint one-input-return bound which is combined by taking a minimum.
    """
    if not (0 <= weight <= width and 0 < y <= 1 and 1 <= tau <= q-1):
        raise ValueError('invalid weight or tilts')
    if isinstance(tau, int):
        tau = Fraction(tau)
    if isinstance(y, int):
        y = Fraction(y)
    gamma = root_factors(tau, q, state)
    zero = y*0
    matrix = [[zero, zero], [zero, zero]]
    if weight == 0:
        matrix[0][0] = y**0
    else:
        matrix[0][1] = y**weight
        if weight >= state+1:
            matrix[0][0] = y**weight/tau**state
    for i in range(weight+1):
        for z in range(min(state-1, width-weight)+1):
            selected = i+z
            term = (comb(weight, i)*comb(width-weight, z)
                    *y**(width-selected)*(1-y)**selected*gamma[z]/tau**i)
            matrix[1][1] += term
            if i < weight:
                matrix[1][0] += term/tau**min(weight-i, state)
    if joint_one and weight == 1:
        # Return implies all state coordinates are nonzero, because the
        # injection column is full-support. Reveal the input diagonal FIRST,
        # then use all state diagonals as alive pivots. No root factor is
        # multiplied into this alternative whole-packet bound.
        matrix[1][0] = min(matrix[1][0], y**width/tau**state)
    return matrix


def transfer(x, y, tau, q=256, width=16, state=4, joint_one=True):
    """Sum all input supports with mark x, for global route/rank bounds."""
    if x <= 0:
        raise ValueError('positive input mark required')
    zero = y*0
    matrix = [[zero, zero], [zero, zero]]
    for weight in range(width+1):
        local = local_rank_transfer(weight, y, tau, q, width, state, joint_one)
        multiplicity = comb(width, weight)*x**weight
        for r in range(2):
            for k in range(2):
                matrix[r][k] += multiplicity*local[r][k]
    return matrix


if __name__ == '__main__':
    import json
    y, tau = Fraction(9, 10), Fraction(14)
    print(json.dumps({'y': str(y), 'tau': str(tau),
                      'local_zero_input': [[str(x) for x in row] for row in local_rank_transfer(0, y, tau)],
                      'local_one_input': [[str(x) for x in row] for row in local_rank_transfer(1, y, tau)],
                      'rank_existence_not_word_count': True}, indent=2))
