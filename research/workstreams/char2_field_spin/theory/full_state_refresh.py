"""Local rank-action backup with a fresh uniform-entry four-state map.

Singular maps are allowed.  These conditional local upper actions are not
a whole-code certificate or a native implementation.  The preferred
construction remains the diagonal-state refresh while its mask proof runs.
"""

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb

from cauchy_joint_return import matrix_rank, packet_forms, projective_strata


def canonical_factor(tau, rank=4):
    """Cost of splitting by the first nonzero independent state form."""
    if not 1 <= tau <= 255 or not 0 <= rank <= 4:
        raise ValueError('require 1<=tau<=255 and state-form rank0..4')
    if isinstance(tau, int):
        tau = Fraction(tau)
    return sum((tau**(-j) for j in range(rank)), tau*0)


@lru_cache(maxsize=None)
def return_terms(phase, inputs):
    """Integer (selected-zero count, input-rank charge, count) profile.

    The projective arrangement is ONLY a feasibility oracle.  Each feasible
    selected-zero set is counted once, irrespective of its number of words
    or incoming state-support patterns.
    """
    injection, output = packet_forms(phase, tuple(inputs))
    feasible = set()
    for _, zero_mask, _ in projective_strata(injection, output):
        subset = zero_mask
        while True:
            feasible.add(subset)
            if not subset:
                break
            subset = (subset-1) & zero_mask
    profile = Counter()
    for mask in feasible:
        rank = matrix_rank([row for j, row in enumerate(output) if mask >> j & 1], len(inputs))
        if rank >= len(inputs):
            raise ArithmeticError('feasible selected-zero system kills all inputs')
        profile[mask.bit_count(), rank] += 1
    return tuple((zeros, rank, count) for (zeros, rank), count in sorted(profile.items()))


@lru_cache(maxsize=None)
def return_profiles(weight, phase=None, average_inputs=False):
    if weight not in (1, 2, 3) or (phase is not None and phase not in range(4)):
        raise ValueError('require weight1..3 and phase0..3 or None')
    profiles = set()
    for current in range(4) if phase is None else (phase,):
        accumulated = Counter()
        for inputs in combinations(range(16), weight):
            profile = return_terms(current, inputs)
            if average_inputs:
                for zeros, rank, count in profile:
                    accumulated[zeros, rank] += count
            else:
                profiles.add(profile)
        if average_inputs:
            profiles.add(tuple((zeros, rank, Fraction(count, comb(16, weight)))
                               for (zeros, rank), count in sorted(accumulated.items())))
    return tuple(sorted(profiles))


def return_rank_action(weight, y, tau, *, state_rank=4, phase=None, average_inputs=False):
    """Nonzero incoming witness state -> zero next state, uniformly in history.

    The input-diagonal zero equations are imposed first.  The state-form
    rank can decrease there, so ``state_rank`` is an upper bound, not a
    claim that the rank is unchanged.  Four fresh matrix rows are then
    consumed after the canonical first-nonzero split.  The safe common
    parameter range is 1<=tau<=255.
    """
    if not 0 < y <= 1 or weight not in (0, 1, 2, 3):
        raise ValueError('require mark0<y<=1 and weight0..3')
    if isinstance(y, int):
        y = Fraction(y)
    if isinstance(tau, int):
        tau = Fraction(tau)
    factor = canonical_factor(tau, state_rank)/tau**4
    if weight == 0:
        return factor  # On a zero-input return every output is zero.
    return factor*max(sum(count*y**(16-z)*(1-y)**z/tau**rank
                          for z, rank, count in profile)
                       for profile in return_profiles(weight, phase, average_inputs))


def marginal_rank_action(weight, y, tau, *, state_rank=4):
    """Uniform positive zero-set marginal, including next-state returns.

    rank([A_Z,E_ZI])=i+min(z,4), where i selected zeros lie on active
    input positions and z on zero-input positions.  Do not multiply this
    action by return_rank_action: they expose the same packet randomness.
    """
    if weight not in range(17) or not 0 < y <= 1:
        raise ValueError('invalid input weight or output mark')
    if isinstance(y, int):
        y = Fraction(y)
    if isinstance(tau, int):
        tau = Fraction(tau)
    return canonical_factor(tau, state_rank)*sum(
        comb(weight, i)*comb(16-weight, z)*y**(16-i-z)*(1-y)**(i+z)
        /tau**(i+min(z, 4))
        for i in range(weight+1) for z in range(17-weight))


def empty_nonreturn_action(y, tau, *, state_rank=4):
    """Empty packet, nonzero incoming and nonzero outgoing witness state."""
    if not 0 < y <= 1:
        raise ValueError('invalid output mark')
    if isinstance(y, int):
        y = Fraction(y)
    if isinstance(tau, int):
        tau = Fraction(tau)
    return canonical_factor(tau, state_rank)*sum(
        comb(16, z)*y**(16-z)*(1-y)**z/tau**z for z in range(4))


if __name__ == '__main__':
    import json
    y, tau = Fraction(9997, 10000), Fraction(35)
    for weight in range(4):
        value = return_rank_action(weight, y, tau, average_inputs=True)
        print(json.dumps(dict(weight=weight, y=str(y), tau=str(tau),
                              action=str(value), approximate=float(value),
                              type='local-rank-upper-action', whole_code=False)))
