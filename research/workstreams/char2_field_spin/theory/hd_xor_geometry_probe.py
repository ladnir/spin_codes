"""Bounded FLOAT diagnostics for fixed constituent-pair geometry.

No certificate or encoder failure is claimed. The small Hölder screen
enumerates only r<=3, with exact group-multiset multiplicities. Repeated-pair
controls use only 28 geometries at any supported r.
"""
from collections import Counter
from itertools import combinations_with_replacement
from math import comb, factorial, log2

import numpy as np

from hd_xor_column_potential import PAIRS, category_costs, multiset_column_matrices
from hd_xor_joint_columns import (joint_log_count,
    occupancy_envelope_caps_float, resource_envelope_float)
from hd_xor_multigroup import LIMIT, Q, REGION_LENGTH, region_operators


def multiset_multiplicity(values):
    result = factorial(len(values))
    for count in Counter(values).values():
        result //= factorial(count)
    return result


def assignment_multiplicities(single, double):
    """Weights for independent uniform singleton rows and unordered pairs."""
    if single < 0 or double < 0 or single+double > 3:
        raise ValueError('This exact diagnostic is capped at three source groups')
    singles = [(rows, multiset_multiplicity(rows))
               for rows in combinations_with_replacement(range(8), single)]
    doubles = [(edges, multiset_multiplicity(edges))
               for edges in combinations_with_replacement(range(28), double)]
    weights = [left*right for rows, left in singles for edges, right in doubles]
    denominator = 8**single*28**double
    if sum(weights) != denominator:
        raise ArithmeticError('The ordered geometry assignments were not counted exactly')
    return weights, denominator


def holder_geometry_screen(groups, z, log_potential):
    """Hölder order 32 charges geometry shared by all 32 source columns."""
    if not 1 <= groups <= 3:
        raise ValueError('Only the bounded one-to-three-group diagnostic is implemented')
    matrices = multiset_column_matrices(z, groups)
    maxima, _ = category_costs(matrices, log_potential)
    moments = {}
    for (single, double), values in matrices.items():
        contractions = np.logaddexp2(values[:, :, 0], values[:, :, 1]+log_potential)
        contractions[:, 1] -= log_potential
        costs = np.max(contractions, axis=1)
        counts, denominator = assignment_multiplicities(single, double)
        logs = np.array([log2(count)-log2(denominator) for count in counts])
        moments[single, double] = float(np.logaddexp2.reduce(32*costs+logs)/32)
    prefix = (log2(comb(REGION_LENGTH, groups))+groups*log2(28)-log2(Q-1)
              -LIMIT*log2(z)-log_potential)
    return dict(active_groups=groups, mark=z, log2_potential=log_potential,
                maximum_category_log2_bound=prefix+joint_log_count(maxima, groups),
                holder_category_log2_bound=prefix+joint_log_count(moments, groups),
                category_maxima={str(key): value for key, value in maxima.items()},
                category_holder_caps={str(key): value for key, value in moments.items()},
                holder_order=32, exact_multiset_multiplicities=True,
                geometry_is_union_index_not_setup_randomness=True,
                directed=False, whole_code_certificate=False)


def repeated_pair_matrices(regions, groups):
    """All groups use one fixed pair; retain every cancellation-count branch."""
    result = []
    for pair in PAIRS:
        mixture = np.full((2, 2), -np.inf)
        for cancelled in range(groups+1):
            word = np.array([[0., -np.inf], [-np.inf, 0.]])
            for index in range(8):
                matrix = regions[groups if pair >> index & 1 else groups-cancelled]
                word = np.logaddexp2(word[:, 0, None]+matrix[None, 0, :],
                                    word[:, 1, None]+matrix[None, 1, :])
            weight = (log2(comb(groups, cancelled))+cancelled*log2(1/(Q-1))
                      +(groups-cancelled)*log2((Q-2)/(Q-1)))
            mixture = np.logaddexp2(mixture, word+weight)
        result.append(mixture)
    return np.array(result)


def repeated_pair_diagnostic(groups=17, z=.996548217, log_potential=-63.5):
    if not 1 <= groups <= 32:
        raise ValueError('Only one to thirty-two groups are supported')
    regions = region_operators(z, groups)
    matrices = repeated_pair_matrices(regions, groups)
    contractions = np.logaddexp2(matrices[:, :, 0], matrices[:, :, 1]+log_potential)
    contractions[:, 1] -= log_potential
    maxima = np.max(contractions, axis=0)
    pairs = [PAIRS[index] for index in np.argmax(contractions, axis=0)]
    caps = occupancy_envelope_caps_float(groups, z, log_potential)
    terms = []
    for cancelled in range(groups+1):
        action = resource_envelope_float(regions, log_potential, groups-cancelled, groups)[2*cancelled]
        weight = (log2(comb(groups, cancelled))+cancelled*log2(1/(Q-1))
                  +(groups-cancelled)*log2((Q-2)/(Q-1)))
        terms.append(action+weight-np.array([0., log_potential]))
    return dict(active_groups=groups, mark=z, log2_potential=log_potential,
                repeated_pair_row_caps=maxima.tolist(), maximizing_pair_masks=pairs,
                resource_envelope_pure_double_cap=caps[0, groups],
                dominant_cancellation_counts=np.argmax(terms, axis=0).tolist(),
                repeated_pair_geometries=28, full_geometry_count=28**groups,
                log2_geometry_count_difference=(groups-1)*log2(28),
                directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--holder', action='store_true')
    parser.add_argument('--groups', type=int, default=17)
    parser.add_argument('--mark', type=float, default=.996548217)
    parser.add_argument('--potential', type=float, default=-63.5)
    args = parser.parse_args()
    result = (holder_geometry_screen(args.groups, args.mark, args.potential) if args.holder else
              repeated_pair_diagnostic(args.groups, args.mark, args.potential))
    print(json.dumps(result), flush=True)
