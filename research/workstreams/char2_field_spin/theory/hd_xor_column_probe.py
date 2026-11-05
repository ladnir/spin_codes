"""Exact support probes for H=I+11^T, for even group width in characteristic two.

Research probe only. No whole-code certificate or native performance claim.
The row encoders are independent copies of one MDS[n,k] encoder. Independent
nonzero diagonals precede H in each column. An optional independent output
diagonal preserves the support formulas and randomizes values given support.
"""
from fractions import Fraction
from math import comb, log2


def mds_prescribed_support(q, n, k, weight):
    if not 0 <= weight <= n or not 0 < k < n:
        raise ValueError('Invalid MDS parameters')
    if weight == 0:
        return 1
    distance = n-k+1
    if weight < distance:
        return 0
    return (q-1)*sum((-1)**j * comb(weight-1, j) * q**(weight-distance-j)
                     for j in range(weight-distance+1))


def xor_column(values):
    if not values or len(values) % 2:
        raise ValueError('Even nonempty column required')
    total = 0
    for value in values:
        total ^= value
    return tuple(value ^ total for value in values)


def column_weight_counts(q, group, active):
    """Counts among (q-1)^active nonzero amplitudes on a fixed input support."""
    if q < 2 or q & (q-1) or group < 2 or group % 2 or not 0 <= active <= group:
        raise ValueError('Power-of-two field and even group required')
    counts = [0]*(group+1)
    if active == 0:
        counts[0] = 1
        return counts
    zero_sum = ((q-1)**active + (q-1)*(-1)**active)//q
    counts[active] += zero_sum
    for equal_sum in range(active+1):
        remaining = active-equal_sum
        if remaining == 0:
            residual_count = equal_sum % 2
        else:
            residual_count = ((q-2)**remaining + (q//2-1)*(-2)**remaining)//q
        counts[group-equal_sum] += comb(active, equal_sum)*(q-1)*residual_count
    assert sum(counts) == (q-1)**active
    return counts


def column_support_counts(q, group, input_mask):
    """Exact output-mask counts for nonzero amplitudes on input_mask."""
    if q < 2 or q & (q-1) or group < 2 or group % 2 or not 0 <= input_mask < 1 << group:
        raise ValueError('Invalid field, group, or mask')
    counts = [0]*(1 << group)
    active = input_mask.bit_count()
    if not active:
        counts[0] = 1
        return counts
    counts[input_mask] += ((q-1)**active+(q-1)*(-1)**active)//q
    equal_mask = input_mask
    while True:
        equal_sum = equal_mask.bit_count()
        remaining = active-equal_sum
        residual_count = (equal_sum % 2 if remaining == 0 else
                          ((q-2)**remaining+(q//2-1)*(-2)**remaining)//q)
        counts[((1 << group)-1) ^ equal_mask] += (q-1)*residual_count
        if not equal_mask:
            break
        equal_mask = (equal_mask-1) & input_mask
    assert sum(counts) == (q-1)**active
    return counts


def uniform_column_support_counts(q, group, active_rows):
    """Counts for H U when U is uniform on the coordinate space F^A."""
    result = [0]*(1 << group)
    support = active_rows
    while True:
        result = [a+b for a, b in zip(result, column_support_counts(q, group, support))]
        if not support:
            break
        support = (support-1) & active_rows
    assert sum(result) == q**active_rows.bit_count()
    return result


def membership_mass(q, n, k, output_columns):
    """Exact expected membership for H D, excluding the zero message."""
    if len(output_columns) != n or not output_columns:
        raise ValueError('Wrong column count')
    group = len(output_columns[0])
    if any(len(column) != group for column in output_columns):
        raise ValueError('Inconsistent group widths')
    unmixed = tuple(xor_column(column) for column in output_columns)
    weights = tuple(sum(bool(column[i]) for column in unmixed) for i in range(group))
    if not any(weights):
        return Fraction(0)
    result = Fraction(1)
    for weight in weights:
        result *= Fraction(mds_prescribed_support(q, n, k, weight), (q-1)**weight)
    return result


def two_active_row_weight_enumerator(q=256, n=32, k=16, group=8):
    """FLOAT expected outer counts for two fixed nonzero constituent rows.

    Both rows must be nonzero. The expectation is over the input diagonals;
    an independent output diagonal does not change this weight enumerator.
    This computes a moment, not the distribution of minimum code distance.
    """
    result = [0.0]*(n*group+1)
    for first in range(n-k+1, n+1):
        first_count = mds_prescribed_support(q, n, k, first)
        for second in range(n-k+1, n+1):
            second_count = mds_prescribed_support(q, n, k, second)
            for overlap in range(max(0, first+second-n), min(first, second)+1):
                support_count = comb(n, first)*comb(first, overlap)*comb(n-first, second-overlap)
                messages = support_count*first_count*second_count
                singleton_weight = (group-1)*(first+second-2*overlap)
                for cancellations in range(overlap+1):
                    weight = singleton_weight + group*(overlap-cancellations) + 2*cancellations
                    probability = (comb(overlap, cancellations)*(q-2)**(overlap-cancellations)
                                   / (q-1)**overlap)
                    result[weight] += messages*probability
    return result


def main():
    q, group = 256, 8
    for active in range(1, group+1):
        counts = column_weight_counts(q, group, active)
        probabilities = [(weight, count/(q-1)**active) for weight, count in enumerate(counts) if count]
        print('column active', active, 'weight probabilities', probabilities)
    enumerator = two_active_row_weight_enumerator()
    print('two fixed active rows, FLOAT expected projective outer counts')
    for threshold in (34, 40, 46, 52, 58, 64, 70, 76, 82, 89, 100, 119, 128):
        count = sum(enumerator[:threshold+1])/(q-1)
        print('weight <=', threshold, 'log2 expectation', log2(count) if count else '-inf')
    print('total log2 messages', log2(sum(enumerator)), 'expected', 2*log2(q**16-1))


if __name__ == '__main__':
    main()
