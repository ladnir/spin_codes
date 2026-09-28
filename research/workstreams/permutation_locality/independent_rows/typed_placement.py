"""Positive normalized placement of two labeled group types.

The caller supplies B[j1,j2], an epoch operator conditional on j1 type-one
groups and j2 type-two groups occupying distinct slots. An epoch contains
`slots` positions. Selected groups may produce zero packets; that randomness
belongs inside B and must not be removed from the type counts.

For E epochs, the returned operator averages over all disjoint placements
of n1 and n2 selected groups in E*slots positions, in chronological order.
Equivalently it is the coefficient of z1^n1*z2^n2 in

  (sum_{j1+j2<=slots} multinomial(slots;j1,j2,slots-j1-j2)
                         z1^j1*z2^j2 B[j1,j2])**E,

divided by multinomial(E*slots;n1,n2,E*slots-n1-n2).

The recurrence normalizes after every epoch. Conditional on a total of
`a,b` groups in the first `e` epochs, the last epoch has type counts j1,j2
with probability

  C(a,j1) C(b,j2) C(e*slots-a-b,slots-j1-j2) / C(e*slots,slots).

These exact hypergeometric weights are nonnegative and sum to one. Thus no
large multinomial coefficient is accumulated inside the matrix recurrence.
Multiplication is previous_prefix @ current_epoch; matrices need not commute.
The default implementation uses outward Arb endpoints. For exact tests use
matrix=fmpq_mat and rounding=lambda value: value.
"""
from math import comb

from flint import arb, arb_mat


def upper_matrix(value):
    """Round every finite entry upward to an exact Arb endpoint."""
    for i in range(value.nrows()):
        for j in range(value.ncols()):
            assert value[i,j].is_finite()
    return arb_mat([[arb(value[i,j].upper()) for j in range(value.ncols())]
                    for i in range(value.nrows())])


def _validate(operators, counts, epochs, slots):
    if (not isinstance(epochs,int) or epochs < 0
            or not isinstance(slots,int) or slots < 1
            or len(counts) != 2 or any(not isinstance(x,int) or x < 0 for x in counts)
            or sum(counts) > epochs*slots):
        raise ValueError('invalid epoch count, slot count, or pair of type counts')
    if not operators:
        raise ValueError('at least one epoch operator is required')
    first = next(iter(operators.values()))
    size = first.nrows()
    if size < 1 or first.ncols() != size:
        raise ValueError('nonempty square operators required')
    for key,value in operators.items():
        if (len(key) != 2 or any(not isinstance(x,int) or x < 0 for x in key)
                or sum(key) > slots):
            raise ValueError('invalid local type counts')
        if value.nrows() != size or value.ncols() != size:
            raise ValueError('all operators must have the same square dimensions')
        if any(not value[i,j] >= 0 for i in range(size) for j in range(size)):
            raise ValueError('nonnegative epoch operators required')
    return size,first[0,0]*0+1


def typed_placement(operators, counts, *, epochs=64, slots=32,
                    matrix=arb_mat, rounding=upper_matrix, return_grid=False):
    """Average two-type placement, using supplied conditional epoch matrices.

    counts=(n1,n2) fixes the desired type occupancy. With return_grid=True,
    return every feasible pair bounded componentwise by counts; otherwise
    prune prefix counts that cannot reach the requested final pair.

    A missing epoch operator raises ValueError only if the recurrence uses
    it. This permits sparse operator dictionaries at forced occupancies.
    """
    counts = tuple(counts)
    size,one = _validate(operators,counts,epochs,slots)
    n1,n2 = counts
    current = {(0,0):matrix([[int(i==j) for j in range(size)] for i in range(size)])}
    first_choices = [[comb(a,j) for j in range(min(a,slots)+1)] for a in range(n1+1)]
    second_choices = [[comb(b,j) for j in range(min(b,slots)+1)] for b in range(n2+1)]
    for e in range(1,epochs+1):
        positions = e*slots
        previous_positions = positions-slots
        future_positions = (epochs-e)*slots
        denominator = comb(positions,slots)
        inverse = one/denominator
        empty_choices = {}
        updated = {}
        for a in range(min(n1,positions)+1):
            for b in range(min(n2,positions-a)+1):
                if not return_grid and n1+n2-a-b > future_positions:
                    continue
                total = a+b
                if total not in empty_choices:
                    empty = positions-total
                    empty_choices[total] = [comb(empty,j) if j <= empty else 0
                                            for j in range(slots+1)]
                choices = empty_choices[total]
                value = matrix(size,size)
                normalization = 0
                for j1,c1 in enumerate(first_choices[a]):
                    for j2 in range(min(b,slots-j1)+1):
                        numerator = c1*second_choices[b][j2]*choices[slots-j1-j2]
                        if not numerator:
                            continue
                        previous = a-j1,b-j2
                        if sum(previous) > previous_positions:
                            raise AssertionError('hypergeometric term exceeds previous capacity')
                        if previous not in current:
                            raise AssertionError('required normalized prefix was pruned')
                        local = j1,j2
                        if local not in operators:
                            raise ValueError(f'missing epoch operator for type counts {local}')
                        # Form the small positive weight before multiplying
                        # matrices, not an unnormalized polynomial coefficient.
                        weight = numerator*inverse
                        value += (current[previous]*operators[local])*weight
                        normalization += numerator
                assert normalization == denominator
                updated[a,b] = rounding(value)
        current = updated
    return current if return_grid else current[counts]
