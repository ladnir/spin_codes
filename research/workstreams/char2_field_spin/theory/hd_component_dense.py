"""Dense tail upper and actual-code first-moment obstruction, not distance failure.

The obstruction is a rare regional routing event carrying a large low-weight
code subspace. It applies to the ideal independent-uniform routing ensemble.
Both placement formulas match hd_affine_two_state.cpp; no native code runs.
"""

from hashlib import sha256
from math import comb
from random import Random

from hd_affine_proof import admissible_tuples, inverse, mul
from hd_two_state_regional_strides import STRIDES

ROWS = 32768
PERMITTED = 19660
WINDOW = 13107
MASK_SHA256 = 'eab8a631796d80d096c309dff05dca4985e3947b732733e264d89126749cd436'
EXPECTED = {
    1: (11759,11766,11744,11815,11797,11940,11764,11812,11746,11773,11806,11808,11749,11754,11788,11772,
        11740,11754,11743,11808,11927,11743,11798,11748,11751,11693,11762,11738,11801,11835,11777,11694,
        11797,11771,11732,11812,11828,11774,11748,11776,11789,11734,11788,11828,11802,11795,11812,11784,
        11825,11812,11835,11812,11792,11872,11783,11810,11775,11736,11808,11773,11801,11867,11800,11778),
    4: (11833,11837,11809,11813,11832,11778,11808,11849,11825,11847,11800,11812,11817,11805,11831,11809,
        11820,11816,11797,11810,11810,11833,11775,11828,11850,11784,11802,11810,11843,11836,11795,11811,
        11809,11799,11820,11812,11867,11810,11811,11812,11842,11812,11834,11832,11805,11803,11824,11851,
        11810,11814,11806,11779,11778,11833,11775,11796,11835,11802,11815,11771,11784,11780,11815,11828),
}


def fixed_mask(rows=ROWS, permitted=PERMITTED):
    mask = bytearray(rows)
    for row in Random(0).sample(range(rows), permitted):
        mask[row] = 1
    if rows == ROWS and permitted == PERMITTED:
        if sha256(mask).hexdigest() != MASK_SHA256:
            raise ArithmeticError('The fixed witness set changed')
    return mask


def placement_order(rows, stride, group):
    if rows <= 0 or rows & (rows-1) or group <= 0 or rows % group or stride % 2 != 1:
        raise ValueError('Need power-of-two row geometry and odd strides')
    return tuple(group*((stride*(j//group)) % (rows//group))+j % group for j in range(rows))


def unit_positions(mask, strides=STRIDES, group=1):
    rows = len(mask)
    selected = []
    for stride in strides:
        order = placement_order(rows, stride, group)
        selected.append(tuple(j for j in range(rows-1) if mask[order[j]] and mask[order[j+1]]))
    return selected


def log_comb_arb(n, k):
    from flint import arb
    return arb(n+1).lgamma()-arb(k+1).lgamma()-arb(n-k+1).lgamma()


def obstruction_receipt(group, precision=160):
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        mask = fixed_mask()
        counts = tuple(map(len, unit_positions(mask, group=group)))
        if counts != EXPECTED[group] or max(counts) > WINDOW:
            raise ArithmeticError('Frozen unit-component receipt changed')
        regions, rho = len(STRIDES), 32
        variables = sum(counts)
        dimension = variables-2*regions-rho*PERMITTED
        if dimension <= 1:
            raise ArithmeticError('Need a positive actual-code dimension')
        log_probability = sum(log_comb_arb(WINDOW, h)-log_comb_arb(ROWS, h) for h in counts)/arb(2).log()
        # (q^d-1)/(q-1) > q^(d-2), using q=2^16 and d>=2.
        log_first_moment_lower = log_probability+16*(dimension-2)
        certified_bits = 420595 if group == 1 else 441457
        if not log_first_moment_lower > certified_bits:
            raise ArithmeticError('The first-moment lower bound was not certified')
        return dict(group=group, permitted_rows=PERMITTED, mask_sha256=MASK_SHA256,
                    unit_counts=counts, unit_variables=variables, post_exit_constraints=2*regions,
                    outer_check_constraints=rho*PERMITTED, actual_code_dimension_lower=dimension,
                    output_weight_upper=regions*WINDOW, forbidden_cutoff=(2*regions*ROWS)//5,
                    route_event_log2=str(log_probability),
                    actual_projective_first_moment_log2_lower=str(log_first_moment_lower),
                    permitted_class_log2_lower=str(log_first_moment_lower-1),
                    first_moment_strictly_exceeds_power_of_two=certified_bits,
                    rules_out_uncapped_first_moment=True, rejects_40_bit_distance=False,
                    precision=precision)
    finally:
        ctx.prec = previous


def dense_tail_upper(start=29943, precision=160):
    """Whole a>=start tail, using an information set and a geometric ratio."""
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        rows, b, rho, Q = ROWS, 64, 32, 65535
        cutoff = (2*b*rows)//5
        dimension = b*start
        if not 1 <= start <= rows or 2*cutoff > dimension:
            raise ValueError('This simple ball bound uses cutoff <= dimension/2')
        ratio = (arb(rows-start)/(start+1)*arb(Q)**(-rho)
                 * (arb(dimension+1)/(dimension-cutoff+1))**b)
        if not ratio < 1:
            raise ArithmeticError('Tail ratio is not below one')
        logarithm = (arb(cutoff).log()+log_comb_arb(rows, start)+log_comb_arb(dimension, cutoff)
                     +(cutoff-1-rho*start)*arb(Q).log()-(1-ratio).log())/arb(2).log()
        return dict(first_class=start, last_class=rows, number_of_classes=rows-start+1,
                    fraction=start/rows, entire_tail_log2_upper=str(logarithm),
                    entire_tail_below_2_to_minus_80=bool(logarithm < -80),
                    applies_to_exact_active_and_positive_permitted_classes=True,
                    independent_of_placement_and_post_setup=True,
                    whole_distance_certificate=False, precision=precision)
    finally:
        ctx.prec = previous


def nullspace(matrix, width):
    """Small GF4 RREF oracle; not a native implementation."""
    work = [list(row) for row in matrix]
    pivots = []
    for column in range(width):
        pivot = next((i for i in range(len(pivots), len(work)) if work[i][column]), None)
        if pivot is None:
            continue
        top = len(pivots)
        work[top], work[pivot] = work[pivot], work[top]
        scale = inverse(work[top][column])
        work[top] = [mul(scale, value) for value in work[top]]
        for i in range(len(work)):
            if i != top and work[i][column]:
                scale = work[i][column]
                work[i] = [left ^ mul(scale, right) for left, right in zip(work[i], work[top])]
        pivots.append(column)
    result = []
    for free in range(width):
        if free not in pivots:
            vector = [0]*width
            vector[free] = 1
            for i, pivot in enumerate(pivots):
                vector[pivot] = work[i][free]
            result.append(vector)
    return result


def small_actual_code_check(group=1, three_product=False):
    """Construct actual GF4 codewords on the routing event, with nonzero prefix weights."""
    rows, b, k, permitted, window = 64, 4, 2, 38, 25
    strides = (1, 3, 5, 7)
    mask = fixed_mask(rows, permitted)
    selected = unit_positions(mask, strides, group)
    assert max(map(len, selected)) <= window
    variables = [(c, j) for c in range(b) for j in selected[c]]
    width = len(variables)
    allowed = [r for r in range(rows) if mask[r]]
    row_index = {r: i for i, r in enumerate(allowed)}
    prefix = [1+(7*i+1) % 3 for i in range(b*rows)]
    choices = admissible_tuples()
    coefficients = [choices[(13*i+7) % len(choices)] for i in range(b*rows)]
    if three_product:
        coefficients = [(d1, d1, d2, f) for _, d1, d2, f in coefficients]
    # Map selected inputs injectively to the beginning of each regional window.
    destinations = {(c, j): c*rows+t for c in range(b) for t, j in enumerate(selected[c])}
    constraints = [[0]*width for _ in range(2*b+2*permitted)]
    for variable, (c, j) in enumerate(variables):
        position = destinations[c, j]
        _, d1, d2, _ = coefficients[position]
        constraints[2*c][variable] = d1
        constraints[2*c+1][variable] = d2
        order = placement_order(rows, strides[c], group)
        for endpoint in (j, j+1):
            row = row_index[order[endpoint]]
            scale = inverse(prefix[c*rows+endpoint])
            constraints[2*b+2*row][variable] ^= scale
            constraints[2*b+2*row+1][variable] ^= mul(c, scale)
    basis = nullspace(constraints, width)
    lower = width-2*b-(b-k)*permitted
    if len(basis) < lower or not basis:
        raise ArithmeticError('Actual-code dimension inequality failed')
    for vector in basis:
        middle = [0]*(b*rows)
        post_input = [0]*(b*rows)
        for value, (c, j) in zip(vector, variables):
            middle[c*rows+j] = value
            post_input[destinations[c, j]] = value
        # Check every actual GRS[4,2] row, after the inverse prefix diagonal.
        source = [[0]*b for _ in range(rows)]
        previous = 0
        for c in range(b):
            order = placement_order(rows, strides[c], group)
            for j, r in enumerate(order):
                position = c*rows+j
                source[r][c] = mul(inverse(prefix[position]), middle[position] ^ previous)
                previous = middle[position]
        for row in source:
            check0 = check1 = 0
            for c, value in enumerate(row):
                check0 ^= value
                check1 ^= mul(c, value)
            if check0 or check1:
                raise ArithmeticError('The constructed vector failed actual MDS checks')
        s1 = s2 = 0
        output = []
        for position, (x, (c, d1, d2, f)) in enumerate(zip(post_input, coefficients)):
            output.append(s1 ^ mul(f, s2) ^ mul(c, x))
            s1 ^= mul(d1, x)
            s2 ^= mul(d2, x)
            if position % rows == window-1 and (s1 or s2):
                raise ArithmeticError('Post state was not zero at the window exit')
        if not any(output) or any(value for i, value in enumerate(output) if i % rows >= window):
            raise ArithmeticError('A nonzero codeword escaped the promised output windows')
    return dict(group=group, three_product=three_product, variables=width,
                dimension_lower=lower, actual_dimension=len(basis),
                checked_basis_words=len(basis), output_weight_upper=b*window)


def main():
    print('dense tail:', dense_tail_upper())
    for group in (1, 4):
        print('actual-code first-moment obstruction:', obstruction_receipt(group))
        for three_product in (False, True):
            print('small actual-code regression:', small_actual_code_check(group, three_product))


if __name__ == '__main__':
    main()
