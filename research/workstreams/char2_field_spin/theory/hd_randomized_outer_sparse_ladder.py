"""Exact regional coefficients for the randomized-outer sparse ladder.

The float path selects witnesses only.  Directed replay uses independently
constructed Arb local operators and positive matrix-polynomial arithmetic.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, inf, log2

from hd_randomized_outer_state_profile import local_profiles, power, total
from hd_randomized_outer_proof import prescribed_mds_count


@lru_cache(None)
def occupied_tile_counts(j):
    """Exact numbers of j-subsets occupying each possible number of FFT8 tiles."""
    return tuple((u, comb(8, u)*sum((-1)**i*comb(u, i)*comb(8*(u-i), j)
                                    for i in range(u+1) if j <= 8*(u-i)))
                 for u in range(9))


def refine_residual_returns(profile, log_z):
    """Charge untouched emissions jointly with the live-to-zero probability."""
    from hd_randomized_outer_state_screen import total
    q, d = 65536, 65536**4-1
    result = list(profile)
    for j in range(1, len(profile)):
        rank = min(j, 4)
        unweighted = -16*rank+log2(1/q+(1-1/q)*(q**rank-1)/d)
        untouched = total(log2(count)+max(0, 64-8*u-3)*float(log_z)
                          for u, count in occupied_tile_counts(j) if count)-log2(comb(64, j))
        rows = [list(row) for row in profile[j]]
        rows[2][0] = min(rows[2][0], unweighted+untouched)
        result[j] = tuple(tuple(row) for row in rows)
    return tuple(result)


def refine_residual_membership(profile, log_z):
    """Average the fixed state's subspace membership before its supremum.

    For j<4 MDS columns, a fixed nonzero vector belongs to their random span
    with probability at most j/64.  I/L preserve its scalar line; J is flat.
    """
    q, d = 65536, 65536**4-1
    result = list(profile)
    for j in range(1, min(4, len(profile))):
        probability = (1/q)*(j/64)+(1-1/q)*(q**j-1)/d
        bound = -16*j+log2(probability)+max(0, 61-8*j)*float(log_z)
        rows = [list(row) for row in result[j]]
        rows[2][0] = min(rows[2][0], bound)
        result[j] = tuple(tuple(row) for row in rows)
    return tuple(result)


def refine_flat_density(profile, log_z):
    """FLOAT-only adapter for the independently developed flat-state moment."""
    from hd_randomized_outer_flat_moment import flat_moments
    from hd_randomized_outer_state_screen import add, subtract, total
    field, d = 65536, 65536**4-1
    log_z = float(log_z)
    new_uniform = flat_moments(F(2.0**log_z))
    t = add(-16, log2(1-1/field)+log_z)
    result = list(profile)
    for j in range(len(profile)):
        gamma = total(log2(count)+min(j*t, max(0, 64-8*u-3)*log_z)
                      for u, count in occupied_tile_counts(j) if count)-log2(comb(64, j))
        gamma0 = float('-inf')
        if j:
            gamma0 = total(log2(count)+min(subtract(j*t, -j*16),
                max(0, 64-8*u-3)*log_z+log2(1-field**(-min(j, 4))))
                for u, count in occupied_tile_counts(j) if count)-log2(comb(64, j))
        rows = [list(row) for row in profile[j]]
        rows[1][1] = min(rows[1][1], add(log2(1-1/field)+new_uniform[j], -16+gamma))
        if j:
            rows[1][0] = min(rows[1][0], add(log2(1-1/field)+new_uniform[j], -16+gamma0)-log2(d))
        result[j] = tuple(tuple(row) for row in rows)
    return tuple(result)


def float_polynomial_product(left, right, degree):
    """Small positive matrix convolution; no coefficient subtraction."""
    import numpy as np
    length = min(degree+1, len(left)+len(right)-1)
    result = np.zeros((length, 3, 3))
    for i in range(3):
        for k in range(3):
            for j in range(3):
                result[:, i, j] += np.convolve(left[:, i, k], right[:, k, j])[:length]
    return result


@lru_cache(maxsize=64)
def float_block_paths(log_z, maximum, steps, provider='base'):
    """All count paths in a short block; products precede fractional powers."""
    import numpy as np
    if steps not in (1, 2, 4, 8, 16) or comb(maximum+steps, steps) > 100000:
        raise ValueError('Need a power-of-two short block with at most 100000 count paths')
    profile = local_profiles(float(log_z))
    supported = ('base', 'return', 'columns', 'columns_return', 'flat',
                 'flat_return', 'flat_columns_return', 'membership', 'refined_membership')
    if provider not in supported:
        raise ValueError('Unknown local-operator provider')
    if provider == 'refined_membership':
        from hd_randomized_outer_refined_profile import refined_profiles
        profile = refined_profiles(float(log_z))
    if 'flat' in provider:
        profile = refine_flat_density(profile, log_z)
    if 'columns' in provider:
        from hd_randomized_outer_column_moment import refine_columns
        profile = refine_columns(profile, float(log_z))
    if 'return' in provider or 'membership' in provider:
        profile = refine_residual_returns(profile, log_z)
    if 'membership' in provider:
        profile = refine_residual_membership(profile, log_z)
    local = np.exp2(np.array(profile[:maximum+1]))
    degrees, weights = np.array([0]), np.array([1.0])
    products = np.eye(3).reshape((1, 3, 3))
    for _ in range(steps):
        next_degrees, next_weights, next_products = [], [], []
        for j in range(maximum+1):
            keep = degrees+j <= maximum
            if not np.any(keep):
                continue
            next_degrees.append(degrees[keep]+j)
            next_weights.append(weights[keep]*comb(64, j))
            next_products.append(products[keep]@local[j])
        degrees = np.concatenate(next_degrees)
        weights = np.concatenate(next_weights)
        products = np.concatenate(next_products)
    return degrees, weights, products


def float_regions(log_z, alpha=1, *, maximum=16, steps=1, provider='base'):
    """All regional coefficients through maximum, for one weight witness."""
    import numpy as np
    degrees, weights, products = float_block_paths(log_z, maximum, steps, provider)
    coefficients = np.zeros((maximum+1, 3, 3))
    weighted = weights[:, None, None]*products**float(alpha)
    for j in range(maximum+1):
        coefficients[j] = weighted[degrees == j].sum(axis=0)
    result = np.eye(3).reshape((1, 3, 3))
    exponent = 128//steps
    while exponent:
        if exponent & 1:
            result = float_polynomial_product(result, coefficients, maximum)
        exponent >>= 1
        if exponent:
            coefficients = float_polynomial_product(coefficients, coefficients, maximum)
    return tuple(result[r]/comb(8192, r) for r in range(maximum+1))


def float_endpoint(rows, region, log_z, alpha=1):
    import numpy as np
    with np.errstate(divide='ignore'):
        logarithms = np.log2(region)
    run = power(tuple(tuple(float(x) for x in row) for row in logarithms), 256)
    return (log2(comb(8192, rows))
            +float(alpha)*(rows*(4096-128*log2(65535))-log2(65535)-838860*float(log_z))
            +total(run[0]))


def screen(maximum=16, *, alphas=(F(1),), steps=1, provider='base'):
    """Bounded screen of EVERY occupancy 1..maximum; not a certificate."""
    if not 1 <= maximum <= 64:
        raise ValueError('This sparse screen is restricted to at most 64 rows')
    marks = sorted(set([F(-i, 512) for i in range(1, 2*maximum+1)]
                       +[F(-1, 4096), F(-1, 2048), F(-1, 1024)]))
    best = [None]*(maximum+1)
    for alpha in alphas:
        for log_z in marks:
            regions = float_regions(log_z, alpha, maximum=maximum, steps=steps, provider=provider)
            for rows in range(1, maximum+1):
                value = float_endpoint(rows, regions[rows], log_z, alpha)
                if best[rows] is None or value < best[rows]['screen_log2']:
                    best[rows] = dict(active_rows=rows, log2_z=str(log_z), alpha=str(alpha),
                                      steps=steps, provider=provider, screen_log2=value, directed=False)
    return tuple(best[1:])


def _minimum(left, right):
    """Both arguments enclose valid upper expressions; either is safe."""
    return left if left < right else right


def arb_local_profiles(log_z, maximum, provider='base'):
    """Directed reconstruction of the unchanged Z/B/R positive operator."""
    from flint import arb
    q, d = 65536, 65536**4-1
    zero, one = arb(0), arb(1)
    epsilon, refresh = one/q, arb(q-1)/q
    log_z = F(log_z)
    z = arb(2)**(arb(log_z.numerator)/log_z.denominator)
    t = one/q+refresh*z
    tile = []
    for a in range(9):
        value = sum((comb(8, w)*prescribed_mds_count(w, 8, a, q)*z**w
                     for w in range(9)), zero)/q**a
        tile.append(value*comb(8, a))
    moments, occupancies = {(0, 0): one}, {(0, 0): 1}
    for _ in range(8):
        following, next_occupancies = {}, {}
        for (j, occupied), value in moments.items():
            for a in range(min(8, maximum-j)+1):
                key = j+a, occupied+(a != 0)
                following[key] = following.get(key, zero)+value*tile[a]
                next_occupancies[key] = next_occupancies.get(key, 0)+occupancies[j, occupied]*comb(8, a)
        moments, occupancies = following, next_occupancies
    state_moments = {}
    for occupied in range(9):
        length = 64-8*occupied
        rank = min(length, 4)
        multiplicity = q**(4-rank)
        value = arb(multiplicity-1)
        value += sum((multiplicity*comb(length, w)
                      *prescribed_mds_count(w, length, rank, q)*z**w
                      for w in range(1, length+1)), zero)
        state_moments[occupied] = value/d
    result = []
    for j in range(maximum+1):
        normalization = comb(64, j)
        plain = live = uniform = gamma = gamma_zero = untouched_average = zero
        for occupied in range(9):
            key = j, occupied
            if key not in moments:
                continue
            value, count = moments[key], occupancies[key]
            untouched = z**max(0, 64-8*occupied-3)
            plain += value
            live += value*untouched
            uniform += value*state_moments[occupied]
            gamma += count*_minimum(t**j, untouched)
            untouched_average += count*untouched
            if j:
                gamma_zero += count*_minimum(t**j-one/q**j,
                                              untouched*(one-one/q**min(j, 4)))
        plain, live, uniform, gamma, gamma_zero, untouched_average = (
            value/normalization for value in (plain, live, uniform, gamma, gamma_zero, untouched_average))
        rank = min(j, 4)
        zero_zero = _minimum(plain, one/q**rank)
        zero_live = plain-one/q**j
        residual_zero = zero
        if j:
            residual_zero = _minimum(arb(q**3)/d*live,
                (epsilon+refresh*arb(q**rank-1)/d)/q**rank)
            if provider in ('return', 'membership'):
                residual_zero = _minimum(residual_zero,
                    untouched_average*(epsilon+refresh*arb(q**rank-1)/d)/q**rank)
            if provider == 'membership' and j < 4:
                residual_zero = _minimum(residual_zero, z**max(0, 61-8*j)
                    *(epsilon*arb(j)/64+refresh*arb(q**j-1)/d)/q**j)
        bulk_zero = (refresh*uniform+epsilon*gamma_zero)/d if j else zero
        result.append(((zero_zero, zero, zero_live),
                       (bulk_zero, refresh*uniform+epsilon*gamma, zero),
                       (residual_zero, refresh*live, epsilon*live)))
    return tuple(result)


def arb_matrix_product(left, right):
    from flint import arb
    return tuple(tuple(sum((left[i][k]*right[k][j] for k in range(3)), arb(0))
                       for j in range(3)) for i in range(3))


def arb_polynomial_product(left, right, degree):
    from flint import arb
    length = min(degree+1, len(left)+len(right)-1)
    result = [[[arb(0) for _ in range(3)] for _ in range(3)] for _ in range(length)]
    for a, left_matrix in enumerate(left):
        for b, right_matrix in enumerate(right[:length-a]):
            product = arb_matrix_product(left_matrix, right_matrix)
            for i in range(3):
                for j in range(3):
                    result[a+b][i][j] += product[i][j]
    return tuple(tuple(tuple(row) for row in matrix) for matrix in result)


def directed_endpoint(rows, log_z, alpha=F(1), *, steps=1, precision=192, provider='base'):
    """One complete occupancy, including all active-row choices and messages."""
    from flint import arb, ctx
    if not 1 <= rows <= 16 or steps not in (1, 2, 4, 8, 16):
        raise ValueError('Outside the bounded sparse ladder')
    if comb(rows+steps, steps) > 100000:
        raise ValueError('Too many block count paths for this directed replay')
    log_z, alpha = F(log_z), F(alpha)
    if not log_z < 0 or not 0 < alpha <= 1:
        raise ValueError('Need a negative tilt and 0<alpha<=1')
    if provider not in ('base', 'return', 'membership'):
        raise ValueError('Unsupported directed local-operator provider')
    previous = ctx.prec
    ctx.prec = precision
    try:
        zero = arb(0)
        identity = tuple(tuple(arb(int(i == j)) for j in range(3)) for i in range(3))
        local = arb_local_profiles(log_z, rows, provider)
        paths = [(0, 1, identity)]
        for _ in range(steps):
            paths = [(degree+j, weight*comb(64, j), arb_matrix_product(matrix, local[j]))
                     for degree, weight, matrix in paths for j in range(rows-degree+1)]
        aa = arb(alpha.numerator)/alpha.denominator
        coefficients = [[[arb(0) for _ in range(3)] for _ in range(3)] for _ in range(rows+1)]
        for degree, weight, matrix in paths:
            for i in range(3):
                for j in range(3):
                    coefficients[degree][i][j] += weight*matrix[i][j]**aa
        result, base, exponent = (identity,), coefficients, 128//steps
        while exponent:
            if exponent & 1:
                result = arb_polynomial_product(result, base, rows)
            exponent >>= 1
            if exponent:
                base = arb_polynomial_product(base, base, rows)
        region = tuple(tuple(value/comb(8192, rows) for value in row) for row in result[rows])
        run, base, exponent = identity, region, 256
        while exponent:
            if exponent & 1:
                run = arb_matrix_product(run, base)
            exponent >>= 1
            if exponent:
                base = arb_matrix_product(base, base)
        logarb = lambda x: x.log()/arb(2).log()
        value = (logarb(arb(comb(8192, rows)))
                 +aa*(rows*(4096-128*logarb(arb(65535)))-logarb(arb(65535))
                       -838860*arb(log_z.numerator)/log_z.denominator)
                 +logarb(sum(run[0], zero)))
        if not value.is_finite():
            raise ArithmeticError('The directed endpoint is not finite')
        return dict(active_rows=rows, log2_z=str(log_z), alpha=str(alpha), steps=steps, provider=provider,
                    failure_log2_upper=str(value), below_2_to_minus_40=bool(value < -40),
                    covers_all_words=True, directed=True, exact_regional_coefficient=True,
                    whole_code_certificate=False)
    finally:
        ctx.prec = previous


def sparse_receipt(provider='membership'):
    """Directed union of all r=1,...,5 classes for the current witnesses."""
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = 192
    try:
        witnesses = ((1, F(-1, 512), F(1), 1), (2, F(-1, 256), F(1), 1),
                     (3, F(-3, 512), F(11, 16) if provider == 'membership' else F(5, 8), 4),
                     (4, F(-1, 128), F(5, 16), 4),
                     (5, F(-5, 512), F(1, 4), 4))
        records = tuple(directed_endpoint(r, z, a, steps=b, provider=provider)
                        for r, z, a, b in witnesses)
        total_bound = sum((arb(2)**arb(record['failure_log2_upper']) for record in records), arb(0))
        logarithm = total_bound.log()/arb(2).log()
        if not logarithm < -40:
            raise ArithmeticError('The joined sparse range does not meet 40 bits')
        return dict(active_rows_first=1, active_rows_last=5, provider=provider, endpoints=records,
                    union_log2_upper=str(logarithm), directed=True, covers_all_words=True,
                    whole_code_certificate=False)
    finally:
        ctx.prec = previous


def main():
    import json
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', action='store_true')
    parser.add_argument('--steps', type=int, default=1)
    parser.add_argument('--fractional', action='store_true')
    parser.add_argument('--provider', choices=('base', 'return', 'membership', 'refined_membership'),
                        default='membership')
    arguments = parser.parse_args()
    if arguments.receipt:
        print(json.dumps(sparse_receipt(arguments.provider), indent=2))
    else:
        alphas = tuple(F(i, 16) for i in range(1, 17)) if arguments.fractional else (F(1),)
        for item in screen(steps=arguments.steps, alphas=alphas, provider=arguments.provider):
            print(json.dumps(item), flush=True)


if __name__ == '__main__':
    main()
