"""Three-component weighted-state envelope; floating screens, not certificates.

Components are zero mass, a cap relative to uniform nonzero-state density,
and arbitrary nonzero residual mass. Births enter the residual component.
"""

from functools import lru_cache
from math import comb, inf, log2

from hd_randomized_outer_proof import prescribed_mds_count
from hd_randomized_outer_state_screen import (
    NEG, add, total, subtract, mds_moment, local_operators as two_state_operators,
)


def multiply(left, right):
    size = len(left)
    return tuple(tuple(total(left[i][k]+right[k][j] for k in range(size))
                       for j in range(size)) for i in range(size))


def power(matrix, exponent):
    size = len(matrix)
    result = tuple(tuple(0.0 if i == j else NEG for j in range(size))
                   for i in range(size))
    while exponent:
        if exponent & 1:
            result = multiply(result, matrix)
        matrix = multiply(matrix, matrix)
        exponent >>= 1
    return result


def punctured_state_moment(length, log_z, q=65536, state=4):
    """Exact-enumerator FLOAT moment for a uniform nonzero MDS state word."""
    rank = min(length, state)
    multiplicity = q**(state-rank)
    terms = []
    if multiplicity > 1:
        terms.append(log2(multiplicity-1))
    for weight in range(1, length+1):
        count = (multiplicity*comb(length, weight)
                 *prescribed_mds_count(weight, length, rank, q))
        if count:
            terms.append(log2(count)+weight*log_z)
    return total(terms)-log2(q**state-1)


@lru_cache(None)
def local_profiles(log_z, mixer='rank_one', q=65536):
    """Positive cone transfer in log2 arithmetic; indices are (Z,B,R)."""
    if not log_z < 0 or mixer not in ('rank_one', 'uniform'):
        raise ValueError('Need a negative weight tilt and a supported mixer')
    log_q, log_d = log2(q), log2(q**4-1)
    epsilon = 1/q if mixer == 'rank_one' else 0
    log_epsilon = -log_q if epsilon else NEG
    log_c = log2(1-epsilon)
    log_t = add(-log_q, log2(1-1/q)+log_z)
    tile = [log2(comb(8, a))+mds_moment(a, log_z, q) for a in range(9)]
    moments, occupancies = {(0, 0):0.0}, {(0, 0):0.0}
    for unused in range(8):
        following, next_occupancies = {}, {}
        for (j, occupied), value in moments.items():
            for a in range(9):
                key = j+a, occupied+(a != 0)
                following[key] = add(following.get(key, NEG), value+tile[a])
                next_occupancies[key] = add(next_occupancies.get(key, NEG),
                    occupancies[j, occupied]+log2(comb(8, a)))
        moments, occupancies = following, next_occupancies
    state_moments = {u:punctured_state_moment(64-8*u, log_z, q) for u in range(9)}
    two_state = two_state_operators(log_z, mixer, q)
    result = []
    for j in range(65):
        normalization = log2(comb(64, j))
        live = two_state[j][1][1]
        uniform = total(value+state_moments[u] for (jj, u), value in moments.items()
                        if jj == j)-normalization
        gamma = total(value+min(j*log_t, max(0, 64-8*u-3)*log_z)
                      for (jj, u), value in occupancies.items() if jj == j)-normalization
        gamma_zero = NEG
        if j:
            nonzero_syndrome = log2(1-q**(-min(j, 4)))
            nonzero_input = subtract(j*log_t, -j*log_q)
            gamma_zero = total(value+min(nonzero_input,
                max(0, 64-8*u-3)*log_z+nonzero_syndrome)
                for (jj, u), value in occupancies.items() if jj == j)-normalization
        bulk_zero = (add(log_c+uniform, log_epsilon+gamma_zero)-log_d
                     if j else NEG)
        bulk_bulk = add(log_c+uniform, log_epsilon+gamma)
        zero_zero, zero_live = two_state[j][0]
        residual_zero = two_state[j][1][0]
        result.append(((zero_zero, NEG, zero_live),
                       (bulk_zero, bulk_bulk, NEG),
                       (residual_zero, log_c+live, log_epsilon+live)))
    return tuple(result)


@lru_cache(None)
def local_pairs(log_z, mixer='rank_one', q=65536):
    local = local_profiles(log_z, mixer, q)
    chooses = [log2(comb(64, j)) for j in range(65)]
    return tuple((j+k, chooses[j]+chooses[k], multiply(local[j], local[k]))
                 for j in range(65) for k in range(65))


@lru_cache(None)
def pair_coefficients(log_z, alpha, mixer='rank_one', q=65536):
    coefficients = [[[NEG]*3 for unused in range(3)] for degree in range(129)]
    for degree, weight, pair in local_pairs(log_z, mixer, q):
        for a in range(3):
            for b in range(3):
                coefficients[degree][a][b] = add(
                    coefficients[degree][a][b], weight+alpha*pair[a][b])
    return tuple(tuple(tuple(row) for row in coefficient) for coefficient in coefficients)


def matrix_add(left, right):
    return tuple(tuple(add(x, y) for x, y in zip(a, b)) for a, b in zip(left, right))


def coefficient_one(zero, one, exponent):
    """Coefficient of x in (zero+x*one)^exponent, with matrix products."""
    size = len(zero)
    identity = tuple(tuple(0.0 if i == j else NEG for j in range(size))
                     for i in range(size))
    empty = tuple((NEG,)*size for unused in range(size))
    result = identity, empty
    base = zero, one

    def product_pair(left, right):
        return (multiply(left[0], right[0]),
                matrix_add(multiply(left[0], right[1]), multiply(left[1], right[0])))

    while exponent:
        if exponent & 1:
            result = product_pair(result, base)
        base = product_pair(base, base)
        exponent >>= 1
    return result[1]


def route_bound_one_row(log_z, alpha, mixer='rank_one', *, steps=1):
    """Exact regional one-mark coefficient; no Bernoulli conditioning cost."""
    local = local_profiles(log_z, mixer)
    fractional = lambda matrix: tuple(tuple(alpha*x for x in row) for row in matrix)
    if steps == 1:
        zero, one = fractional(local[0]), fractional(local[1])
    elif steps == 2:
        zero = fractional(multiply(local[0], local[0]))
        one = matrix_add(fractional(multiply(local[0], local[1])),
                         fractional(multiply(local[1], local[0])))
        one = tuple(tuple(x-1 for x in row) for row in one)
    else:
        raise ValueError('Only one- and two-packet blocks are screened')
    blocks = 128//steps
    region = coefficient_one(zero, one, blocks)
    region = tuple(tuple(x-log2(blocks) for x in row) for row in region)
    run = power(region, 256)
    log_beta = 4096-128*log2(65535)
    return (13+alpha*(log_beta-log2(65535)-838860*log_z)+total(run[0]))


def one_row_receipt(precision=192):
    """Directed actual-rank-one bound for every word with one active row.

    Only j=0,1 packet operators occur. The mark is the positive real number
    z=2^(-1/512); all occurrences are enclosed with Arb, not binary floats.
    """
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        q, d = 65536, 65536**4-1
        epsilon, refresh = arb(1)/q, arb(q-1)/q
        z = arb(2)**(-arb(1)/512)

        def state_moment(length):
            return sum((arb(comb(length, weight)
                            *prescribed_mds_count(weight, length, 4, q))*z**weight
                        for weight in range(length-3, length+1)), arb(0))/d

        live0 = z**61
        uniform0 = state_moment(64)
        input1 = arb(1)/q+refresh*z**8
        live1 = input1*z**53
        uniform1 = input1*state_moment(56)
        gamma1 = z**53
        gamma_zero1 = refresh*z**53
        residual_zero1 = (epsilon+refresh*arb(q-1)/d)/q
        if not residual_zero1 < arb(q**3)/d*live1:
            raise ArithmeticError('The fixed residual-return branch changed')
        zero, one = arb(0), arb(1)
        packet0 = ((one, zero, zero),
                   (zero, refresh*uniform0+epsilon*live0, zero),
                   (zero, refresh*live0, epsilon*live0))
        packet1 = ((arb(1)/q, zero, refresh*z**8),
                   ((refresh*uniform1+epsilon*gamma_zero1)/d,
                    refresh*uniform1+epsilon*gamma1, zero),
                   (residual_zero1, refresh*live1, epsilon*live1))
        identity = tuple(tuple(arb(int(i == j)) for j in range(3)) for i in range(3))
        empty = tuple((zero,)*3 for unused in range(3))

        def matmul(left, right):
            return tuple(tuple(sum((left[i][k]*right[k][j] for k in range(3)), zero)
                               for j in range(3)) for i in range(3))

        def matadd(left, right):
            return tuple(tuple(x+y for x, y in zip(a, b)) for a, b in zip(left, right))

        def pairmul(left, right):
            return (matmul(left[0], right[0]),
                    matadd(matmul(left[0], right[1]), matmul(left[1], right[0])))

        result, base, exponent = (identity, empty), (packet0, packet1), 128
        while exponent:
            if exponent & 1:
                result = pairmul(result, base)
            base = pairmul(base, base)
            exponent >>= 1
        region = tuple(tuple(x/128 for x in row) for row in result[1])
        run, base, exponent = identity, region, 256
        while exponent:
            if exponent & 1:
                run = matmul(run, base)
            base = matmul(base, base)
            exponent >>= 1
        bound = (arb(8192)*q**256/(arb(q-1)**129)*z**(-838860)*sum(run[0], zero))
        logarithm = bound.log()/arb(2).log()
        if not logarithm < -40:
            raise ArithmeticError('The one-row class did not close')
        return dict(active_rows=1, covers_all_words=True, mixer='rank_one',
                    weight_limit=838860, mark='2^(-1/512)', alpha='1',
                    failure_log2_upper=str(logarithm),
                    exact_regional_coefficient=True, directed=True,
                    whole_code_certificate=False)
    finally:
        ctx.prec = previous


def route_bound(active_rows, log_z, alpha, marker, mixer='rank_one', *, steps=1):
    if not 0 < active_rows < 8192 or not 0 < marker < 1 or not 0 < alpha <= 1:
        raise ValueError('Invalid route screen parameters')
    if steps == 1:
        local = local_profiles(log_z, mixer)
        weights = [log2(comb(64, j))+j*log2(marker)+(64-j)*log2(1-marker)
                   for j in range(65)]
        transfer = tuple(tuple(total(weights[j]+alpha*local[j][a][b]
                                    for j in range(65)) for b in range(3))
                         for a in range(3))
    elif steps == 2:
        coefficients = pair_coefficients(log_z, alpha, mixer)
        weights = [j*log2(marker)+(128-j)*log2(1-marker) for j in range(129)]
        transfer = tuple(tuple(total(weights[j]+coefficients[j][a][b]
                                    for j in range(129)) for b in range(3))
                         for a in range(3))
    else:
        raise ValueError('Only one- and two-packet blocks are screened')
    run = power(transfer, 32768//steps)
    choose = log2(comb(8192, active_rows))
    normalizer = choose+active_rows*log2(marker)+(8192-active_rows)*log2(1-marker)
    log_beta = 4096-128*log2(65535)
    return (choose+alpha*(active_rows*log_beta-log2(65535)-838860*log_z)
            +total(run[0])-256*normalizer)


def screen(active_rows, mixer='rank_one', steps=1, *, exact_region=False):
    """A bounded diagnostic grid; no optimization or directed acceptance."""
    best = (inf, None)
    if exact_region and active_rows != 1:
        raise ValueError('Exact regional coefficient is implemented only for one row')
    for log_z in (-1/4096, -1/2048, -1/1024, -1/512, -1/256, -1/128,
                  -1/64, -.0625, -.25, -1, -4, -8, -12, -16):
        for alpha in (1/32, 1/16, 1/8, 1/4, 1/2, 1):
            if exact_region:
                value = route_bound_one_row(log_z, alpha, mixer, steps=steps)
                if value < best[0]:
                    best = value, dict(log_z=log_z, alpha=alpha)
                continue
            for factor in (.25, .5, 1, 2, 4, 8):
                marker = active_rows/8192*factor
                if not 0 < marker < 1:
                    continue
                value = route_bound(active_rows, log_z, alpha, marker, mixer, steps=steps)
                if value < best[0]:
                    best = value, dict(log_z=log_z, alpha=alpha, marker=marker)
    return dict(active_rows=active_rows, mixer=mixer, steps=steps,
                exact_regional_coefficient=exact_region,
                log2_upper_diagnostic=best[0], witness=best[1],
                directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rows', nargs='+', type=int, default=[1, 16, 128, 512, 1024, 4096])
    parser.add_argument('--mixer', choices=['rank_one', 'uniform'], default='rank_one')
    parser.add_argument('--steps', type=int, choices=[1, 2], default=1)
    parser.add_argument('--exact-r1', action='store_true')
    parser.add_argument('--receipt', action='store_true')
    arguments = parser.parse_args()
    if arguments.receipt:
        print(json.dumps(one_row_receipt(), indent=2))
        raise SystemExit(0)
    for active_rows in arguments.rows:
        print(json.dumps(screen(active_rows, arguments.mixer, arguments.steps,
                                exact_region=arguments.exact_r1)), flush=True)
