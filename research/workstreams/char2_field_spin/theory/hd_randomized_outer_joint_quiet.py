"""Joint quiet/output moment using tile roots and residual feedback rank.

Only selected occupancies are screened.  This module does not issue a range
or whole-code certificate.  The numerical screen is replayed with Arb at its
chosen rational witnesses before an occupancy endpoint is reported.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, lgamma, log, log2

from hd_randomized_outer_proof import Parameters


@lru_cache(None)
def tile_residual_terms(mark, field=65536):
    """Exact coefficients, scaled by q^a; keyed by (a,min(4,a-v)).

    a is the potential input count, v the selected output-zero count.  Each
    actual mask has output constraint rank min(a,v).  Residual feedback rank
    is min(4,max(a-v,0)) by polynomial degree filtration.
    """
    mark = F(mark)
    if not 0 < mark <= 1:
        raise ValueError('Need 0<z<=1')
    result = [[F(0) for _ in range(5)] for _ in range(9)]
    for a in range(9):
        for v in range(9):
            dimension = max(a-v, 0)
            result[a][min(4, dimension)] += (
                comb(8, a)*comb(8, v)*mark**(8-v)*(1-mark)**v
                *field**dimension)
    return tuple(tuple(row) for row in result)


def _joint_moments(mark, field=65536, *, directed=False, distinct_lines=True):
    """K_j upper bounds and exact unconstrained T_j, before fractional powers.

    DP states keep total potential slots and a feedback-rank lower bound.
    Two different tiles' one-dimensional residual images are distinct.
    Scaling by q^j avoids underflow in the float screen.
    """
    if directed:
        from flint import arb
        convert = lambda value: arb(value.numerator)/value.denominator
        zero, one = arb(0), arb(1)
    else:
        convert = float
        zero, one = 0.0, 1.0
    if mark == 1:
        if directed:
            return (tuple(one/field**min(j, 4) for j in range(65)), (one,)*65)
        return (tuple(-min(j, 4)*log2(field) for j in range(65)), (0.0,)*65)
    local = tuple(tuple(convert(value) for value in row)
                  for row in tile_residual_terms(mark, field))
    dp = [[zero for _ in range(5)]]
    dp[0][0] = one
    for _ in range(8):
        updated = [[zero for _ in range(5)] for _ in range(len(dp)+8)]
        for j, row in enumerate(dp):
            for state, previous in enumerate(row):
                if not previous:
                    continue
                for a, local_row in enumerate(local):
                    for residual, coefficient in enumerate(local_row):
                        if coefficient:
                            successor = (2 if distinct_lines and state == residual == 1
                                         else max(state, residual))
                            updated[j+a][successor] += previous*coefficient
        dp = updated
    quiet, free = [], []
    for j, row in enumerate(dp):
        normalizer = field**j*comb(64, j)
        weighted = sum(row[m]/field**m for m in range(5))
        unweighted = sum(row)
        if directed:
            bound = weighted/normalizer
            cap = arb(1)/field**min(j, 4)
            quiet.append(bound if bound < cap else cap)
            free.append(unweighted/normalizer)
        else:
            # Return logarithms instead of potentially subnormal probabilities.
            denominator = j*log2(field)+log2(comb(64, j))
            quiet.append(min(log2(weighted)-denominator, -min(j, 4)*log2(field)))
            free.append(log2(unweighted)-denominator)
    return tuple(quiet), tuple(free)


def joint_moment_logs(mark, field=65536, *, distinct_lines=True):
    return _joint_moments(F(mark), field, directed=False, distinct_lines=distinct_lines)


def joint_moment_enclosures(mark, field=65536, *, distinct_lines=True):
    return _joint_moments(F(mark), field, directed=True, distinct_lines=distinct_lines)


def _logsum2(values):
    largest = max(values)
    return largest+log2(sum(2**(value-largest) for value in values))


def screen_selected(occupancies=(954, 1024, 1536, 2048, 4096, 5120, 5632, 6144, 6565)):
    """Bounded floating search only; no sweep over the middle occupancy range."""
    p = Parameters()
    marks = tuple(F(1, 1 << shift) for shift in (0, 1, 2, 4, 6, 8, 10, 12, 14, 16))
    marks += (F(2, 3*(p.field-1)), F(133, 128*(p.field-1)))
    alphas = tuple(F(a, 32) for a in (1, 2, 3, 4, 6, 8, 12, 16, 24, 32))
    best = {r: None for r in occupancies}
    log_combinations = tuple(log2(comb(64, j)) for j in range(65))
    beta_log = 256*log2(p.field)-128*log2(p.field-1)
    for mark in marks:
        quiet, free = joint_moment_logs(mark)
        for rows in occupancies:
            choose = (lgamma(p.rows+1)-lgamma(rows+1)-lgamma(p.rows-rows+1))/log(2)
            center = rows/p.rows
            markers = sorted(set([F(min(127, max(1, round(128*center*scale))), 128)
                                  for scale in (0.5, 0.75, 1, 1.25, 1.5, 2, 3)]
                                 +[F(1, 8)]))
            for marker in markers:
                m = float(marker)
                log_weights = tuple(log_combinations[j]+j*log2(m)+(64-j)*log2(1-m)
                                    for j in range(65))
                correction = 256*(choose+rows*log2(m)+(p.rows-rows)*log2(1-m))
                for alpha in alphas:
                    a = float(alpha)
                    ordinary = _logsum2(tuple(log_weights[j]+a*quiet[j] for j in range(65)))
                    last = _logsum2(tuple(log_weights[j]+a*free[j] for j in range(65)))
                    value = (choose+a*(rows*beta_log-log2(p.field-1)-p.limit*log2(float(mark)))
                             +(p.packets-1)*ordinary+last-correction)
                    if best[rows] is None or value < best[rows]['screen_log2']:
                        best[rows] = dict(active_rows=rows, mark=str(mark), marker=str(marker),
                                          alpha=str(alpha), screen_log2=value)
    return tuple(best[r] for r in occupancies)


def directed_endpoint(rows, mark, marker, alpha):
    """Directed endpoint for exactly one occupancy, unioned over row subsets."""
    from flint import arb
    p = Parameters()
    mark, marker, alpha = map(F, (mark, marker, alpha))
    quiet, free = joint_moment_enclosures(mark)
    convert = lambda value: arb(value.numerator)/value.denominator
    z, m, a = map(convert, (mark, marker, alpha))
    ordinary, last = arb(0), arb(0)
    for j in range(65):
        weight = comb(64, j)*m**j*(1-m)**(64-j)
        ordinary += weight*quiet[j]**a
        last += weight*free[j]**a
    logarb = lambda value: value.log()/arb(2).log()
    choose = logarb(arb(comb(p.rows, rows)))
    correction = 256*(choose+rows*logarb(m)+(p.rows-rows)*logarb(1-m))
    beta_log = 256*logarb(arb(p.field))-128*logarb(arb(p.field-1))
    result = (choose+a*(rows*beta_log-logarb(arb(p.field-1))-p.limit*logarb(z))
              +(p.packets-1)*logarb(ordinary)+logarb(last)-correction)
    return dict(active_rows=rows, mark=str(mark), marker=str(marker), alpha=str(alpha),
                failure_log2_upper=str(result), below_2_to_minus_40=bool(result < -40),
                covers_only_this_occupancy=True, whole_code_certificate=False)


def directed_interval(first=5940, last=8192, mark=F(1, 16384),
                      marker=F(127, 128), alpha=F(3, 8)):
    """Certify an entire interval with two endpoints, using convexity in r.

    With fixed witnesses, log bound = affine(r)-255*log2(binomial(8192,r)).
    Binomial coefficients are log-concave, hence the bound is convex.  Its
    maximum on the interval occurs at an endpoint.  This is not interpolation
    from a heuristic screen; convexity proves the interior bound.
    """
    from flint import arb
    if not 1 <= first < last <= 8192:
        raise ValueError('Need a nontrivial valid occupancy interval')
    lower_endpoint = directed_endpoint(first, mark, marker, alpha)
    upper_endpoint = directed_endpoint(last, mark, marker, alpha)
    lower = arb(lower_endpoint['failure_log2_upper'])
    upper = arb(upper_endpoint['failure_log2_upper'])
    if lower > upper:
        maximum = lower
    elif upper > lower:
        maximum = upper
    else:
        raise ArithmeticError('Increase precision to order the endpoint enclosures')
    union = maximum+arb(last-first+1).log()/arb(2).log()
    return dict(class_kind='quiet_low_weight', active_rows_first=first,
                active_rows_last=last, mark=str(mark), marker=str(marker),
                alpha=str(alpha), first_endpoint_log2_upper=str(lower),
                last_endpoint_log2_upper=str(upper), union_log2_upper=str(union),
                below_2_to_minus_40=bool(union < -40),
                interior_justification='convexity of affine(r)-255*log2(binomial(8192,r))',
                occupancy_endpoints_evaluated=2, whole_code_certificate=False)


def tile_ranks(tile, input_mask, zero_mask):
    """Actual beta64 output rank and residual feedback rank over GF256."""
    from hd_blockwise_shared_maps import matrices, rank
    feedthrough, _, feedback, _, _ = matrices(64)
    columns = tuple(8*tile+i for i in range(8) if input_mask >> i & 1)
    rows = tuple(8*tile+i for i in range(8) if zero_mask >> i & 1)
    output = tuple(tuple(feedthrough[i][j] for j in columns) for i in rows)
    combined = output+tuple(tuple(row[j] for j in columns) for row in feedback)
    output_rank = rank(output)
    return output_rank, rank(combined)-output_rank


def main():
    import json
    from flint import ctx
    ctx.prec = 160
    selected = screen_selected()
    for witness in selected:
        print(json.dumps(directed_endpoint(witness['active_rows'], witness['mark'],
                                            witness['marker'], witness['alpha'])))
    print(json.dumps(directed_interval()))


if __name__ == '__main__':
    main()
