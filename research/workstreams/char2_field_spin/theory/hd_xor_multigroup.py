"""Nonzero-amplitude operators and multigroup XOR2 screens and receipts.

Conditional on the output supports, D2 H D1 supplies independent uniform
nonzero amplitudes. These operators retain feedback cancellation at every
occupancy. The scalar potential comparison retains exact MDS support floors
and covers every constituent-row choice in the stated homogeneous class.
Saved directed witnesses cover the one-constituent class for 1..8192 groups.
No whole-code certificate is claimed.
"""
from functools import lru_cache
from fractions import Fraction
from math import ceil, comb, log2
from time import perf_counter

import numpy as np
from scipy.optimize import minimize_scalar

from hd_bucket_feedback import profiles
from hd_common_refresh_screen import local_return_cap
from hd_grouped_quiet import quiet_screen
from hd_grouped_subfield_screen import density_bits
from hd_randomized_outer_state_profile import total
from hd_xor_column_probe import (mds_prescribed_support,
                                 two_active_row_weight_enumerator)

Q = 256
REGION_LENGTH = 8192
REGION_PACKETS = 128
LIMIT = (1 << 21)*7//20


def nonzero_kernel_probability(occupancy, rank=4, q=Q):
    """Exact probability C_U X=0 for independent X_i uniform in F*.

    Every rank columns of C are independent. Its kernel on l>rank selected
    columns is an MDS[l,l-rank] code; require its entire support to be nonzero.
    """
    if occupancy < 0 or rank < 1 or q < 2:
        raise ValueError('Invalid MDS kernel parameters')
    if occupancy == 0:
        return Fraction(1)
    if occupancy <= rank:
        return Fraction(0)
    return Fraction(mds_prescribed_support(q, occupancy, occupancy-rank, occupancy),
                    (q-1)**occupancy)


@lru_cache(None)
def profile_arrays():
    items = profiles(16, 64, 4)
    return np.array([row for row, _ in items]), np.array([count for _, count in items], float)


def nonzero_matrices(z, q=Q):
    """Profile-wise positive two-state caps, before any fractional power.

    The first row is exact before the refresh density b. The flat-state
    row uses conditioning and the existing return/live caps. Empty inputs
    retain the established inverse-refresh bound.
    """
    if not 0 < z <= 1:
        raise ValueError('Weight mark must lie in (0,1]')
    length, rank, state = 16, 4, 16
    occ, counts = profile_arrays()
    degrees = occ.sum(axis=1)
    kernel = np.array([float(nonzero_kernel_probability(j, rank, q)) for j in range(length+1)])
    quiet = z**np.arange(length+1)*kernel
    uniform_quiet = np.array([sum(comb(j, c)*z**(j-c)*(1-z)**c*q**(-c-min(rank, j-c))
                                  for c in range(j+1)) for j in range(length+1)])
    free = (z+(1-z)/q)**np.arange(length+1)
    condition = (q/(q-1))**np.arange(length+1)
    ret = np.minimum(condition*np.maximum(0, free-uniform_quiet),
                     (1-kernel)*z**np.maximum(1, length-rank+1-np.arange(length+1)))
    ret = np.minimum(ret, condition*np.array([local_return_cap(j, z, rank, length, float(q))
                                              for j in range(length+1)]))
    base = float(q)**(-rank)
    excess = np.array([sum(comb(n, c)*z**(n-c)*(1-z)**c*(q**(-c)-base)
                           for c in range(min(rank, n+1))) for n in range(length+1)])
    packet_ret, packet_total = np.zeros(len(occ)), np.ones(len(occ))
    live_nonzero, live_total = np.zeros(len(occ)), np.ones(len(occ))
    for block in range(4):
        l = occ[:, block]
        packet_ret = packet_ret*quiet[l]+packet_total*ret[l]
        packet_total *= quiet[l]+ret[l]
        live_nonzero = live_nonzero*base+live_total*excess[length-l]
        live_total *= base+excess[length-l]
    live = np.prod(free[occ], axis=1)*live_nonzero/(1-float(q)**(-state))
    live *= (q/(q-1))**degrees
    d = q**state-1
    packet_ret = np.minimum(packet_ret, d*live)
    feedback_quiet = np.prod(kernel[occ], axis=1)
    emitted = z**degrees
    b = 2**density_bits(q, state)
    values = np.stack((emitted*feedback_quiet,
                       b*emitted*(1-feedback_quiet), packet_ret/d, b*live), axis=1).reshape((-1, 2, 2))
    values[degrees == 0, 1, 1] = np.minimum(values[degrees == 0, 1, 1],
        quiet_screen(z, q=q, state=state, width=64, h=rank))
    return degrees, counts, values


def nonzero_coefficients(z, alpha=1., maximum=64):
    if not 0 < alpha <= 1 or not 0 <= maximum <= 64:
        raise ValueError('Invalid fractional exponent or degree')
    degrees, counts, values = nonzero_matrices(z)
    result = np.zeros((maximum+1, 2, 2))
    keep = degrees <= maximum
    for row in range(2):
        for column in range(2):
            result[:, row, column] = np.bincount(degrees[keep],
                weights=counts[keep]*values[keep, row, column]**alpha, minlength=maximum+1)
    with np.errstate(divide='ignore'):
        return np.log2(result)


def log_polynomial_product(left, right, maximum):
    """Truncated positive matrix-polynomial product in base-two logarithms."""
    size = min(maximum+1, len(left)+len(right)-1)
    result = np.full((size, 2, 2), -np.inf)
    for index, matrix in enumerate(left[:size]):
        count = min(len(right), size-index)
        product = np.logaddexp2(matrix[:, 0, None]+right[:count, 0, None, :],
                               matrix[:, 1, None]+right[:count, 1, None, :])
        result[index:index+count] = np.logaddexp2(result[index:index+count], product)
    return result


def region_operators(z, maximum, alpha=1.):
    """Average over all subsets of m regional positions, for 0<=m<=maximum.

    For alpha<1, packet caps are raised entrywise before route averaging;
    callers must separately justify their fractional-moment aggregation.
    """
    if not 0 <= maximum <= REGION_LENGTH:
        raise ValueError('Invalid regional occupancy')
    if maximum >= 256 and alpha == 1.:
        return region_operators_arb_float(z, maximum)
    base = nonzero_coefficients(z, alpha, min(maximum, 64))
    result = np.array([[[0., -np.inf], [-np.inf, 0.]]])
    exponent = REGION_PACKETS
    while exponent:
        if exponent & 1:
            result = log_polynomial_product(result, base, maximum)
        exponent >>= 1
        if exponent:
            base = log_polynomial_product(base, base, maximum)
    normalizers = np.array([log2(comb(REGION_LENGTH, m)) for m in range(maximum+1)])
    return result-normalizers[:, None, None]


def potential_envelope(regions, log_potential):
    """Least log(u) with R_m v <= u^m v, v=(1,2^log_potential).

    m=0 must be a contraction in this potential. Return the worst positive
    occupancy and the two-state witness; this is an entrywise inequality.
    """
    potential = np.array([0., log_potential])
    contractions = np.logaddexp2(regions[:, :, 0]+potential[0],
                                regions[:, :, 1]+potential[1])-potential
    if np.max(contractions[0]) > 1e-12:
        raise ArithmeticError('The empty region is not a contraction')
    per_coordinate = np.max(contractions[1:], axis=1)/np.arange(1, len(regions))
    worst = int(np.argmax(per_coordinate))+1
    return float(per_coordinate[worst-1]), worst


@lru_cache(None)
def row_weight_counts():
    return tuple((weight, log2(comb(32, weight)*mds_prescribed_support(Q, 32, 16, weight)))
                 for weight in range(17, 33))


def one_row_group_moment(log_support_mark):
    # A nonzero one-row column spreads to exactly seven nonzero outputs.
    return total(count+7*weight*log_support_mark for weight, count in row_weight_counts())


@lru_cache(None)
def two_row_weight_counts():
    return tuple((weight, log2(count)) for weight, count in
                 enumerate(two_active_row_weight_enumerator()) if count)


def exact_group_moment(active, log_support_mark):
    if active == 1:
        return one_row_group_moment(log_support_mark)
    if active == 2:
        return total(count+weight*log_support_mark for weight, count in two_row_weight_counts())
    raise ValueError('Exact MDS enumerators are available for one or two active rows')


def homogeneous_bound(groups, active, z, log_potential, regions=None):
    """First moment for r groups with a fixed active-constituent count.

    All constituent-row choices per group are included. Scalar-potential
    domination uses only the number of occupied regional positions; it does
    not assume independence of output masks or replace a packet by one input.
    """
    if not 1 <= groups <= REGION_LENGTH:
        raise ValueError('Invalid active-group count')
    if regions is None:
        regions = region_operators(z, groups)
    mark, worst = potential_envelope(regions, log_potential)
    value = (log2(comb(REGION_LENGTH, groups))+log2(comb(8, active))*groups
             -log2(Q-1)-LIMIT*log2(z)+max(0., -log_potential)
             +groups*exact_group_moment(active, mark))
    return value, mark, worst


def proposal(groups, iterations=45, active=1):
    start = perf_counter()
    best = None
    def objective(ell):
        nonlocal best
        z = 2**(-ell)
        regions = region_operators(z, groups)
        selected = minimize_scalar(lambda p: homogeneous_bound(groups, active, z, p, regions)[0],
                                   bounds=(-150., 0.), method='bounded', options={'xatol': 1e-7})
        log_potential = float(selected.x)
        value, support_mark, worst = homogeneous_bound(groups, active, z, log_potential, regions)
        if best is None or value < best['log2_upper_screen']:
            best = dict(active_groups=groups, active_rows_per_group=active, mark=z,
                        log2_potential=log_potential, log2_support_mark=support_mark,
                        worst_regional_occupancy=worst, log2_upper_screen=value,
                        exact_mds_support=True, independent_nonzero_amplitudes=True,
                        all_row_set_choices=True, alpha=1., directed=False,
                        whole_code_certificate=False)
        return value
    minimize_scalar(objective, bounds=(1e-6, max(.04, groups/1500+.015)), method='bounded',
                    options={'maxiter': iterations, 'xatol': 1e-8})
    best['elapsed_seconds'] = perf_counter()-start
    return best


def interval_screen(upper, iterations=45, target=-40.):
    """Replay one FLOAT potential for every integer group count up to upper.

    The operator family for upper dominates occupancies at every smaller
    group count. This evaluates each integer count; it does not interpolate.
    """
    witness = proposal(upper, iterations, 1)
    z, log_potential = witness['mark'], witness['log2_potential']
    regions = region_operators(z, upper)
    values = [homogeneous_bound(groups, 1, z, log_potential, regions)[0]
              for groups in range(1, upper+1)]
    lower = upper
    while lower > 1 and values[lower-2] < target:
        lower -= 1
    if values[upper-1] >= target:
        lower = None
    witness.update(covered_integer_interval=[lower, upper] if lower else None,
                   interval_log2_upper_screen=max(values[lower-1:]) if lower else None,
                   requested_log2_target=target)
    return witness


def nonzero_arb_coefficients(z, maximum):
    """Directed version of the nonzero packet caps, with exact support counts."""
    from flint import arb
    from hd_common_refresh_screen import principal_defects
    from hd_grouped_quiet import quiet_cap
    q, h, length, state = Q, 4, 16, 16
    zero, one = arb(0), arb(1)
    zz = arb(z.numerator)/z.denominator
    kernel = []
    for occupancy in range(17):
        probability = nonzero_kernel_probability(occupancy)
        kernel.append(arb(probability.numerator)/probability.denominator)
    quiet = [zz**j*kernel[j] for j in range(17)]
    free = [(zz+(one-zz)/q)**j for j in range(17)]
    uniform_quiet = [sum((comb(j, c)*zz**(j-c)*(one-zz)**c/q**(c+min(h, j-c))
                         for c in range(j+1)), zero) for j in range(17)]
    condition = [(arb(q)/(q-1))**j for j in range(17)]
    ret = [zero]
    defects = principal_defects(h)
    for j in range(1, 17):
        rank_cap = zero
        for c in range(j+1):
            for e in range(min(h, length-j+1)):
                first = arb(q)**(-max(0, c-min(defects[c], max(0, h-j+c))))
                second = arb(q)**(max(0, j-c-h)-j)*(q**min(min(h, j), h-e)-1)
                caps = (one-arb(q)**(-min(h, j)), first, second)
                selected = caps[0]
                for candidate in caps[1:]:
                    if candidate < selected:
                        selected = candidate
                rank_cap += comb(j, c)*comb(length-j, e)*zz**(length-c-e)*(one-zz)**(c+e)*selected
        caps = (condition[j]*(free[j]-uniform_quiet[j]),
                (one-kernel[j])*zz**max(1, length-h+1-j), condition[j]*rank_cap)
        selected = caps[0]
        for candidate in caps[1:]:
            if candidate < selected:
                selected = candidate
        ret.append(selected)
    base = one/q**h
    excess = [sum((comb(n, c)*zz**(n-c)*(one-zz)**c*(one/q**c-base)
                   for c in range(min(h, n+1))), zero) for n in range(17)]
    d = q**state-1
    b = arb(d)/(q-1)**state
    empty_cap = quiet_cap(z, q=q, state=state, width=64, h=h)
    inverse = arb(empty_cap.numerator)/empty_cap.denominator
    result = [[zero for _ in range(maximum+1)] for _ in range(4)]
    for profile, count in profiles(16, 64, 4):
        degree = sum(profile)
        if degree > maximum:
            continue
        packet_ret, packet_total = zero, one
        live_nonzero, live_total, feedback_quiet, free_product = zero, one, one, one
        for occupancy in profile:
            packet_ret = packet_ret*quiet[occupancy]+packet_total*ret[occupancy]
            packet_total *= quiet[occupancy]+ret[occupancy]
            live_nonzero = live_nonzero*base+live_total*excess[length-occupancy]
            live_total *= base+excess[length-occupancy]
            feedback_quiet *= kernel[occupancy]
            free_product *= free[occupancy]
        live = free_product*live_nonzero/(one-one/q**state)*(arb(q)/(q-1))**degree
        if d*live < packet_ret:
            packet_ret = d*live
        actual_live = b*live
        if degree == 0 and inverse < actual_live:
            actual_live = inverse
        entries = (zz**degree*feedback_quiet, b*zz**degree*(one-feedback_quiet),
                   packet_ret/d, actual_live)
        for index, entry in enumerate(entries):
            result[index][degree] += count*entry
    return result


def region_operators_arb_float(z, maximum, precision=256):
    """Fast FLOAT proposal using Arb polynomial multiplication as its engine.

    Returned logarithms are ordinary floats, not directed output. The saved
    witness checker separately replays all coefficients and inequalities.
    """
    from flint import arb, arb_poly, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        coeff = nonzero_arb_coefficients(Fraction(z), min(maximum, 64))
        tilt = arb(maximum)/REGION_LENGTH
        base = [arb_poly([value*tilt**degree for degree, value in enumerate(entry)]) for entry in coeff]
        result = [arb_poly([1]), arb_poly([]), arb_poly([]), arb_poly([1])]
        def multiply_polynomials(left, right):
            return [(left[2*i]*right[j]+left[2*i+1]*right[2+j]).truncate(maximum+1)
                    for i in range(2) for j in range(2)]
        exponent = REGION_PACKETS
        while exponent:
            if exponent & 1:
                result = multiply_polynomials(result, base)
            exponent >>= 1
            if exponent:
                base = multiply_polynomials(base, base)
        output = np.full((maximum+1, 2, 2), -np.inf)
        logarithm_two = arb(2).log()
        log_tilt = log2(maximum/REGION_LENGTH)
        normalizer = 0.
        for occupancy in range(maximum+1):
            if occupancy:
                normalizer += log2((REGION_LENGTH-occupancy+1)/occupancy)+log_tilt
            for index, entry in enumerate(result):
                coefficient = entry[occupancy]
                if coefficient == 0:
                    continue
                if not coefficient > 0:
                    raise ArithmeticError('Arb FLOAT proposal lost a positive coefficient')
                output[occupancy, index//2, index % 2] = float(coefficient.log()/logarithm_two)-normalizer
        return output
    finally:
        ctx.prec = previous


def directed_one_row_receipt(lower, upper, z, potential_bits, support_bits,
                             precision=512):
    """Validate a fixed potential and sum the a=1 class over lower<=r<=upper.

    z is rational, v=(1,2^potential_bits), and u=2^support_bits; exponents are
    integers or fractions. Every 0<=m<=upper is explicitly checked with Arb.
    This certifies only the homogeneous one-active-constituent class.
    """
    from flint import arb, arb_poly, ctx
    if not 1 <= lower <= upper <= REGION_LENGTH:
        raise ValueError('Invalid integer interval')
    z, potential_bits, support_bits = Fraction(z), Fraction(potential_bits), Fraction(support_bits)
    if not 0 < z < 1 or potential_bits > 0:
        raise ValueError('Invalid rational witness')
    previous = ctx.prec
    ctx.prec = precision
    start = perf_counter()
    try:
        zero, one = arb(0), arb(1)
        zz = arb(z.numerator)/z.denominator
        v = arb(2)**(arb(potential_bits.numerator)/potential_bits.denominator)
        u = arb(2)**(arb(support_bits.numerator)/support_bits.denominator)
        coeff = nonzero_arb_coefficients(z, min(upper, 64))
        # Scale the polynomial variable; its coefficient cancels in normalization.
        tilt = arb(upper)/REGION_LENGTH
        base = [arb_poly([value*tilt**degree for degree, value in enumerate(entry)]) for entry in coeff]
        result = [arb_poly([one]), arb_poly([]), arb_poly([]), arb_poly([one])]
        def multiply_polynomials(left, right):
            return [(left[2*i]*right[j]+left[2*i+1]*right[2+j]).truncate(upper+1)
                    for i in range(2) for j in range(2)]
        exponent = REGION_PACKETS
        while exponent:
            if exponent & 1:
                result = multiply_polynomials(result, base)
            exponent >>= 1
            if exponent:
                base = multiply_polynomials(base, base)
        normalizer, support_bound = one, one
        for occupancy in range(upper+1):
            if occupancy:
                normalizer = normalizer*(REGION_LENGTH-occupancy+1)*tilt/occupancy
                support_bound *= u
            entries = [entry[occupancy]/normalizer for entry in result]
            zero_row, live_row = entries[0]+entries[1]*v, entries[2]+entries[3]*v
            # The empty zero row equals one exactly; do not demand strictness.
            if occupancy == 0:
                if not zero_row == one or not live_row < v:
                    raise ArithmeticError('Directed empty-region potential check failed')
            elif not zero_row < support_bound or not live_row < support_bound*v:
                raise ArithmeticError(f'Directed regional potential check failed at occupancy {occupancy}')
        moment = sum((comb(32, weight)*mds_prescribed_support(Q, 32, 16, weight)*u**(7*weight)
                      for weight in range(17, 33)), zero)
        term = comb(REGION_LENGTH, lower)*(8*moment)**lower
        union = term
        for groups in range(lower+1, upper+1):
            term = term*(REGION_LENGTH-groups+1)*(8*moment)/groups
            union += term
        bound = zz**(-LIMIT)*union/(Q-1)/v
        logarithm = bound.log()/arb(2).log()
        if not logarithm < -40:
            raise ArithmeticError(f'Directed one-row class exceeds target: {logarithm}')
        return dict(active_group_interval=[lower, upper], active_rows_per_group=1,
                    mark=str(z), log2_potential=str(potential_bits), log2_support_mark=str(support_bits),
                    failure_log2_upper=str(logarithm), all_row_set_choices=True,
                    exact_mds_support=True, directed=True, whole_code_certificate=False,
                    precision_bits=precision, elapsed_seconds=perf_counter()-start)
    finally:
        ctx.prec = previous


def directed_coverage(maximum=128, iterations=30, target=-80.):
    """Generate contiguous directed receipts for the homogeneous a=1 class.

    Floating proposals only select rational witnesses and integer intervals.
    Every accepted interval is independently replayed by the Arb evaluator.
    """
    upper = maximum
    while upper:
        proposed = proposal(upper, iterations, 1)
        z = Fraction(round(proposed['mark']*10**9), 10**9)
        potential = Fraction(round(proposed['log2_potential']*64), 64)
        regions = region_operators(float(z), upper)
        support, _ = potential_envelope(regions, float(potential))
        support = Fraction(ceil(support*4096)+1, 4096)
        def estimate(groups):
            return (log2(comb(REGION_LENGTH, groups))+3*groups-log2(Q-1)-LIMIT*log2(float(z))
                    -float(potential)+groups*one_row_group_moment(float(support)))
        lower = upper
        while lower > 1 and estimate(lower-1) < target:
            lower -= 1
        if estimate(upper) >= target:
            raise ArithmeticError('Proposal did not close the next integer interval')
        yield directed_one_row_receipt(lower, upper, z, potential, support)
        upper = lower-1


def directed_remaining_one_row(target=-160.):
    """Generate 129..8192 witnesses from bounded interpolated FLOAT proposals.

    Every proposed interval is replayed. Interpolation only selects z; it
    never certifies an omitted count or a regional potential inequality.
    """
    anchors = ((128, .967753099), (256, .947219946), (512, .910213757),
               (1024, .833836287), (2048, .668565828), (4096, .341165697),
               (8192, .022463076))
    positions = np.array([log2(r) for r, z in anchors])
    tilts = np.array([log2(-log2(z)) for r, z in anchors])
    upper = REGION_LENGTH
    while upper >= 129:
        ell = 2**float(np.interp(log2(upper), positions, tilts))
        z = Fraction(round(2**(-ell)*10**9), 10**9)
        regions = region_operators(float(z), upper)
        chosen = minimize_scalar(lambda p: homogeneous_bound(upper, 1, float(z), p, regions)[0],
                                 bounds=(-150., 0.), method='bounded', options={'xatol': 1e-7})
        potential = Fraction(round(float(chosen.x)*64), 64)
        support, _ = potential_envelope(regions, float(potential))
        support = Fraction(ceil(support*4096)+1, 4096)
        slope = 3+one_row_group_moment(float(support))
        prefix = -log2(Q-1)-LIMIT*log2(float(z))-float(potential)
        choose = log2(comb(REGION_LENGTH, upper))
        if prefix+choose+upper*slope >= target:
            raise ArithmeticError(f'Interpolated proposal did not close at {upper} groups')
        lower = upper
        while lower > 129:
            previous_choose = choose+log2(lower/(REGION_LENGTH-lower+1))
            if prefix+previous_choose+(lower-1)*slope >= target:
                break
            choose = previous_choose
            lower -= 1
        yield directed_one_row_receipt(lower, upper, z, potential, support)
        upper = lower-1


SAVED_ONE_ROW_128 = (
    (1, 1, '62078423/62500000', '-4011/64', '-139041/2048'),
    (2, 2, '248213319/250000000', '-4013/64', '-37191/1024'),
    (3, 3, '992487709/1000000000', '-2007/32', '-105081/4096'),
    (4, 4, '992146531/1000000000', '-4015/64', '-81/4'),
    (5, 5, '993179513/1000000000', '-1703/32', '-56503/4096'),
    (6, 6, '496185589/500000000', '-3757/64', '-53625/4096'),
    (7, 7, '123944813/125000000', '-4025/64', '-25795/2048'),
    (8, 8, '247736031/250000000', '-2011/32', '-48771/4096'),
    (9, 10, '30949121/31250000', '-1005/16', '-41783/4096'),
    (11, 12, '247460353/250000000', '-2011/32', '-4621/512'),
    (13, 14, '989330251/1000000000', '-4023/64', '-16729/2048'),
    (15, 16, '197766979/200000000', '-4025/64', '-15391/2048'),
    (17, 19, '494058103/500000000', '-4027/64', '-27761/4096'),
    (20, 22, '197484299/200000000', '-4029/64', '-25511/4096'),
    (23, 25, '986746343/1000000000', '-4031/64', '-23765/4096'),
    (26, 29, '246467873/250000000', '-4033/64', '-21957/4096'),
    (30, 33, '246255149/250000000', '-4035/64', '-20553/4096'),
    (34, 37, '984189903/1000000000', '-4037/64', '-9713/2048'),
    (38, 42, '491587701/500000000', '-4039/64', '-18293/4096'),
    (43, 48, '981987363/1000000000', '-2021/32', '-269/64'),
    (49, 54, '490413149/500000000', '-4045/64', '-4089/1024'),
    (55, 61, '979500323/1000000000', '-253/4', '-15545/4096'),
    (62, 68, '244550111/250000000', '-4051/64', '-7441/2048'),
    (69, 76, '976742249/1000000000', '-2027/32', '-14259/4096'),
    (77, 85, '975131317/1000000000', '-4057/64', '-13681/4096'),
    (86, 95, '486686529/500000000', '-4061/64', '-6577/2048'),
    (96, 105, '485821569/500000000', '-4065/64', '-12715/4096'),
    (106, 116, '969768633/1000000000', '-4069/64', '-3077/1024'),
    (117, 128, '967753099/1000000000', '-4073/64', '-11935/4096'),
)


SAVED_ONE_ROW_REMAINING = (
    (129, 133, '60428651/62500000', '-2043/32', '-5907/2048'),
    (134, 146, '964577841/1000000000', '-2053/32', '-11523/4096'),
    (147, 160, '962188979/1000000000', '-2059/32', '-5623/2048'),
    (161, 175, '191939699/200000000', '-4125/64', '-10985/4096'),
    (176, 191, '957112443/1000000000', '-2065/32', '-10739/4096'),
    (192, 207, '954591843/1000000000', '-4133/64', '-10521/4096'),
    (208, 225, '951826753/1000000000', '-517/8', '-5151/2048'),
    (226, 244, '474490701/500000000', '-517/8', '-10095/4096'),
    (245, 264, '118243623/125000000', '-1035/16', '-4957/2048'),
    (265, 285, '942657427/1000000000', '-4149/64', '-2441/1024'),
    (286, 307, '58704629/62500000', '-1039/16', '-4811/2048'),
    (308, 330, '93580237/100000000', '-2081/32', '-593/256'),
    (331, 355, '37283953/40000000', '-521/8', '-2339/1024'),
    (356, 381, '464159337/500000000', '-4173/64', '-9231/4096'),
    (382, 408, '462232307/500000000', '-2089/32', '-9113/4096'),
    (409, 437, '36816013/40000000', '-2091/32', '-2249/1024'),
    (438, 468, '458067891/500000000', '-2091/32', '-4441/2048'),
    (469, 500, '911814331/1000000000', '-261/4', '-8773/4096'),
    (501, 534, '906729459/1000000000', '-525/8', '-8711/4096'),
    (535, 569, '225307081/250000000', '-1055/16', '-4337/2048'),
    (570, 606, '895467381/1000000000', '-529/8', '-8639/4096'),
    (607, 645, '889453781/1000000000', '-4243/64', '-2151/1024'),
    (646, 685, '441673337/500000000', '-4253/64', '-8573/4096'),
    (686, 727, '21924959/25000000', '-533/8', '-8543/4096'),
    (728, 771, '870416163/1000000000', '-1069/16', '-8515/4096'),
    (772, 818, '863460277/1000000000', '-67', '-1061/512'),
    (819, 867, '171257727/200000000', '-4301/64', '-8463/4096'),
    (868, 919, '106095621/125000000', '-2157/32', '-4219/2048'),
    (920, 974, '840902057/1000000000', '-541/8', '-8415/4096'),
    (975, 1033, '83230823/100000000', '-543/8', '-2101/1024'),
    (1034, 1095, '205451357/250000000', '-1091/16', '-8457/4096'),
    (1096, 1160, '810843861/1000000000', '-4385/64', '-133/64'),
    (1161, 1230, '399551437/500000000', '-551/8', '-4285/2048'),
    (1231, 1304, '196692619/250000000', '-277/4', '-8629/4096'),
    (1305, 1383, '1208911/1562500', '-4459/64', '-4345/2048'),
    (1384, 1468, '18994127/25000000', '-561/8', '-8755/4096'),
    (1469, 1559, '74499257/100000000', '-4519/64', '-4411/2048'),
    (1560, 1658, '729106721/1000000000', '-2277/32', '-4447/2048'),
    (1659, 1765, '356082519/500000000', '-287/4', '-8969/4096'),
    (1766, 1882, '693922437/1000000000', '-4633/64', '-9049/4096'),
    (1883, 2011, '42135193/62500000', '-4679/64', '-4567/2048'),
    (2012, 2152, '649272129/1000000000', '-2369/32', '-4661/2048'),
    (2153, 2306, '310521181/500000000', '-4807/64', '-4779/2048'),
    (2307, 2476, '590441297/1000000000', '-2441/32', '-2453/1024'),
    (2477, 2664, '55739929/100000000', '-1241/16', '-2521/1024'),
    (2665, 2874, '10432207/20000000', '-5051/64', '-10381/4096'),
    (2875, 3111, '241390251/500000000', '-5149/64', '-2677/1024'),
    (3112, 3384, '2201287/5000000', '-5249/64', '-5539/2048'),
    (3385, 3704, '393574657/1000000000', '-1337/16', '-719/256'),
    (3705, 4087, '171154843/500000000', '-2723/32', '-6003/2048'),
    (4088, 4542, '273105249/1000000000', '-5565/64', '-13025/4096'),
    (4543, 5106, '50173871/250000000', '-5679/64', '-1791/512'),
    (5107, 5848, '128004237/1000000000', '-1445/16', '-16107/4096'),
    (5849, 6894, '62457379/1000000000', '-5839/64', '-18705/4096'),
    (6895, 8192, '5615769/250000000', '0', '-21473/4096'),
)


SAVED_LARGE_ONE_ROW_INTERVALS = (
    (972, 1024, '833836287/1000000000', '-2171/32', '-41/20'),
    (1924, 2048, '167141457/250000000', '-1173/16', '-11181/5000'),
)


def replay_one_row_receipts(witnesses, maximum, precision=512):
    """Replay and sum contiguous saved witnesses; consult no FLOAT optimizer."""
    from flint import arb, ctx
    previous = ctx.prec
    ctx.prec = precision
    try:
        receipts, next_group, union = [], 1, arb(0)
        for lower, upper, z, potential, support in witnesses:
            if lower != next_group:
                raise ArithmeticError('Saved interval coverage has a gap or overlap')
            receipt = directed_one_row_receipt(lower, upper, z, potential, support, precision)
            union += arb(2)**arb(receipt['failure_log2_upper'])
            receipts.append(receipt)
            next_group = upper+1
        if next_group != maximum+1:
            raise ArithmeticError('Saved interval coverage is incomplete')
        logarithm = union.log()/arb(2).log()
        if not logarithm < -97:
            raise ArithmeticError('Saved homogeneous class receipt exceeded its bound')
        return dict(active_group_interval=[1, maximum], active_rows_per_group=1,
                    interval_count=len(receipts), failure_log2_upper=str(logarithm),
                    scalar_route=True, field_order=Q, weight_limit=LIMIT,
                    all_row_set_choices=True, exact_mds_support=True,
                    directed=True, whole_code_certificate=False,
                    precision_bits=precision)
    finally:
        ctx.prec = previous


def replay_one_row_128(precision=512):
    return replay_one_row_receipts(SAVED_ONE_ROW_128, 128, precision)


def replay_one_row_all(precision=512):
    return replay_one_row_receipts(SAVED_ONE_ROW_128+SAVED_ONE_ROW_REMAINING,
                                   REGION_LENGTH, precision)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups', type=int, nargs='+', default=(2, 4, 16, 128))
    parser.add_argument('--active', type=int, nargs='+', choices=(1, 2), default=(1, 2))
    parser.add_argument('--iterations', type=int, default=45)
    parser.add_argument('--intervals', action='store_true')
    parser.add_argument('--directed-through', type=int)
    parser.add_argument('--replay-one-row-128', action='store_true')
    parser.add_argument('--replay-one-row-all', action='store_true')
    parser.add_argument('--replay-large-intervals', action='store_true')
    args = parser.parse_args()
    if args.replay_one_row_all:
        print(json.dumps(replay_one_row_all()), flush=True)
        raise SystemExit(0)
    if args.replay_one_row_128:
        print(json.dumps(replay_one_row_128()), flush=True)
        raise SystemExit(0)
    if args.replay_large_intervals:
        for witness in SAVED_LARGE_ONE_ROW_INTERVALS:
            print(json.dumps(directed_one_row_receipt(*witness)), flush=True)
        raise SystemExit(0)
    if args.directed_through:
        for receipt in directed_coverage(args.directed_through, args.iterations):
            print(json.dumps(receipt), flush=True)
        raise SystemExit(0)
    for active in args.active:
        for groups in args.groups:
            if args.intervals and active != 1:
                continue
            result = interval_screen(groups, args.iterations) if args.intervals else proposal(groups, args.iterations, active)
            print(json.dumps(result), flush=True)
