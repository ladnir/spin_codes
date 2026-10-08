"""Positive outer-freshness transfer and a profile-free packet envelope.

The exact helpers accept Fractions.  ``floating_log_bound`` is only a screen,
not a directed certificate.  Packet input matrices must already satisfy the
conditional entrywise contract in OUTER_FRESHNESS_TRANSFER.md.
"""

from fractions import Fraction
from itertools import product
from math import comb, exp, expm1, factorial, log, log1p, log2


def row_polynomial(fresh_marks, dependent_marks, parity):
    """Sum valid row supports; scan backwards, accepting empty or >p active."""
    if len(fresh_marks) != len(dependent_marks) or parity < 0:
        raise ValueError("equal mark lengths and nonnegative parity required")
    states = [0] * (parity + 2)
    states[0] = 1
    for fresh, dependent in zip(reversed(fresh_marks), reversed(dependent_marks)):
        nxt = states.copy()  # inactive coordinate
        for count in range(parity):
            nxt[count + 1] += states[count] * dependent
        nxt[parity + 1] += (states[parity] + states[parity + 1]) * fresh
        states = nxt
    return states[0] + states[parity + 1]


def valid_row_patterns(length, parity):
    """Exhaustive reference only: 0=inactive, 1=fresh, 2=dependent."""
    yield (0,) * length
    for mask in product((0, 1), repeat=length):
        active = [j for j, bit in enumerate(mask) if bit]
        if len(active) <= parity:
            continue
        split = len(active) - parity
        result = [0] * length
        for j in active[:split]:
            result[j] = 1
        for j in active[split:]:
            result[j] = 2
        yield tuple(result)


def ternary_count(size, fresh, dependent):
    if fresh < 0 or dependent < 0 or fresh + dependent > size:
        return 0
    return comb(size, fresh) * comb(size - fresh, dependent)


def identity(size):
    return [[int(i == j) for j in range(size)] for i in range(size)]


def matmul(left, right):
    size = len(left)
    return [[sum(left[i][k] * right[k][j] for k in range(size))
             for j in range(size)] for i in range(size)]


def matpow(matrix, power):
    if power < 0:
        raise ValueError("nonnegative matrix exponent required")
    result = identity(len(matrix))
    while power:
        if power & 1:
            result = matmul(result, matrix)
        power >>= 1
        if power:
            matrix = matmul(matrix, matrix)
    return result


def packet_envelope(packet_sums, width, fresh_mark, dependent_mark):
    """Entrywise max K_fd / (multinomial(width;f,d,*) v**f w**d).

    ``packet_sums`` maps (f,d) to a nonnegative square matrix.  Integer
    normalizers are converted to Fraction so exact rational inputs stay exact.
    Missing (f,d) pairs are interpreted as zero, not as unknown matrices.
    """
    if fresh_mark <= 0 or dependent_mark <= 0 or not packet_sums:
        raise ValueError("positive marks and nonempty packet table required")
    size = len(next(iter(packet_sums.values())))
    result = [[0] * size for _ in range(size)]
    for (fresh, dependent), matrix in packet_sums.items():
        count = ternary_count(width, fresh, dependent)
        if not count:
            raise ValueError("invalid packet counts")
        denominator = (Fraction(count) * fresh_mark**fresh
                       * dependent_mark**dependent)
        for i in range(size):
            for j in range(size):
                if matrix[i][j] < 0:
                    raise ValueError("packet operator must be nonnegative")
                result[i][j] = max(result[i][j], matrix[i][j] / denominator)
    return result


def packet_scalar_envelope(packet_sums, width, fresh_mark, dependent_mark,
                           continuation):
    """Smallest rho with N_fd V <= rho v**f w**d V, for fixed V>0."""
    if min(continuation) <= 0 or fresh_mark <= 0 or dependent_mark <= 0:
        raise ValueError("positive continuation and marks required")
    result = 0
    for (fresh, dependent), matrix in packet_sums.items():
        denominator = (Fraction(ternary_count(width, fresh, dependent))
                       * fresh_mark**fresh * dependent_mark**dependent)
        for i, row in enumerate(matrix):
            result = max(result, sum(entry*v for entry, v in zip(row, continuation))
                         / (denominator * continuation[i]))
    return result


def worst_modal_counts(rows, active_max):
    """Counts minimizing multinomial mass at its mean, with f+d<=active_max."""
    if rows <= 0 or not 0 <= active_max <= rows:
        raise ValueError("invalid active-row cap")
    zero = max(rows-active_max, rows//3)
    fresh = (rows-zero)//2
    return fresh, rows-zero-fresh, zero


def exact_conditioning_cost(rows, active_max):
    counts = worst_modal_counts(rows, active_max)
    numerator = rows**rows
    denominator = factorial(rows)
    for count in counts:
        numerator *= factorial(count)
        denominator *= count**count  # Python's 0**0=1 is the desired convention.
    return Fraction(numerator, denominator)


def thinned_coefficients(coefficients, width, density):
    """Exact Bernstein coefficients after restricting u+d<=density.

    coefficients[(f,d)] is a scalar full-simplex Bernstein coefficient.
    The new triangle has vertices (0,0), (density,0), (0,density).
    """
    density = Fraction(density)
    if not 0 <= density <= 1:
        raise ValueError("density must lie in [0,1]")
    result = {}
    for i in range(width+1):
        for j in range(width-i+1):
            result[i, j] = sum(
                comb(i, f)*comb(j, d)*density**(f+d)
                *(1-density)**(i+j-f-d)*coefficients[f, d]
                for f in range(i+1) for d in range(j+1))
    return result


def exact_bound(packet_sums, *, rows, length, parity, width, tau, x, y,
                fresh_marks, dependent_marks, input_max, output_max,
                input_min=1, continuation=None):
    """Exact un-clipped bound from the theorem, for rational input data.

    The packet table is evaluated at tau,y before this call.  No validity of
    that table is inferred here.  Nonempty support is excluded BEFORE the
    matrix envelope is multiplied, using H**rows-1.
    """
    if rows <= 0 or rows % width or length <= parity or parity < 0:
        raise ValueError("invalid outer shape or packet alignment")
    if len(fresh_marks) != length or len(dependent_marks) != length:
        raise ValueError("one pair of marks per region required")
    if tau < 1 or x <= 0 or not 0 < y <= 1:
        raise ValueError("require tau>=1, x>0, and 0<y<=1")
    if not 1 <= input_min <= input_max <= rows * length:
        raise ValueError("invalid nonempty input interval")
    if not 0 <= output_max <= rows * length:
        raise ValueError("invalid output cutoff")
    tau, x, y = Fraction(tau), Fraction(x), Fraction(y)
    h = row_polynomial([tau*x*v for v in fresh_marks],
                       [x*w for w in dependent_marks], parity)
    if continuation is None:
        transfer = identity(len(next(iter(packet_sums.values()))))
        for v, w in zip(fresh_marks, dependent_marks):
            envelope = packet_envelope(packet_sums, width, v, w)
            transfer = matmul(transfer, matpow(envelope, rows // width))
        inner_bound = sum(transfer[0])
    else:
        inner_bound = Fraction(continuation[0]) / min(continuation)
        for v, w in zip(fresh_marks, dependent_marks):
            rho = packet_scalar_envelope(packet_sums, width, v, w, continuation)
            inner_bound *= rho**(rows // width)
    input_endpoint = input_max if x <= 1 else input_min
    return ((h**rows - 1) * inner_bound
            / (tau * x**input_endpoint * y**output_max))


def exact_rank_packet_table(tau, y, q=256, width=16, state=4, law='unscaled'):
    """Adapter to the separate, exact integer packet coefficient table."""
    from rank_packet_operator import coefficient_table
    tau, y = Fraction(tau), Fraction(y)
    if not 1 <= tau <= q - 1 or not 0 < y <= 1:
        raise ValueError("invalid rank/output tilts")
    result = {}
    for (f, d, v, r, k, charge, exceptional), multiplicity in coefficient_table(
            width, state, law).items():
        matrix = result.setdefault((f, d), [[Fraction(0) for _ in range(state + 1)]
                                           for _ in range(state + 1)])
        matrix[r][k] += (multiplicity * y**v / tau**charge
                         / (q - 1)**exceptional)
    return result


def _logadd(a, b):
    if a == float('-inf'):
        return b
    if b == float('-inf'):
        return a
    maximum = max(a, b)
    return maximum + log1p(exp(min(a, b) - maximum))


def log_row_polynomial(fresh_logs, dependent_logs, parity):
    states = [float('-inf')] * (parity + 2)
    states[0] = 0.0
    for fresh, dependent in zip(reversed(fresh_logs), reversed(dependent_logs)):
        nxt = states.copy()
        for count in range(parity):
            nxt[count + 1] = _logadd(nxt[count + 1], states[count] + dependent)
        nxt[parity + 1] = _logadd(nxt[parity + 1],
                                _logadd(states[parity], states[parity + 1]) + fresh)
        states = nxt
    return _logadd(states[0], states[parity + 1])


def floating_log_bound(*, rows, length=32, parity=10, width=16, state=4,
                       q=256, tau, x, y, fresh_marks, dependent_marks,
                       input_max, output_max, input_min=1, law='unscaled',
                       continuation=None, bellman=False):
    """Exploratory log2 upper expression; no directed rounding guarantee."""
    import numpy as np
    from rank_packet_operator import packet_action, packet_matrices
    if rows % width or len(fresh_marks) != length or len(dependent_marks) != length:
        raise ValueError("invalid alignment or region marks")
    if bellman and continuation is None:
        raise ValueError("Bellman action requires a positive continuation")
    if not bellman:
        normalized = packet_matrices(tau, y, q=q, t=width, s=state,
                                     normalized=True, law=law)
    logs = log_row_polynomial([log(tau*x*v) for v in fresh_marks],
                             [log(x*w) for w in dependent_marks], parity)
    exponent = rows * logs
    if exponent == 0:
        # Underflow of the nonempty row term is not a valid zero probability.
        return float('inf')
    outer_log = exponent + log(-expm1(-exponent))
    if continuation is not None:
        potential = np.asarray(continuation, dtype=float)
        if min(potential) <= 0:
            raise ValueError("positive continuation required")
        if bellman:
            actions = packet_action(tau, y, potential, q=q, t=width, s=state,
                                    normalized=True, law=law) / potential
        else:
            actions = np.einsum('fdrk,k->fdr', normalized, potential) / potential
        transfer_log = log(potential[0] / min(potential))
        for v, w in zip(fresh_marks, dependent_marks):
            rho = max(np.max(actions[f, d]) / (v**f * w**d)
                      for f in range(width + 1) for d in range(width - f + 1))
            transfer_log += rows / width * log(rho)
        endpoint = input_max if x <= 1 else input_min
        return ((outer_log + transfer_log - log(tau) - endpoint*log(x)
                 - output_max*log(y)) / log(2))
    vector = np.zeros(state + 1)
    vector[0] = 1.0
    transfer_log = 0.0
    for v, w in zip(fresh_marks, dependent_marks):
        matrix = np.zeros((state + 1, state + 1))
        for f in range(width + 1):
            for d in range(width - f + 1):
                matrix = np.maximum(matrix, normalized[f, d] / (v**f * w**d))
        scale = np.max(matrix)
        if not scale > 0 or not np.isfinite(scale):
            return float('inf')
        matrix /= scale
        transfer_log += rows / width * log(scale)
        exponent = rows // width
        while exponent:
            if exponent & 1:
                vector = vector @ matrix
                scale_v = np.max(vector)
                if not scale_v > 0 or not np.isfinite(scale_v):
                    return float('inf')
                vector /= scale_v
                transfer_log += log(scale_v)
            exponent >>= 1
            if exponent:
                matrix = matrix @ matrix
                scale_m = np.max(matrix)
                if not scale_m > 0 or not np.isfinite(scale_m):
                    return float('inf')
                matrix /= scale_m
                # matrix now stands for 2**bit original packets.  Its removed
                # scale must be multiplied by every remaining use.
                transfer_log += exponent * log(scale_m)
    endpoint = input_max if x <= 1 else input_min
    return ((outer_log + transfer_log + log(sum(vector)) - log(tau)
             - endpoint*log(x) - output_max*log(y)) / log(2))


def replay_screen():
    """Replay fixed FLOAT witnesses, including the entire remaining interval.

    These are upper-envelope evaluations, all positive.  The bounded search
    producing the witnesses was not a global optimization.
    """
    import math
    import numpy as np
    from scipy.special import logsumexp
    from rank_packet_operator import packet_matrices
    from return_cover import return_screen_objective
    rows, length = 47680, 32
    total = rows * length
    witnesses = (
        (False, [1.8106895908903058, -2.0559615454011655, -0.7602128371775881,
                 -0.11335589191212625, 0.018454490115716423,
                 0.7889324057065451, 0.37491928958505377,
                 0.09170883409136563, 0.010368516023640244]),
        (True, [2.0740253570659575, -2.0838225724959565, -0.7594901176774869,
                -0.36805644550911715, 0.04871528337710086,
                0.5518347171952381, 0.08438030873985185,
                0.019267177984175904, 0.00013850427068583498]))
    result = {'rigorous_numeric_certificate': False, 'law': 'eta-scaled',
              'rows': rows, 'N': total, 'input_interval': [1, 289894],
              'output_max': 152575, 'whole_interval': [], 'diagnostics': []}
    for structured, logarithms in witnesses:
        tau, x, y, fresh, dependent = np.exp(logarithms[:5])
        continuation = np.exp([0] + logarithms[5:])
        fs, ds = [fresh]*length, [dependent]*length
        if structured:
            fs[22:] = [math.exp(12)]*10
            ds[0] = math.exp(12)
        value = floating_log_bound(rows=rows, tau=tau, x=x, y=y,
            fresh_marks=fs, dependent_marks=ds, input_min=1,
            input_max=289894, output_max=152575, law='eta-scaled',
            continuation=continuation)
        result['whole_interval'].append(dict(endpoint_aware_marks=structured,
            log2_upper_per_N=value/total, log_parameters=logarithms))
    root_witnesses = (
        (.01, [-3.0258570984074526, -.5857953682093398, 2.5381883850844997,
               2.014679343059534, 0]),
        (.10, [.70708872410493, -1.9100611711286484, 2.6167343763275186,
               -.2683299405157173, 0]),
        (.18, [2.1865170506895493, -3.3950117622093154, 3.492359740395451,
               -1.3396605412597715, 0]))
    for alpha, values in root_witnesses:
        lx, ly, lt, ll, le = values
        old = return_screen_objective(values, alpha)
        row_log = logsumexp([32*math.log1p(-alpha), -10*lt-11*le
            +32*np.logaddexp(math.log1p(-alpha), lt+ll+le+math.log(alpha))])
        old_packet = 16*(old*math.log(2)+alpha*lx+.1*ly+alpha*ll-row_log/32)
        diagnostic = dict(alpha=alpha, return_only_per_N=float(old),
                          finite_boundary_factors_included=False)
        for law in ('unscaled', 'eta-scaled'):
            packet = packet_matrices(math.exp(lt), math.exp(ly), law=law)
            matrix = sum(packet[0, d]*math.exp(d*lx) for d in range(17))
            radius_log = math.log(max(abs(np.linalg.eigvals(matrix))))
            delta = (radius_log-old_packet)/(16*math.log(2))
            diagnostic[law] = dict(packet_delta_per_N=float(delta),
                                   combined_per_N=float(old+delta))
        result['diagnostics'].append(diagnostic)
    return result


def active_row_log_bound(*, rows, active_min, active_max, tau, x, y,
                         fresh_marks, dependent_marks, continuation,
                         input_max, output_max, input_min=1,
                         length=32, parity=10, width=16, state=4, q=256,
                         components=False):
    """Count-aware FLOAT bound for an interval of active outer row counts.

    Uses eta Bellman action, exact canonical row polynomial, worst modal
    conditioning probability, and a Bernstein upper bound on the iid simplex.
    No optimization over iid probabilities is hidden in the evaluation.
    """
    import numpy as np
    from scipy.special import gammaln, logsumexp
    from rank_packet_operator import packet_action
    if not 1 <= active_min <= active_max <= rows or rows % width:
        raise ValueError("invalid active interval or alignment")
    if len(fresh_marks) != length or len(dependent_marks) != length:
        raise ValueError("one pair of marks per region required")
    low_weight = max(input_min, (parity+1)*active_min)
    high_weight = min(input_max, length*active_max)
    if low_weight > high_weight:
        return float('-inf')
    potential = np.asarray(continuation, dtype=float)
    if min(potential) <= 0:
        raise ValueError("positive continuation required")
    log_h = log_row_polynomial([log(tau*x*v) for v in fresh_marks],
                              [log(x*w) for w in dependent_marks], parity)
    if log_h == 0:
        return float('inf')
    log_nonempty = log_h+log(-expm1(-log_h))
    active = np.arange(active_min, active_max+1)
    log_source = logsumexp(gammaln(rows+1)-gammaln(active+1)
                          -gammaln(rows-active+1)+active*log_nonempty)
    counts = worst_modal_counts(rows, active_max)
    log_condition = rows*log(rows)-gammaln(rows+1)
    for count in counts:
        log_condition += gammaln(count+1)-(count*log(count) if count else 0)
    normalized = packet_action(tau, y, potential, q=q, t=width, s=state,
                               normalized=True, law='eta-scaled') / potential
    density = active_max/rows
    thin = np.zeros((width+1, width+1))
    for i in range(width+1):
        for f in range(i+1):
            thin[i, f] = comb(i, f)*density**f*(1-density)**(i-f)
    packet_log = 0.0
    for region, (v, w) in enumerate(zip(fresh_marks, dependent_marks)):
        marks = np.array([[v**f*w**d for d in range(width+1)]
                          for f in range(width+1)])
        coefficients = normalized / marks[:, :, None]
        thinned = np.empty_like(coefficients)
        for r in range(state+1):
            thinned[:, :, r] = thin @ coefficients[:, :, r] @ thin.T
        choices = [np.max(thinned[f, d])
                   for f in range(width+1) for d in range(width-f+1)
                   if not (region >= length-parity and f)
                   and not (region == 0 and d)]
        rho = max(choices)
        packet_log += rows/width*log(rho)
    endpoint = high_weight if x <= 1 else low_weight
    terms = dict(source=log_source, conditioning=length*log_condition,
                 packets=packet_log,
                 boundary=log(potential[0]/min(potential))-log(tau),
                 input_cut=-endpoint*log(x), output_cut=-output_max*log(y))
    value = sum(terms.values())/log(2)
    if components:
        return value, {key: term/log(2) for key, term in terms.items()}
    return value


def replay_active_screen():
    """Replay all five active-count bins and three fixed-weight diagnostics."""
    import numpy as np
    rows, total = 47680, 1525760
    bins = (
        ((1, 16), [2.092246605548231, -1.7106558992926244, -.4970239882342654,
                   .002884954145168526, -.0023481499537360215,
                   -.7348581050241048, -2.5839350128635656,
                   -4.877478149009899, -3.932713761790888]),
        ((17, 256), [1.533327466150359, -2.284738084588092, -.5008016466438501,
                     -.00010063481977929439, .0008658097456542128,
                     -1.0720021453795918, -2.5916759564740666,
                     -2.5584199729937893, -4.724937188307415]),
        ((257, 4096), [3.312645070448044, -.00005065327999609653,
                      -1.2530688449140495, -3.088177969470618, -2.5052720059691076,
                      -7.508203838371078, -9.658149126739865,
                      -11.995655427674256, -14.60719179800419]),
        ((4097, 16384), [.5629192509188128, -.3081108708447371, -.8963963131537054,
                        .0016940058290952718, .0013275313355073804,
                        -1.7928418338324343, -2.6686716755368654,
                        -3.2020587435728167, -3.5864610269917865]),
        ((16385, 26354), [1.0532798985809442, -1.8891910941501653,
                         -.647137485207861, .001612018755553517, .003943287905601504,
                         -1.4842149734692098, -2.041552873528861,
                         -2.780424104926073, -2.0311646954519187]))
    fixed = (
        (4768, [1.769015009722092, 9.999558253322135, -1.9099350911744848,
                -2.31641928268981, -2.314045047096304, -8.649632394063095,
                -8.449662604912499, -10.867221022261223, -10.408518475532857]),
        (9536, [1.6631808911146955, -.07099525084031973, -1.9100738118427347,
                -1.6163954162044143, -1.615340853492173, -6.184886777082708,
                -7.042103323809984, -10.867212502556427, -9.865738811885764]),
        (13870, [1.5690410671017223, -8.639779304282339, -1.910013739717808,
                 -1.2348796626768155, -1.2347167092257358, -6.192243964123771,
                 -6.997424758224934, -10.867211504291438, -9.726215527442118]))
    result = {'rigorous_numeric_certificate': False, 'law': 'eta-scaled',
              'rows': rows, 'N': total, 'input_interval': [1, 289894],
              'output_max': 152575, 'active_bins': [], 'fixed_weight': []}
    for section, entries in (('active_bins', bins), ('fixed_weight', fixed)):
        for selected, logarithms in entries:
            lo, hi = selected if section == 'active_bins' else (selected, selected)
            wmin, wmax = (1, 289894) if section == 'active_bins' else (152576, 152576)
            tau, x, y, fresh, dependent = np.exp(logarithms[:5])
            value, parts = active_row_log_bound(rows=rows, active_min=lo,
                active_max=hi, tau=tau, x=x, y=y, fresh_marks=[fresh]*32,
                dependent_marks=[dependent]*32,
                continuation=np.exp([0]+logarithms[5:]),
                input_min=wmin, input_max=wmax, output_max=152575, components=True)
            result[section].append(dict(active_interval=[lo, hi],
                input_interval=[wmin, wmax], log2_upper_per_N=value/total,
                components_per_N={key: part/total for key, part in parts.items()},
                log_parameters=logarithms))
    return result


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--screen', action='store_true',
                        help='replay fixed exploratory witnesses; no search')
    parser.add_argument('--active-screen', action='store_true',
                        help='replay count-aware active-row witnesses; no search')
    args = parser.parse_args()
    if args.screen:
        print(json.dumps(replay_screen(), indent=2))
    if args.active_screen:
        print(json.dumps(replay_active_screen(), indent=2))
