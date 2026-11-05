"""Bounded FLOAT screens for the literal extension-field paired-route lift.

The fixed GF256 maps embed in GF(2**m), with m a multiple of eight. Outer
invertible Toeplitz4 matrices and both refresh diagonals are sampled over
that entire field. Routing, state persistence, and all dimensions stay fixed.
Only the requested active-group weights are screened. No output from this
module is directed arithmetic or a whole-code distance certificate.

Local and polynomial entries remain in log2 form, including before alpha.
In particular q**16 is never converted to a binary float. Existing helpers
are imported for field-independent geometry and are never monkeypatched.
"""

from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb, expm1, isfinite, lgamma, log, log2
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize

from hd_common_refresh_packets import principal_defects
from hd_grouped_packet_screen import OCC
from hd_grouped_quiet import block_support_counts, weight_transition
from hd_randomized_outer_state_profile import multiply, power, total
from hd_shifted_good_pair_screen import good_geometry
from hd_shifted_pairs_count import profile_count_upper
from hd_shifted_swap_sparse import sparse_geometry


FIELD_ORDER = 1 << 128
GROUPS, PAIRS, REGIONS, STEPS = 8192, 4096, 256, 128
DIMENSION, OUTPUT_LENGTH = 1 << 20, 1 << 21
MAXIMUM = 64
NEG = -float('inf')
LN2 = log(2)
OUTER_MEASURES = ('inherited-diagonal', 'invertible-toeplitz')


@dataclass
class Budget:
    """Cooperative wall-time guard, checked between bounded numeric kernels."""

    max_seconds: float = 60.

    def __post_init__(self):
        if not isfinite(self.max_seconds) or not 0 < self.max_seconds <= 3600:
            raise ValueError('max_seconds must lie in (0, 3600]')
        self.started = monotonic()

    def check(self):
        if monotonic()-self.started >= self.max_seconds:
            raise TimeoutError('Bounded FLOAT screen exhausted max_seconds')


def checked_q(q, *, construction=False):
    if not isinstance(q, int) or q < 2 or q & (q-1):
        raise ValueError('q must be an exact characteristic-two field order')
    bits = q.bit_length()-1
    if bits > 128:
        raise ValueError('This bounded implementation supports field bits <= 128')
    if construction and bits % 8:
        raise ValueError('The fixed GF256 maps require field bits divisible by eight')
    return bits


def log1mexp2(value):
    """log2(1-2**value), retaining accuracy near zero and at -2048."""
    if value > 0:
        raise ValueError('log1mexp2 requires a nonpositive argument')
    return NEG if value == 0 else log(-expm1(value*LN2))/LN2


def log_fraction(value):
    return NEG if not value else log2(value.numerator)-log2(value.denominator)


def log_multinomial(counts):
    return (lgamma(sum(counts)+1)-sum(lgamma(n+1) for n in counts))/LN2


def outer_log_density(q=FIELD_ORDER, mode='invertible-toeplitz'):
    """Log density relative to the uniform 256-symbol output measure."""
    bits = checked_q(q, construction=True)
    if mode not in OUTER_MEASURES:
        raise ValueError('Unknown established outer-measure bound')
    # Writing log2(q-1)=bits+log2(1-1/q) avoids cancellation of large terms.
    correction = log1mexp2(-bits)
    return 128*bits-(128 if mode == 'inherited-diagonal' else 32)*correction


@lru_cache(8)
def quiet_transition_logs(q):
    return np.array([[log_fraction(p) for p in row] for row in weight_transition(q, 16)])


def quiet_log(z, q=FIELD_ORDER):
    """Inverse-refresh quiet bound without converting tiny probabilities."""
    checked_q(q)
    if q < 32:
        raise ValueError('The H16 interpolation map needs at least 32 field elements')
    lz = log2(z)
    values = []
    for b, row in enumerate(block_support_counts(16, 4)):
        values.append(np.logaddexp2.reduce([log2(n)+13*a*lz
                                           for a, n in enumerate(row) if n])-log2(comb(16, b)))
    return float(np.max(np.logaddexp2.reduce(quiet_transition_logs(q)+values[1:], axis=1)))


@lru_cache(32)
def local_moment_logs(z, q=FIELD_ORDER):
    """Unrounded h4 local moments; q16 is a formal algebraic replay only."""
    bits = checked_q(q)
    if not isfinite(z) or not 0 < z < 1:
        raise ValueError('Expected 0 < z < 1')
    lz, lw = log2(z), log2(1-z)
    tilt = float(np.logaddexp2(lz, lw-bits))
    quiet, excess, returns, births = [], [], [], []
    nu = principal_defects(4)
    for length in range(17):
        quiet.append(float(np.logaddexp2.reduce([
            log2(comb(length, c))+(length-c)*lz+c*lw-(c+min(4, length-c))*bits
            for c in range(length+1)])))
        # Expand t_l-k_l positively. Subtracting rounded logarithms would
        # erase this mass when z is much smaller than 1/q.
        births.append(float(np.logaddexp2.reduce([
            log2(comb(length, c))+(length-c)*lz+c*lw-c*bits
            + log1mexp2(-min(4, length-c)*bits)
            for c in range(length)])) if length else NEG)
        excess.append(float(np.logaddexp2.reduce([
            log2(comb(length, c))+(length-c)*lz+c*lw-c*bits+log1mexp2(-(4-c)*bits)
            for c in range(min(4, length+1))])))
        rho = min(4, length)
        terms = []
        for c in range(length+1):
            for e in range(min(4, 17-length)):
                rank = max(0, c-min(nu[c], max(0, 4-length+c)))
                dimension = min(rho, 4-e)
                count = ((dimension-length+max(0, length-c-4))*bits
                         + log1mexp2(-dimension*bits))
                probability = min(log1mexp2(-rho*bits), -rank*bits, count)
                terms.append(log2(comb(length, c)*comb(16-length, e))
                             +(16-c-e)*lz+(c+e)*lw+probability)
        returns.append(min(births[-1],
                           log1mexp2(-rho*bits)+max(0, 13-length)*lz,
                           float(np.logaddexp2.reduce(terms))))
    return np.array(quiet), np.array(excess), np.array(returns), np.array(births), tilt


def local_logs(z, q=FIELD_ORDER, occupancies=None, *, tighten_quiet=True):
    """Four-entry packet operator in logs, with exact-zero entries at -inf.

    q16 may test the algebra with tighten_quiet=False. It does not instantiate
    the fixed GF256 construction. Public source screens enforce the embedding.
    """
    bits = checked_q(q)
    occ = OCC if occupancies is None else np.asarray(occupancies, dtype=int)
    if occ.ndim != 2 or occ.shape[1] != 4 or np.any(occ < 0) or np.any(occ > 16):
        raise ValueError('Expected four occupancies in 0..16')
    quiet, excess, returns, births, tilt = local_moment_logs(float(z), q)
    degree = occ.sum(axis=1)
    k = quiet[occ].sum(axis=1)
    t = degree*tilt
    running, birth = np.zeros(len(occ)), np.full(len(occ), NEG)
    for block in range(4):
        positions = occ[:, block]
        birth = np.logaddexp2(birth+quiet[positions], running+births[positions])
        running += positions*tilt
    base, running, nonzero = -4*bits, np.zeros(len(occ)), np.full(len(occ), NEG)
    for block in range(4):
        e = excess[16-occ[:, block]]
        nonzero = np.logaddexp2(nonzero+base, running+e)
        running += np.logaddexp2(base, e)
    log_nonzero_fraction = log1mexp2(-16*bits)
    v = t+nonzero-log_nonzero_fraction
    running, ret = np.zeros(len(occ)), np.full(len(occ), NEG)
    for block in range(4):
        kb, jb = quiet[occ[:, block]], returns[occ[:, block]]
        ret = np.logaddexp2(ret+kb, running+jb)
        running += np.logaddexp2(kb, jb)
    log_d = 16*bits+log_nonzero_fraction
    ret = np.minimum(ret, log_d+v)
    density = log_nonzero_fraction-16*log1mexp2(-bits)
    entries = np.stack((k, density+birth, ret-log_d, density+v), axis=1)
    empty = degree == 0
    if tighten_quiet and np.any(empty):
        entries[empty, 3] = np.minimum(entries[empty, 3], quiet_log(z, q))
    if np.any(np.isnan(entries)) or np.any(np.isposinf(entries)):
        raise ArithmeticError('Nonfinite packet logarithm')
    return entries.reshape((-1, 2, 2))


def _empty(maximum):
    return np.full((4, maximum+1, maximum//2+1), NEG)


def single_coefficients(z, alpha, maximum, q, budget):
    budget.check()
    geometry = sparse_geometry(maximum)
    budget.check()
    occupancies = sorted({row[3] for row in geometry})
    logs = local_logs(z, q, occupancies).reshape((-1, 4))
    lookup = dict(zip(occupancies, logs))
    result = _empty(maximum)
    for index, (s, f, count, occ) in enumerate(geometry):
        if index % 128 == 0:
            budget.check()
        result[:, s, f] = np.logaddexp2(result[:, s, f], log2(count)-s+alpha*lookup[occ])
    return result


def good_pair_coefficients(z, alpha, maximum, q, budget):
    budget.check()
    s, f, counts, indices = good_geometry()
    keep = s+2*f <= maximum
    s, f, counts, indices = s[keep], f[keep], counts[keep], indices[keep]
    logs = local_logs(z, q, OCC[indices])
    weights = np.log2(counts)
    result = _empty(maximum)
    for i in range(len(s)):
        budget.check()
        for j in range(len(s)):
            ss, ff = int(s[i]+s[j]), int(f[i]+f[j])
            if ss+2*ff > maximum:
                continue
            composed = np.logaddexp2(logs[i, :, 0, None]+logs[j, None, 0, :],
                                     logs[i, :, 1, None]+logs[j, None, 1, :])
            value = weights[i]+weights[j]+alpha*composed.ravel()
            result[:, ss, ff] = np.logaddexp2(result[:, ss, ff], value)
    return result


def polynomial_product(left, right, maximum, budget):
    """Positive truncated matrix convolution, entirely in log2 arithmetic."""
    result = _empty(maximum)
    for s in range(min(maximum+1, left.shape[1])):
        budget.check()
        for f in range(min((maximum-s)//2+1, left.shape[2])):
            if np.all(np.isneginf(left[:, s, f])):
                continue
            ns = min(right.shape[1], maximum-s+1)
            nf = min(right.shape[2], (maximum-s-2*f)//2+1)
            if ns <= 0 or nf <= 0:
                continue
            for row in range(2):
                for col in range(2):
                    value = np.logaddexp2(left[2*row, s, f]+right[col, :ns, :nf],
                                         left[2*row+1, s, f]+right[2+col, :ns, :nf])
                    view = result[2*row+col, s:s+ns, f:f+nf]
                    np.logaddexp2(view, value, out=view)
    weight = np.arange(maximum+1)[:, None]+2*np.arange(maximum//2+1)
    result[:, weight > maximum] = NEG
    return result


def polynomial_power(base, exponent, maximum, budget):
    result = np.full((4, 1, 1), NEG)
    result[(0, 3), 0, 0] = 0.
    while exponent:
        budget.check()
        if exponent & 1:
            result = polynomial_product(result, base, maximum, budget)
        exponent >>= 1
        if exponent:
            base = polynomial_product(base, base, maximum, budget)
    return result


def regional_screen(z, alpha, maximum, *, q=FIELD_ORDER, paired=True, max_seconds=60., _budget=None):
    """Regional matrices for all coefficients s+2*f<=maximum, maximum<=64."""
    checked_q(q, construction=True)
    if (not isinstance(maximum, int) or not 1 <= maximum <= MAXIMUM
            or not 0 < z < 1 or not 0 < alpha <= 1):
        raise ValueError('Expected maximum in 1..64 and a valid witness')
    budget = _budget or Budget(max_seconds)
    single = single_coefficients(z, alpha, maximum, q, budget)
    if paired and alpha != 1:
        good, dirty = single.copy(), single.copy()
        good[:, 2:, :] = NEG
        dirty[:, :2, :] = NEG
        base = good_pair_coefficients(z, alpha, maximum, q, budget)
        for left, right in ((good, dirty), (dirty, good), (dirty, dirty)):
            base = np.logaddexp2(base, polynomial_product(left, right, maximum, budget))
        run = polynomial_power(base, 64, maximum, budget)
    else:
        run = polynomial_power(single, 128, maximum, budget)
    regions = {}
    for r in range(maximum+1):
        for f in range(r//2+1):
            s = r-2*f
            values = run[:, s, f]-log_multinomial((PAIRS-s-f, s, f))
            regions[r, f] = tuple(map(tuple, values.reshape((2, 2))))
    budget.check()
    return regions


def weight_screen(r, z, alpha, *, q=FIELD_ORDER, distance=F(7, 20), paired=True,
                  outer_measure='invertible-toeplitz', max_seconds=60., _budget=None):
    """FLOAT union over every matching profile at the selected source weight r."""
    bits = checked_q(q, construction=True)
    distance = F(distance)
    if not 0 < distance < F(1, 2):
        raise ValueError('Expected a target distance strictly between zero and one half')
    budget = _budget or Budget(max_seconds)
    forbidden = int(distance*OUTPUT_LENGTH)
    density = outer_log_density(q, outer_measure)
    regions = regional_screen(z, alpha, r, q=q, paired=paired, _budget=budget)
    prefix = alpha*(r*density-(bits+log1mexp2(-bits))-forbidden*log2(z))
    profiles = []
    for a in range(r//2+1):
        budget.check()
        for b in range(r//2+1):
            count = profile_count_upper(PAIRS, r, a, b)
            if count:
                run = power(multiply(regions[r, a], regions[r, b]), REGIONS//2)
                profiles.append((float(log2(count)+prefix+total(run[0])), a, b))
    value = float(total(row[0] for row in profiles))
    if not isfinite(value):
        raise ArithmeticError('Selected-class FLOAT result is not finite')
    budget.check()
    return dict(active_groups=r, z=z, alpha=alpha, field_order=q, field_bits=bits,
                target_distance=str(distance), forbidden_weight=forbidden,
                outer_groups=GROUPS, source_dimension=DIMENSION, output_length=OUTPUT_LENGTH,
                outer_measure=outer_measure, log2_outer_density=density,
                composition_steps=2 if paired else 1, pair_swaps=True,
                log2_upper_screen=value, worst_profile=max(profiles), profile_count=len(profiles),
                coefficient_method='positive-truncated-bivariate-log2-FLOAT',
                covers_all_supports_at_selected_weight=True, covers_all_even_odd_splits=True,
                arithmetic='FLOAT', directed=False, whole_code_certificate=False,
                scope='only the selected active-group weight; no claims for other weights')


def proposal(r, *, q=FIELD_ORDER, distance=F(7, 20), paired=True,
             outer_measure='invertible-toeplitz', max_seconds=60., max_evaluations=24,
             initial_z=None, initial_alpha=None):
    """Bounded witness proposal; returning a negative value does not certify it."""
    bits = checked_q(q, construction=True)
    if not isinstance(max_evaluations, int) or not 1 <= max_evaluations <= 1000:
        raise ValueError('max_evaluations must lie in 1..1000')
    if not isinstance(r, int) or not 1 <= r <= MAXIMUM:
        raise ValueError('Selected source weight must lie in 1..64')
    budget = Budget(max_seconds)
    old_z = 1-r/(1200 if r <= 17 else 1400)
    old_alpha = 1 if r <= 8 else .4 if r <= 17 else .2 if r <= 32 else .1
    z = old_z**(bits/8) if initial_z is None else float(initial_z)
    alpha = min(1., old_alpha*8/bits) if initial_alpha is None else float(initial_alpha)
    if not 0 < z < 1 or not 0 < alpha <= 1:
        raise ValueError('Invalid initial witness')
    best, evaluations, exhausted = None, 0, False

    def objective(theta):
        nonlocal best, evaluations
        budget.check()
        ell, a = theta
        if not 1e-8 <= ell <= 64 or not 1e-6 <= a <= 1:
            return 1e100
        row = weight_screen(r, 2.**-ell, float(a), q=q, distance=distance, paired=paired,
                            outer_measure=outer_measure, _budget=budget)
        evaluations += 1
        if best is None or row['log2_upper_screen'] < best['log2_upper_screen']:
            best = row
        return row['log2_upper_screen']

    try:
        minimize(objective, [-log2(z), alpha], method='Nelder-Mead',
                 options=dict(maxfev=max_evaluations, xatol=1e-7, fatol=1e-5))
    except TimeoutError:
        exhausted = True
    if best is None:
        raise TimeoutError('Budget ended before a complete selected-weight screen')
    best.update(proposal_evaluations=evaluations, resource_limit_reached=exhausted,
                max_seconds=max_seconds, max_evaluations=max_evaluations,
                elapsed_seconds=monotonic()-budget.started)
    return best


def source_pins():
    """Pin this screen and directly reused numeric/geometry source files."""
    directory = Path(__file__).resolve().parent
    names = ('hd_field128_lift_screen.py', 'hd_common_refresh_packets.py',
             'hd_common_refresh_principal.py', 'hd_grouped_packet_screen.py',
             'hd_grouped_quiet.py', 'hd_randomized_outer_state_profile.py',
             'hd_randomized_outer_state_screen.py', 'hd_shifted_good_pair_screen.py',
             'hd_shifted_pairs_count.py', 'hd_shifted_swap_sparse.py')
    return {name: sha256((directory/name).read_bytes()).hexdigest() for name in names}


def main():
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--field-bits', type=int, default=128)
    parser.add_argument('--r', type=int, nargs='+', default=[1, 2, 4, 8, 16])
    parser.add_argument('--distance', default='7/20')
    parser.add_argument('--z', help='Fixed mark; omit to use the field-scaled initial witness')
    parser.add_argument('--alpha', help='Fixed exponent; omit to use the field-scaled initial witness')
    parser.add_argument('--optimize', action='store_true')
    parser.add_argument('--max-seconds', type=float, default=60., help='Total wall-time budget for this invocation')
    parser.add_argument('--max-evaluations', type=int, default=24, help='Maximum optimizer evaluations per weight')
    parser.add_argument('--outer-measure', choices=OUTER_MEASURES, default='invertible-toeplitz')
    parser.add_argument('--unpaired', action='store_true')
    parser.add_argument('--output', help='New JSON file; existing files are never overwritten')
    args = parser.parse_args()
    if not 8 <= args.field_bits <= 128 or args.field_bits % 8:
        parser.error('--field-bits must be a multiple of eight in 8..128')
    if (args.z is None) != (args.alpha is None):
        parser.error('--z and --alpha must appear together')
    if any(not 1 <= r <= MAXIMUM for r in args.r):
        parser.error('--r values must lie in 1..64')
    if args.output and Path(args.output).exists():
        parser.error('--output already exists')
    budget = Budget(args.max_seconds)
    rows, failure = [], None
    q = 1 << args.field_bits
    for r in args.r:
        remaining = args.max_seconds-(monotonic()-budget.started)
        try:
            budget.check()
            if args.optimize:
                row = proposal(r, q=q, distance=F(args.distance), paired=not args.unpaired,
                               outer_measure=args.outer_measure, max_seconds=remaining,
                               max_evaluations=args.max_evaluations,
                               initial_z=None if args.z is None else float(F(args.z)),
                               initial_alpha=None if args.alpha is None else float(F(args.alpha)))
            else:
                z = ((1-r/(1200 if r <= 17 else 1400))**(args.field_bits/8)
                     if args.z is None else float(F(args.z)))
                alpha = ((1 if r <= 8 else .4 if r <= 17 else .2 if r <= 32 else .1)*8/args.field_bits
                         if args.alpha is None else float(F(args.alpha)))
                row = weight_screen(r, z, alpha, q=q, distance=F(args.distance), paired=not args.unpaired,
                                    outer_measure=args.outer_measure, _budget=budget)
            rows.append(row)
        except (TimeoutError, ArithmeticError) as exc:
            failure = str(exc)
            break
    result = dict(requested_weights=args.r, selected_weight_screens=rows,
                  completed_requested_weights=len(rows) == len(args.r), failure=failure,
                  max_seconds=args.max_seconds, elapsed_seconds=monotonic()-budget.started,
                  source_sha256=source_pins(), pin_scope='listed source files; not all transitive imports',
                  directed=False, whole_code_certificate=False,
                  scope='FLOAT screens for the listed weights only; no complete distance claim')
    encoded = json.dumps(result, indent=2, allow_nan=False)
    if args.output:
        with open(args.output, 'x', encoding='utf-8') as stream:
            stream.write(encoded+'\n')
    print(encoded)


if __name__ == '__main__':
    main()
