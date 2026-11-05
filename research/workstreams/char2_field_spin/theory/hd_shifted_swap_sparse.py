"""Truncated sparse-support transfer for the independently swapped pair route.

Only nonnegative bidegrees s,f occur, with active-group weight s+2*f.
Discarding larger weights preserves every requested coefficient exactly.
FLOAT screens are proposals, never directed or whole-code certificates.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb, log2

import numpy as np
from scipy.signal import convolve2d
from flint import arb, arb_poly, ctx

from hd_grouped_packet_screen import LIMIT, OCC, local_logs, log_multinomial
from hd_randomized_outer_state_profile import multiply, power, total
from hd_shifted_good_pair_screen import good_geometry
from hd_shifted_pairs_count import profile_count_upper
from hd_common_refresh_packets import local_moments, matrix
from hd_grouped_quiet import quiet_cap
from hd_shifted_alias_certificate import rational
from hd_shifted_coverage import rational_ceiling
from hd_shifted_cover_backend import profile_count
from hd_shifted_swap_certificate import SWAP_CONSTRUCTION
from hd_sumrank_paired_certificate import round_up, scalar_power


# Separate from the frozen per-weight 1..64 family. Each disjoint band uses
# maximum=128, alpha=1/10, and the literal rational witness shown here.
BATCH_SCHEDULE = ((65, 72, F(333, 350), F(1, 10)),
                  (73, 90, F(659, 700), F(1, 10)),
                  (91, 110, F(13, 14), F(1, 10)),
                  (111, 128, F(32, 35), F(1, 10)))

# Use maximum=256. The boundary at 225 is intentional: the last witness
# alone does not certify weight 225 at the 40-bit target.
BATCH_SCHEDULE_256 = ((129, 160, F(157, 175), F(1, 10)),
                      (161, 192, F(153, 175), F(1, 10)),
                      (193, 225, F(149, 175), F(1, 10)),
                      (226, 256, F(29, 35), F(1, 10)))


@lru_cache(4)
def sparse_geometry(maximum):
    """Exact labeled arrangements, aggregated before dividing by 2**s."""
    halves = []
    for a in range(17):
        for b in range(17-a):
            for f in range(17-a-b):
                if a+b+2*f<=maximum:
                    count = comb(16, a)*comb(16-a, b)*comb(16-a-b, f)
                    halves.append((a, b, f, count))
    counts = {}
    for i, (a, b, f, count) in enumerate(halves):
        for j in range(i, len(halves)):
            u, v, g, other = halves[j]
            s, full = a+b+u+v, f+g
            if s+2*full>maximum:
                continue
            occupancy = tuple(sorted((a+f, b+f, u+g, v+g)))
            key = s, full, occupancy
            counts[key] = counts.get(key, 0)+count*other*(1 if i==j else 2)
    return tuple((s, f, n, occ) for (s, f, occ), n in counts.items())


@lru_cache(4)
def packet_logs_screen(z):
    return local_logs(z).reshape(-1, 4)


_OCC_INDEX = {tuple(row): i for i, row in enumerate(OCC)}


@lru_cache(4)
def swap_single_screen(z, alpha, maximum):
    """Unmarked swap average, after each packet's entrywise alpha power."""
    logs = packet_logs_screen(z)
    result = np.zeros((4, maximum+1, maximum//2+1))
    for s, f, count, occ in sparse_geometry(maximum):
        result[:, s, f] += np.exp2(log2(count)-s+alpha*logs[_OCC_INDEX[occ]])
    return result


def good_pair_screen(z, alpha, maximum):
    """Retain only good-good geometries of target-relevant weight."""
    s, f, counts, indices = good_geometry()
    keep = s+2*f<=maximum
    s, f, counts, indices = s[keep], f[keep], counts[keep], indices[keep]
    logs = packet_logs_screen(z)[indices].reshape(-1, 2, 2)
    composed = np.logaddexp2(logs[:, None, :, 0, None]+logs[None, :, None, 0, :],
                            logs[:, None, :, 1, None]+logs[None, :, None, 1, :])
    count_logs = np.log2(counts)
    weighted = np.exp2(count_logs[:, None, None, None]+count_logs[None, :, None, None]+alpha*composed)
    result = np.zeros((4, maximum+1, maximum//2+1))
    for i in range(len(s)):
        for j in range(len(s)):
            ss, ff = s[i]+s[j], f[i]+f[j]
            if ss+2*ff<=maximum:
                result[:, ss, ff] += weighted[i, j].ravel()
    return result


def truncate_screen(coefficients, maximum):
    """Keep the down-set s+2*f<=maximum in a four-entry matrix array."""
    result = np.zeros((4, maximum+1, maximum//2+1))
    ns = min(result.shape[1], coefficients.shape[1])
    nf = min(result.shape[2], coefficients.shape[2])
    result[:, :ns, :nf] = coefficients[:, :ns, :nf]
    weight = np.arange(maximum+1)[:, None]+2*np.arange(maximum//2+1)
    result[:, weight>maximum] = 0.
    return result


def product_screen(left, right, maximum):
    """Ordinary ordered matrix product with positive direct convolutions."""
    shape = (left.shape[1]+right.shape[1]-1, left.shape[2]+right.shape[2]-1)
    result = np.empty((4,)+shape)
    for row in range(2):
        for col in range(2):
            result[2*row+col] = (convolve2d(left[2*row], right[col])
                                + convolve2d(left[2*row+1], right[2+col]))
    return truncate_screen(result, maximum)


def power_screen(base, exponent, maximum):
    """Coefficient-exact truncation; arithmetic here remains FLOAT."""
    result = np.zeros((4, 1, 1))
    result[(0, 3), 0, 0] = 1.
    while exponent:
        if exponent & 1:
            result = product_screen(result, base, maximum)
        exponent >>= 1
        if exponent:
            base = product_screen(base, base, maximum)
    return result


def regional_screen(z, alpha, maximum, *, paired=True):
    """Return every regional log2 matrix with s+2*f<=maximum.

    Equal positive markers reduce coefficient magnitudes. They cancel
    explicitly through each target normalizer and are not approximations.
    """
    if not isinstance(maximum, int) or not 1<=maximum<=64 or not 0<z<1 or not 0<alpha<=1:
        raise ValueError('Expected 1<=maximum<=64 and a valid witness')
    marker = maximum/4096
    degree = np.arange(maximum+1)[:, None]+np.arange(maximum//2+1)
    single = swap_single_screen(z, alpha, maximum).copy()
    single *= marker**degree
    if paired:
        good, dirty = single.copy(), single.copy()
        good[:, 2:, :] = 0.
        dirty[:, :2, :] = 0.
        composed = good_pair_screen(z, alpha, maximum)
        composed *= marker**degree
        base = composed+product_screen(good, dirty, maximum)+product_screen(dirty, good, maximum)+product_screen(dirty, dirty, maximum)
        run = power_screen(base, 64, maximum)
    else:
        run = power_screen(single, 128, maximum)
    regions = {}
    for r in range(maximum+1):
        for f in range(r//2+1):
            s = r-2*f
            values = run[:, s, f]
            if np.any(values<0) or not np.all(np.isfinite(values)):
                raise ArithmeticError('Invalid positive sparse coefficient')
            with np.errstate(divide='ignore'):
                logs = np.log2(values).reshape(2, 2)
            logs -= log_multinomial((4096-s-f, s, f))+(s+f)*log2(marker)
            regions[r, f] = tuple(map(tuple, logs))
    return regions


def weight_screen(r, z, alpha, *, paired=True):
    """Union all matching profiles at this r, including all lane balances."""
    regions = regional_screen(z, alpha, r, paired=paired)
    prefix = alpha*(r*(2048-128*log2(255))-log2(255)-LIMIT*log2(z))
    profiles = []
    for a in range(r//2+1):
        for b in range(r//2+1):
            count = profile_count_upper(4096, r, a, b)
            if not count:
                continue
            run = power(multiply(regions[r, a], regions[r, b]), 128)
            value = log2(count)+prefix+total(run[0])
            profiles.append((value, a, b))
    return dict(active_groups=r, z=z, alpha=alpha, composition_steps=2 if paired else 1,
                log2_upper_screen=total(row[0] for row in profiles),
                worst_profile=max(profiles), profile_count=len(profiles),
                coefficient_method='positive-truncated-bivariate',
                directed=False, whole_code_certificate=False)


@lru_cache(4)
def sparse_packet_caps(z, maximum, bits=60):
    """The existing exact packet cap, restricted to relevant occupancies."""
    z = F(z)
    local = local_moments(z, q=256, h=4)
    density, quiet = F(256**16-1, 255**16), quiet_cap(z, q=256, h=4)
    result = {}
    for _, _, _, occ in sparse_geometry(maximum):
        if occ in result:
            continue
        k, birth, ret, live = matrix(occ, z, local, q=256, h=4)
        entries = (k, density*birth, ret,
                   min(density*live, quiet) if not sum(occ) else density*live)
        result[occ] = tuple(round_up(x, bits) for x in entries)
    return result


def valid_degrees(maximum):
    return tuple((s, f) for f in range(maximum//2+1) for s in range(maximum-2*f+1))


def encode_arb(coefficients, maximum, marker=F(1)):
    """Carry-free encoding s+(2*maximum+1)*f with positive rational tilt."""
    stride = 2*maximum+1
    degree = stride*(maximum//2)+(maximum%2)
    result = [[arb(0) for _ in range(degree+1)] for _ in range(4)]
    marker = rational(marker)
    for (s, f), values in coefficients.items():
        if s+2*f>maximum:
            continue
        tilt = marker**(s+f)
        for entry in range(4):
            result[entry][s+stride*f] = values[entry]*tilt
    return tuple(arb_poly(row) for row in result)


def product_arb(left, right, maximum):
    """Positive coefficient majorant of a truncated ordered matrix product.

    Both inputs must already have weight s+2*f<=maximum. The radix then
    exceeds every singleton degree in one multiplication, so
    discarded singleton overflow cannot masquerade as a full-pair degree.
    Arb encloses each product; upper endpoints preserve coefficient order.
    """
    stride = 2*maximum+1
    degree = stride*(maximum//2)+(maximum%2)
    indices = tuple(s+stride*f for s, f in valid_degrees(maximum))
    result = []
    for row in range(2):
        for col in range(2):
            value = (left[2*row]*right[col]+left[2*row+1]*right[2+col]).truncate(degree+1)
            kept = [arb(0) for _ in range(degree+1)]
            for index in indices:
                upper = value[index].upper()
                if upper<0 or not upper.is_finite():
                    raise ArithmeticError('Invalid directed positive coefficient')
                kept[index] = upper
            result.append(arb_poly(kept))
    return tuple(result)


def power_arb(base, exponent, maximum):
    result = tuple(arb_poly([int(i in (0, 3))]) for i in range(4))
    while exponent:
        if exponent & 1:
            result = product_arb(result, base, maximum)
        exponent >>= 1
        if exponent:
            base = product_arb(base, base, maximum)
    return result


def directed_parts(z, alpha, maximum, bits=60):
    """Swap-average single steps and good-good products before alpha."""
    table = sparse_packet_caps(F(z), maximum, bits)
    exponent = rational(alpha)
    powered = {occ: tuple(rational(x)**exponent if x else arb(0) for x in values)
               for occ, values in table.items()}
    single = {}
    for s, f, n, occ in sparse_geometry(maximum):
        values = single.setdefault((s, f), [arb(0) for _ in range(4)])
        coefficient = rational(F(n, 2**s))
        for entry in range(4):
            values[entry] += coefficient*powered[occ][entry]
    if F(alpha)==1:
        return single, None
    s, f, counts, indices = good_geometry()
    groups = [(int(ss), int(ff), int(n), tuple(rational(x) for x in table[tuple(OCC[index])]))
              for ss, ff, n, index in zip(s, f, counts, indices) if ss+2*ff<=maximum]
    composed = {}
    for s, f, n, left in groups:
        for t, g, m, right in groups:
            if s+t+2*(f+g)>maximum:
                continue
            values = composed.setdefault((s+t, f+g), [arb(0) for _ in range(4)])
            for row in range(2):
                for col in range(2):
                    value = left[2*row]*right[col]+left[2*row+1]*right[2+col]
                    if not value.is_zero():
                        values[2*row+col] += n*m*value**exponent
    return single, composed


def regional_directed(z, alpha, maximum, bits=60):
    """Every regional coefficient through maximum, at caller precision.

    A larger maximum changes the positive marker and retained down-set,
    not a requested exact coefficient after normalization. The resource
    guard is deliberately bounded; it is not a mathematical cutoff.
    """
    if not isinstance(maximum, int) or not 1<=maximum<=256:
        raise ValueError('Directed sparse maximum must lie in 1..256')
    marker = F(maximum, 4096)
    single, composed = directed_parts(z, alpha, maximum, bits)
    if F(alpha)==1:
        run = power_arb(encode_arb(single, maximum, marker), 128, maximum)
    else:
        good = encode_arb({key: value for key, value in single.items() if key[0]<=1}, maximum, marker)
        dirty = encode_arb({key: value for key, value in single.items() if key[0]>=2}, maximum, marker)
        parts = (encode_arb(composed, maximum, marker), product_arb(good, dirty, maximum),
                 product_arb(dirty, good, maximum), product_arb(dirty, dirty, maximum))
        base = tuple(sum((part[i] for part in parts), arb_poly()) for i in range(4))
        run = power_arb(base, 64, maximum)
    regions = {}
    stride = 2*maximum+1
    for s, f in valid_degrees(maximum):
        count = comb(4096, s)*comb(4096-s, f)
        normalizer = count*rational(marker)**(s+f)
        regions[s+2*f, f] = tuple(run[i][s+stride*f].upper()/normalizer for i in range(4))
    return regions


def _weight_from_regions(r, z, alpha, regions, precision, bits):
    """All-profile union; caller retains the regional precision context."""
    union, profiles = arb(0), 0
    for a in range(r//2+1):
        for b in range(r//2+1):
            count = profile_count_upper(4096, r, a, b)
            if not count:
                continue
            left, right = regions[r, a], regions[r, b]
            pair = tuple(left[2*i]*right[j]+left[2*i+1]*right[2+j] for i in range(2) for j in range(2))
            run = scalar_power(pair, 128)
            union += count*(run[0]+run[1])
            profiles += 1
    log = lambda x: x.log()/arb(2).log()
    prefix = rational(alpha)*(r*(2048-128*log(arb(255)))-log(arb(255))-LIMIT*log(rational(z)))
    value = prefix+log(union)
    ceiling = rational_ceiling(value)
    enclosure = str(value)
    ranges = ((r, r), (0, 4096), (0, 4096))
    return dict(active_groups=r, profiles=profiles, ranges=ranges,
                integer_profile_triples=profile_count(ranges),
                covers_all_supports_at_weight=True, covers_all_even_odd_splits=True,
                log2_upper_rational=str(ceiling), log2_upper_enclosure=enclosure,
                witness=dict(z=str(z), alpha=str(alpha), precision=precision, packet_mantissa_bits=bits),
                coefficient_method='carry-free-Arb-truncated-bivariate',
                construction=SWAP_CONSTRUCTION,
                pair_swaps=True, field_order=256, outer_groups=8192, source_dimension=1048576,
                output_length=2097152, forbidden_weight=LIMIT,
                directed=True, whole_code_certificate=False)


def weight_receipt(r, z, alpha, *, precision=256, bits=60):
    """Frozen 1..64 API: every nonempty support with exactly r groups."""
    z, alpha = F(z), F(alpha)
    if not isinstance(r, int) or not 1<=r<=64 or not 0<z<1 or not 0<alpha<=1 or precision<128 or bits<16:
        raise ValueError('Invalid sparse scope or rational witness')
    with ctx.workprec(precision):
        regions = regional_directed(z, alpha, r, bits)
        return _weight_from_regions(r, z, alpha, regions, precision, bits)


def batch_receipt(first, last, z, alpha, *, maximum=128, precision=256, bits=60, progress=None):
    """Directed all-support union on a contiguous, nonempty weight band.

    Compute regional_directed exactly once with a common rational witness.
    No interpolation in r is used: each weight retains every matching
    profile and its own outer/projective factor. The shared marker cancels
    separately in each regional coefficient. Only first..last is claimed,
    even if maximum is larger. This API does not extend family_witness.
    """
    z, alpha = F(z), F(alpha)
    if (not all(isinstance(x, int) for x in (first, last, maximum))
            or not 1<=first<=last<=maximum<=256 or not 0<z<1 or not 0<alpha<=1
            or precision<128 or bits<16):
        raise ValueError('Invalid sparse batch scope or rational witness')
    rows = []
    with ctx.workprec(precision):
        regions = regional_directed(z, alpha, maximum, bits)
        for r in range(first, last+1):
            row = _weight_from_regions(r, z, alpha, regions, precision, bits)
            row['witness'].update(coefficient_maximum=maximum, coefficient_marker=str(F(maximum, 4096)))
            rows.append(row)
            if progress:
                progress(row)
        union = sum((arb(2)**rational(F(row['log2_upper_rational'])) for row in rows), arb(0))
        ceiling = rational_ceiling(union.log()/arb(2).log())
    ranges = ((first, last), (0, 4096), (0, 4096))
    return dict(active_group_range=[first, last], ranges=ranges,
                integer_profile_triples=profile_count(ranges), covers_all_even_odd_splits=True,
                covers_all_supports_in_range=True, log2_upper_rational=str(ceiling),
                weight_receipts=rows,
                witness=dict(z=str(z), alpha=str(alpha), precision=precision, packet_mantissa_bits=bits,
                             coefficient_maximum=maximum, coefficient_marker=str(F(maximum, 4096))),
                coefficient_method='carry-free-Arb-truncated-bivariate',
                directed=True, whole_code_certificate=False, construction=SWAP_CONSTRUCTION,
                pair_swaps=True, field_order=256, outer_groups=8192, source_dimension=1048576,
                output_length=2097152, forbidden_weight=LIMIT,
                scope='all indexed nonempty support sets in the stated weight range')


def family_witness(r):
    """Fixed rational proposal; selecting it does not assert acceptance."""
    if not isinstance(r, int) or not 1<=r<=64:
        raise ValueError('Sparse witness family is scoped to 1..64')
    if r<=17:
        return F(1)-F(r, 1200), F(1) if r<=8 else F(2, 5)
    return F(1)-F(r, 1400), F(1, 5) if r<=32 else F(1, 10)


def family_receipt(last=64, *, precision=256, bits=60, progress=None):
    """Directed disjoint union of all support weights 1 through last."""
    if not isinstance(last, int) or not 1<=last<=64:
        raise ValueError('Sparse family is scoped to 1..64')
    rows = []
    for r in range(1, last+1):
        z, alpha = family_witness(r)
        row = weight_receipt(r, z, alpha, precision=precision, bits=bits)
        rows.append(row)
        if progress:
            progress(row)
    with ctx.workprec(precision):
        union = sum((arb(2)**rational(F(row['log2_upper_rational'])) for row in rows), arb(0))
        ceiling = rational_ceiling(union.log()/arb(2).log())
    ranges = ((1, last), (0, 4096), (0, 4096))
    return dict(active_group_range=[1, last], ranges=ranges,
                integer_profile_triples=profile_count(ranges), covers_all_even_odd_splits=True,
                log2_upper_rational=str(ceiling),
                weight_receipts=rows, directed=True, whole_code_certificate=False,
                construction=SWAP_CONSTRUCTION,
                pair_swaps=True, field_order=256, forbidden_weight=LIMIT,
                scope='all indexed nonempty support sets in the stated weight range')


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r', type=int, default=1)
    parser.add_argument('--z', default='0.999', help='Exact rational/decimal witness for directed mode')
    parser.add_argument('--alpha', default='0.2', help='Exact rational/decimal witness for directed mode')
    parser.add_argument('--unpaired', action='store_true')
    parser.add_argument('--directed', action='store_true')
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--family', type=int, help='Replay every weight 1..N, N<=64')
    parser.add_argument('--batch', type=int, nargs=2, metavar=('FIRST', 'LAST'),
                        help='Replay all weights in FIRST..LAST with one witness, LAST<=256')
    parser.add_argument('--maximum', type=int, default=128, help='Batch coefficient maximum, at most 256')
    args = parser.parse_args()
    if args.family is not None and args.batch is not None:
        parser.error('--family and --batch are mutually exclusive')
    if args.batch is not None:
        result = batch_receipt(*args.batch, F(args.z), F(args.alpha), maximum=args.maximum,
                               precision=args.precision,
                               progress=lambda row: print(json.dumps(row), flush=True))
        result.pop('weight_receipts')
        print(json.dumps(result))
    elif args.family is not None:
        result = family_receipt(args.family, precision=args.precision,
                                progress=lambda row: print(json.dumps(row), flush=True))
        result.pop('weight_receipts')
        print(json.dumps(result))
    elif args.directed:
        print(json.dumps(weight_receipt(args.r, F(args.z), F(args.alpha), precision=args.precision)))
    else:
        print(json.dumps(weight_screen(args.r, float(F(args.z)), float(F(args.alpha)), paired=not args.unpaired)))
