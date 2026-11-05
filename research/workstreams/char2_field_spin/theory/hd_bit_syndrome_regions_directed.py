"""Directed original-one-state NZ regional caps for the H* selected gate.

The production interface rebuilds the exact-rational mark993/1000 at
256 or384 Arb bits. It retains every histogram s+2f<=18, f<=1:36 matrices.
Packet amplitudes are independent uniform nonzero values conditional on
the complete outer mask. Pair permutations and fair swaps are averaged
once; one4096-pair multinomial normalizes the chronological128-packet
power. No old outer receipt, FLOAT cap, or post-refresh bank is imported.

Profile formulas reproduce the reviewed hd_xor_pair_route cap algebra,
without calling or changing its maximum17 interface. The new coefficient
cutoff is a positive down-set. The radix exceeds every singleton degree
sum before truncation, so discarded terms cannot carry into kept terms.
These regional bounds are not a distance certificate or a whole-code claim.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
from math import comb
from pathlib import Path

from flint import arb, arb_poly, ctx

from hd_bit_syndrome_outer import CONSTRUCTION
from hd_common_refresh_screen import principal_defects
from hd_grouped_quiet import quiet_cap
from hd_shifted_swap_sparse import sparse_geometry
from hd_xor_group_geometry import matching_max_counts
from hd_xor_multigroup import nonzero_kernel_probability

DIRECTORY = Path(__file__).resolve().parent
MARK, MAXIMUM, FULL_MAX = F(993, 1000), 18, 1
PACKETS, PAIRS_PER_PACKET, PAIRS, SOURCE_SIZE = 128, 32, 4096, 8192
SOURCE_PINS = {
    'hd_xor_pair_route.py': 'dec1f2e2ba064f60529c31675a3779ef53dc6021035bcd024366e5614b3861bc',
    'hd_xor_multigroup.py': 'e64a98c2db2a04d34ee61a7814825a85fba6724d281e0963215c6e59442d1bd0',
    'hd_xor_column_probe.py': '06b3ecbd39c2aa917851745fa4f197f5c5641ea266d6700362a39fb8ead1b563',
    'hd_bucket_feedback.py': '66a1e97608046ae9e1fa6e37f6b6e5634de586fd03ae819278014d6a8ae16d38',
    'hd_common_refresh_screen.py': 'e1af79f18c7634dd3c2489185bb6c95d0b54a207d3a360fcfeff17b94093922f',
    'hd_common_refresh_principal.py': '1ecd8f4243b81ca095f045d6c6b9c28529900bcfff1c5fa4382a2cd8ca259cb5',
    'hd_grouped_quiet.py': '4d84222bcaaeab1d134bb3ab8235d4082f39ebb9970e47efa110f52f9bf5a369',
    'hd_shifted_swap_sparse.py': '9b8875ec1433a2bc79dbd1e42d3b90815528a70ce0234069cfd03d61b030c58c',
    'hd_xor_group_geometry.py': '880b51f316341ed4d77d5020694627805cfd708cdf6278a95e142b8ddf5aff17',
    'hd_shifted_general_count.py': '74db770b5918a3dfbfce167fad2f136c81d4cfe7665bc411261d51e0bff56209',
    'hd_shifted_pairs_count.py': 'bb1700278361343325007cf558751af2db6df19dbf6c2090502941f1a95397be',
    'hd_systematic_mds.py': 'fe04df5613156b79215784e960e8a85e3cf6d4e17f8dc0c268c0eb911e66c5e2',
    'hd_sumrank_maps.py': '704713311f17fe327aeeece946326ed53e2487e752f86fe96aeb82e135d2577b',
    'hd_blockwise_shared_transfer.py': '0245af670491b5aaab89cb6d9e620a2ebc490907eef4dea8dc1ee804085a260b',
    'hd_blockwise_maps.py': '1f80f1f8c6f89509c69043b4542ba9598d56d6e1afe2805159c2c50848becab4',
    'hd_bit_syndrome_outer.py': 'a2c2d0fee734fad2c4b74da7850589f8e84180798dcf5f79beab30e06939292d',
}


def checked_dependencies():
    for name, expected in SOURCE_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Reviewed dependency changed: '+name)


def rational(value):
    value = F(value)
    return arb(value.numerator)/value.denominator


def _mark(value):
    if not isinstance(value, (int, F)) or isinstance(value, bool) or not 0 < value <= 1:
        raise ValueError('Exact rational mark in (0,1] required; FLOAT marks are not accepted')
    return F(value)


def _scope(maximum, full_max):
    if (type(maximum) is not int or not 0 <= maximum <= MAXIMUM
            or type(full_max) is not int or not 0 <= full_max <= maximum//2):
        raise ValueError('Explicit maximum0..18 and downward full-pair cutoff required')
    return tuple((s, f) for f in range(full_max+1) for s in range(maximum-2*f+1))


def _upper(value):
    upper = value.upper()
    if not upper.is_finite() or upper < 0:
        raise ArithmeticError('A finite nonnegative directed upper endpoint is required')
    return upper


def _upper_min(*values):
    return min(_upper(value) for value in values)


def packet_caps_arb(mark, maximum=MAXIMUM, *, selected_profiles=None, check_time=lambda: None):
    """Profile-wise original NZ caps at the caller's current Arb precision.

    selected_profiles is a small exact-test interface. Production uses all
    sorted four-block occupancy profiles of total degree<=maximum.
    All choices between overlapping balls use explicit upper endpoints.
    """
    mark = _mark(mark)
    _scope(maximum, 0)
    check_time()
    profiles = (tuple(profile for profile in combinations_with_replacement(range(17), 4)
                      if sum(profile) <= maximum) if selected_profiles is None else tuple(selected_profiles))
    if (len(set(profiles)) != len(profiles) or any(len(p) != 4 or tuple(sorted(p)) != p
            or any(type(x) is not int or not 0 <= x <= 16 for x in p)
            or sum(p) > maximum for p in profiles)):
        raise ValueError('Distinct sorted four-block profiles within the declared maximum required')
    q, h, length, state = 256, 4, 16, 16
    zero, one, zz = arb(0), arb(1), rational(mark)
    kernel = [rational(nonzero_kernel_probability(j)) for j in range(17)]
    quiet = [zz**j*kernel[j] for j in range(17)]
    free = [(zz+(one-zz)/q)**j for j in range(17)]
    uniform_quiet = [sum((comb(j, c)*zz**(j-c)*(one-zz)**c/q**(c+min(h, j-c))
                         for c in range(j+1)), zero) for j in range(17)]
    condition = [(arb(q)/(q-1))**j for j in range(17)]
    ret = [zero]
    defects = principal_defects(h)  # Exact GF256 ranks, despite the legacy module name.
    for j in range(1, 17):
        rank_cap = zero
        for c in range(j+1):
            for e in range(min(h, length-j+1)):
                first = arb(q)**(-max(0, c-min(defects[c], max(0, h-j+c))))
                second = arb(q)**(max(0, j-c-h)-j)*(q**min(min(h, j), h-e)-1)
                selected = _upper_min(one-arb(q)**(-min(h, j)), first, second)
                rank_cap += comb(j, c)*comb(length-j, e)*zz**(length-c-e)*(one-zz)**(c+e)*selected
        ret.append(_upper_min(condition[j]*(free[j]-uniform_quiet[j]),
                              (one-kernel[j])*zz**max(1, length-h+1-j), condition[j]*rank_cap))
    base = one/q**h
    excess = [sum((comb(n, c)*zz**(n-c)*(one-zz)**c*(one/q**c-base)
                   for c in range(min(h, n+1))), zero) for n in range(17)]
    dim = q**state-1
    refresh_density = arb(dim)/(q-1)**state
    inverse = rational(quiet_cap(mark, q=q, state=state, width=64, h=h))
    result = {}
    for profile in profiles:
        check_time()
        degree = sum(profile)
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
        packet_ret = _upper_min(packet_ret, dim*live)
        actual_live = refresh_density*live
        if degree == 0:
            actual_live = _upper_min(actual_live, inverse)
        result[profile] = tuple(_upper(value) for value in (
            zz**degree*feedback_quiet, refresh_density*zz**degree*(one-feedback_quiet),
            packet_ret/dim, actual_live))
    return result


def packet_from_table(table, maximum, full_max, *, geometry=None):
    """Unmarked B[s,f]=sum n(s,f,o)*2^-s M_o; multiplicities stay exact."""
    valid = set(_scope(maximum, full_max))
    rows = sparse_geometry(maximum) if geometry is None else geometry
    result = {key: [arb(0) for _ in range(4)] for key in valid}
    for s, f, count, profile in rows:
        if (any(type(x) is not int or x < 0 for x in (s, f, count))
                or s+2*f > maximum):
            raise ValueError('Nonnegative exact in-scope geometry rows required')
        if f > full_max:
            continue
        values = table[profile]
        if len(values) != 4:
            raise ValueError('Four ordered packet entries required')
        weight = rational(F(count, 2**s))
        for entry in range(4):
            result[s, f][entry] += weight*_upper(arb(values[entry]))
    return {key: tuple(_upper(value) for value in values) for key, values in result.items()}


def encode_arb(coefficients, maximum, full_max, marker=F(1)):
    """Tilt and carry-free encode only the declared positive down-set."""
    degrees = _scope(maximum, full_max)
    marker = _mark(marker)
    stride = 2*maximum+1
    last = maximum+(stride-2)*full_max
    result = [[arb(0) for _ in range(last+1)] for _ in range(4)]
    zz = rational(marker)
    valid = set(degrees)
    for (s, f), values in coefficients.items():
        if (s, f) not in valid or len(values) != 4:
            raise ValueError('Packet coefficient exceeds the declared down-set')
        weight = zz**(s+f)
        for entry in range(4):
            result[entry][s+stride*f] = _upper(arb(values[entry])*weight)
    return tuple(arb_poly(row) for row in result)


def _checked_polys(values, maximum, full_max):
    degrees = _scope(maximum, full_max)
    stride = 2*maximum+1
    valid = {s+stride*f for s, f in degrees}
    if len(values) != 4 or any(not isinstance(value, arb_poly) for value in values):
        raise ValueError('Four Arb polynomials required')
    for polynomial in values:
        for index in range(len(polynomial)):
            coefficient = polynomial[index]
            if (not coefficient.is_finite() or not coefficient >= 0
                    or (index not in valid and not coefficient.is_zero())):
                raise ValueError('Each input must already be a nonnegative truncated polynomial')
    return degrees


def product_arb(left, right, maximum, full_max):
    """Ordered positive product; no singleton carry, no discarded feedback."""
    degrees = _checked_polys(left, maximum, full_max)
    _checked_polys(right, maximum, full_max)
    stride = 2*maximum+1
    last = maximum+(stride-2)*full_max
    indices = tuple(s+stride*f for s, f in degrees)
    result = []
    for row in range(2):
        for column in range(2):
            value = (left[2*row]*right[column]+left[2*row+1]*right[2+column]).truncate(last+1)
            kept = [arb(0) for _ in range(last+1)]
            for index in indices:
                kept[index] = _upper(value[index])
            result.append(arb_poly(kept))
    return tuple(result)


def power_arb(base, exponent, maximum, full_max, *, check_time=lambda: None):
    _checked_polys(base, maximum, full_max)
    if type(exponent) is not int or exponent < 0:
        raise ValueError('Nonnegative exact packet exponent required')
    result = tuple(arb_poly([int(i in (0, 3))]) for i in range(4))
    while exponent:
        check_time()
        if exponent & 1:
            result = product_arb(result, base, maximum, full_max)
        exponent >>= 1
        if exponent:
            base = product_arb(base, base, maximum, full_max)
    return result


def normalize_arb(powered, maximum, full_max, *, pairs, marker=F(1)):
    """One GLOBAL pair normalizer; no orientation or half normalization here."""
    degrees = _checked_polys(powered, maximum, full_max)
    marker = _mark(marker)
    if type(pairs) is not int or pairs < maximum:
        raise ValueError('Enough global pair slots for every requested histogram required')
    stride, zz = 2*maximum+1, rational(marker)
    regions, normalizers = {}, {}
    for s, f in degrees:
        count = comb(pairs, s)*comb(pairs-s, f)
        normalizers[s+2*f, f] = count
        divisor = count*zz**(s+f)
        regions[s+2*f, f] = tuple(_upper(entry[s+stride*f]/divisor) for entry in powered)
    return regions, normalizers


def source_counts():
    """Exact indexed r18 support partition; no saved FLOAT artifact is read."""
    counts = matching_max_counts(SOURCE_SIZE, MAXIMUM, max_reduced_weight=MAXIMUM)
    total = comb(SOURCE_SIZE, MAXIMUM)
    n0_numerator = SOURCE_SIZE*comb(SOURCE_SIZE-MAXIMUM, MAXIMUM)
    n0, remainder = divmod(n0_numerator, SOURCE_SIZE-MAXIMUM)
    if remainder or counts[0] != n0 or sum(counts.values()) != total:
        raise ArithmeticError('Exact matching partition / independent-set identity failed')
    return dict(F0=counts[0], F1=counts[1], F_le1=counts[0]+counts[1],
                F_ge2=total-counts[0]-counts[1], all_sources=total,
                all_strata=counts)


def regional_bank(precision=256, *, mark=MARK, check_time=lambda: None):
    """Fresh selected bank; caller owns serialization and later contraction.

    Returned entries are flat row-major Arb upper endpoints, keyed(m,f).
    normalizers contains exact integer multinomials (without the cancelling
    dyadic coefficient marker). The source count is NOT multiplied into R.
    """
    if type(precision) is not int or precision not in (256, 384) or _mark(mark) != MARK:
        raise ValueError('This frozen adapter uses mark993/1000 and precision256/384 only')
    checked_dependencies()
    check_time()
    marker = F(MAXIMUM, PAIRS)
    counts = source_counts()
    with ctx.workprec(precision):
        table = packet_caps_arb(mark, MAXIMUM, check_time=check_time)
        check_time()
        coefficients = packet_from_table(table, MAXIMUM, FULL_MAX)
        base = encode_arb(coefficients, MAXIMUM, FULL_MAX, marker)
        powered = power_arb(base, PACKETS, MAXIMUM, FULL_MAX, check_time=check_time)
        regions, normalizers = normalize_arb(powered, MAXIMUM, FULL_MAX, pairs=PAIRS, marker=marker)
        check_time()
    if len(regions) != 36:
        raise AssertionError('Incomplete selected regional histogram bank')
    return dict(regions=regions, normalizers=normalizers, source_counts=counts,
        metadata=dict(outer_construction=CONSTRUCTION,
            intended_scope='homogeneous arity8/r18/Fmax<=1 scalar class; regional caps only',
            inner='original-one-state h4; four RS[16,4] blocks; state16; packet64',
            recurrence="Y=X+A Q; Qprime=M(Q+C X)", state_order=['zero', 'normalized-uniform-nonzero'],
            matrix_entry_order=['00', '01', '10', '11'],
            field_order=256, mark=str(MARK), regional_alpha=1, precision_bits=precision,
            amplitudes='independent-uniform-nonzero conditional on complete outer mask',
            coefficient_maximum=MAXIMUM, coefficient_full_max=FULL_MAX,
            coefficient_marker=str(marker), carry_free_stride=2*MAXIMUM+1,
            coefficient_method='positive-Arb-ordered-truncated-bivariate',
            packet_cap_method='original-NZ-profile-formulas-with-upper-endpoint-minima',
            total_packets=PACKETS, pairs_per_packet=PAIRS_PER_PACKET, regional_pairs=PAIRS,
            regional_histograms=36, regional_length=SOURCE_SIZE,
            independent_regional_pair_permutations=True, independent_fair_swaps=True,
            alternating_matchings=True, singleton_orientation_denominator='2^s once inside packet',
            global_normalization_count=1, source_geometry_resampled=False,
            source_count_in_regional_entries=False, source_counts_recomputed_exactly=True,
            post_refresh_bank_used=False, uniform_amplitude_provider_used=False,
            maximum17_guard_modified=False, saved_FLOAT_caps_used=False,
            source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
            reviewed_dependency_sha256=dict(SOURCE_PINS), directed=True,
            regional_upper_bounds=True, distance_certificate=False, whole_code_certificate=False))
