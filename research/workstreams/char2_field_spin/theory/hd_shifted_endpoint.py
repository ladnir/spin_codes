"""Fixed-lane comparator obstruction and a fresh-swap endpoint certificate.

The receipt covers only the two alternating supports. Fresh independent
within-pair swaps define a DIFFERENT ensemble from the fixed-lane route.
Fractional powers precede the average over the new swap randomness.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb

from flint import arb, ctx

from hd_common_refresh_packets import local_moments, matrix
from hd_randomized_outer_proof import prescribed_mds_count
from hd_shifted_alias_certificate import LIMIT, rational
from hd_sumrank_paired_certificate import round_up, scalar_power


LENGTH = 1 << 21
STEPS = LENGTH // 64
ACTIVE_GROUPS = 4096
DEFAULT_Z = F(1, 32)
DEFAULT_ALPHA = F(1, 4)


def _exact(value):
    if isinstance(value, float):
        raise TypeError('Use a rational witness, not a floating-point proposal')
    return F(value)


def _ceil(value):
    return F(int((value.upper()*10**12).upper().ceil().fmpz()), 10**12)


def _floor(value):
    return F(int((value.lower()*10**12).lower().floor().fmpz()), 10**12)


@lru_cache(None)
def swap_occupancy_counts(half=16):
    """Exact multiplicities of sorted lane occupancies of 2*half fair swaps."""
    if not isinstance(half, int) or half < 1:
        raise ValueError('A positive integer half-packet size is required')
    counts = {}
    for a in range(half+1):
        for b in range(half+1):
            occupancy = tuple(sorted((a, half-a, b, half-b)))
            counts[occupancy] = counts.get(occupancy, 0)+comb(half,a)*comb(half,b)
    assert sum(counts.values()) == 2**(2*half)
    return tuple(sorted(counts.items()))


@lru_cache(4)
def endpoint_packet_caps(z=DEFAULT_Z, bits=60):
    """Exact rational h4 caps, restricted to occupancies used by fair swaps."""
    z = _exact(z)
    if not 0 < z < 1 or not isinstance(bits,int) or bits < 16:
        raise ValueError('Invalid packet mark or dyadic mantissa width')
    local = local_moments(z, q=256, h=4)
    density = F(256**16-1, 255**16)
    result = []
    for occupancy, count in swap_occupancy_counts():
        quiet, birth, ret, live = matrix(occupancy,z,local,q=256,h=4)
        # Every occupancy has weight32, so the empty-packet quiet cap is unused.
        entries = quiet, density*birth, ret, density*live
        result.append((occupancy,count,tuple(round_up(x,bits) for x in entries)))
    return tuple(result)


def fractional_average(weighted_caps, alpha):
    """Average entrywise cap**alpha; weights are exact nonnegative rationals."""
    alpha = _exact(alpha)
    if not 0 < alpha <= 1:
        raise ValueError('The fractional moment must lie in (0,1]')
    result, total = [arb(0) for _ in range(4)], F(0)
    exponent = rational(alpha)
    for weight, caps in weighted_caps:
        weight, caps = _exact(weight), tuple(map(_exact,caps))
        if weight < 0 or len(caps) != 4 or min(caps) < 0:
            raise ValueError('Expected a nonnegative weight and four positive-operator entries')
        total += weight
        for entry, cap in enumerate(caps):
            if cap:
                result[entry] += rational(weight)*rational(cap)**exponent
    if total != 1:
        raise ValueError('The exact orientation probabilities must sum to one')
    return tuple(result)


def swap_endpoint_receipt(z=DEFAULT_Z, alpha=DEFAULT_ALPHA, *, precision=256, bits=60):
    """Directed union bound for the two alternating supports in the swap ensemble."""
    z, alpha = _exact(z), _exact(alpha)
    if not 0 < z < 1 or not 0 < alpha <= 1 or precision < 128:
        raise ValueError('Invalid endpoint witness or arithmetic precision')
    with ctx.workprec(precision):
        caps = endpoint_packet_caps(z,bits)
        averaged = fractional_average(((F(n,2**32),row) for _,n,row in caps),alpha)
        powered = scalar_power(averaged,STEPS)
        log2 = lambda value:value.log()/arb(2).log()
        # The projective factor belongs INSIDE the fractional power. The
        # number of distinct support classes belongs OUTSIDE that power.
        prefix = rational(alpha)*(ACTIVE_GROUPS*(2048-128*log2(arb(255)))
                                   -log2(arb(255))-LIMIT*log2(rational(z)))
        bound = arb(1)+prefix+log2(powered[0]+powered[1])
        ceiling = _ceil(bound)
        return dict(construction='alternating-pairs-with-fresh-independent-pair-swaps',
                    profile=[4096,0,0],exact_support_count=2,field_order=256,
                    message_length=LENGTH//2,output_length=LENGTH,forbidden_weight=LIMIT,
                    local_rank=4,packet_width=64,packet_steps=STEPS,
                    orientation_average='E_swap[entrywise(packet_cap**alpha)]',
                    witness=dict(z=str(z),alpha=str(alpha),precision=precision,
                                 packet_mantissa_bits=bits),
                    log2_upper_rational=str(ceiling),log2_upper_enclosure=str(bound),
                    directed=True,whole_code_certificate=False)


def fixed_lane_quiet_coefficient(weight=22):
    """Exact [Z**weight] W(Z)**2 / 256**32 for the two occupied h4 blocks."""
    if not isinstance(weight,int) or not 0 <= weight <= 32:
        raise ValueError('Expected a weight in 0..32')
    coefficients = tuple(comb(16,w)*prescribed_mds_count(w,16,12,256)
                         for w in range(17))
    return F(sum(coefficients[w]*coefficients[weight-w]
                 for w in range(max(0,weight-16),min(16,weight)+1)),256**32)


def fixed_lane_obstruction(*, precision=256):
    """Uniform-in-z LOWER bound on a numerical comparator, not on failure."""
    if precision < 128:
        raise ValueError('At least 128 bits are required')
    coefficient = fixed_lane_quiet_coefficient(22)
    assert 22*STEPS <= LIMIT
    with ctx.workprec(precision):
        log2 = lambda value:value.log()/arb(2).log()
        bound = (ACTIVE_GROUPS*(2048-128*log2(arb(255)))-log2(arb(255))
                 +STEPS*log2(rational(coefficient)))
        return dict(kind='lower-bound-on-fixed-lane-comparator-not-on-failure',
                    profile=[4096,0,0],packet_weight=22,total_selected_weight=22*STEPS,
                    forbidden_weight=LIMIT,packet_steps=STEPS,
                    quiet_coefficient=str(coefficient),precision=precision,
                    log2_lower_rational=str(_floor(bound)),log2_enclosure=str(bound),
                    uniform_mark_domain='0<z<1',directed=True,
                    actual_distance_refutation=False,whole_code_certificate=False)


if __name__ == '__main__':
    import argparse,json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision',type=int,default=256)
    args = parser.parse_args()
    print(json.dumps(dict(obstruction=fixed_lane_obstruction(precision=args.precision),
                          swap=swap_endpoint_receipt(precision=args.precision)),indent=2))
