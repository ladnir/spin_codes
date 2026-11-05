"""Exact low-nonquiet-column setup event for the inverse Hstar schedule.

For fixed A and quiet ranks r_j, the enlarged array count with at most h
nonquiet columns is q^sum(a-r_j) [Y^<=h] prod(1+(q^r_j-1)Y)-1.
Only the global zero array is removed. The RS counting-measure domination,
projective divisor and source union are then applied once. This is a setup
event, not a distance certificate or a conditional independence assertion.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path

from hd_bit_syndrome_scheduled_expurgation import (
    CONSTRUCTION, DIMENSIONS_RECEIPT, DIRECTORY, INCIDENCE_RECEIPT,
    PINS as BASE_PINS, checked_histogram, deterministic_low_exclusions,
    log_ceiling,
)

MAXIMUM_NONQUIET = {4:4, 5:4, 6:3, 7:2, 8:1, 9:0}
DEPENDENCY_PINS = {**BASE_PINS,
    'hd_bit_syndrome_scheduled_expurgation.py':
    '7132f9f9019c576f46104188deb67dbef3fb4b57d8b65fb0c43082b7072ed4bd'}


def nonquiet_array_count(ranks, arity, maximum_nonquiet, *, q=256):
    """Exact enlarged-array count; ranks may describe arbitrary column kernels."""
    ranks = tuple(ranks)
    if (not ranks or type(arity) is not int or arity < 1
            or any(type(r) is not int or not 0 <= r <= arity for r in ranks)
            or type(maximum_nonquiet) is not int or not 0 <= maximum_nonquiet <= len(ranks)
            or type(q) is not int or q < 2):
        raise ValueError('Exact ranks, arity, field order and nonquiet cutoff required')
    coefficients = [1]+[0]*maximum_nonquiet
    for rank in ranks:
        weight = q**rank-1
        for degree in range(maximum_nonquiet, 0, -1):
            coefficients[degree] += weight*coefficients[degree-1]
    quiet_dimension = sum(arity-r for r in ranks)
    return q**quiet_dimension*sum(coefficients)-1


def rank_profile(arity, quiet_dimension):
    """Recover the rank multiset for a4..9 from its frozen dimension sum.

    Four labels have rank3 or4; five through nine labels rank4 or5.
    This unweighted event depends only on the multiset, not column order.
    """
    if type(arity) is not int or not 4 <= arity <= 9 or type(quiet_dimension) is not int:
        raise ValueError('The proved rank ranges cover only arities4..9')
    low = 3 if arity == 4 else 4
    low_count = quiet_dimension-16*(arity-low-1)
    if not 0 <= low_count <= 16:
        raise ValueError('Dimension incompatible with the proved two-rank range')
    return (low,)*low_count+(low+1,)*(16-low_count)


def event_bound(histogram, arity, maximum_nonquiet):
    """Union over8192 source positions and every fixed labelled row set."""
    if arity not in histogram or sum(histogram[arity].values()) != comb(16, arity):
        raise ValueError('Complete labelled row-set histogram required')
    count = 0
    for dimension, multiplicity in histogram[arity].items():
        if type(multiplicity) is not int or multiplicity < 0:
            raise ValueError('Nonnegative exact row-set multiplicities required')
        count += multiplicity*nonquiet_array_count(
            rank_profile(arity, dimension), arity, maximum_nonquiet)
    return Fraction(8192*count, 255**(8*arity+1))


def run():
    for name, digest in DEPENDENCY_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen dependency changed: '+name)
    histogram = checked_histogram(json.loads((DIRECTORY/DIMENSIONS_RECEIPT).read_text()))
    exclusions = deterministic_low_exclusions(json.loads((DIRECTORY/INCIDENCE_RECEIPT).read_text()))
    if 4 not in exclusions:
        raise ValueError('The separate deterministic four-row exclusion is missing')
    per = {a:event_bound(histogram,a,h) for a,h in MAXIMUM_NONQUIET.items()}
    cells = []
    for a,h in MAXIMUM_NONQUIET.items():
        checks = {str(p):str(log_ceiling(per[a],p)) for p in (256,384)}
        if len(set(checks.values())) != 1:
            raise ArithmeticError('Per-arity precision ceilings disagree')
        cells.append(dict(arity=a, maximum_nonquiet=h, complement_nonquiet_floor=h+1,
                          exact_probability_upper=str(per[a]), log2_upper=checks,
                          deterministic_event_empty=(a==4)))
    unions = []
    for maximum_arity in (8,9):
        for remove_four in (False,True):
            value = sum((v for a,v in per.items() if a <= maximum_arity and (a!=4 or not remove_four)),Fraction(0))
            checks = {str(p):str(log_ceiling(value,p)) for p in (256,384,512)}
            if len(set(checks.values())) != 1:
                raise ArithmeticError('Union precision ceilings disagree')
            unions.append(dict(maximum_arity=maximum_arity,
                variant='deterministic-four-row-event-removed' if remove_four else 'raw-dimension-count',
                exact_probability_upper=str(value), log2_upper=checks,
                strictly_below_2_to_minus45=value<Fraction(1,2**45)))
    return dict(construction=CONSTRUCTION,
        event='Some group has a nonzero RS word with arity a and at most h_a nonquiet columns',
        probability_space='Independent uniform nonzero D1 entries sampled once at setup',
        exact_enlarged_array_count=True, source_groups=8192, field_order=256,
        rs_length=16, rs_dimension=8, denominator_exponent='8*a+1', projective_divisor=255,
        arities=cells, unions=unions,
        deterministic_nonquiet_floors={'1':9,'2':9,'3':9,'4':5},
        arity9_unrestricted_in_maximum8_variant=True,
        conditioning_preserves_independence_claim=False,
        setup_event_only=True, whole_code_certificate=False, regional_provider_runs=0,
        dependency_sha256=DEPENDENCY_PINS,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == '__main__':
    # Serialization is owned by the caller; this command prints exact evidence.
    print(json.dumps(run(),indent=2))
