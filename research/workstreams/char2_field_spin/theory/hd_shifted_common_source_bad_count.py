"""Exact rare-route-setup lemma for r18 sources with Fmax exactly one.

This checks rational coupon probabilities and an integer factorial-moment
union bound. It builds no regional operator and proves no distance bound.
Regional routes are independent for each fixed G; sources need not be.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path

from hd_shifted_common_source_partition import good_probability


DIRECTORY = Path(__file__).resolve().parent
PROBABILITIES = 'hd_shifted_common_source_analytic_gate_float.json'
PROBABILITIES_SHA = '85642b69e127d3a02429019336c9546f14cb738b33230da2f11a3e42e6910780'
PARTITION = 'hd_two_sum_gf16_source_strata_anchors_float.json'
PARTITION_SHA = '5050033694bb6a1ef7852aa309ee3abca252c9740fc713a83aa8945a380fb1cf'
COUPON_SOURCE_SHA = 'fd4bfc68dbb14f53851e00be6d162705f3c9f10cbbc1edffd6417b96918f48f3'
CONSTRUCTION = 'two-sum-gf16-coset-g16-rs16x8-h4-swapped-pairs-k20-rate-half'
SOURCE_COUNT = 154355220469939768130077702840442839490778956275433472
REGIONS, THRESHOLD = 256, 13


def _probability(value):
    value = F(value)
    if not 0 <= value <= 1:
        raise ValueError('A probability in [0,1] is required')
    return value


def count_distribution(probabilities):
    """Exact Poisson-binomial coefficients; intended for small reference tests."""
    answer = [F(1)]
    for probability in probabilities:
        p = _probability(probability)
        following = [F(0)]*(len(answer)+1)
        for degree, value in enumerate(answer):
            following[degree] += (1-p)*value
            following[degree+1] += p*value
        answer = following
    return tuple(answer)


def two_parity_factorial_bound(p_even, p_odd, per_parity, threshold):
    """[Y^b](1+p_even Y)^n(1+p_odd Y)^n bounds Pr[B>=b]."""
    p_even, p_odd = _probability(p_even), _probability(p_odd)
    if (type(per_parity) is not int or per_parity < 1
            or type(threshold) is not int or not 1 <= threshold <= 2*per_parity):
        raise ValueError('Positive parity size and feasible threshold required')
    return sum((F(comb(per_parity, left)*comb(per_parity, threshold-left))
                *p_even**left*p_odd**(threshold-left)
                for left in range(max(0, threshold-per_parity), min(per_parity, threshold)+1)), F(0))


def union_bound(source_count, probability_max, regions, threshold):
    """Uncapped rational bound N*C(regions,b)*pmax^b, counting sources once."""
    p = _probability(probability_max)
    if (type(source_count) is not int or source_count < 1
            or type(regions) is not int or regions < 1
            or type(threshold) is not int or not 1 <= threshold <= regions):
        raise ValueError('Positive exact counts and feasible threshold required')
    return source_count*comb(regions, threshold)*p**threshold


def load_inputs():
    probability_bytes = (DIRECTORY/PROBABILITIES).read_bytes()
    partition_bytes = (DIRECTORY/PARTITION).read_bytes()
    coupon_bytes = (DIRECTORY/'hd_shifted_common_source_partition.py').read_bytes()
    if (sha256(probability_bytes).hexdigest() != PROBABILITIES_SHA
            or sha256(partition_bytes).hexdigest() != PARTITION_SHA
            or sha256(coupon_bytes).hexdigest() != COUPON_SOURCE_SHA):
        raise ValueError('Reviewed exact-input byte pin changed')
    data, partition = json.loads(probability_bytes), json.loads(partition_bytes)
    if (data.get('construction') != CONSTRUCTION or data.get('active_groups') != 18
            or data.get('source_good_event_common_to_all_masks') is not True
            or partition.get('source_size') != 8192):
        raise ValueError('The fixed common-source construction is required')
    cases = data.get('cases', [])
    if [(row.get('potential_singletons'), row.get('potential_full_pairs')) for row in cases] != [(18, 0), (16, 1)]:
        raise ValueError('Exactly the two r18/Fmax1 potential histograms are required')
    probabilities = []
    for row in cases:
        good = good_probability(row['potential_singletons'], row['potential_full_pairs'])
        bad = F(row['bad_probability'])
        if good != F(row['good_probability']) or bad != 1-good:
            raise ArithmeticError('Saved rational differs from exact coupon coefficient')
        probabilities.append(bad)
    banks = [row for row in partition.get('banks', []) if row.get('active_groups') == 18]
    if len(banks) != 1:
        raise ValueError('One exact source-count partition is required')
    rows = banks[0].get('source_strata', [])
    counts = {row['full_max']: row['count'] for row in rows}
    if (len(rows) != 10 or set(counts) != set(range(10))
            or any(type(value) is not int or value <= 0 for value in counts.values())
            or sum(counts.values()) != comb(8192, 18) or counts[1] != SOURCE_COUNT):
        raise ArithmeticError('The exact Fmax1 source count changed')
    return tuple(probabilities)


def exact_lemma():
    probabilities = load_inputs()
    maximum = max(probabilities)
    bound = union_bound(SOURCE_COUNT, maximum, REGIONS, THRESHOLD)
    # This integer comparison, not the FLOAT logarithm, proves the margin.
    if (bound.numerator << 63) >= bound.denominator:
        raise ArithmeticError('The requested strict 63-bit rare-setup bound failed')
    return dict(construction=CONSTRUCTION, active_groups=18, source_size=8192,
        source_full_pair_max_exactly=1, source_count=str(SOURCE_COUNT), regions=REGIONS,
        even_regions=128, odd_regions=128, bad_region_threshold=THRESHOLD,
        retained_bad_regions_at_most=THRESHOLD-1,
        event='some size18 Fmax1 source has at least13 common-source bad regions',
        good_block_capacity=4, emission_block_length=16,
        bad_probabilities_by_full_pairs=list(map(str, probabilities)),
        probability_max=str(maximum), exact_union_bound=str(bound),
        exact_strict_bound='2^-63', exact_integer_check=True,
        diagnostic_log2_FLOAT=log2(bound.numerator)-log2(bound.denominator),
        independence='regional routes for each fixed source; not between sources',
        source_count_once=True, rowset_or_message_count_charged=False,
        conditioning_on_global_good_event=False, distance_bound=False,
        whole_code_certificate=False, regional_provider_calls=0,
        probability_artifact=PROBABILITIES, probability_artifact_sha256=PROBABILITIES_SHA,
        source_partition=PARTITION, source_partition_sha256=PARTITION_SHA,
        exact_coupon_source_sha256=COUPON_SOURCE_SHA,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite an exact lemma receipt')
    answer = exact_lemma()
    with args.output.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({key: answer[key] for key in ('exact_strict_bound', 'diagnostic_log2_FLOAT',
        'retained_bad_regions_at_most', 'distance_bound')}))
