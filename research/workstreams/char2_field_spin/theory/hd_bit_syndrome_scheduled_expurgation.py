"""Exact setup-event unions for the fixed inverse-column Hstar schedule.

For a fixed original RS row set A, the sixteen quiet spaces have total
dimension D_A.  The ideal independent-D1 union is
8192 * sum_A (256**D_A - 1) / 255**(8*len(A) + 1).
The scalar quotient is projective, not a row-set or column multiplicity.
This helper consumes authenticated dimension counts; it builds no provider.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import perf_counter


DIRECTORY = Path(__file__).resolve().parent
CONSTRUCTION = 'G16-Hstar-column-inverse-input-permutations-research'
DIMENSIONS_RECEIPT = 'hd_bit_syndrome_column_dimensions_exact.json'
INCIDENCE_RECEIPT = 'hd_bit_syndrome_column_permutations_exact.json'
PINS = {
    'hd_bit_syndrome_column_dimensions.py': '765e27e8da66a6c2946eafb888653815f9f50b70e6042e1ce99de091981bf3fc',
    DIMENSIONS_RECEIPT: 'c78f3c8caa755e6d9fcbbfdd7cd01df56341a4191294373ca75e0ea907e16e1c',
    'hd_bit_syndrome_column_permutations.py': 'fecd7c914a221901aad426a7d93f93b603805642377288431a26c38a9cc8163f',
    INCIDENCE_RECEIPT: '1d7b0e5b86c378e32f28ce48730bee51444ece43b5792ce495a5939ba512d6bd',
    'hd_bit_syndrome_outer.py': 'a2c2d0fee734fad2c4b74da7850589f8e84180798dcf5f79beab30e06939292d',
}
EXPECTED = {
    4: {0: 120, 1: 1200, 2: 480, 4: 20},
    5: {5: 1488, 6: 1440, 7: 960, 8: 240, 9: 240},
    6: {16: 480, 17: 2640, 18: 4080, 19: 640, 20: 120, 21: 48},
    7: {32: 7600, 33: 3840},
    8: {48: 12390, 49: 480},
    9: {64: 11440},
}


def checked_histogram(summary):
    """Validate the complete labelled arity4..9 dimension partition."""
    if (summary.get('construction') != CONSTRUCTION
            or summary.get('field_order') != 256
            or summary.get('columns') != 16
            or summary.get('all_requested_arities_complete') is not True
            or summary.get('source_rows_fixed_across_columns') is not True):
        raise ValueError('Complete fixed-row inverse-schedule dimensions required')
    histogram = {}
    for record in summary.get('records', []):
        a = record.get('arity')
        if type(a) is not int or a not in EXPECTED or a in histogram:
            raise ValueError('Each requested arity must occur exactly once')
        counts = {int(d): n for d, n in record['quiet_dimension_sum_histogram'].items()}
        if (record.get('complete') is not True
                or record.get('labelled_sets_checked') != comb(16, a)
                or record.get('expected_labelled_sets') != comb(16, a)
                or any(type(n) is not int or n <= 0 for n in counts.values())
                or counts != EXPECTED[a]
                or sum(counts.values()) != comb(16, a)):
            raise ValueError('Incomplete or changed labelled dimension histogram')
        histogram[a] = counts
    if set(histogram) != set(EXPECTED):
        raise ValueError('Missing arity in the dimension partition')
    # Any <=3 distinct augmented binary labels are independent over GF256.
    return {1: {0: 16}, 2: {0: 120}, 3: {0: 560}, **histogram}


def deterministic_low_exclusions(receipt):
    """Require the separately authenticated RS support-incidence exclusion."""
    four = receipt.get('four_row_summary', {})
    incidence = receipt.get('activity_incidence', {})
    five = [r for r in incidence.get('results', []) if r.get('arity') == 5]
    if (four.get('construction') != CONSTRUCTION
            or four.get('fixed_source_sets_checked') != comb(16, 4)
            or four.get('maximum_quiet_columns') != 4
            or incidence.get('construction') != CONSTRUCTION
            or incidence.get('all_requested_arities_complete') is not True
            or len(five) != 1
            or five[0].get('complete') is not True
            or five[0].get('labelled_sets_checked') != comb(16, 5)
            or five[0].get('excluded_all_quiet_sets') != comb(16, 5)
            or five[0].get('not_excluded_sets') != 0
            or five[0].get('universal_nonquiet_column_floor') != 2):
        raise ValueError('The deterministic arity4/5 exclusion was not established')
    return (1, 2, 3, 4, 5)


def event_union(histogram, maximum, *, excluded_arities=(), groups=8192,
                q=256, row_dimension=8):
    """Pure rational count; histogram dimensions already sum over columns."""
    if (type(maximum) is not int or maximum < 1
            or type(groups) is not int or groups < 1
            or type(q) is not int or q < 2
            or type(row_dimension) is not int or row_dimension < 1):
        raise ValueError('Positive exact event parameters required')
    excluded = set(excluded_arities)
    if any(type(a) is not int or not 1 <= a <= maximum for a in excluded):
        raise ValueError('Excluded arities must belong to the requested event')
    per = {}
    for a in range(1, maximum + 1):
        if a not in histogram:
            raise ValueError('Missing arity histogram')
        value = Fraction(0)
        for dimension, multiplicity in histogram[a].items():
            if (type(dimension) is not int or dimension < 0
                    or type(multiplicity) is not int or multiplicity < 0):
                raise ValueError('Nonnegative exact dimensions and counts required')
            if a not in excluded:
                value += Fraction(groups * multiplicity * (q**dimension - 1),
                                  (q - 1)**(row_dimension*a + 1))
        per[a] = value
    return per, sum(per.values(), Fraction(0))


def log_ceiling(value, precision):
    """Directed Arb upper endpoint, rounded upward to 10^-12 bits."""
    from flint import arb, ctx
    if not isinstance(value, Fraction) or value <= 0 or precision < 64:
        raise ValueError('Positive rational and adequate Arb precision required')
    with ctx.workprec(precision):
        log_value = (arb(value.numerator) / value.denominator).log() / arb(2).log()
        endpoint = log_value.upper().fmpq()
        upper = Fraction(int(endpoint.p), int(endpoint.q))
        scale = 10**12
        return Fraction(-((-upper.numerator * scale) // upper.denominator), scale)


def run():
    started = perf_counter()
    if not PINS or DIMENSIONS_RECEIPT not in PINS or INCIDENCE_RECEIPT not in PINS:
        raise ValueError('Final reviewed dependencies must be pinned before execution')
    for name, digest in PINS.items():
        if sha256((DIRECTORY / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Frozen dependency changed: {name}')
    dimensions = json.loads((DIRECTORY / DIMENSIONS_RECEIPT).read_text(encoding='utf-8'))
    summary = dimensions.get('summary', dimensions)
    histogram = checked_histogram(summary)
    incidence = json.loads((DIRECTORY / INCIDENCE_RECEIPT).read_text(encoding='utf-8'))
    exclusions = deterministic_low_exclusions(incidence)
    results = []
    for maximum in (8, 9):
        for exclusion in ((), exclusions):
            per, total = event_union(histogram, maximum, excluded_arities=exclusion)
            logs = [dict(precision_bits=p, log2_ceiling=str(log_ceiling(total, p)))
                    for p in (256, 384)]
            if logs[0]['log2_ceiling'] != logs[1]['log2_ceiling']:
                raise ArithmeticError('The independently rounded ceilings disagree')
            results.append(dict(maximum_arity=maximum,
                variant='dimension-only' if not exclusion else 'deterministic-arity1-through5-excluded',
                deterministically_excluded_arities=list(exclusion),
                exact_per_arity_probability_upper={str(a): str(v) for a, v in per.items()},
                exact_probability_upper=str(total), directed_log_checks=logs,
                strictly_below_2_to_minus40=total < Fraction(1, 2**40)))
    return dict(construction=CONSTRUCTION,
        event='Some source group has a nonzero RS word whose scaled permuted columns are all quiet',
        probability_space='Ideal independent uniform NZ D1 entries sampled once at setup',
        source_groups=8192, field_order=256, rs_length=16, rs_dimension=8,
        rowsets_fixed_across_columns=True, denominator_exponent='8*a+1',
        projective_divisor=255, output_weight_mark_used=False,
        setup_event_only=True, whole_code_certificate=False, native_implementation_claim=False,
        histogram={str(a): {str(d): n for d, n in sorted(counts.items())}
                   for a, counts in sorted(histogram.items())}, results=results,
        dependency_sha256=PINS, source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        regional_provider_runs=0, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite an exact setup-event receipt')
    result = run()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: result[k] for k in ('results', 'elapsed_seconds')}, indent=2))
