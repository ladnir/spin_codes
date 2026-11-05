"""Exact global setup-event bound for all-quiet groups of arity at most7.

This counts nonzero RS row words in the ideal independent-D1 setup law.
It uses no output-weight mark, regional provider, or distance gate.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import perf_counter

from hd_bit_syndrome_outer import affine_rank

DIRECTORY = Path(__file__).resolve().parent
SOURCE = 'hd_bit_syndrome_outer.py'
SOURCE_SHA = 'a2c2d0fee734fad2c4b74da7850589f8e84180798dcf5f79beab30e06939292d'


def quiet_histogram(maximum=7):
    if type(maximum) is not int or not 1 <= maximum <= 7:
        raise ValueError('This setup-event calculation covers arities1..7')
    counts = Counter()
    for mask in range(1, 1 << 16):
        a = mask.bit_count()
        if a <= maximum:
            counts[a, a-affine_rank(mask)] += 1
    for a in range(1, maximum+1):
        if sum(n for (weight, _), n in counts.items() if weight == a) != comb(16, a):
            raise ArithmeticError('Labelled source row sets do not partition')
    return counts


def event_upper(counts, groups=8192, q=256, columns=16, rs_dimension=8):
    if type(groups) is not int or groups < 1:
        raise ValueError('Positive source-group count required')
    per_arity = {}
    for (a, d), multiplicity in counts.items():
        if not 0 <= d <= a or multiplicity < 0:
            raise ValueError('Valid exact support/dimension histogram required')
        value = F(groups*multiplicity*(q**(columns*d)-1), (q-1)**(rs_dimension*a+1))
        per_arity[a] = per_arity.get(a, F(0))+value
    return per_arity, sum(per_arity.values(), F(0))


def run():
    from flint import arb, ctx
    start = perf_counter()
    if sha256((DIRECTORY/SOURCE).read_bytes()).hexdigest() != SOURCE_SHA:
        raise ValueError('Frozen outer algebra changed')
    histogram = quiet_histogram()
    per_arity, total = event_upper(histogram)
    if not total < F(1, 2**40):
        raise ArithmeticError('The exact rational40-bit setup-event test failed')
    logs = []
    for precision in (256, 384):
        with ctx.workprec(precision):
            value = (arb(total.numerator)/total.denominator).log()/arb(2).log()
            endpoint = value.upper().fmpq()
            upper = F(int(endpoint.p), int(endpoint.q))
            ceiling = F(-((-upper.numerator*10**12)//upper.denominator), 10**12)
            logs.append(dict(precision_bits=precision, log2_ceiling=str(ceiling)))
    return dict(construction='G16-RS16x8-D2-bit-syndrome-D1-research',
        event='Some one of8192 groups has a nonzero all-quiet RS word with1..7 active rows',
        sampling='Ideal independent uniform nonzero D1 entries sampled once at setup',
        arities=list(range(1, 8)), source_groups=8192, columns=16,
        rowsets_fixed_across_columns=True, output_weight_mark_used=False,
        histogram=[dict(arity=a, quiet_dimension=d, labelled_rowsets=n)
                   for (a, d), n in sorted(histogram.items())],
        exact_per_arity_probabilities={str(a): str(v) for a, v in per_arity.items()},
        exact_probability_upper=str(total), strictly_below_2_to_minus40=True,
        directed_log_checks=logs, regional_provider_runs=0, whole_code_certificate=False,
        dependency_sha256={SOURCE: SOURCE_SHA},
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        elapsed_seconds=perf_counter()-start)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite the setup-event receipt')
    result = run()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(result, indent=2, allow_nan=False))
