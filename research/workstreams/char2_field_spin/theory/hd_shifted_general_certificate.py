"""Directed single-class receipts for arbitrary shifted-pair histograms.

The two regional operators use three independently chosen rational markers.
This certifies the displayed support class only; it is not whole-code
coverage. The interval FFT and exact packet caps are shared with the first
shifted-pair certificate.
"""
from fractions import Fraction as F

from flint import arb, ctx

from hd_shifted_alias_certificate import (
    LIMIT, _multinomial, odd_coefficients, positive_alias_power, rational,
)
from hd_shifted_general_count import (
    exact_joint_count, joint_count_upper, joint_run_count_upper,
    matching_histograms,
)
from hd_sumrank_paired_certificate import scalar_power


PAIRS = 4096
DEFAULT_PROFILE = (512, 512, 400, 400)
DEFAULT_MARK = F('0.4826158925972776')
DEFAULT_ALPHA = F('0.1026151436509744')
DEFAULT_MARKERS = (F('41044534/928270931'), F('12934738/292533347'),
                   F('195385295/968824422'))


def general_coefficients(z, alpha, markers, *, bits=60, progress=None):
    """Arb coefficients tilted independently along the 10, 01, 11 axes."""
    markers = tuple(map(F, markers))
    if len(markers) != 3 or min(markers) <= 0:
        raise ValueError('Three positive rational markers are required')
    return odd_coefficients(z, alpha, markers[0], markers[2],
                            singleton01=markers[1], bits=bits, progress=progress)


def truncate_coefficients(coefficients, target):
    """Drop input monomials that cannot contribute to the target power."""
    if len(target) != 3 or any(not isinstance(value, int) or value < 0 for value in target):
        raise ValueError('Three nonnegative integer target degrees are required')
    width = round(len(coefficients[0]) ** (1/3))
    if len(coefficients) != 4 or any(len(row) != width**3 for row in coefficients):
        raise ValueError('Expected four cubic coefficient arrays')
    result = [list(row) for row in coefficients]
    for a in range(width):
        for b in range(width):
            for c in range(width):
                if a > target[0] or b > target[1] or c > target[2]:
                    degree = a*width**2+b*width+c
                    for entry in range(4):
                        result[entry][degree] = arb(0)
    return result


def regional(z, alpha, histogram, markers, *, grid=64, bits=60, progress=None):
    """Enclose one regional operator for a fixed (00,10,01,11) histogram."""
    histogram, markers = tuple(histogram), tuple(map(F, markers))
    if (len(histogram) != 4 or any(not isinstance(n, int) or n < 0 for n in histogram)
            or sum(histogram) != PAIRS or len(markers) != 3 or min(markers) <= 0):
        raise ValueError('Invalid regional histogram or rational markers')
    target = histogram[1:]
    coeff = general_coefficients(z, alpha, markers, bits=bits, progress=progress)
    if any(n < 32 for n in target):
        coeff = truncate_coefficients(coeff, target)
    if not any(target):
        # The constant coefficient requires no similarity normalization.
        values = scalar_power(tuple(row[0] for row in coeff), 128)
        diagnostic = dict(method='constant-coefficient', directed=True)
    else:
        values, diagnostic = positive_alias_power(coeff, target, grid, progress=progress)
    normalizer = arb(_multinomial(histogram))
    for marker, degree in zip(markers, target):
        normalizer *= rational(marker)**degree
    values = tuple(entry/normalizer for entry in values)
    return values, dict(histogram=list(histogram), markers=list(map(str, markers)),
                        target=list(target), composition_steps=1, **diagnostic)


def receipt(profile=DEFAULT_PROFILE, z=DEFAULT_MARK, alpha=DEFAULT_ALPHA,
            markers=(DEFAULT_MARKERS, DEFAULT_MARKERS), *, grid=64, precision=256,
            bits=60, threshold=-40, support_count='exact', progress=None):
    """Bound one oriented profile (even,odd,even_full,odd_full)."""
    profile = tuple(profile)
    z, alpha = F(z), F(alpha)
    markers = tuple(tuple(map(F, row)) for row in markers)
    counters = dict(exact=exact_joint_count, upper=joint_count_upper,
                    run_upper=joint_run_count_upper)
    if (len(profile) != 4 or any(not isinstance(n, int) for n in profile)
            or support_count not in counters or len(markers) != 2
            or any(len(row) != 3 or min(row) <= 0 for row in markers)
            or not 0 < z < 1 or not 0 < alpha <= 1 or precision < 128):
        raise ValueError('Invalid profile, rational witness, precision, or support count mode')
    count = counters[support_count](PAIRS, *profile)
    if count <= 0 or sum(profile[:2]) == 0:
        raise ValueError('Empty or zero-message support class')
    histograms = matching_histograms(PAIRS, *profile)
    with ctx.workprec(precision):
        cache = {}
        regions, diagnostics = [], []
        for parity in range(2):
            key = (histograms[parity], markers[parity])
            if key not in cache:
                if progress:
                    progress(dict(stage='regional-operator', parity=parity, histogram=list(key[0])))
                cache[key] = regional(z, alpha, *key, grid=grid, bits=bits, progress=progress)
            region, diagnostic = cache[key]
            regions.append(region)
            diagnostics.append(dict(parity=parity, **diagnostic))
        first, second = regions
        paired = tuple(first[2*r]*second[s]+first[2*r+1]*second[2+s]
                       for r in range(2) for s in range(2))
        product = scalar_power(paired, 128)
        log2 = lambda value: value.log()/arb(2).log()
        active = sum(profile[:2])
        exponent = rational(alpha)*(active*(2048-128*log2(arb(255)))
                                   -log2(arb(255))-LIMIT*log2(rational(z)))
        bound = log2(arb(count))+exponent+log2(product[0]+product[1])
        ceiling = F(int((bound.upper()*10**12).upper().ceil().fmpz()), 10**12)
        return dict(construction='shifted-pairs-alternating-general-h4',
            profile=list(profile), active_groups=active, matching_histograms=[list(row) for row in histograms],
            mark=str(z), alpha=str(alpha), markers=[list(map(str, row)) for row in markers],
            grid=grid, precision=precision, packet_mantissa_bits=bits,
            support_count_method=support_count, support_count=str(count),
            log2_support_count=str(log2(arb(count))), forbidden_weight=LIMIT, field_order=256,
            log2_upper_enclosure=str(bound), log2_upper_rational=str(ceiling),
            accepted=bool(bound < threshold), threshold=threshold,
            regional_entry_log2_upper=[[str(log2(entry)) for entry in row] for row in regions],
            regional_diagnostics=diagnostics, distinct_regional_replays=len(cache),
            full_under_either_matching=any(hist[1] == hist[2] == 0 for hist in histograms),
            directed=True, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', type=int, nargs=4, default=DEFAULT_PROFILE)
    parser.add_argument('--mark', type=F, default=DEFAULT_MARK)
    parser.add_argument('--alpha', type=F, default=DEFAULT_ALPHA)
    parser.add_argument('--even-markers', type=F, nargs=3, default=DEFAULT_MARKERS)
    parser.add_argument('--odd-markers', type=F, nargs=3, default=DEFAULT_MARKERS)
    parser.add_argument('--grid', type=int, default=64)
    parser.add_argument('--precision', type=int, default=256)
    parser.add_argument('--bits', type=int, default=60)
    parser.add_argument('--threshold', type=int, default=-40)
    parser.add_argument('--support-count', choices=('exact', 'upper', 'run_upper'), default='exact')
    args = parser.parse_args()
    print(json.dumps(receipt(args.profile, args.mark, args.alpha,
        (args.even_markers, args.odd_markers), grid=args.grid, precision=args.precision,
        bits=args.bits, threshold=args.threshold, support_count=args.support_count,
        progress=lambda item: print(json.dumps(item), flush=True))), flush=True)
