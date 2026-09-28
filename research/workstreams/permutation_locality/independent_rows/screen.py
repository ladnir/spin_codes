"""Exact rational outer-support screen; displayed logs are diagnostics only."""
import argparse
from fractions import Fraction as Q
from itertools import accumulate
from math import log2
from pathlib import Path
import sys

from support import union_cdf_upper, union_shells, fixed_weight_union, expected_union, weighted_cdf_upper


def bits(value):
    return log2(value.numerator)-log2(value.denominator) if value else float('-inf')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compare-shared', action='store_true')
    parser.add_argument('--penalty', action='append', default=[])
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from bch_joint_support import authenticated_caps, support_caps
    spectrum = authenticated_caps()
    upper = union_cdf_upper(spectrum, 1 << 128)
    print('NEW ENSEMBLE: independent row coordinate permutations; exact rational CDF caps')
    print('support,log2_independent_CDF_cap')
    for u in (38,57,64,80,96,128,144,160,176,192,224,256):
        print(u, f'{bits(upper[u]):.9f}', sep=',', flush=True)
    print('Known total check:', upper[-1] == (1 << 512)-1)
    # Shellwise caps are a different object from CDF differences.
    shells = union_shells([spectrum]*4)
    shells[0] -= 1
    capped_cdf = list(accumulate(shells))
    assert all(a <= b for a,b in zip(upper, capped_cdf))
    print('Known-total CDF tightens componentwise-shell-cap CDF at every support')
    print('four_equal_row_weights,expected_union,log2_probability_union_equals_row_weight')
    for w in (38,64,96,128,192,256):
        distribution = fixed_weight_union(256, [w]*4)
        print(w, float(expected_union(256,[w]*4)), bits(distribution[w]), sep=',')
    if args.compare_shared:
        from shortened_bound import dimension_caps
        from basis_lattice import improve_caps
        from shortening_moments import improve
        dimensions = dimension_caps()
        old, _ = improve(improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum),dimensions)
        print('support,log2_shared_CDF_cap,log2_independent_CDF_cap')
        for u in (38,57,64,80,96,128,144,160,176,192,224,256):
            shared = sum(row[u] for row in old)
            print(u, log2(shared) if shared else float('-inf'), bits(upper[u]), sep=',')
        print('Shared comparison uses basis-lattice and shortening-moment caps, not later dual refinements.')
    for text in args.penalty:
        rho = Q(text)
        if not 0 < rho <= 1:
            parser.error('penalty must be in (0,1]')
        weighted = weighted_cdf_upper(spectrum,1 << 128,full_weight=1/rho)
        print('All-one penalty',text,'support,log2_weighted_CDF_cap')
        for u in (38,80,96,128,144,160,176,192,224,256):
            print(u,bits(weighted[u]),sep=',',flush=True)
    print('Outer counts only. No inner transfer, full-code margin, or certificate follows.')


if __name__ == '__main__':
    main()
