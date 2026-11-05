"""Directed unsigned operators for fresh independent swaps inside routed pairs.

The oriented C coefficients already contain entrywise fractional packet
powers. Average these coefficients with weight 2**(-singleton_count);
do not take a fractional power after averaging the raw packet caps.
The good-good polynomial is unchanged because good-step caps are exactly
orientation-invariant. See HD_SHIFTED_PAIR_SWAPS.md for the conditioning
argument and the limits of inheritance from the fixed-lane construction.
"""
from fractions import Fraction as F
from functools import lru_cache

from flint import arb,ctx

from hd_shifted_cover_backend import CachedBound
from hd_shifted_general_certificate import general_coefficients
from hd_shifted_good_pair_certificate import good_coefficients


SWAP_CONSTRUCTION = 'toeplitz-gf256-h4-alternating-swapped-pairs-k20-rate-half'
BOUNDARY_BOX = ((4064,4128),(0,32),(0,32))
BOUNDARY_MARKERS = ((F(512),F(1)),)*2


def average_coefficients(oriented, *, width=33):
    """Return B[s,f] = 2**(-s) sum_u C[u,s-u,f], as Arb upper endpoints.

    The four flattened cubic arrays use axes (10,01,11). The width argument
    permits exact small fixtures; production packet arrays have width33.
    Inputs must be nonnegative coefficient enclosures, not frequency samples.
    """
    if (not isinstance(width,int) or width < 1 or len(oriented)!=4
            or any(len(row)!=width**3 for row in oriented)):
        raise ValueError('Expected four cubic oriented coefficient arrays')
    result = {}
    square = width*width
    for singleton in range(width):
        divisor = 1 << singleton
        for full in range(width-singleton):
            indices = tuple(first*square+(singleton-first)*width+full
                            for first in range(singleton+1))
            values = tuple((sum((row[index] for index in indices),arb(0))/divisor).upper()
                           for row in oriented)
            if any(not value.is_zero() for value in values):
                result[singleton,full] = values
    return result


@lru_cache(8)
def swap_unmarked(z,alpha,precision=256,bits=60):
    """Cache swap coefficients separately by exact z, alpha, precision, bits."""
    z,alpha = F(z),F(alpha)
    if not 0<z<1 or not 0<alpha<=1 or precision<128 or bits<16:
        raise ValueError('Invalid directed coefficient witness')
    with ctx.workprec(precision):
        oriented = general_coefficients(z,alpha,(F(1),)*3,bits=bits)
        single = average_coefficients(oriented)
        paired = good_coefficients(z,alpha,bits)
    return single,paired


class SwapCachedBound(CachedBound):
    """Positive, local-alias, or uniform-alias box bounds for the swap ensemble."""
    construction = SWAP_CONSTRUCTION

    @classmethod
    def unmarked_coefficients(cls,z,alpha,precision,bits):
        return swap_unmarked(z,alpha,precision,bits)

    def witness(self):
        return dict(super().witness(),pair_swaps=True,
                    orientation_aggregation='independent-fair-swaps')


def boundary_receipt(*, precision=256):
    """Replay a 47,377-profile box near the alternating-support boundary."""
    bound = SwapCachedBound(F(1,32),F(1,4),BOUNDARY_MARKERS,precision=precision)
    return dict(bound.receipt(BOUNDARY_BOX),witness=bound.witness())


if __name__=='__main__':
    import argparse,json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--grid',type=int,default=None)
    parser.add_argument('--z',default='1/32')
    parser.add_argument('--alpha',default='1/4')
    parser.add_argument('--profile',type=int,nargs=3,help='Default: the saved boundary box')
    parser.add_argument('--markers',nargs=4,default=('512','1','512','1'))
    args = parser.parse_args()
    markers = (args.markers[:2],args.markers[2:])
    bound = SwapCachedBound(F(args.z),F(args.alpha),markers,precision=args.precision,
                            grid=args.grid,progress=lambda item:print(json.dumps(item),flush=True))
    ranges = BOUNDARY_BOX if args.profile is None else tuple((n,n) for n in args.profile)
    print(json.dumps(dict(receipt=bound.receipt(ranges),
                          witness=bound.witness()),indent=2))
