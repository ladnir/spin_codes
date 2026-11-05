"""FLOAT construction search at 35% distance; no result is a certificate.

Keep N=2^21 and rate 1/2. Changing outer_length changes both the number
of independent RS rows and the number of independently permuted regions.
The old 40% construction and its pinned witnesses remain unchanged.
"""
from fractions import Fraction
from time import perf_counter

from hd_systematic_mds_screen import optimize

LENGTH = 1 << 21
DISTANCE = Fraction(7, 20)
FORBIDDEN_WEIGHT = (LENGTH * DISTANCE).__floor__()


def proposal(rows, *, outer_length=256, state=16, width=64,
             family='mds', iterations=80, regional=False, local_rank=4,
             forbidden_weight=FORBIDDEN_WEIGHT):
    options = {}
    if family == 'bucket':
        from hd_bucket_feedback_screen import paired_coefficients
        options['pair_source'] = lambda z, a, s, w: paired_coefficients(
            z, a, s, w, local_rank=local_rank)
    elif family == 'common-single':
        from hd_common_refresh_screen import matrices
        import numpy as np
        def single(z,a,s,w):
            degrees,weight,logs=matrices(z,s,w,h=local_rank)
            result=np.zeros((w+1,2,2))
            for r in range(2):
                for c in range(2):
                    result[:,r,c]=np.bincount(degrees,weights=np.exp2(weight+a*logs[:,r,c]),minlength=w+1)
            with np.errstate(divide='ignore'):return np.log2(result)
        options['pair_source']=single
        options['coefficient_width']=width
    elif family == 'common-average':
        from hd_common_refresh_screen import matrices
        from math import comb, log2
        import numpy as np
        def averaged(z, s, w):
            degrees, weight, logs = matrices(z,s,w,h=local_rank)
            result=np.zeros((w+1,2,2))
            for a in range(2):
                for b in range(2):
                    result[:,a,b]=np.bincount(degrees,weights=np.exp2(weight+logs[:,a,b]),minlength=w+1)
            with np.errstate(divide='ignore'):
                return np.log2(result)-np.array([log2(comb(w,j)) for j in range(w+1)])[:,None,None]
        options['packet_source'] = averaged
    elif family == 'common':
        from hd_common_refresh_screen import coefficients
        options['pair_source'] = lambda z, a, s, w: coefficients(
            z, a, s, w, h=local_rank)
    elif family == 'identity':
        options['mixer'] = False
    elif family == 'rank1':
        options['feedback'] = 'rank1_cone'
    elif family == 'paired-mds':
        from hd_systematic_mds_screen import coefficients, packet_logs
        options['pair_source'] = lambda z, a, s, w: coefficients(
            packet_logs(z, s, width=w), a)
    elif family != 'mds':
        raise ValueError('Unknown construction family')
    start = perf_counter()
    result = optimize(rows, state, iterations, regional, width=width,
                      outer_length=outer_length, forbidden_weight=forbidden_weight,
                      **options)
    result.update(candidate_family=family, local_rank=local_rank,
                  elapsed_seconds=perf_counter()-start)
    if family.startswith('common'):
        result.update(construction='common-diagonal-Cauchy-refresh',state_mixer=False)
    return result


def main():
    import argparse
    import json
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rows', nargs='+', type=int, default=(128, 512, 2048))
    p.add_argument('--outer-lengths', nargs='+', type=int, choices=(64, 128, 256), default=(256,))
    p.add_argument('--states', nargs='+', type=int, default=(16,))
    p.add_argument('--width', type=int, default=64)
    p.add_argument('--family', choices=('mds', 'paired-mds', 'identity', 'rank1', 'bucket', 'common', 'common-single', 'common-average'), default='mds')
    p.add_argument('--local-rank', type=int, default=4)
    p.add_argument('--regional', action='store_true')
    p.add_argument('--iterations', type=int, default=80)
    args = p.parse_args()
    for length in args.outer_lengths:
        for state in args.states:
            for rows in args.rows:
                print(json.dumps(proposal(rows, outer_length=length, state=state,
                    width=args.width, family=args.family, iterations=args.iterations,
                    regional=args.regional, local_rank=args.local_rank)), flush=True)


if __name__ == '__main__':
    main()
