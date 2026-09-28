"""High-rank two-group diagnostic retaining the number of all-one columns.

Selected points only; exact integer flag counts but binary64 probabilities.
"""
import argparse
from itertools import accumulate
from math import comb,log

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from group_rank_one_verify import TILTS
from two_column_moment import census
from two_group_moment import collision_census,worst_regions
from two_group_memory import raw_regions
from two_group_weight import optimize_symmetric
from rank_four_flags import self_test


def flag_cap(u,a4,caps,dimensions,cdf):
    v = u-a4
    if not 0 <= v <= u <= 256 or not caps[2][v] or dimensions[v] < 3:
        return 0
    spaces = caps[2][v]//2520
    extensions = (comb(256-v,a4)*(1 << dimensions[v])-(8 if a4==0 else 0))//8
    extensions = min(extensions,cdf[(2*u-v)//2])
    return min(caps[3][u],1344*spaces*extensions)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--supports',type=int,nargs='+',default=[216])
    parser.add_argument('--ones',type=int,nargs='+',default=[0,16,32,64,96,128])
    parser.add_argument('--memory',action='store_true')
    args = parser.parse_args()
    self_test()
    spectrum = authenticated_caps()
    dimensions = dimension_caps()
    caps = support_caps(spectrum,g=4,dimensions=dimensions)
    cdf = list(accumulate([0]+spectrum[1:]))
    single,pair = census(),collision_census()
    counts = {(u,a):flag_cap(u,a,caps,dimensions,cdf) for u in args.supports for a in args.ones}
    best = {key:(0.,0.,1.) for key,count in counts.items() if count}
    for text in TILTS:
        tilt = float(text)
        raw = (raw_regions(single,pair,tilt) if args.memory
               else worst_regions(single,pair,tilt,return_shapes=True))
        for u,a in best:
            value,rho,_ = optimize_symmetric(raw,u,a,tilt,'all-ones')
            if value < best[u,a][0]:
                best[u,a] = value,tilt,rho
        print('BINARY64 all-one-count tilt',text,flush=True)
    for (u,a),(value,tilt,rho) in best.items():
        score = (value+2*log(counts[u,a])+log(comb(2048,2)))/log(2)
        print('point u,a4',u,a,'log2 probability',round(value/log(2),4),
              'log2 flag-count contribution',round(score,4),'tilt,rho',tilt,round(rho,6),flush=True)
    print('Selected symmetric points only; no interval or asymmetric coverage; no certificate.')


if __name__=='__main__':
    main()
