"""Outward replay of the complete one-active-group rank-three bound.

Uses the fixed sparse split at support 120 and the grid screened in
rank_three_flags.py. Does not certify rank four or multiple active groups.
"""
import argparse
from itertools import accumulate,combinations_with_replacement
from math import comb
from flint import arb,arb_mat,ctx

from bch_joint_support import authenticated_caps,support_caps
from basis_lattice import improve_caps,self_test as lattice_test
from shortened_bound import dimension_caps
from rank_two_types import subspace_cap,self_test as type_test
from rank_three_flags import TILTS,flag_test
from random_group_verify import epoch_transfers,transfer_tests
from random_rank_one import regions,conditional_numerators,self_test as transfer_self_test
from rank_two_moment import census
from group_rank_one_verify import up

PENALTIES=('1','0.875','0.75','0.625','0.5','0.375','0.25','0.125','0')
MAXIMUM=120


def penalty_regions(shapes,region,rho):
    weighted=[region[k+1]*rho**shape.count(4) for k,shape in enumerate(shapes)]
    result=[region[0]]
    for b in range(1,5):
        choices=[weighted[k] for k,shape in enumerate(shapes) if sum(w!=0 for w in shape)==b]
        result.append(arb_mat([[max(up(row[i,j]) for row in choices) for j in range(7)] for i in range(7)]))
    return result


def conditional_bounds(data):
    best=[[[[arb(1) for _ in range(257)] for _ in range(4)] for _ in range(64)] for _ in PENALTIES]
    for tilt in TILTS:
        shapes,epoch=epoch_transfers(data,tilt)
        region=regions(epoch)
        factor=(arb(tilt)*209715).exp()
        for r,text in enumerate(PENALTIES):
            rho=arb(text)
            last=256 if text in ('0','1') else MAXIMUM
            for l,values in conditional_numerators(penalty_regions(shapes,region,rho)):
                for b in range(1,5):
                    for u in range(67,min(last,4*l+b)+1):
                        value=up(values[u-b,b-1]*factor/comb(4*l,u-b))
                        best[r][l][b-1][u]=min(best[r][l][b-1][u],value)
        print('Outward',ctx.prec,'bits: completed flag tilt',tilt,flush=True)
    return best


def decreasing(values):
    for u in range(len(values)-2,-1,-1):
        values[u]=max(values[u],values[u+1])
    assert all(a>=b for a,b in zip(values,values[1:]))
    return values


def simple_probabilities(best):
    result=[arb(1)]*67
    for u in range(67,257):
        result.append(up(sum((comb(4,b)*comb(4*l,u-b)*best[l][b-1][u]
                              for l in range(64) for b in range(1,5) if 0<=u-b<=4*l),arb(0))/comb(256,u)))
    return decreasing(result)


def flag_probabilities(best,v,factors):
    result=[arb(1)]*(max(67,v))
    for u in range(max(67,v),MAXIMUM+1):
        total=arb(0)
        correction=[row[u-v] for row in factors]
        for l in range(64):
            for b in range(1,5):
                if 0<=u-b<=4*l:
                    p=min([arb(1)]+[up(best[r][l][b-1][u]*correction[r]) for r in range(8)])
                    total+=comb(4,b)*comb(4*l,u-b)*p
        result.append(up(total/comb(256,u)))
    return decreasing(result)


def cdf_sum(counts,weights):
    assert len(counts)==len(weights)
    assert all(a<=b for a,b in zip(counts,counts[1:]))
    assert all(a>=b for a,b in zip(weights,weights[1:]))
    return up(2048*(counts[-1]*weights[-1]+sum((counts[u]*(weights[u]-weights[u+1])
                                               for u in range(len(counts)-1)),arb(0))))


def partition_test():
    # Check sparse/dense CDF partition and decreasing-majorant summation with
    # exact integers before the final run. No difference of caps is a shell.
    counts=[0,2,3,5,7,11]
    weights=[arb(1)/2**u for u in range(6)]
    full=cdf_sum(list(accumulate(counts)),weights)
    for cut in range(5):
        low=cdf_sum(list(accumulate(counts[:cut+1])),weights[:cut+1])
        high_counts=list(accumulate([0 if u<=cut else c for u,c in enumerate(counts)]))
        high=cdf_sum(high_counts,weights)
        assert low+high==full
    print('Sparse/tail partition checked with exact dyadic weights',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision',type=int,default=192)
    args=parser.parse_args()
    assert args.precision>=128
    lattice_test()
    type_test()
    flag_test()
    transfer_self_test()
    partition_test()
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)[2]
    shells=[0]*(MAXIMUM+1)
    allowed=[w for w in range(1,257) if spectrum[w]]
    for weights in combinations_with_replacement(allowed,3):
        if sum(weights)%2==0 and sum(weights)//2<=MAXIMUM:
            shells[sum(weights)//2]+=subspace_cap(weights,spectrum)
    data=census()
    ctx.prec=args.precision
    transfer_tests(data)
    best=conditional_bounds(data)
    factors=[[up(arb(text)**-m) for m in range(MAXIMUM+1)] for text in PENALTIES[:-1]]
    outer_cdf=list(accumulate(spectrum))
    sparse=arb(0)
    for v in range(57,MAXIMUM+1):
        counts=[0 if u<max(67,v) else min(caps[u],168*shells[v]*outer_cdf[(2*u-v)//2])
                for u in range(MAXIMUM+1)]
        sparse+=cdf_sum(counts,flag_probabilities(best,v,factors))
    sparse=up(sparse)
    absent=cdf_sum([(cap*8)//15 for cap in caps],simple_probabilities(best[-1]))
    present=[0 if u<=MAXIMUM else (cap*7)//15 for u,cap in enumerate(caps)]
    dense=cdf_sum(present,simple_probabilities(best[0]))
    total=up(sparse+absent+dense)
    for name,value in (('sparse flags',sparse),('all no-all-one tuples',absent),('dense remainder',dense),('rank three total',total)):
        print(name,'outward upper',value,'margin',-value.log()/arb(2).log(),flush=True)
    assert 0<total<arb(2)**-40
    print('VERIFIED: all rank-three messages in one group, union over 2048 groups, <2^-40',flush=True)
    print('Full distance certificate: NO. Rank four and multiple active groups remain.')


if __name__=='__main__':
    main()
