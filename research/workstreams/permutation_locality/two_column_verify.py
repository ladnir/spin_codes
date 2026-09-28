"""Outward one-active-group bound for uniform shared g=4,c=1 or c=2 routing.

K=2^20, N=2^21, IMT(128,19), output weight <=209715. All four ranks and
2048 group locations are included; multiple active groups are not included.
No total-weight buckets, OA refinement, or orbit-memory extension is needed.
"""
import argparse
from math import comb
from flint import arb,arb_mat,ctx

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from group_rank_one_verify import TILTS,up
from random_rank_one import regions,conditional_numerators,self_test
from random_group_verify import epoch_transfers
from two_column_moment import census
from rank_two_moment import region_transfers


def worst_regions(data,tilt):
    width=len(next(iter(data[1])))
    shapes,epoch=epoch_transfers(data,tilt,windows=32//width)
    region=regions(epoch,epochs=64*width)
    result=[region[0]]
    for b in range(1,width+1):
        candidates=[region[j+1] for j,s in enumerate(shapes) if sum(w!=0 for w in s)==b]
        result.append(arb_mat([[max(row[i,j] for row in candidates) for j in range(7)] for i in range(7)]))
    return result


def transfer_tests(data):
    width=len(next(iter(data[1])))
    shapes,zero=epoch_transfers(data,'0',windows=32//width)
    assert all(not data[3][s][0] for s in shapes)
    for matrix in zero:
        for i in range(7):
            assert up(sum((matrix[i,j] for j in range(7)),arb(0)))>=1
    # Independently compare the epoch-averaging implementation with binary64.
    # This is a consistency test; the verifier below uses only outward values.
    for tilt in ('0','0.0004','0.0025'):
        _,epoch=epoch_transfers(data,tilt,windows=32//width)
        exact=regions(epoch,epochs=64*width)
        other_shapes,empty,active=region_transfers(data,float(tilt),windows=32//width,epochs=64*width)
        assert shapes==other_shapes
        for a,matrix in enumerate(exact):
            approximate=empty if a==0 else active[a-1]
            for i in range(7):
                for j in range(7):
                    assert abs(float(matrix[i,j])-approximate[i,j])<1e-12*max(1,abs(approximate[i,j]))
    print('Columns',width,'zero-tilt and independent region-transfer checks passed',flush=True)


def probabilities(data):
    width=len(next(iter(data[1])))
    length=256//width
    best=[[[arb(1) for _ in range(257)] for _ in range(width)] for _ in range(length)]
    for tilt in TILTS:
        factor=up((arb(tilt)*209715).exp())
        for l,values in conditional_numerators(worst_regions(data,tilt),length=length,width=width):
            for b in range(1,width+1):
                for u in range(b,min(257,width*l+b+1)):
                    value=up(values[u-b,b-1]*factor/comb(width*l,u-b))
                    best[l][b-1][u]=min(best[l][b-1][u],value)
        print('Outward',ctx.prec,'bits: columns',width,'tilt',tilt,flush=True)
    result=[arb(1)]
    for u in range(1,257):
        numerator=sum((comb(width,b)*comb(width*l,u-b)*best[l][b-1][u]
                       for l in range(length) for b in range(1,width+1) if 0<=u-b<=width*l),arb(0))
        result.append(min(arb(1),up(numerator/comb(256,u))))
    for u in range(255,-1,-1):
        result[u]=max(result[u],result[u+1])
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--columns',type=int,choices=(1,2),default=2)
    args=parser.parse_args()
    assert args.precision>=128
    self_test()
    spectrum=authenticated_caps()
    caps=support_caps(spectrum,g=4,dimensions=dimension_caps())
    data=census(args.columns)
    # Historical authentication imports change the global FLINT precision.
    ctx.prec=args.precision
    transfer_tests(data)
    weights=probabilities(data)
    assert all(a>=b for a,b in zip(weights,weights[1:]))
    totals=[]
    for h,row in enumerate(caps,1):
        total=up(2048*(row[-1]*weights[-1]+sum((row[u]*(weights[u]-weights[u+1]) for u in range(256)),arb(0))))
        totals.append(total)
        print('rank',h,'outward upper',total,'margin',-total.log()/arb(2).log(),flush=True)
    total=up(sum(totals,arb(0)))
    print('Combined single-group upper',total,'margin',-total.log()/arb(2).log(),flush=True)
    assert 0<total<arb(2)**-40
    print('VERIFIED: all single-active-group ranks, union over all 2048 groups, <2^-40',flush=True)
    print('Full certificate: NO. Multiple active groups remain uncovered.')


if __name__=='__main__':
    main()
