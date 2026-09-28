"""Outward one-active-group bounds with worst-shape transfer envelopes.

Checks rank two against 2^-50. Reports other ranks without certifying them.
Multiple active groups and full SPIN distance remain outside this bound.
"""
import argparse
from math import comb
from flint import arb,arb_mat,ctx

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from group_rank_one_verify import TILTS,up
from random_rank_one import regions,conditional_numerators,transfers as rank_one_transfers
from rank_two_moment import census


def epoch_transfers(data,tilt,windows=8):
    spectrum,allowed,moments,atoms,cancel=data
    levels=sorted(spectrum)
    m=(1<<19)-1
    powers=[up((-arb(tilt)*w).exp()) for w in range(145)]
    empty=[[arb(0) for _ in range(7)] for _ in range(7)]
    empty[0][0]=arb(1)
    empty[1][1]=powers[48]/2
    for j,w in enumerate(levels):
        empty[1][j+2]=powers[48]*spectrum[w]/(2*m)
    for i,v in enumerate(levels):
        empty[i+2][i+2]=powers[v]/2
        for j,w in enumerate(levels):
            empty[i+2][j+2]+=powers[v]*spectrum[w]/(2*m)
    shapes=sorted(allowed)
    result=[arb_mat([[up(v) for v in row] for row in empty])]
    for shape in shapes:
        choices=windows*len(allowed[shape])
        zero=atoms[shape][0]
        active=[[arb(0) for _ in range(7)] for _ in range(7)]
        active[0][0]=powers[sum(shape)]*zero/choices
        active[0][1]=powers[sum(shape)]*(choices-zero)/choices
        f=powers[48-sum(shape)]
        maxatom=max(count for syndrome,count in atoms[shape].items() if syndrome)
        least=min(w for row in cancel[shape].values() for w in row)
        c=min(up(f),up(powers[least]*maxatom/choices))
        for i in range(1,7):
            if i>=2:
                v=levels[i-2]
                f=up(sum((count*powers[w] for w,count in moments[shape][v].items()),arb(0))/(choices*spectrum[v]))
                c=up(sum((count*powers[w] for w,count in cancel[shape][v].items()),arb(0))/(choices*spectrum[v]))
            active[i][0]=c/2+f/(2*m)
            active[i][1]=f/2
            for j,w in enumerate(levels):
                active[i][j+2]=f*spectrum[w]/(2*m)
        result.append(arb_mat([[up(v) for v in row] for row in active]))
    return shapes,result


def worst_regions(data,tilt):
    shapes,epoch=epoch_transfers(data,tilt)
    region=regions(epoch)
    result=[region[0]]
    for b in range(1,5):
        candidates=[region[k+1] for k,shape in enumerate(shapes) if sum(w!=0 for w in shape)==b]
        result.append(arb_mat([[max(row[i,j] for row in candidates) for j in range(7)] for i in range(7)]))
    return result


def transfer_tests(data):
    spectrum,allowed,moments,atoms,cancel=data
    keyed_moments,keyed_max,keyed_cancel={},{},{}
    for a in range(1,5):
        for b in range(1,5):
            shape=(0,)*(4-b)+(a,)*b
            assert atoms[shape][0]==0
            keyed_moments[a,b]=moments[shape]
            keyed_max[a,b]=max(atoms[shape].values())
            keyed_cancel[a,b]=cancel[shape]
    for tilt in ('0','0.0004','0.0025'):
        shapes,epoch=epoch_transfers(data,tilt)
        for a in range(1,5):
            reference=rank_one_transfers((spectrum,keyed_moments,keyed_max,keyed_cancel),tilt,a)
            for b in range(5):
                actual=epoch[0] if b==0 else epoch[1+shapes.index((0,)*(4-b)+(a,)*b)]
                for i in range(7):
                    for j in range(7):
                        assert abs(actual[i,j]-reference[b][i,j])<arb(2)**(-ctx.prec//2)
        if tilt=='0':
            for matrix in epoch:
                for i in range(7):
                    assert up(sum((matrix[i,j] for j in range(7)),arb(0)))>=1
    print('All rank-one orbit transfers match the prior implementation; zero-tilt mass checks passed',flush=True)


def probabilities(data):
    best=[[[arb(1) for _ in range(257)] for _ in range(4)] for _ in range(64)]
    for tilt in TILTS:
        factor=(arb(tilt)*209715).exp()
        for l,values in conditional_numerators(worst_regions(data,tilt)):
            for b in range(1,5):
                for u in range(b,min(257,4*l+b+1)):
                    value=up(values[u-b,b-1]*factor/comb(4*l,u-b))
                    best[l][b-1][u]=min(best[l][b-1][u],value)
        print('Outward',ctx.prec,'bits: completed worst-shape tilt',tilt,flush=True)
    result=[arb(1)]
    for u in range(1,257):
        value=sum((comb(4,b)*comb(4*l,u-b)*best[l][b-1][u]
                   for l in range(64) for b in range(1,5) if 0<=u-b<=4*l),arb(0))/comb(256,u)
        result.append(min(arb(1),up(value)))
    for u in range(255,-1,-1):
        result[u]=max(result[u],result[u+1])
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision',type=int,default=192)
    args=parser.parse_args()
    assert args.precision>=128
    spectrum=authenticated_caps()
    caps=support_caps(spectrum,g=4,dimensions=dimension_caps())
    data=census()
    ctx.prec=args.precision
    transfer_tests(data)
    weights=probabilities(data)
    assert all(a>=b for a,b in zip(weights,weights[1:]))
    for h,row in enumerate(caps,1):
        # Differences are of exact dyadic probability majorants, not CDF caps.
        total=up(2048*(row[-1]*weights[-1]+sum((row[u]*(weights[u]-weights[u+1]) for u in range(256)),arb(0))))
        print('rank',h,'outward upper',total,'margin',-total.log()/arb(2).log(),flush=True)
        if h==2:
            assert 0<total<arb(2)**-50
            print('VERIFIED: all rank-two messages in one group, union over 2048 groups, <2^-50',flush=True)
    print('Full distance certificate: NO. Higher ranks and multiple groups remain.')


if __name__=='__main__':
    main()
