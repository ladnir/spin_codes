"""Outward two-active-group bound for selected pairs of row-span ranks.

Uniform g=4,c=1 or c=2 routing at K=2^20, IMT(128,19), weight <=209715.
Includes same-epoch cancellation and every pair of group locations. Other
rank pairs and three or more active groups require separate bounds.
"""
import argparse
from math import comb

import numpy as np
from flint import arb, arb_mat, ctx

from bch_joint_support import authenticated_caps, support_caps
from shortened_bound import dimension_caps
from group_rank_one_verify import TILTS, up
from random_group_verify import epoch_transfers
from two_column_moment import census
from two_group_moment import collision_census, placement_regions, self_test
from two_group_screen import optimized_bound, log_binomial_mass, log_power_moment


def collision_transfers(data, tilt):
    spectrum, counts = data
    levels = sorted(spectrum)
    m = (1 << 19)-1
    powers = [up((-arb(tilt)*w).exp()) for w in range(145)]
    result = {}
    for (a,b), (choices,zero,maxatom,cancel) in counts.items():
        weight = sum(a)+sum(b)
        matrix = [[arb(0) for _ in range(7)] for _ in range(7)]
        matrix[0][0] = powers[weight]*zero/choices
        matrix[0][1] = powers[weight]*(choices-zero)/choices
        least = min(w for row in cancel.values() for w in row)
        for i in range(1,7):
            v = 48 if i == 1 else levels[i-2]
            f = powers[v-weight]
            if i == 1:
                c = min(up(f), up(powers[least]*maxatom/choices))
            else:
                c = sum((count*powers[w] for w,count in cancel[v].items()),arb(0))/(choices*spectrum[v])
            matrix[i][0] = c/2+f/(2*m)
            matrix[i][1] = f/2
            for j,w in enumerate(levels):
                matrix[i][j+2] = f*spectrum[w]/(2*m)
        result[a,b] = arb_mat([[up(x) for x in row] for row in matrix])
    return result


def maximum(matrices):
    matrices = list(matrices)
    return arb_mat([[max(a[i,j] for a in matrices) for j in range(7)] for i in range(7)])


def rounded(matrix):
    return arb_mat([[up(matrix[i,j]) for j in range(7)] for i in range(7)])


def regions(single_data,pair_data,tilt):
    width=len(next(iter(single_data[1])))
    epochs,windows=64*width,32//width
    occupancies=range(1,width+1)
    shapes, epoch = epoch_transfers(single_data,tilt,windows=windows)
    empty = epoch[0]
    single = {b: maximum(epoch[i+1] for i,s in enumerate(shapes) if sum(w!=0 for w in s)==b)
              for b in occupancies}
    pairs = collision_transfers(pair_data,tilt)
    collision = {(a,b): maximum(matrix for (s,t),matrix in pairs.items()
                               if sum(w!=0 for w in s)==a and sum(w!=0 for w in t)==b)
                 for a in occupancies for b in occupancies}
    one = {b: arb_mat(7,7) for b in occupancies}
    two = {key: arb_mat(7,7) for key in collision}
    same = {key: arb_mat(7,7) for key in collision}
    power = arb_mat([[int(i==j) for j in range(7)] for i in range(7)])
    for _ in range(epochs):
        two = {(a,b): rounded(two[a,b]*empty + one[a]*single[b] + one[b]*single[a]) for a,b in two}
        same = {key: rounded(same[key]*empty + power*collision[key]) for key in same}
        one = {b: rounded(one[b]*empty + power*single[b]) for b in one}
        power = rounded(power*empty)
    result = {(0,0): power}
    for b in occupancies:
        result[b,0] = result[0,b] = rounded(one[b]/epochs)
    for key in two:
        result[key] = rounded((windows*two[key]+(windows-1)*same[key])/(epochs*2047))
    # Independent binary64 placement implementation, consistency test only.
    def floating(matrix):
        return np.array([[float(matrix[i,j]) for j in range(7)] for i in range(7)])
    check = placement_regions(floating(empty),[floating(single[b]) for b in occupancies],
                              {key:floating(matrix) for key,matrix in collision.items()},epochs,windows)
    for key in result:
        assert np.allclose(floating(result[key]),check[key],rtol=2e-12,atol=1e-280)
    return result,{key:floating(matrix) for key,matrix in result.items()}


def buckets(spectrum,caps,rank,step):
    row = [sum(values) for values in zip(*caps)] if rank==0 else caps[rank-1]
    minimum = next(u for u,count in enumerate(row) if count)
    result = []
    boundaries=sorted(set(list(range(minimum,257,step))+[257]+[
        next(u for u,count in enumerate(values) if count) for values in caps if rank==0]))
    for lo,end in zip(boundaries,boundaries[1:]):
        hi=end-1
        count = row[hi]
        if rank == 1:
            occupied = [w for w in range(lo,hi+1) if spectrum[w]]
            if not occupied:
                continue
            lo,hi = occupied[0],occupied[-1]
            count = min(count,15*sum(spectrum[w] for w in occupied))
        result.append((lo,hi,count))
    return result


def select_witness(operators,first,second):
    lo,hi,_ = first
    low,high,_ = second
    best = None
    for tilt,(_,approximate) in operators.items():
        width=max(a for a,b in approximate)
        _,(p,q) = optimized_bound(approximate,(lo+hi)/2,(low+high)/2,float(tilt),True)
        # Replay the rounded rational witness, not an optimizer promise.
        denominator = 10**9
        pa = denominator if lo == 256 else max(1,min(denominator-1,round(p*denominator)))
        pb = denominator if low == 256 else max(1,min(denominator-1,round(q*denominator)))
        p,q = pa/denominator,pb/denominator
        matrix = sum(comb(width,a)*p**a*(1-p)**(width-a)*comb(width,b)*q**b*(1-q)**(width-b)*approximate[a,b]
                     for a in range(width+1) for b in range(width+1))
        value = log_power_moment(matrix,256//width)+float(tilt)*209715
        value -= min(log_binomial_mass(256,u,p) for u in (lo,hi))
        value -= min(log_binomial_mass(256,v,q) for v in (low,high))
        candidate = (value,tilt,pa,pb,denominator)
        if best is None or candidate[0] < best[0]:
            best = candidate
    return best[1:]


def rectangle_bound(region,tilt,pa,pb,denominator,first,second):
    p,q = arb(pa)/denominator,arb(pb)/denominator
    matrix = arb_mat(7,7)
    width=max(a for a,b in region)
    for a in range(width+1):
        for b in range(width+1):
            matrix += region[a,b]*(comb(width,a)*p**a*(1-p)**(width-a)*comb(width,b)*q**b*(1-q)**(width-b))
    powered = matrix**(256//width)
    value = sum((powered[0,j] for j in range(7)),arb(0))*(arb(tilt)*209715).exp()
    # Binomial mass is log-concave in its integer index. Its minimum on an
    # interval is attained at an endpoint; use a lower endpoint enclosure.
    for probability,(lo,hi,_) in ((p,first),(q,second)):
        masses = [arb(comb(256,u))*probability**u*(1-probability)**(256-u) for u in (lo,hi)]
        lower = min(arb(x.lower()) for x in masses)
        assert lower > 0
        value /= lower
    return min(arb(1),up(value))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--ranks',type=int,nargs=2,default=[1,1],help='Use 0 0 for all rank pairs')
    parser.add_argument('--step',type=int,default=4)
    parser.add_argument('--columns',type=int,choices=(1,2),default=2)
    args = parser.parse_args()
    assert args.precision >= 128 and args.step > 0
    h,k = sorted(args.ranks)
    assert (h==k==0) or 1 <= h <= k <= 4
    print('Geometry columns',args.columns,'ranks',h,k,'(0 means all ranks)',flush=True)
    self_test()
    spectrum = authenticated_caps()
    caps = support_caps(spectrum,g=4,dimensions=dimension_caps())
    single,pair = census(args.columns),collision_census(args.columns)
    ctx.prec = args.precision
    operators = {tilt: regions(single,pair,tilt) for tilt in TILTS}
    first,second = buckets(spectrum,caps,h,args.step),buckets(spectrum,caps,k,args.step)
    total = arb(0)
    worst = []
    for i,a in enumerate(first):
        for j,b in enumerate(second):
            if h == k and j < i:
                continue
            tilt,pa,pb,denominator = select_witness(operators,a,b)
            probability = rectangle_bound(operators[tilt][0],tilt,pa,pb,denominator,a,b)
            symmetry = 2 if h == k and i != j else 1
            term = up(symmetry*a[2]*b[2]*probability)
            total = up(total+term)
            worst.append((term,a[:2],b[:2]))
        if i%8 == 0:
            print('Outward',ctx.prec,'rank pair',h,k,'support bucket',i+1,'/',len(first),flush=True)
    total = up(total*comb(2048,2)*(2 if h != k else 1))
    print('Rank pair',h,k,'all group locations: upper',total,'margin',-total.log()/arb(2).log(),flush=True)
    print('Dominant support rectangles:',[(a,b) for _,a,b in sorted(worst,key=lambda item:float(item[0]),reverse=True)[:5]],flush=True)
    if 0 < total < arb(2)**-40:
        print('VERIFIED <2^-40 for '+('all two-group rank pairs' if h==k==0 else 'this rank pair only'),flush=True)
    else:
        print('This rank pair does NOT close at 40 bits with this bound',flush=True)
    print('Full certificate: NO. '+('Three or more active groups remain.' if h==k==0
                                    else 'Uncovered rank pairs and higher occupancies remain.'))


if __name__ == '__main__':
    main()
