"""Nine-coordinate envelope retaining fresh feedback and mature-state density.

Coordinates: zero; fresh mass; mature mass; mature pointwise density cap;
five uniform E-weight classes. The density coordinate is NOT mass and is
excluded from terminal evaluation. No full-code certificate is asserted.
"""
import argparse
from collections import Counter
from itertools import product
from math import comb,log,prod

import numpy as np
from flint import arb,arb_mat,ctx

from group_rank_one_verify import up
from random_group_verify import epoch_transfers
from two_group_verify import collision_transfers
from occupancy_model import local_data,atom_cap,placement
from feedback_convolution import counts as convolution_counts
from occupancy_screen import optimize

Z,F,M,C,U=0,1,2,3,4
TERMINAL=np.array([1.,1.,1.,0.,1.,1.,1.,1.,1.])


def rounded(matrix):
    return arb_mat([[up(matrix[i,j]) for j in range(matrix.ncols())] for i in range(matrix.nrows())])


def maximum(matrices):
    matrices=list(matrices)
    return arb_mat([[max(t[i,j] for t in matrices) for j in range(9)] for i in range(9)])


def prepare(data):
    atoms={s:Counter({x:c for x,c in values.items() if c}) for s,values in data[0][3].items()}
    return data,convolution_counts(atoms)


def rational(value):
    return arb(value.numerator)/value.denominator


def shape_classes(q,with_full=False,with_odd=False):
    """Exact (total weight, largest binom(4,w)) classes of nonzero shapes."""
    result=set()
    for n1 in range(q+1):
        for n2 in range(q-n1+1):
            for n3 in range(q-n1-n2+1):
                n4=q-n1-n2-n3
                value=(n1+2*n2+3*n3+4*n4,6 if n2 else 4 if n1+n3 else 1)
                result.add(value+(n4,n1+n3) if with_odd else value+(n4,) if with_full else value)
    return result


def coarse_epoch(spectrum,q,tilt,pair_bound=None,input_penalty=1,full_penalty=1,odd_penalty=1):
    """All weight shapes, without a high-occupancy kernel enumeration.

    A feedback atom has probability <=1/((33-q)*max_i binom(4,w_i)).
    Only its upper bound is known; complements therefore use the lower
    bound zero, never one minus that upper bound.
    """
    assert 1<=q<=32
    levels=sorted(spectrum)
    states=(1<<19)-1
    powers=[up((-arb(tilt)*w).exp()) for w in range(129)]
    candidates=[]
    shapes=(shape_classes(q,True,True) if odd_penalty!=1 else
            {(w,c,f,0) for w,c,f in shape_classes(q,True)} if full_penalty!=1 else
            {(w,c,0,0) for w,c in shape_classes(q)})
    for weight,choices,full,odd in sorted(shapes):
        atom=up(arb(1)/((33-q)*choices))
        if pair_bound is not None:
            # Fix q-2 groups. At most the global pair census's peak number
            # of completions remains among the (34-q)(33-q) available
            # ordered windows. The distinct-window pair has no zero atom.
            atom=min(atom,up(rational(pair_bound)*32*31/((34-q)*(33-q))))
        f=powers[max(0,48-weight)]
        row=[[arb(0) for _ in range(9)] for _ in range(9)]
        row[Z][Z]=powers[weight]*atom
        row[Z][M]=powers[weight]
        row[Z][C]=powers[weight]*atom
        fresh=f*min(atom,arb(1)/32)
        row[F][Z]=fresh/2+f/(2*states)
        row[F][M]=f/2
        row[F][C]=fresh/2
        row[M][Z]=f/(2*states)
        row[M][M]=f/2
        row[C][Z]=f/2
        row[C][C]=f/2
        for i in (F,M):
            for j,w in enumerate(levels):
                row[i][U+j]=f*spectrum[w]/(2*states)
        for i,v in enumerate(levels):
            peak=powers[max(0,v-weight)]
            row[U+i][Z]=peak/(2*spectrum[v])+peak/(2*states)
            row[U+i][M]=peak/2
            row[U+i][C]=peak/(2*spectrum[v])
            for j,w in enumerate(levels):
                row[U+i][U+j]=peak*spectrum[w]/(2*states)
        candidate=arb_mat(row)
        if input_penalty!=1:
            candidate*=arb(input_penalty)**weight
        if full_penalty!=1:
            candidate*=arb(full_penalty)**full
        if odd_penalty!=1:
            candidate*=arb(odd_penalty)**odd
        candidates.append(rounded(candidate))
    return maximum(candidates)


def epoch_operators(prepared,tilt,detailed=False,maximum_groups=None,pair_conditioned=False,input_penalty=1,full_penalty=1,odd_penalty=1):
    (single,pair,zeros),convolutions=prepared
    spectrum,allowed,_,_,cancel=single
    levels=sorted(spectrum)
    states=(1<<19)-1
    powers=[up((-arb(tilt)*w).exp()) for w in range(145)]
    shapes,one=epoch_transfers(single,tilt,windows=32)
    pairs=collision_transfers(pair,tilt)
    empty=[[arb(0) for _ in range(9)] for _ in range(9)]
    empty[Z][Z]=arb(1)
    for i in (F,M):
        f=powers[48]
        empty[i][i]=f/2
        for j,v in enumerate(levels):
            empty[i][U+j]=f*spectrum[v]/(2*states)
    empty[C][C]=powers[48]/2
    for i,v in enumerate(levels):
        f=powers[v]
        empty[U+i][U+i]+=f/2
        for j,w in enumerate(levels):
            empty[U+i][U+j]+=f*spectrum[w]/(2*states)
    result=[rounded(arb_mat(empty))]
    known_groups=max(zeros,default=2)
    maximum_groups=known_groups if maximum_groups is None else maximum_groups
    assert known_groups<=maximum_groups<=32
    for q in range(1,known_groups+1):
        candidates={}
        for weights in product(range(1,5),repeat=q):
            weight=sum(weights)
            f=powers[max(0,48-weight)]
            if q==1:
                b=(weights[0],)
                choices=32*len(allowed[b])
                zero=arb(0)
                feedback_cap=arb(1)/choices
                old=one[1+shapes.index(b)]
                # Equal feedback requires identical window and lane mask.
                # Every fresh shape except b has disjoint feedback support.
                fresh_cancel=sum((count*powers[w] for row in cancel[b].values()
                                  for w,count in row.items()),arb(0))/(choices*choices)
                fresh_cap=f*max(up(rational(convolutions[tuple(sorted((a,b)))][1])) for a in allowed)
            elif q==2:
                a,b=((weights[0],),(weights[1],))
                choices,z,peak,_=pair[1][a,b]
                zero=arb(z)/choices
                feedback_cap=arb(peak)/choices
                old=pairs[a,b]
                fresh_cancel=f*feedback_cap
                fresh_cap=f*feedback_cap
            else:
                choices=prod(range(32-q+1,33))*prod(comb(4,w) for w in weights)
                zero=arb(zeros[q][weights])/choices
                feedback_cap=atom_cap(weights)
                old=None
                fresh_cancel=f*feedback_cap
                fresh_cap=f*feedback_cap
            row=[[arb(0) for _ in range(9)] for _ in range(9)]
            row[Z][Z]=powers[weight]*zero
            if q==1:
                row[Z][F]=powers[weight]
            else:
                row[Z][M]=powers[weight]*(1-zero)
                row[Z][C]=powers[weight]*feedback_cap
            row[F][Z]=fresh_cancel/2+f/(2*states)
            row[F][M]=f/2
            row[F][C]=fresh_cap/2
            row[M][Z]=f/(2*states)
            row[M][M]=f/2
            row[C][Z]=f*(1-zero)/2
            row[C][C]=f/2
            for i in (F,M):
                for j,w in enumerate(levels):
                    row[i][U+j]=f*spectrum[w]/(2*states)
            for i,v in enumerate(levels):
                peak=powers[max(0,v-weight)]
                if old is None:
                    moment=peak
                    cancellation=peak*(1-zero)/spectrum[v]
                    row[U+i][Z]=cancellation/2+moment/(2*states)
                else:
                    moment=2*old[i+2,1]
                    row[U+i][Z]=old[i+2,0]
                row[U+i][M]=moment/2
                row[U+i][C]=peak/(2*spectrum[v])
                for j,w in enumerate(levels):
                    row[U+i][U+j]=moment*spectrum[w]/(2*states)
            candidate=arb_mat(row)
            if input_penalty!=1:
                candidate*=arb(input_penalty)**weight
            if full_penalty!=1:
                candidate*=arb(full_penalty)**weights.count(4)
            if odd_penalty!=1:
                candidate*=arb(odd_penalty)**sum(w%2 for w in weights)
            candidates[weights]=rounded(candidate)
        result.append(candidates if detailed else maximum(candidates.values()))
    assert not detailed or maximum_groups==known_groups
    from fractions import Fraction
    pair_bound=max(Fraction(row[2],row[0]) for row in pair[1].values()) if pair_conditioned else None
    result.extend(coarse_epoch(spectrum,q,tilt,pair_bound,input_penalty,full_penalty,odd_penalty) for q in range(known_groups+1,maximum_groups+1))
    return result


def build(prepared,tilt,maximum_groups=None,pair_conditioned=False):
    exact=placement(epoch_operators(prepared,tilt,maximum_groups=maximum_groups,pair_conditioned=pair_conditioned),rounding=rounded)
    arrays=[np.array([[float(t[i,j]) for j in range(9)] for i in range(9)]) for t in exact]
    return exact,arrays


def self_test(prepared):
    matrices=epoch_operators(prepared,'0')
    terminal=arb_mat([[int(x)] for x in TERMINAL])
    for t in matrices:
        mass=t*terminal
        assert all(mass[i,0].upper()>=terminal[i,0] for i in range(9))
    # The auxiliary cap must not be counted as probability mass.
    assert TERMINAL[C]==0 and all(TERMINAL[i]==1 for i in range(9) if i!=C)
    assert matrices[0][C,Z]==0 and matrices[1][Z,C]==0
    assert all(matrices[1][i,F]==0 for i in range(1,9))
    print('Memory-envelope zero-tilt mass and coordinate checks passed',flush=True)


def main():
    from bch_joint_support import authenticated_caps,support_caps
    from shortened_bound import dimension_caps
    from basis_lattice import improve_caps
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,choices=range(3,33),default=3)
    parser.add_argument('--pair-conditioned',action='store_true')
    parser.add_argument('--tilts',nargs='+',default=['.00064','.001','.00125','.001375','.0016','.002','.0025','.0032','.004','.005','.0064','.008','.01'])
    args=parser.parse_args()
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)
    prepared=prepare(local_data(min(args.groups,4)))
    ctx.prec=192
    self_test(prepared)
    points=[(u,)*args.groups for u in (38,57,72,76,80,83,84,96,128,176,216,256)]
    best={s:(0.,None) for s in points}
    for tilt in args.tilts:
        _,region=build(prepared,tilt,args.groups,args.pair_conditioned)
        for supports in points:
            value,_=optimize(region,supports,float(tilt),TERMINAL)
            if value<best[supports][0]:
                best[supports]=value,tilt
        print('BINARY64 memory screen tilt',tilt,flush=True)
    results=[]
    for supports,(value,tilt) in best.items():
        score=(value+sum(log(sum(row[u] for row in caps)) for u in supports)+log(comb(2048,args.groups)))/log(2)
        results.append((score,supports,value/log(2),tilt))
    for row in sorted(results,reverse=True):
        print(row,flush=True)
    print('Selected points only; no complete support coverage or outward union.')


if __name__=='__main__':
    main()
