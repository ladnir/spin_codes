"""One-column occupancy operators for several active four-row groups.

Exact without-replacement placement counts; outward local transfer envelopes.
No full-code certificate is asserted by this module.
"""
from itertools import combinations,product
from math import comb,prod

import numpy as np
from flint import arb,arb_mat,ctx,fmpq_mat

from group_moment import maps
from group_rank_one_verify import up
from random_group_verify import epoch_transfers
from two_column_moment import census as single_census
from two_group_moment import collision_census
from two_group_verify import collision_transfers,maximum,rounded,regions as pair_regions
from window_feedback import census as kernel_census


def local_data(maximum_groups):
    assert 2<=maximum_groups<=4
    single,pair=single_census(1),collision_census(1)
    # Every two four-bit window maps have trivial intersection. This is the
    # hypothesis for the atom bound used at occupancies three and four.
    assert all(row[1]==0 for row in pair[1].values())
    _,columns,_=maps()
    zeros={q:kernel_census(columns,q)[1] for q in range(3,maximum_groups+1)}
    return single,pair,zeros


def atom_cap(weights,windows=32):
    # Condition on the other q-1 groups. A target syndrome is represented
    # in at most one available window and by at most one lane mask there.
    return arb(1)/((windows-len(weights)+1)*max(comb(4,w) for w in weights))


def epoch_operators(data,tilt):
    single,pair,zeros=data
    spectrum=single[0]
    levels=sorted(spectrum)
    states=(1<<19)-1
    _,one=epoch_transfers(single,tilt,windows=32)
    result=[one[0],maximum(one[1:]),maximum(collision_transfers(pair,tilt).values())]
    powers=[up((-arb(tilt)*w).exp()) for w in range(145)]
    for q,zero_counts in sorted(zeros.items()):
        candidates=[]
        slots=prod(range(32-q+1,33))
        for weights in product(range(1,5),repeat=q):
            weight=sum(weights)
            choices=slots*prod(comb(4,w) for w in weights)
            zero=arb(zero_counts[weights])/choices
            atom=atom_cap(weights)
            assert zero<=atom
            row=[[arb(0) for _ in range(7)] for _ in range(7)]
            row[0][0]=powers[weight]*zero
            row[0][1]=powers[weight]*(1-zero)
            for i in range(1,7):
                v=48 if i==1 else levels[i-2]
                f=powers[max(0,v-weight)]
                # Weighted cancellation <= f times the feedback atom bound
                # for an arbitrary state. For a uniform E-weight class it
                # is <= f*(1-zero)/class_size, without any uniformity claim
                # about the feedback distribution itself.
                c=f*atom if i==1 else f*(1-zero)/spectrum[v]
                row[i][0]=c/2+f/(2*states)
                row[i][1]=f/2
                for j,w in enumerate(levels):
                    row[i][j+2]=f*spectrum[w]/(2*states)
            candidates.append(arb_mat([[up(v) for v in r] for r in row]))
        result.append(maximum(candidates))
    return result


def placement(operators,epochs=64,windows=32,matrix=arb_mat,rounding=rounded,maximum_groups=None):
    """Coefficient of (sum binom(W,k) x^k T_k)^E / binom(EW,r)."""
    local_degree=min(windows,len(operators)-1)
    degree=local_degree if maximum_groups is None else maximum_groups
    assert 0<=degree<=epochs*windows
    assert local_degree>=min(degree,windows), 'missing local occupancy operators'
    n=operators[0].nrows()
    coefficients=[matrix([[int(i==j) for j in range(n)] for i in range(n)])]
    for e in range(epochs):
        next_coefficients=[]
        for r in range(min(degree,windows*(e+1))+1):
            value=matrix(n,n)
            for k in range(min(r,local_degree)+1):
                if r-k<len(coefficients):
                    value+=coefficients[r-k]*operators[k]*comb(windows,k)
            next_coefficients.append(rounding(value))
        coefficients=next_coefficients
    return [rounding(value/comb(epochs*windows,r)) for r,value in enumerate(coefficients)]


def build(data,tilt):
    exact=placement(epoch_operators(data,tilt))
    arrays=[np.array([[float(row[i,j]) for j in range(7)] for i in range(7)]) for row in exact]
    return exact,arrays


def self_test(data):
    # Direct slot-subset enumeration, including multiple occupied slots in
    # one epoch. Noncommuting matrices check chronological multiplication.
    operators=[fmpq_mat([[1,1],[0,1]]),fmpq_mat([[1,0],[1,1]]),
               fmpq_mat([[2,1],[0,1]]),fmpq_mat([[1,2],[1,1]])]
    for epochs in (1,2,3):
        for windows in (2,3):
            actual=placement(operators[:windows+1],epochs,windows,fmpq_mat,lambda x:x,
                             maximum_groups=epochs*windows)
            for r,value in enumerate(actual):
                total=fmpq_mat(2,2)
                for chosen in combinations(range(epochs*windows),r):
                    result=fmpq_mat([[1,0],[0,1]])
                    for e in range(epochs):
                        occupancy=sum(x//windows==e for x in chosen)
                        result=result*operators[occupancy]
                    total+=result
                assert value==total/comb(epochs*windows,r)
    for tilt in ('0','.0004','.0025'):
        new,_=build(data,tilt)
        old,_=pair_regions(data[0],data[1],tilt)
        # Rounded equivalent formulas need not coincide in their last bits.
        for a,b in ((0,0),(1,0),(1,1)):
            for i in range(7):
                for j in range(7):
                    assert abs(new[a+b][i,j]-old[a,b][i,j])<arb(2)**(-ctx.prec//2)
        if tilt=='0':
            for row in epoch_operators(data,tilt):
                assert all(sum((row[i,j] for j in range(7)),arb(0)).upper()>=1 for i in range(7))
    print('Exact occupancy polynomial and prior two-group regression checks passed',flush=True)


if __name__=='__main__':
    ctx.prec=192
    self_test(local_data(4))
