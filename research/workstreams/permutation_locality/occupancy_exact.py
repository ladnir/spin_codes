"""Direct trivariate positive coefficients for the difficult three-group range.

Binary64 diagnostic only. Covers all support vectors in a truncated cube,
but neither outward rounding nor support vectors outside that cube.
"""
import argparse
from itertools import combinations,product
from math import comb,log

import numpy as np
from scipy.special import logsumexp
from flint import ctx

from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from basis_lattice import improve_caps
from occupancy_model import local_data,build
from occupancy_screen import optimize


def coefficients(region,maximum,length=256):
    q=len(region)-1
    n=len(region[0])
    current=np.zeros((1,)*q+(n,))
    current[(0,)*q+(0,)]=1.
    for step in range(length):
        width=current.shape[0]
        next_width=min(maximum+1,width+1)
        updated=np.zeros((next_width,)*q+(n,))
        for r in range(q+1):
            value=(current.reshape(-1,n)@region[r]).reshape(current.shape)
            for active in combinations(range(q),r):
                shift=tuple(int(i in active) for i in range(q))
                sizes=tuple(min(width,next_width-d) for d in shift)
                target=tuple(slice(d,d+s) for d,s in zip(shift,sizes))
                source=tuple(slice(0,s) for s in sizes)
                updated[target]+=value[source]
        current=updated
        if length==256 and (step+1)%64==0:
            print('Direct coefficients: regions',step+1,'/',length,'axis',next_width,flush=True)
    return current.sum(axis=-1)


def self_test():
    region=[np.array([[1.,1.],[0.,1.]]),np.array([[1.,0.],[1.,1.]]),
            np.array([[2.,1.],[0.,1.]]),np.array([[1.,2.],[1.,1.]])]
    for length in (1,2,3):
        for maximum in range(length+1):
            actual=coefficients(region,maximum,length)
            expected=np.zeros((maximum+1,)*3)
            for supports in product(range(1<<length),repeat=3):
                counts=tuple(s.bit_count() for s in supports)
                if max(counts)>maximum:
                    continue
                row=np.array([1.,0.])
                for j in range(length):
                    row=row@region[sum((s>>j)&1 for s in supports)]
                expected[counts]+=row.sum()
            assert np.array_equal(actual,expected)
    print('Direct support coefficients agree with exhaustive noncommuting toy enumeration',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maximum',type=int,default=100)
    parser.add_argument('--tilts',nargs='+',default=['.00064','.001','.00125','.0016'])
    parser.add_argument('--test-only',action='store_true')
    args=parser.parse_args()
    self_test()
    if args.test_only:
        return
    assert 38<=args.maximum<=128
    spectrum=authenticated_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimension_caps()),spectrum)
    counts=np.array([sum(row[u] for row in caps) for u in range(args.maximum+1)],dtype=object)
    data=local_data(3)
    ctx.prec=192
    shape=(args.maximum+1,)*3
    normalization=np.array([log(comb(256,u)) for u in range(args.maximum+1)])
    denominator=normalization[:,None,None]+normalization[None,:,None]+normalization[None,None,:]
    best=np.zeros(shape)
    for tilt in args.tilts:
        _,region=build(data,tilt)
        numerators=coefficients(region,args.maximum)
        assert np.all(numerators[38:,38:,38:]>0),'binary64 underflow: cannot use this diagnostic'
        with np.errstate(divide='ignore'):
            bounds=np.log(numerators)+float(tilt)*209715-denominator
        for u in (38,57,72,76,80,96,100):
            if u<=args.maximum:
                cauchy,_=optimize(region,(u,)*3,float(tilt))
                assert bounds[u,u,u]<=cauchy+1e-9,('coefficient exceeds Cauchy',u,tilt)
        best=np.minimum(best,bounds)
        print('BINARY64 direct-coefficient tilt',tilt,'point log2',
              [(u,round(best[u,u,u]/log(2),4)) for u in (38,57,72,76,80,96) if u<=args.maximum],flush=True)
    log_counts=np.array([log(c) if c else -np.inf for c in counts])
    terms=best+log_counts[:,None,None]+log_counts[None,:,None]+log_counts[None,None,:]+log(comb(2048,3))
    print('BINARY64 complete cube 38..',args.maximum,'union log2',logsumexp(terms)/log(2),flush=True)
    flat=np.argsort(terms.ravel())[-10:][::-1]
    print('Dominant supports:',[(np.unravel_index(int(i),shape),float(terms.ravel()[i]/log(2))) for i in flat],flush=True)
    from occupancy_cdf import majorant,differences,self_test as cdf_test
    cdf_test()
    # Differences are a dominating measure ONLY for decreasing test functions;
    # they are not upper bounds on unknown shell counts.
    monotone=majorant(best)
    increments=differences(counts)
    log_increments=np.array([log(c) if c else -np.inf for c in increments])
    cdf_terms=monotone+log_increments[:,None,None]+log_increments[None,:,None]+log_increments[None,None,:]+log(comb(2048,3))
    print('BINARY64 decreasing-majorant CDF cube union log2',logsumexp(cdf_terms)/log(2),flush=True)
    flat=np.argsort(cdf_terms.ravel())[-10:][::-1]
    print('Dominant CDF supports:',[(np.unravel_index(int(i),shape),float(cdf_terms.ravel()[i]/log(2))) for i in flat],flush=True)
    print('No certificate: binary64 arithmetic and omitted supports above maximum.')


if __name__=='__main__':
    main()
