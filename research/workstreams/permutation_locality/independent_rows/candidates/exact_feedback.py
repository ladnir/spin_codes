"""Exact feedback peaks beyond six windows using guarded limb inversion.

This opt-in candidate leaves the shared-route census and production cache
unchanged. Two signed-int64 Walsh transforms replace a potentially wider
single transform. Reconstruction divides by the state count BEFORE shifting
the high limb, and checks the remaining integer bounds before arithmetic.
"""
import argparse
from fractions import Fraction as Q

import numpy as np

from mass_density_screen import baseline
from feedback_character_census import character_polynomials
from pair_tail import hadamard
from group_moment import maps


def inverse_counts(coefficients, denominator):
    values=np.asarray(coefficients)
    if (values.ndim!=1 or values.dtype!=np.dtype(np.int64) or not len(values)
            or len(values)&(len(values)-1) or len(values)>1<<30
            or not isinstance(denominator,int) or isinstance(denominator,bool)
            or not 0<denominator<1<<63 or int(values[0])!=denominator):
        raise ValueError('signed-int64 transform, power-of-two size, and exact denominator required')
    size=len(values)
    magnitude=max(abs(int(values.min())),abs(int(values.max())))
    if magnitude>denominator:
        raise ValueError('character coefficient exceeds the total count')
    if magnitude*size<1<<63:
        numerator=hadamard(values)
        if np.any(numerator%size):
            raise ValueError('nonintegral inverse transform')
        result=numerator//size
    else:
        # Floor division is intentional, including for negative coefficients.
        # values = high * 2^32 + low, with 0 <= low < 2^32.
        low=values&np.int64((1<<32)-1)
        high=values>>32
        assert size*((1<<32)-1)<1<<63
        assert size*max(abs(int(high.min())),abs(int(high.max())))<1<<63
        low=hadamard(low); high=hadamard(high)
        # size divides 2^32, so integrality requires size to divide W(low).
        if np.any(low%size):
            raise ValueError('nonintegral low-limb inverse transform')
        low//=size
        scale=(1<<32)//size
        bound=scale*max(abs(int(high.min())),abs(int(high.max())))
        bound+=max(abs(int(low.min())),abs(int(low.max())))
        if bound>=1<<63:
            raise ValueError('reconstruction exceeds its signed-int64 guard')
        result=high*scale+low
    if int(result.min())<0 or int(result.max())>denominator:
        raise ValueError('inverse is not a nonnegative counting measure')
    # Exact inversion and values[0] already imply this identity. Python
    # integers provide an independent check without an unchecked int64 sum.
    if sum(map(int,result))!=denominator:
        raise ValueError('inverse does not preserve total count')
    return result


def build(maximum=10,known=None):
    if not isinstance(maximum,int) or isinstance(maximum,bool) or not 1<=maximum<=10:
        raise ValueError('one through ten occupied windows supported')
    shapes,rows,denominators,multiplicities,inverse=character_polynomials(maximum)
    images,columns,spectrum=maps()
    size=len(images)
    assert all(c.bit_count()%2==1 for c in columns)
    states=np.arange(size,dtype=np.uint32)
    parity=np.bitwise_count(states)&1
    expansion=np.array([image.bit_count() for image in images],dtype=np.uint8)
    signs={a:1-2*(np.bitwise_count(states&a)&1).astype(np.int64) for a in (1,17,size-1)}
    indices={v:expansion==v for v in spectrum}
    result={}; checked=0
    for number,(shape,row) in enumerate(zip(shapes,rows)):
        weights=tuple(b for b in range(1,5) for _ in range(shape[b-1]))
        if not weights:
            continue
        denominator=int(denominators[shape])
        counts=inverse_counts(row[inverse],denominator)
        assert not np.any(counts[parity!=(sum(weights)%2)])
        zero=int(counts[0]); peak=int(counts[1:].max())
        classes={v:int(counts[index].sum()) for v,index in indices.items()}
        assert zero+sum(classes.values())==denominator
        assert zero*size==sum(int(a)*int(b) for a,b in zip(row,multiplicities))
        # Since counts are nonnegative with total < 2^63, every partial
        # signed sum is bounded in absolute value by that total.
        for a,sign in signs.items():
            assert int(counts@sign)==int(row[inverse[a]])
        record=zero,peak,denominator,classes
        if known is not None and weights in known:
            z,p,d,c=known[weights]
            assert z*denominator==zero*d and p*denominator==peak*d
            assert all(c[v]*denominator==classes[v]*d for v in c)
            checked+=1
        result[weights]=record
        if number%100==0:
            print('EXACT LIMB FEEDBACK shapes',number,'of',len(shapes)-1,flush=True)
    for j in range(1,maximum+1):
        worst=max((Q(max(z,p),d),shape) for shape,(z,p,d,c) in result.items() if len(shape)==j)
        print('EXACT FEEDBACK all-target peak',j,worst,flush=True)
    print('EXACT LIMB FEEDBACK independent census comparisons',checked,
          '; all totals, parity, zero and selected character checks pass',flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maximum',type=int,default=10)
    args=parser.parse_args()
    from full_feedback_census import census
    build(args.maximum,census(min(2,args.maximum)))
