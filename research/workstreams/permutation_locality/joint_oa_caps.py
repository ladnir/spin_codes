"""Exact Christoffel shell caps for joint supports under OA strength 29.

This component checks its algebra; it does not replay the BCH dual-distance
certificate. The caller must establish the indicated orthogonal-array strength.
"""
from fractions import Fraction as F
from itertools import product
from math import comb,log2


def krawtchouk(n,q,u,degree):
    out=[1]
    if degree:
        out.append((q-1)*n-q*u)
    for j in range(1,degree):
        numerator=((q-1)*n-q*u-(q-2)*j)*out[-1]-(q-1)*(n-j+1)*out[-2]
        value,remainder=divmod(numerator,j+1)
        assert remainder==0
        out.append(value)
    return out


def shell_caps(n,k,h,strength=29):
    q=1<<h;degree=min(n,strength//2)
    norms=[comb(n,j)*(q-1)**j for j in range(degree+1)]
    result=[]
    for u in range(n+1):
        values=krawtchouk(n,q,u,degree)
        kernel=sum((F(v*v,norm) for v,norm in zip(values,norms)),F(0))
        bound=F(1<<(k*h))/kernel
        result.append(bound.numerator//bound.denominator)
    return result


def self_test():
    for n in (4,8,12):
        for q in (2,4,8,16):
            degree=min(n,5)
            rows=[krawtchouk(n,q,u,degree) for u in range(n+1)]
            for u,row in enumerate(rows):
                for j,value in enumerate(row):
                    direct=sum((-1)**l*(q-1)**(j-l)*comb(u,l)*comb(n-u,j-l)
                               for l in range(max(0,j-(n-u)),min(u,j)+1))
                    assert value==direct
            for i in range(degree+1):
                for j in range(degree+1):
                    inner=sum(comb(n,u)*(q-1)**u*rows[u][i]*rows[u][j] for u in range(n+1))
                    assert inner==(q**n*comb(n,i)*(q-1)**i if i==j else 0)
    words=[w for w in range(16) if w.bit_count()%2==0]
    for h in (1,2,3,4):
        actual=[0]*5
        for tup in product(words,repeat=h):
            union=0
            for w in tup:
                union|=w
            actual[union.bit_count()]+=1
        assert all(a<=b for a,b in zip(actual,shell_caps(4,3,h,strength=3)))
    print('Joint OA polynomial recurrence, orthogonality, and small-code shell bounds pass exact checks',flush=True)


if __name__=='__main__':
    self_test()
    for h in (3,4):
        caps=shell_caps(256,128,h)
        print('rank',h,[(u,log2(caps[u])) for u in (72,96,128,160,176,192,208,224,240,256)],flush=True)
