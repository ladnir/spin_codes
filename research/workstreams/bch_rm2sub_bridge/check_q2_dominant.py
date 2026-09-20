"""Independent 512-bit Arb check of the dominant weight pair, read-only.

Uses polynomial binary powering for regions and normalized conditional
support probabilities for the outer pair (not the certificate count loop).
"""
import math
import sys
from fractions import Fraction as F

import numpy as np
from flint import arb,ctx
import bridge as base


def regions(t,s,spectrum,z):
    den=(1<<s)-1; kappa=arb(den)/(den-1); d=min(spectrum)
    epoch=[]
    for j in range(3):
        moment=arb(0)
        for w,count in spectrum.items():
            for v in range(j+1):
                if v<=w and j-v<=t-w:
                    moment+=arb(count*math.comb(w,v)*math.comb(t-w,j-v))*z**(w+j-2*v)/(den*math.comb(t,j))
        if j==0:
            matrix=[[arb(1),arb(0),arb(0)],[arb(0),arb(0),z**d],[arb(0),arb(0),kappa*moment]]
        else:
            matrix=[[arb(0),z**j,arb(0)],[z**(d-j)/den,arb(0),z**(d-j)],
                    [kappa*moment/den,arb(0),kappa*moment]]
        epoch.append([[v*math.comb(t,j) for v in row] for row in matrix])
    def mul(a,b):
        return [[sum((a[i][k]*b[k][j] for k in range(3)),arb(0)) for j in range(3)] for i in range(3)]
    def polymul(a,b):
        output=[]
        for degree in range(3):
            terms=[mul(a[j],b[degree-j]) for j in range(degree+1)]
            output.append([[sum((m[i][k] for m in terms),arb(0)) for k in range(3)] for i in range(3)])
        return output
    identity=[[arb(int(i==j)) for j in range(3)] for i in range(3)]
    empty=[[arb(0) for _ in range(3)] for _ in range(3)]
    result=[identity,empty,empty]; power=epoch; exponent=8192//t
    while exponent:
        if exponent&1: result=polymul(result,power)
        exponent>>=1
        if exponent: power=polymul(power,power)
    return [[[v/math.comb(8192,j) for v in row] for row in result[j]] for j in range(3)]


def pair(region,a,b):
    current=[[(arb(1),arb(0),arb(0))]]
    for n in range(1,257):
        output=[]
        for i in range(min(n,a)+1):
            row=[]
            for j in range(min(n,b)+1):
                out=[arb(0),arb(0),arb(0)]
                for x,y in ((0,0),(1,0),(0,1),(1,1)):
                    old_i,old_j=i-x,j-y
                    if not (0<=old_i<len(current) and 0<=old_j<len(current[old_i])):continue
                    factor=(i if x else n-i)*(j if y else n-j)
                    if not factor:continue
                    v=current[old_i][old_j]; matrix=region[x+y]
                    for k in range(3):
                        out[k]+=(v[0]*matrix[0][k]+v[1]*matrix[1][k]+v[2]*matrix[2][k])*factor/(n*n)
                row.append(tuple(out))
            output.append(row)
        current=output
    return sum(current[a][b],arb(0))


if __name__=='__main__':
    ctx.prec=512; name='t128_s15'
    t,s,spectrum=base.load_map(name)
    receipt=base.read(base.HERE/'generated'/f'{name}_q2_j-75.json')
    a,b=receipt['dominant_pairs'][0]['weights']
    lam=(arb(-75)/10).exp(); z=(-lam).exp()
    value=pair(regions(t,s,spectrum,z),a,b)*(209716*lam).exp()
    with np.load(base.HERE/'generated'/f'{name}_q2_j-75.npz') as saved:
        stored=float(saved['conditional_pair_upper'][a,b])
    assert value<arb(stored)
    ratio=arb(stored)/value
    assert ratio<arb(1)+arb(1)/10**10
    print('Independent normalized 512-bit Arb check passed:',(a,b),'stored/exact-envelope',ratio,flush=True)
