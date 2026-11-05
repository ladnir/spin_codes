"""Exact local polynomial and conserved-route checks; no certificate."""
from fractions import Fraction as F
from itertools import combinations
import math
import shutil
import tempfile
import time
from test_random_pairing_birth_rows import matchings


def product_polynomial(factors):
    out=[F(1)]
    for a,b in factors:
        new=[F(0)]*(len(out)+1)
        for i,c in enumerate(out):new[i]+=c*a;new[i+1]+=c*b
        out=new
    return out


def main():
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic();polynomials=0;route_checks=0
    for n,w in ((4,2),(6,2),(6,3),(8,4)):
        H=(1<<w)-1;pairs=list(matchings(tuple(range(n))))
        next_masks=[sum(1<<j for j in c) for c in combinations(range(n),w)]
        for x,h,u in ((F(3,2),F(2,3),F(5,4)),(F(1,4),F(1,8),F(1,3))):
            direct=formula=envelope=F(0)
            for matching in pairs:
                live=0
                for a,b in matching:
                    pair=(1<<a)|(1<<b)
                    if H&pair:live|=pair
                k=(live&~H).bit_count()
                poly=product_polynomial([(x,F(1))]*w+[(h,F(1))]*k
                                        +[(F(1),F(1))]*(n-w-k))
                formula+=poly[w]
                envelope+=sum(c*u**j for j,c in enumerate(poly))/u**w
                for next_mask in next_masks:
                    J=(live&~(H|next_mask)).bit_count();D=(H&~next_mask).bit_count()
                    direct+=x**D*h**J
            assert direct==formula and formula<=envelope
            polynomials+=1
    M=4
    for s in range(1,M):
        p=F(s,M)
        probability=math.comb(M,s)*p**s*(1-p)**(M-s)
        assert probability>=F(1,M+1)
        for V in range(M+1):
            for v,t in ((F(1,4),F(2)),(F(3,4),F(1,3))):
                actual=sum(F(math.comb(V,j)*math.comb(M-V,s-j),math.comb(M,s))*v**j
                           for j in range(max(0,s-(M-V)),min(V,s)+1))
                a0=1-p+p*t;a1=1-p+p*t*v
                upper=t**(-s)*a0**(M-V)*a1**V/probability
                assert actual<=upper
                route_checks+=1
    assert time.monotonic()-began<25
    print(dict(joint_polynomial_checks=polynomials,conserved_target_checks=route_checks,
               seconds=time.monotonic()-began))


if __name__=='__main__':main()
