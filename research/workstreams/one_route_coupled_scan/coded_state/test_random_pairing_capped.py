"""Exact checks for capped geometry and the Bernoulli comparison.

Standalone, serial <=30s, no data files; can run through stdin remotely.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import product
import math
import time
import shutil
import tempfile


def matchings(vertices):
    if not vertices:
        yield ();return
    a=vertices[0]
    for i in range(1,len(vertices)):
        for rest in matchings(vertices[1:i]+vertices[i+1:]):yield ((a,vertices[i]),)+rest


def run():
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic();domination=0;capped=0;identities=0
    for n in (2,4,6,8):
        all_pairs=list(matchings(tuple(range(n))))
        for w in range(1,n+1):
            H=(1<<w)-1;outside=((1<<n)-1)^H
            frequencies=Counter()
            for pairs in all_pairs:
                external=0
                for a,b in pairs:
                    if ((H>>a)&1)!=((H>>b)&1):external|=(1<<a)|(1<<b)
                frequencies[external&outside]+=1
            rho=F(max(frequencies.values())*2**(n-w),len(all_pairs))
            assert rho>=1
            for births in range(n-w+1):
                birthmask=((1<<births)-1)<<w;U=H|birthmask;allowed=outside&~birthmask
                counts=Counter()
                for mask,count in frequencies.items():counts[mask&allowed]+=count
                coordinates=[i for i in range(n) if allowed>>i&1]
                coin_masks=[sum(1<<i for j,i in enumerate(coordinates) if x>>j&1)
                            for x in range(1<<len(coordinates))]
                for mask in coin_masks:
                    assert F(counts[mask],len(all_pairs))<=rho/F(2**len(coordinates))
                    domination+=1
                for cost in range(4):
                    def integrand(mask):
                        cover=U|mask
                        dimension=sum(max(0,((cover>>(2*i))&3).bit_count()-1) for i in range(n//2))
                        return min(F(1),F(4)**(dimension-mask.bit_count()-cost))
                    actual=sum((F(count,len(all_pairs))*integrand(mask) for mask,count in counts.items()),F(0))
                    independent=sum((integrand(mask) for mask in coin_masks),F(0))/len(coin_masks)
                    assert actual<=rho*independent;capped+=1
    for R in range(2,7):
        for parity in range(1,R):
            for b in range(R+1):
                for j in range(R-b+1):
                    assert max(b+j-parity,0)-j==b-parity+max(parity-b-j,0)
                    identities+=1
    conditioning=0
    for sizes in ((0,1,3),(1,2,3),(2,2,2),(0,0,4),(4,4,4)):
        columns=[];probability=F(1)
        for a in sizes:
            p=F(a,4);probability*=math.comb(4,a)*p**a*(1-p)**(4-a)
            law={mask:F(math.comb(mask.bit_count(),a),math.comb(4,a)*2**(4-a))
                 for mask in range(16) if mask.bit_count()>=a}
            assert sum(law.values())==1;columns.append(law)
        coverage=(1+F(sum(sizes),12))/2
        for parity in (1,2,3):
            for r in (F(3,2),F(2)):
                actual=F(0)
                for masks in product(*(list(law) for law in columns)):
                    counts=[sum((mask>>i)&1 for mask in masks) for i in range(4)]
                    mass=math.prod(law[mask] for law,mask in zip(columns,masks))
                    actual+=mass*r**sum(max(0,parity-b) for b in counts)
                moment=sum((math.comb(3,b)*coverage**b*(1-coverage)**(3-b)
                            *r**max(0,parity-b) for b in range(4)),F(0))
                assert actual<=moment**4/probability
                conditioning+=1
    assert time.monotonic()-began<30
    print(dict(endpoint_domination_checks=domination,capped_average_checks=capped,
               deficiency_identities=identities,route_coverage_checks=conditioning,
               seconds=time.monotonic()-began))


if __name__=='__main__':run()
