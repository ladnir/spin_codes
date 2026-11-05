"""Small exact tests of active-row geometry and conserved route moments.

Serial, <=30s, standard library only; no output files.
"""
from fractions import Fraction as F
from itertools import product
import math
import shutil
import tempfile
import time


def main():
    assert shutil.disk_usage(tempfile.gettempdir()).free >= 20*(1<<30)
    began=time.monotonic();identities=0;row_checks=0;route_checks=0
    for R in range(2,7):
        for parity in range(1,R):
            for b in range(R+1):
                for j in range(R-b+1):
                    if b+j>parity:
                        assert b+j-parity-j==b-parity
                        identities+=1
    # One row: sum over its U indicators and independent outside-U coins.
    ps=(F(1,4),F(1,2),F(3,4));R=len(ps);parity=1
    for r,t in ((F(2),F(3)),(F(4),F(1,2)),(F(3,2),F(2))):
        active=inactive=F(0)
        for states in product(range(3),repeat=R):
            # 0: neither U nor J; 1: U; 2: outside-U crossing J.
            mass=math.prod(p if v==1 else (1-p)/2 for p,v in zip(ps,states))
            b=states.count(1);j=states.count(2)
            active+=mass*(r/t)**b
            if not j:inactive+=mass*t**(-b)
        assert active==math.prod(1-p+p*r/t for p in ps)
        assert inactive==math.prod((1-p)/2+p/t for p in ps)
        mean=sum(ps)/R
        assert active <= (1-mean+mean*r/t)**R
        assert inactive <= ((1-mean)/2+mean/t)**R
        row_checks+=4
    # Actual uniform regional subsets, then independent crossing coins.
    M=2;full=(1<<M)-1
    for sizes in ((0,1,2),(1,1,1),(2,2,2),(0,0,1)):
        columns=[[u for u in range(1<<M) if u.bit_count()==a] for a in sizes]
        p=F(sum(sizes),M*R)
        conditioning=math.prod(math.comb(M,a)*F(a,M)**a*(1-F(a,M))**(M-a)
                               for a in sizes)
        assert conditioning>=F(1,(M+1)**R)
        for r,t in ((F(2),F(3)),(F(4),F(1,2)),(F(3,2),F(2))):
            total=F(0)
            for us in product(*columns):
                jchoices=[[j for j in range(1<<M) if not j&u] for u in us]
                mass=F(1,math.prod(map(len,columns))*math.prod(map(len,jchoices)))
                for js in product(*jchoices):
                    # Every exact active-row set; permit impossible rows as overcount.
                    terms=[]
                    for S in range(1<<M):
                        if any(j & (full^S) for j in js):continue
                        bs=sum((u&S).bit_count() for u in us)
                        terms.append(r**(bs-parity*S.bit_count()))
                    total+=mass*sum(terms,F(0))
            active=(1-p+p*r/t)**R/r**parity
            inactive=((1-p)/2+p/t)**R
            bound=t**sum(sizes)*(active+inactive)**M/conditioning
            assert total<=bound,(sizes,r,t,total,bound)
            route_checks+=1
    assert time.monotonic()-began<30
    print(dict(active_dimension_identities=identities,row_moment_checks=row_checks,
               conserved_route_checks=route_checks,seconds=time.monotonic()-began))


if __name__=='__main__':main()
