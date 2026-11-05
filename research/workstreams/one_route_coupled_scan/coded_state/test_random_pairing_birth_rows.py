"""Exact birth-row moments and a tiny uniform-image forward transfer.

The uniform-image transfer is a comparison model, not the accepted SL2 law.
Serial <=30s, no output files.
"""
from fractions import Fraction as F
from itertools import product
import math
import shutil
import tempfile
import time


def matchings(vertices):
    if not vertices:
        yield ();return
    a=vertices[0]
    for i in range(1,len(vertices)):
        for tail in matchings(vertices[1:i]+vertices[i+1:]):
            yield ((a,vertices[i]),)+tail


def main():
    assert shutil.disk_usage(tempfile.gettempdir()).free>=20*(1<<30),'disk guard'
    began=time.monotonic();row_checks=0;forward_checks=0
    regions=((F(1,2),F(1,8)),(F(3,4),F(1,4)),(F(1,4),F(0)))
    for r,t,v in ((F(2),F(3),F(1,2)),(F(4),F(1,2),F(3)),(F(3,2),F(2),F(1))):
        active=inactive=F(0)
        # Types: ordinary U, birth without/with crossing, outside without/with crossing.
        for states in product(range(5),repeat=3):
            mass=F(1);ucount=births=0;allowed=True
            for state,(p,beta) in zip(states,regions):
                mass*=p-beta if state==0 else beta/2 if state<3 else (1-p)/2
                ucount+=state<3;births+=state in (1,2)
                allowed &= state not in (1,4)
            active+=mass*(r/t)**ucount*v**(-births)
            if allowed:inactive+=mass*t**(-ucount)*v**(-births)
        aa=math.prod((p-b)*r/t+b*r/(t*v)+1-p for p,b in regions)
        bb=math.prod((p-b)/t+b/(2*t*v)+(1-p)/2 for p,b in regions)
        assert active==aa and inactive==bb
        p=sum(p for p,b in regions)/3;b=sum(b for p,b in regions)/3
        assert aa<=((p-b)*r/t+b*r/(t*v)+1-p)**3
        assert bb<=((p-b)/t+b/(2*t*v)+(1-p)/2)**3
        row_checks+=4
    n=4;pairs=list(matchings(tuple(range(n))))
    a0=F(3,4);a1=F(5,3);y=F(1,3);z=F(2);x=F(3,2)
    for q in (2,4):
        for w in range(n+1):
            H=(1<<w)-1;direct=F(0);formula=F(0)
            for matching in pairs:
                live=0
                for a,b in matching:
                    if H & ((1<<a)|(1<<b)):live|=(1<<a)|(1<<b)
                live_indices=[j for j in range(n) if live>>j&1]
                l=len(live_indices);k=l-w
                formula+=((a0+a1)/q)**l*(x*y+(1+(q-2)*y)*z)**w \
                    *(1+(q-1)*y*z)**k*(a0+a1*y*z)**(n-l)
                for X in range(1<<n):
                    source=a0**(n-X.bit_count())*a1**X.bit_count()/q**l
                    for values in product(range(q),repeat=l):
                        image=dict(zip(live_indices,values));next_weight=out_weight=deaths=0
                        for j in range(n):
                            next_value=image.get(j,0)^((X>>j)&1)
                            next_weight+=next_value!=0
                            out_weight+=(next_value^((H>>j)&1))!=0
                            deaths+=bool((H>>j)&1) and next_value==0
                        direct+=source*z**next_weight*y**out_weight*x**deaths
            assert direct==formula,(q,w,direct,formula)
            forward_checks+=1
            assert time.monotonic()-began<25
    print(dict(birth_moment_checks=row_checks,uniform_image_transfer_checks=forward_checks,
               seconds=time.monotonic()-began))


if __name__=='__main__':main()
