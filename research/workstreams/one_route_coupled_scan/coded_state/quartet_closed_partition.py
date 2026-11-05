"""Check population-conditioned quartet composition against labelled routing.

Width8 exact identities; serial Peach90s,<256MiB,no generated files.
"""
from collections import defaultdict
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
from quartet_closed_transfer import local,advance,live,terminal


def composed(hist):
    out=defaultdict(int)
    for first in itertools.combinations(range(8),4):
        second=tuple(j for j in range(8) if j not in first)
        for (d1,u1,v1,a1,e1),m1 in local(tuple(sorted(hist[j] for j in first))):
            for (d2,u2,v2,a2,e2),m2 in local(tuple(sorted(hist[j] for j in second))):
                key=tuple(sorted(d1+d2)),u1+u2,v1+v2,a1+a2,e1+e2
                out[key]+=m1*m2
    return dict(out)


def labelled(hist):
    out=defaultdict(int);old=[(t//3,t%3!=0,t%3==2) for t in hist]
    assignments=tuple(itertools.product(range(3),repeat=8))
    for rest in itertools.combinations(range(1,8),3):
        first=(0,)+rest;second=tuple(j for j in range(8) if j not in first)
        parts=(first,second);touched=[any(old[j][2] for j in part) for part in parts]
        for cs in assignments:
            ages=[-1]*8;cover=[False]*8
            for part,touch in zip(parts,touched):
                full=all(cs[j]!=0 for j in part)
                for j in part:
                    a,u,s=old[j];ages[j]=advance(a,u and full,s);cover[j]=touch or cs[j]==2
            if any(a<0 or not live(a,c) for a,c in zip(ages,cs)):continue
            assert all(not v or c!=0 for v,c in zip(cover,cs))
            dest=tuple(sorted(3*a+c for a,c in zip(ages,cs)))
            key=dest,sum(c==0 for c in cs),sum(cover),sum(c==2 for c in cs),sum(s and c!=2 for (a,u,s),c in zip(old,cs))
            out[key]+=1
    return dict(out)


def tests():
    coefficients=0
    for hist in ((0,)*8,(1,)*8,(2,)*4+(1,)*4,(0,1,2,4,7,11,14,17)):
        direct=labelled(hist);formula=composed(hist)
        assert formula=={k:2*v for k,v in direct.items()};coefficients+=len(direct)
        print(dict(incoming_types=hist,exact_partition_coefficients=len(direct),
                   ordered_groups=70,unlabelled_partitions=35),flush=True)
        if hist==(0,)*8:
            totals=defaultdict(int)
            for (dest,unpaid,covered,activity,deaths),multiplicity in formula.items():
                if not terminal(dest):continue
                assert deaths==0 and covered==activity
                for c in range(covered+1):totals[unpaid+c,c,activity]+=multiplicity*math.comb(covered,c)
            expected={(r,c,8-r+c):70*math.comb(8,r)*math.comb(r,c) for r in range(9) for c in range(r+1)}
            assert dict(totals)==expected
            print(dict(exact_one_step_free_terminal_coefficients=len(expected)),flush=True)
    print(dict(exact_population_conditioned_coefficients=coefficients),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:tests();print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
