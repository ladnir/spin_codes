"""Exact local automaton for canonical support and used payments.

Six run states, twelve live pending-input types. Exact quartet/local and
two-region width4 tests. Serial Peach90s,<256MiB,no generated files.
"""
from collections import defaultdict
from functools import lru_cache
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
from quartet_closed_paid_sets import tests as closed_tests


def advance(age,z,selected):
    # 0: zero;1,2: unselected run ages;3,4,5: selected run ages1,2,>=3.
    if not z:return 0 if not selected and age in (0,1,2,5) else -1
    if not selected:return {0:1,1:2}.get(age,-1)
    return {0:3,3:4,4:5,5:5}.get(age,-1)


def live(age,choice):
    # choice0=(U,S)=(0,0);1=(1,0);2=(1,1).
    if age in (1,2):return choice!=2
    if age in (3,4):return choice==2
    return True


TYPES=tuple(3*a+c for a in range(6) for c in range(3) if live(a,c))
assert len(TYPES)==12


def terminal(hist):
    return all(advance(t//3,t%3!=0,t%3==2) in (0,3,4,5) for t in hist)


@lru_cache(maxsize=1500)
def local(hist):
    assert len(hist)==4 and all(t in TYPES for t in hist)
    old=[(t//3,t%3!=0,t%3==2) for t in hist]
    touched=any(s for a,u,s in old);out=defaultdict(int)
    for choices in itertools.product(range(3),repeat=4):
        full=all(c!=0 for c in choices)
        ages=[advance(a,u and full,s) for a,u,s in old]
        if any(a<0 or not live(a,c) for a,c in zip(ages,choices)):continue
        cover=[touched or c==2 for c in choices]
        assert all(not v or c!=0 for v,c in zip(cover,choices))
        destination=tuple(sorted(3*a+c for a,c in zip(ages,choices)))
        unpaid=sum(c==0 for c in choices);covered=sum(cover)
        activity=sum(c==2 for c in choices)
        deaths=sum(s and c!=2 for (a,u,s),c in zip(old,choices))
        out[destination,unpaid,covered,activity,deaths]+=1
    return tuple(out.items())


def automaton_tests():
    checked=0
    for length in range(1,11):
        for zs in itertools.product((0,1),repeat=length):
            paths={(0,()):1}
            for z in zs:
                after={}
                for (a,support),multiplicity in paths.items():
                    for s in (0,1):
                        b=advance(a,z,s)
                        if b>=0:after[b,support+(s,)]=multiplicity
                paths=after
            accepted=[s for (a,s) in paths if a in (0,3,4,5)]
            expected=[0]*length;t=0
            while t<length:
                if not zs[t]:t+=1;continue
                end=t
                while end<length and zs[end]:end+=1
                if end==length or end-t>=3:expected[t:end]=[1]*(end-t)
                t=end
            assert accepted==[tuple(expected)];checked+=1
    print(dict(exact_unique_automaton_streams=checked,live_pending_types=len(TYPES)),flush=True)


def regional_test():
    expected=closed_tests();state={((0,0,0,0),0,0,0):1};peak=1
    for step in range(4):
        after=defaultdict(int)
        for (hist,running,paid,score),multiplicity in state.items():
            for (dest,unpaid,covered,activity,deaths),ways in local(hist):
                for c in range(min(covered,2-running-unpaid)+1):
                    count=running+unpaid+c
                    if step%2==1 and count!=2:continue
                    after[dest,(0 if step%2 else count),paid+c,score+2*activity-5*deaths]+=multiplicity*ways*math.comb(covered,c)
        state=dict(after);peak=max(peak,len(state));assert peak<250000
    result=defaultdict(int)
    for (hist,running,paid,score),ways in state.items():
        assert running==0
        if terminal(hist):result[paid,score]+=ways
    assert dict(result)==expected
    assert sum(v for (paid,score),v in result.items() if paid==0)==math.comb(8,2)**2
    print(dict(exact_regional_transfer_coefficients=len(result),maximum_states=peak,
               quartet_kernels=local.cache_info().currsize,terminal_checked_once=True),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        automaton_tests();regional_test();print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
