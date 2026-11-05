"""Exact small audit of a quota-retaining five-type temporal transfer.

No full-width transfer or distance certificate. Serial,90s,<64MiB,no files.
Toy moment weights use lambda=1,tau=5/2,b=4, so all scores are dyadic.
"""
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time


@lru_cache(maxsize=4)
def gate_counts(width):
    current={(0,0):1}
    for _ in range(width//4):
        after=defaultdict(int)
        for (active,full),count in current.items():
            for k in range(5):after[active+k,full+(4 if k==4 else 0)]+=count*math.comb(4,k)
        current=dict(after)
    for active in range(width+1):assert sum(v for (a,k),v in current.items() if a==active)==math.comb(width,active)
    return current


def subcounts(counts,total):
    if len(counts)==1:
        if 0<=total<=counts[0]:yield (total,),math.comb(counts[0],total)
        return
    for first in range(max(0,total-sum(counts[1:])),min(counts[0],total)+1):
        for rest,weight in subcounts(counts[1:],total-first):
            yield (first,)+rest,math.comb(counts[0],first)*weight


@lru_cache(maxsize=600)
def kernel(counts):
    width=sum(counts);result=defaultdict(F)
    for (active,full),multiplicity in gate_counts(width).items():
        denominator=math.comb(width,full)
        for selected,ways in subcounts(counts,full):
            outgoing=(width-active,active-full+selected[0],selected[1],selected[2],selected[3]+selected[4])
            assert min(outgoing)>=0 and sum(outgoing)==width
            exponent=selected[3]+2*selected[4]
            result[outgoing,active,exponent]+=F(multiplicity*ways,denominator)
    assert sum(result.values())==1<<width
    return tuple((key,value) for key,value in result.items())


def coefficients(width,steps,weighted=True):
    state={((width-a,a,0,0,0),a):F(math.comb(width,a)) for a in range(width+1)}
    largest=len(state)
    for step in range(1,steps):
        after=defaultdict(F)
        for (counts,active),value in state.items():
            for (outgoing,new_active,exponent),weight in kernel(counts):
                after[outgoing,active+new_active]+=value*weight*(1<<exponent if weighted else 1)
        state=dict(after);largest=max(largest,len(state))
        assert largest<50000
    result=[F(0)]*(width*steps+1)
    for (counts,active),value in state.items():
        exponent=2*counts[1]+4*counts[2]+6*counts[3]+7*counts[4]
        result[active]+=value*(1<<exponent if weighted else 1)
    if not weighted:assert result==[math.comb(width*steps,a) for a in range(width*steps+1)]
    print(dict(transfer_width=width,steps=steps,weighted=weighted,maximum_exact_states=largest),flush=True)
    return result


def partitions8():
    for rest in itertools.combinations(range(1,8),3):
        first=(0,)+rest;yield first,tuple(j for j in range(8) if j not in first)


def kernel_tests():
    checked=0
    for counts in ((2,2,2,1,1),(0,0,0,0,8),(8,0,0,0,0)):
        types=tuple(i for i,n in enumerate(counts) for _ in range(n));direct=defaultdict(F)
        for parts in partitions8():
            for mask in range(256):
                full=0
                for part in parts:
                    if all((mask>>j)&1 for j in part):full|=sum(1<<j for j in part)
                outgoing=[0]*5;exponent=0
                for j,old in enumerate(types):
                    if (full>>j)&1:
                        new=1 if old==0 else min(4,old+1)
                        exponent+=(1 if old==3 else 2 if old==4 else 0)
                    else:new=1 if (mask>>j)&1 else 0
                    outgoing[new]+=1
                direct[tuple(outgoing),mask.bit_count(),exponent]+=F(1,35);checked+=1
        assert direct==dict(kernel(counts))
    print(dict(exact_partition_input_kernel_cases=checked),flush=True)


def temporal_tests():
    width=4;steps=4;direct=[0]*(width*steps+1)
    for mask in range(1<<(width*steps)):
        inputs=[(mask>>(width*t))&15 for t in range(steps)];score=0
        for j in range(width):
            length=0
            for t in range(steps):
                free=((inputs[t]>>j)&1) and (t+1==steps or inputs[t+1]==15)
                if free:length+=1
                else:score+=max(0,2*length-5);length=0
            score+=2*length
        direct[mask.bit_count()]+=1<<score
    actual=coefficients(width,steps)
    assert actual==direct
    coefficients(width,steps,False)
    print(dict(exact_temporal_input_masks=1<<(width*steps),free_terminal_checked=True,
               exact_fixed_count_coefficients_checked=len(actual)),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        kernel_tests();temporal_tests()
        values=coefficients(8,3);coefficients(8,3,False)
        for active in (12,16,20,24):
            moment=values[active]/math.comb(24,active)
            assert 1<=moment<=4**active
            print(dict(width=8,steps=3,active_inputs=active,exact_fixed_count_moment=str(moment),
                       temporal_overlap_retained=True,uniform_quartet_partitions=True,
                       toy_weights_only=True,whole_code_distance=False),flush=True)
        assert values[-1]==4**24
        longer=coefficients(8,5);coefficients(8,5,False)
        print(dict(width=8,steps=5,active_inputs=28,
                   exact_fixed_count_moment=str(longer[28]/math.comb(40,28)),
                   multiple_overlapping_starts=True,toy_weights_only=True,whole_code_distance=False),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
