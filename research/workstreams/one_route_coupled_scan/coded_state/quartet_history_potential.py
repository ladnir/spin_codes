"""Exact toy audit of count-conditioned potentials and full-update memory.

No full-width certificate. Serial,90s watchdog,<=50000 states,no files.
Toy lambda=1,tau=5/2,b=4 gives h=(1,1,4,16,32).
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
from quartet_window_transfer import kernel,gate_counts
from quartet_regional_transfer import run_exponent


def compositions(total,parts):
    if parts==1:
        yield (total,);return
    for first in range(total+1):
        for rest in compositions(total-first,parts-1):yield (first,)+rest


def potential(counts):
    return math.prod(h**n for h,n in zip((1,1,4,16,32),counts))


def feasible(counts,active,full_run):
    width=sum(counts)
    if active<width:return full_run==0 and counts[0]==width-active
    return 1<=full_run<=4 and not any(counts[:full_run])


def next_label(width,old,l):
    return (l,min(4,old[1]+1)) if l==width else (l,0)


@lru_cache(maxsize=128)
def gammas(width,active,full_run):
    maxima=[F(0)]*(width+1);tested=0
    for counts in compositions(width,5):
        if not feasible(counts,active,full_run):continue
        sums=[F(0)]*(width+1)
        for (outgoing,l,exponent),weight in kernel(counts):
            assert feasible(outgoing,*next_label(width,(active,full_run),l))
            sums[l]+=weight*(1<<exponent)*potential(outgoing)
        for l,value in enumerate(sums):maxima[l]=max(maxima[l],value/potential(counts))
        tested+=1
    assert tested
    return tuple(maxima)


def count_sequence_upper(counts,width=4):
    label=(0,0);value=F(1)
    for l in counts:
        value*=gammas(width,*label)[l]
        label=next_label(width,label,l)
    return value*4**counts[-1]


def stationary_floor_tests():
    width=4;histograms=list(compositions(width,5));checked=0
    for h in ((1,1,4,16,32),(1,2,3,5,7),(1,1,1,1,1)):
        maxima=[F(0)]*(width+1)
        for counts in histograms:
            phi=math.prod(hi**n for hi,n in zip(h,counts));sums=[F(0)]*(width+1)
            for (outgoing,l,exponent),weight in kernel(counts):
                sums[l]+=weight*(1<<exponent)*math.prod(hi**n for hi,n in zip(h,outgoing))
            for l,value in enumerate(sums):maxima[l]=max(maxima[l],value/phi)
        assert all(value>=math.comb(width,l) for l,value in enumerate(maxima))
        assert maxima[-1]>=4**width
        checked+=len(maxima)
    print(dict(stationary_coefficient_floor_checks=checked,
               universal_proof_is_in_note=True),flush=True)


def coefficient_tests():
    h=(1,1,4,16,32);ds=(1,4,16,64,128);checked=0
    cases=[(4,a,0) for a in range(4)]+[(4,4,f) for f in range(1,5)]+[(8,4,0),(8,8,4)]
    for width,a,f in cases:
        permitted=range(1,5) if a<width else range(f,5)
        upper=[F(0)]*(width+1)
        for (l,k),multiplicity in gate_counts(width).items():
            candidates=[]
            for z in (F(1,128),F(1,32),F(1,8),F(1,2),F(2)):
                base=F(1+ds[0]*z,h[0])**(width-a)
                active=max(F(1+ds[i]*z,h[i]) for i in permitted)**a
                candidates.append(base*active/z**k)
            upper[l]+=F(multiplicity*h[1]**(l-k),math.comb(width,k))*min(candidates)
        exact=gammas(width,a,f)
        assert all(a<=b for a,b in zip(exact,upper));checked+=len(upper)
    print(dict(exact_type_separated_cauchy_coefficients_checked=checked,
               full_and_partial_history_labels_checked=True),flush=True)


def tests():
    width=4;steps=4;direct=defaultdict(int)
    for mask in range(1<<(width*steps)):
        counts=tuple(((mask>>(width*t))&15).bit_count() for t in range(steps))
        direct[counts]+=1<<run_exponent(mask,width,steps)
    for counts,value in direct.items():assert value<=count_sequence_upper(counts,width)
    print(dict(exact_partial_input_streams_checked=1<<(width*steps),
               exact_input_count_sequences_checked=len(direct),
               quota_counters_not_dropped=True),flush=True)
    # Every full/empty input sequence is exact, including a free final state.
    steps=10
    for bits in range(1<<steps):
        counts=tuple(width*((bits>>t)&1) for t in range(steps))
        mask=sum(15<<(width*t) for t,l in enumerate(counts) if l)
        assert count_sequence_upper(counts,width)==1<<run_exponent(mask,width,steps)
    print(dict(exact_full_empty_sequences_checked=1<<steps,
               free_terminal_included=True),flush=True)
    burst=(0,4,4,4,4,4,4,0)
    actual=count_sequence_upper(burst)
    assert actual==1<<20
    # Stationary factors at these h are gamma_0=1,gamma_4=256.
    # Remembering only the last active count removes the initial but not
    # the terminal run-history loss: it gives 256^5 instead of 256^6.
    previous_count_only=F(1)
    previous=0
    for l in burst:
        if previous<width:factor=gammas(width,previous,0)[l]
        else:factor=max(gammas(width,width,f)[l] for f in range(1,5))
        previous_count_only*=factor;previous=l
    assert previous_count_only==1<<40
    print(dict(burst=list(burst),exact_and_history_bound=str(actual),
               last_count_only_bound=str(previous_count_only),
               stationary_bound_before_terminal_factor=str(1<<48),
               full_width_distance_claim=False),flush=True)


def freshness_tests():
    h=(1,1,4,16,32);ds=(1,4,16,64,128);checked=0
    for width,a,previous_full in ((4,3,0),(8,4,0),(8,4,4),(8,6,4),(8,8,8)):
        for counts in compositions(width,5):
            if counts[0]!=width-a or counts[1]<a-previous_full:continue
            elementary=[F(1)]+[F(0)]*width;used=0
            for d,n in zip(ds,counts):
                for _ in range(n):
                    used+=1
                    for k in range(used,0,-1):elementary[k]+=d*elementary[k-1]
            phi=potential(counts)
            for k in range(width+1):
                upper=[]
                for z in (F(1,128),F(1,32),F(1,8),F(1,2),F(2)):
                    base=F(1+ds[0]*z,h[0])**(width-a)
                    fresh=F(1+ds[1]*z,h[1])**(a-previous_full)
                    old=max(F(1+ds[i]*z,h[i]) for i in range(1,5))**previous_full
                    upper.append(base*fresh*old/z**k)
                assert elementary[k]/phi<=min(upper);checked+=1
    print(dict(exact_freshness_constrained_coefficients_checked=checked,
               previous_full_quartet_count_retained=True,whole_code_distance=False),flush=True)


def fresh_count_tests():
    h=(1,1,4,16,32);ds=(1,4,16,64,128);checked=0
    for width in (4,8):
        # Audit the entire outgoing fresh-count coefficient, not just its sum.
        for counts in compositions(width,5):
            a=width-counts[0];r=counts[1];rlo=2*(r//2);rhi=min(a,rlo+1)
            upper=defaultdict(F)
            for (l,k),multiplicity in gate_counts(width).items():
                for v0 in range(max(0,k-a),min(k,width-a)+1):
                    candidates=[]
                    for z in (F(1,128),F(1,32),F(1,8),F(1,2),F(2)):
                        fresh=F(1+ds[1]*z,h[1])
                        old=max(F(1+ds[i]*z,h[i]) for i in (2,3,4))
                        candidates.append(max(fresh**rr*old**(a-rr) for rr in (rlo,rhi))/z**(k-v0))
                    factor=F(multiplicity*math.comb(width-a,v0)*ds[0]**v0*h[1]**(l-k),math.comb(width,k))
                    upper[l,l-k+v0]+=factor*min(candidates)
            actual=defaultdict(F)
            for (outgoing,l,exponent),weight in kernel(counts):
                actual[l,outgoing[1]]+=weight*(1<<exponent)*potential(outgoing)/potential(counts)
            for key,value in actual.items():
                assert value<=upper[key];checked+=1
    print(dict(exact_outgoing_fresh_count_coefficients_checked=checked,
               all_histograms_width4_and8=True,whole_code_distance=False),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        stationary_floor_tests();coefficient_tests();tests();freshness_tests();fresh_count_tests()
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
