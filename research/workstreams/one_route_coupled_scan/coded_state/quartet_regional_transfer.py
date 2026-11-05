"""Exact toy audit retaining every regional quota and the total paid count.

Same onset-three weights as quartet_window_transfer; no large-code cert.
Serial,90s watchdog,at most50000 states,no files; numerical work on Peach.
"""
from collections import defaultdict
from fractions import Fraction as F
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
from quartet_window_transfer import kernel,coefficients,gate_counts


def terminal_exponent(counts):
    return 2*counts[1]+4*counts[2]+6*counts[3]+7*counts[4]


def regional_coefficients(width,steps_per_region,regions,inactive,weighted=True):
    m=width*steps_per_region
    assert width%4==0 and 0<=inactive<=m
    # (incoming five-type histogram, this-region active count, paid count)
    state={((width,0,0,0,0),0,0):F(1)};largest=1
    denominator=math.comb(m,inactive)
    for region in range(regions):
        for step in range(steps_per_region):
            after=defaultdict(F)
            remaining=width*(steps_per_region-step-1)
            for (counts,active,paid),value in state.items():
                for (outgoing,l,exponent),weight in kernel(counts):
                    a=active+l
                    if a+remaining<m-inactive:continue
                    contribution=value*weight*(1<<exponent if weighted else 1)
                    if remaining:
                        after[outgoing,a,paid]+=contribution
                    else:
                        c=a-(m-inactive)
                        assert 0<=c<=inactive
                        # For each allowed mask, choose which active positions
                        # were originally inactive and then paid. Normalize the
                        # original quota, not the allowed-mask count.
                        after[outgoing,0,paid+c]+=contribution*F(math.comb(a,c),denominator)
            state=dict(after);largest=max(largest,len(state))
            assert largest<50000
    result=[F(0)]*(regions*inactive+1)
    for (counts,active,paid),value in state.items():
        assert active==0
        result[paid]+=value*(1<<terminal_exponent(counts) if weighted else 1)
    if not weighted:assert result==[math.comb(regions*inactive,c) for c in range(regions*inactive+1)]
    print(dict(regional_width=width,steps_per_region=steps_per_region,regions=regions,
               inactive_per_region=inactive,weighted=weighted,maximum_exact_states=largest,
               state_carried_across_regions=True,whole_code_distance=False),flush=True)
    return result


def run_exponent(mask,width,steps):
    assert width==4  # The quartet partition is unique in this direct audit.
    inputs=[(mask>>(width*t))&15 for t in range(steps)];score=0
    for j in range(width):
        length=0
        for t in range(steps):
            free=((inputs[t]>>j)&1) and (t+1==steps or inputs[t+1]==15)
            if free:length+=1
            else:score+=max(0,2*length-5);length=0
        score+=2*length
    return score


def direct_tests():
    width=4;steps_per_region=2;regions=2;m=8;r=2;N=m*regions
    direct=[F(0)]*(r*regions+1);cases=0
    choices=[tuple(itertools.combinations(range(j*m,(j+1)*m),r)) for j in range(regions)]
    denominator=math.comb(m,r)**regions
    for parts in itertools.product(*choices):
        inactive=tuple(itertools.chain.from_iterable(parts))
        for c in range(len(inactive)+1):
            for paid in itertools.combinations(inactive,c):
                unpaid=set(inactive).difference(paid)
                mask=((1<<N)-1)^sum(1<<j for j in unpaid)
                direct[c]+=F(1<<run_exponent(mask,width,steps_per_region*regions),denominator)
                cases+=1
    actual=regional_coefficients(width,steps_per_region,regions,r)
    assert actual==direct
    regional_coefficients(width,steps_per_region,regions,r,False)
    one=regional_coefficients(width,steps_per_region,1,r)
    block=[sum(one[a]*one[c-a] for a in range(max(0,c-r),min(r,c)+1))
           for c in range(2*r+1)]
    assert all(exact<=upper for exact,upper in zip(actual,block))
    # A deterministic audit of the boundary relaxation, independent of
    # the original quota and the paid-set interpretation. Here scores
    # are measured in half-units: lambda=1, tau=5/2.
    maximum_gap=0
    for mask in range(1<<N):
        full=run_exponent(mask,width,steps_per_region*regions)
        split=sum(run_exponent((mask>>(j*m))&((1<<m)-1),width,steps_per_region)
                  for j in range(regions))
        assert 0<=split-full<=7*width*(regions-1)
        maximum_gap=max(maximum_gap,split-full)
    print(dict(exact_original_quota_and_paid_set_cases=cases,
               all_paid_coefficients_checked=True,free_terminal_checked=True,
               block_moment_coefficientwise_upper_checked=True,
               exact_boundary_masks_checked=1<<N,
               maximum_boundary_gap_in_half_units=maximum_gap),flush=True)


def potential_tests():
    """Exact coefficientwise spatial envelope; it does not drop quotas."""
    cases=0
    for counts in ((2,2,2,1,1),(0,0,0,0,8),(8,0,0,0,0),(1,1,1,2,3)):
        width=sum(counts)
        for h in ((1,2,3,5,7),(1,1,1,1,1),(3,2,4,2,3)):
            # At the toy weights the two transition rewards are 2 and 4.
            d=(h[1],h[2],h[3],2*h[4],4*h[4])
            elementary=[F(1)]+[F(0)]*width
            used=0
            for di,count in zip(d,counts):
                for _ in range(count):
                    used+=1
                    for k in range(used,0,-1):elementary[k]+=di*elementary[k-1]
            mean=F(sum(di*count for di,count in zip(d,counts)),width)
            for k in range(width+1):
                assert elementary[k]/math.comb(width,k)<=mean**k
            direct=[F(0)]*(width+1)
            for (outgoing,l,exponent),weight in kernel(counts):
                direct[l]+=weight*(1<<exponent)*math.prod(hi**n for hi,n in zip(h,outgoing))
            exact=[F(0)]*(width+1)
            upper=[F(0)]*(width+1)
            for (l,k),multiplicity in gate_counts(width).items():
                common=F(multiplicity)*h[0]**(width-l)*h[1]**(l-k)
                exact[l]+=common*elementary[k]/math.comb(width,k)
                upper[l]+=common*mean**k
            # Independently build ((h0+x*h1)^4-(x*h1)^4+(x*mean)^4)^(w/4).
            quartet=[F(math.comb(4,j)*h[0]**(4-j)*h[1]**j) for j in range(4)]+[mean**4]
            polynomial=[F(1)]
            for _ in range(width//4):
                after=[F(0)]*(len(polynomial)+4)
                for i,a in enumerate(polynomial):
                    for j,b in enumerate(quartet):after[i+j]+=a*b
                polynomial=after
            assert direct==exact and upper==polynomial
            assert all(a<=b for a,b in zip(exact,upper))
            cases+=width+1
    print(dict(exact_product_potential_coefficients_checked=cases,
               coefficientwise_maclaurin_envelope_checked=True,
               input_count_retained=True,whole_code_distance=False),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        direct_tests();potential_tests()
        # One region agrees with raw active-count coefficients, including
        # nontrivial random quartet partitions at width eight.
        raw=coefficients(8,2)
        actual=regional_coefficients(8,2,1,4)
        assert actual==[raw[12+c]*F(math.comb(12+c,c),math.comb(16,4)) for c in range(5)]
        actual=regional_coefficients(8,2,2,4)
        regional_coefficients(8,2,2,4,False)
        one=regional_coefficients(8,2,1,4)
        block=[sum(one[a]*one[c-a] for a in range(max(0,c-4),min(4,c)+1))
               for c in range(9)]
        assert all(exact<=upper for exact,upper in zip(actual,block))
        for c in (0,4,8):
            print(dict(paid=c,exact_paid_coefficient=str(actual[c]),
                       toy_lambda='1',toy_tau='5/2',toy_moment_base=4,
                       actual_regional_quotas=True,whole_code_distance=False),flush=True)
        assert actual[-1]==4**32
        # Sharp deterministic length of an all-active suffix after paying c.
        m=8192;r=2100;c=111572
        complete,remainder=divmod(c,r)
        longest=complete*m+(m-r+remainder)
        assert longest==440540 and longest<825856
        print(dict(maximum_all_active_suffix_positions=longest,
                   covered=c,inactive_per_region=r,sharp_geometric_bound=True),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
