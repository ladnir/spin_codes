"""Directed routing bound for full-width bursts, not arbitrary supports.

Every state is either all512 positions or empty. Keep each complete
regional quota exactly. Serial,45s watchdog, tiny scalar arrays, >=20GiB
free, no files. Parameters and probabilities refer to geometric cuts.
"""
import argparse
from fractions import Fraction as F
from itertools import product
import math
import os
from pathlib import Path
import shutil
import threading
import time

from quartet_two_update import log_upper


def up(x,bits=96):
    """Outward dyadic rounding with a bounded significant numerator."""
    x=F(x);assert x>=0
    if not x:return x
    exponent=x.numerator.bit_length()-x.denominator.bit_length()
    unit=F(2)**(exponent-bits)
    ratio=x/unit
    return ((ratio.numerator+ratio.denominator-1)//ratio.denominator)*unit


def power_up(x,n):
    result=F(1)
    while n:
        if n&1:result=up(result*x)
        n//=2
        if n:x=up(x*x)
    return result


def inverse_base(numerator,bits=80):
    """Upper bound on 2^(-numerator/128), checked by integer powers."""
    assert 0<numerator<128*bits
    target=1<<(128*bits-numerator);lo=0;hi=1<<bits
    while hi-lo>1:
        middle=(lo+hi)//2
        if middle**128>=target:hi=middle
        else:lo=middle
    assert hi**128>=target and (hi-1)**128<target
    return F(hi,1<<bits)


def hypergeom_moment(size,inactive,covered,u):
    """E[u^C] for the exact uniform inactive subset, rounded upward."""
    assert 0<=inactive<=size and 0<=covered<=size and 0<u<=1
    low=max(0,inactive-(size-covered));high=min(inactive,covered)
    term=up(F(math.comb(covered,low)*math.comb(size-covered,inactive-low),math.comb(size,inactive)))
    term=up(term*power_up(u,low));total=term
    for c in range(low,high):
        ratio=F((covered-c)*(inactive-c),(c+1)*(size-covered-inactive+c+1))*u
        term=up(term*ratio);total=up(total+term)
    return min(F(1),total)


def region_matrix(steps,activity_reward,death_reward,moments):
    """Keep covered-update count until multiplying its exact quota moment."""
    result=[[F(0),F(0)] for _ in range(2)]
    for incoming in range(2):
        state=[[F(0)]*(steps+1) for _ in range(2)];state[incoming][0]=1
        for t in range(steps):
            nxt=[[F(0)]*(steps+1) for _ in range(2)]
            for old in range(2):
                for new in range(2):
                    h=int(old or new)
                    factor=activity_reward**new*death_reward**int(old and not new)
                    for count in range(t+1):
                        nxt[new][count+h]+=state[old][count]*factor
            state=nxt
        for outgoing in range(2):
            result[incoming][outgoing]=up(sum(state[outgoing][h]*moments[h] for h in range(steps+1)))
    return result


def tests():
    checks=0;matrix_checks=0
    for size in (4,6,8):
        for r in range(1,size):
            moments=[]
            for covered in range(0,size+1,2):
                exact=sum(F(math.comb(covered,c)*math.comb(size-covered,r-c),math.comb(size,r))*F(1,3)**c
                          for c in range(max(0,r-size+covered),min(r,covered)+1))
                bound=hypergeom_moment(size,r,covered,F(1,3))
                assert exact<=bound<=exact*(1+F(1,1<<70));checks+=1
                moments.append(bound)
            steps=size//2;reward=F(2);death=F(1,4)
            matrix=region_matrix(steps,reward,death,moments)
            for incoming in range(2):
                exact=[F(0),F(0)]
                for suffix in product((0,1),repeat=steps):
                    path=(incoming,)+suffix
                    h=sum(int(u or v) for u,v in zip(path,path[1:]))
                    d=sum(int(u and not v) for u,v in zip(path,path[1:]))
                    exact[suffix[-1]]+=reward**sum(suffix)*death**d*moments[h]
                for outgoing in range(2):
                    assert exact[outgoing]<=matrix[incoming][outgoing]
                    assert matrix[incoming][outgoing]<=exact[outgoing]*(1+F(1,1<<70))
                    matrix_checks+=1
    for x in (F(1,3),F(137,17),F(1,1<<300),F(1<<300)):
        assert x<=up(x)<=x*(1+F(1,1<<90))
    defect_checks=0
    for weights in product(range(5),repeat=3):
        deaths=sum(max(0,u-v) for u,v in zip((0,)+weights,weights));best=None
        for bits in product((0,1),repeat=3):
            plus=sum(w for w,b in zip(weights,bits) if not b)
            minus=sum(4-w for w,b in zip(weights,bits) if b)
            starts=sum(4-weights[t] for t in range(3) if bits[t] and (not t or not bits[t-1]))
            ends=sum(4-weights[t] for t in range(2) if bits[t] and not bits[t+1])
            template_deaths=4*sum(bits[t] and not bits[t+1] for t in range(2))
            cost=F(1,4)*plus+F(3,4)*minus+F(4,5)*starts+F(1,5)*ends-F(37,50)*(deaths-template_deaths)
            best=cost if best is None else min(best,cost);defect_checks+=1
        assert defect_charge(list(weights),deaths,width=4)==best
    print(dict(exact_hypergeometric_checks=checks,exact_burst_matrix_checks=matrix_checks,
               exact_defect_template_candidates=defect_checks,
               initial_zero=True,free_terminal=True,arithmetic='directed rational',files_written=0),flush=True)


def defect_charge(weights,deaths,width=512,slope=F(1,4),tau=F(37,50),greedy_threshold=None):
    """Exact deterministic upper charge relative to full/empty templates.

    Return the minimum over templates, or the specified threshold template.
    Counts describe a real support history only when supplied by its caller.
    No independence or support-weight distribution is assumed.
    """
    assert weights and all(0<=w<=width for w in weights) and deaths>=0
    previous=[F(0),None];old_weight=0
    for weight in weights:
        costs=[slope*weight,(1-slope)*(width-weight)]
        target=(0,1) if greedy_threshold is None else (int(weight>=greedy_threshold),)
        nxt=[None,None]
        for bit in target:
            choices=[]
            for old in range(2):
                if previous[old] is None:continue
                transition=F(0)
                if old==0 and bit==1:transition=F(4,5)*(width-weight)
                if old==1 and bit==0:transition=tau*width+F(1,5)*(width-old_weight)
                choices.append(previous[old]+costs[bit]+transition)
            nxt[bit]=min(choices)
        previous=nxt;old_weight=weight
    return min(x for x in previous if x is not None)-tau*deaths


def transition_diagnostic():
    # Rounded weight traces from the failed geometry moment. These are
    # illustrations, not samples of code states or a full path distribution.
    down=[509,509,509,509,509,508,507,504,487,504,507,507,504,489,423,223,5]
    upward=[5,0,0,0,0,0,34,221,486,503,507,508,509,509,509,509,509]
    weights=(down[1:]+upward[1:])*128
    minimum_deaths=sum(max(0,a-b) for a,b in zip([0]+weights,weights))
    best=defect_charge(weights,minimum_deaths)
    greedy=defect_charge(weights,minimum_deaths,greedy_threshold=256)
    assert best<=greedy
    needed=max(F(0),(best-14469)/F(37,50))
    print(dict(illustrative_rounded_trace=True,actual_encoder_distribution=False,
               updates=len(weights),activity=sum(weights),minimum_possible_deaths=minimum_deaths,
               minimum_template_charge_at_minimum_deaths=str(best),
               majority_template_charge_at_minimum_deaths=str(greedy),
               extra_deaths_sufficient_for_14469_charge=str(needed),
               whole_code_distance=False,files_written=0),flush=True)


def certificate(numerator=250,rows=2100):
    # log2(b)=numerator/128. At lambda1/4 and tau37/50, complete
    # 512-position updates give the EXACT integer exponents below.
    assert numerator%25==0 and 25<=numerator<=750 and 1<=rows<8192
    u=inverse_base(numerator)
    moments=[hypergeom_moment(8192,rows,512*h,u) for h in range(17)]
    matrix=region_matrix(16,F(2)**numerator,F(2)**(-74*numerator//25),moments)
    # Work backward from the free terminal state. Round outward each step.
    vector=[F(1),F(1)]
    for _ in range(256):
        vector=[up(sum(matrix[i][j]*vector[j] for j in range(2))) for i in range(2)]
    union_log=log_upper(1,math.comb(8192,rows))+log_upper(1,vector[0])
    logbase=F(numerator,128)
    required=(union_log+60)/logbase
    offset=(required.numerator+required.denominator-1)//required.denominator
    tail=union_log-offset*logbase;assert tail<=-60
    result=dict(inactive_rows=rows,log2_base=str(logbase),slope='1/4',death_weight='37/50',
                all_full_or_empty_histories=True,all_row_sets_of_this_size=True,
                arbitrary_quartet_partitions=True,free_terminal=True,
                exact_regional_quota=True,offset=offset,log2_tail_upper=str(tail),
                union_log_upper=str(union_log),arithmetic='directed rational',
                arbitrary_state_supports=False,whole_code_distance=False,files_written=0)
    print(result,flush=True)
    return offset,tail


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(45,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        parser=argparse.ArgumentParser()
        parser.add_argument('--tests-only',action='store_true')
        parser.add_argument('--numerators',nargs='+',type=int,default=[250])
        args=parser.parse_args();tests()
        assert 1<=len(args.numerators)<=3
        if not args.tests_only:
            for numerator in args.numerators:certificate(numerator)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
