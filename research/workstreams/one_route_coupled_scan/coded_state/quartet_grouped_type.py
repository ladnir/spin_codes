"""Directed subcount: fix incoming types, sum every outgoing assignment.

Serial,90s watchdog,<64MiB numerical storage,>=20GiB free,no files.
Partial-state regional continuation audit, not whole-code distance.
"""
from fractions import Fraction as F
import math
import os
from pathlib import Path
import shutil
import threading
import time
import argparse
import itertools
import numpy as np
from quartet_two_update import log_upper
from quartet_run_components import marked_local
from quartet_marked_type import propose,optimize_histogram

FROZEN={
    131: ([[0,0,10],[0,1,21],[0,2,10],[0,3,1],[1,0,2],[2,0,21],[3,0,40],[4,0,23]],
          ['356968/867687','61454/291555','72933/827572','187475/88384']),
    132: ([[0,0,10],[0,1,19],[0,2,11],[0,3,1],[1,0,3],[2,0,21],[3,0,41],[4,0,22]],
          ['28070/67783','206203/977969','82958/946763','1653067/704607'])}


def polynomials():
    groups={}
    for (a,e,v,d,l,c),multiplicity in marked_local(True,True).items():
        factor=math.comb(4,a)*math.comb(4-a,e);assert multiplicity%factor==0
        assert (a==0 and c==0 and d==0) or (a>0 and c==l)
        key=(v-a+d,d,c,l-c);group=groups.setdefault((a,e),{})
        group[key]=group.get(key,0)+multiplicity//factor
    return groups


def incoming_counts(labels,coverage=27):
    types=marked_local(True,True);keys=list(types)
    ns=optimize_histogram(types,keys,propose(types,labels,coverage));counts={}
    for key,n in zip(keys,ns):counts[key[:2]]=counts.get(key[:2],0)+n
    return {key:n for key,n in counts.items() if n}


def tilt(groups,counts,target):
    packed=[]
    for key,n in counts.items():
        features=np.array(list(groups[key]),float);logs=np.log(np.array(list(groups[key].values()),float))
        packed.append((n,features,logs))
    target=np.array(target,float);parameters=np.zeros(4)
    def evaluate(parameters):
        objective=-parameters@target;gradient=-target.copy();covariance=np.zeros((4,4))
        for n,features,logs in packed:
            raw=logs+features@parameters;maximum=max(raw);p=np.exp(raw-maximum);z=sum(p);p/=z
            mean=p@features;centered=features-mean
            objective+=n*(maximum+math.log(z));gradient+=n*mean;covariance+=n*(centered.T*p)@centered
        return objective,gradient,covariance
    for _ in range(160):
        objective,gradient,covariance=evaluate(parameters)
        if max(abs(gradient))<1e-9:break
        step=np.linalg.solve(covariance+np.eye(4)*1e-12,gradient);step/=max(1.,max(abs(step))/2.)
        scale=1.
        for _ in range(30):
            candidate=parameters-scale*step
            if evaluate(candidate)[0]<=objective+1e-12:break
            scale/=2
        parameters=candidate
    assert max(abs(gradient))<1e-5
    return tuple(F(math.exp(value)).limit_denominator(1<<20) for value in parameters)


def directed_dp(factors,shape):
    coefficient_bits=24;state_bits=38;scale=1<<coefficient_bits
    lower=np.zeros(shape,dtype=np.uint64);upper=np.zeros_like(lower)
    lower.flat[0]=upper.flat[0]=1<<state_bits
    for terms,repeats in factors:
        assert sum(terms.values())==1
        probabilities=[]
        for powers,probability in terms.items():
            lo=probability.numerator*scale//probability.denominator
            hi=-((-probability.numerator*scale)//probability.denominator)
            probabilities.append((powers,lo,hi))
        for _ in range(repeats):
            lo_next=np.zeros_like(lower);hi_next=np.zeros_like(upper)
            assert int(upper.max())<1<<39
            for powers,lo,hi in probabilities:
                if any(power>=bound for power,bound in zip(powers,shape)):continue
                source=tuple(slice(0,bound-power) for power,bound in zip(powers,shape))
                dest=tuple(slice(power,bound) for power,bound in zip(powers,shape))
                lo_next[dest]+=(lower[source]*np.uint64(lo))>>np.uint64(coefficient_bits)
                hi_next[dest]+=(upper[source]*np.uint64(hi)+np.uint64(scale-1))>>np.uint64(coefficient_bits)
            lower,upper=lo_next,hi_next
    return lower,upper,state_bits


def step(labels,coverage=27,proposal=False):
    groups=polynomials()
    if proposal:counts=incoming_counts(labels,coverage)
    else:
        assert coverage==27
        counts={(a,e):n for a,e,n in FROZEN[labels][0]}
    assert sum(counts.values())==128
    assert sum(a*n for (a,e),n in counts.items())==256
    assert sum(e*n for (a,e),n in counts.items())==44
    target=(44,44,coverage,labels-coverage)
    markers=tilt(groups,counts,target) if proposal else tuple(map(F,FROZEN[labels][1]))
    factors=[[],[]];normalizers=[]
    for key,n in counts.items():
        weights={powers:F(m)*math.prod(marker**power for marker,power in zip(markers,powers)) for powers,m in groups[key].items()}
        z=sum(weights.values());normalizers.append((n,z));terms={}
        for powers,weight in weights.items():
            index=powers[:3] if key[0] else (powers[0],powers[3])
            terms[index]=terms.get(index,F(0))+weight/z
        factors[0 if key[0] else 1].append((terms,n))
    loA,hiA,bits=directed_dp(factors[0],(45,45,coverage+1))
    loE,hiE,_=directed_dp(factors[1],(45,labels-coverage+1))
    low=sum(int(loA[g,44,coverage])*int(loE[44-g,labels-coverage]) for g in range(45))
    high=sum(int(hiA[g,44,coverage])*int(hiE[44-g,labels-coverage]) for g in range(45))
    assert low>0
    probability_lower=F(low,1<<(2*bits));probability_upper=F(high,1<<(2*bits))
    factor=math.factorial(128)//math.prod(math.factorial(n) for n in counts.values())
    factor*=math.prod((math.comb(4,a)*math.comb(4-a,e))**n for (a,e),n in counts.items())
    ingredients=[(1,F(factor,math.comb(512,256)*math.comb(256,44))),*normalizers]
    ingredients.extend((-power,marker) for power,marker in zip(target,markers))
    lower=-log_upper(-1,probability_lower)-sum(log_upper(-power,value) for power,value in ingredients)
    upper=log_upper(1,probability_upper)+sum(log_upper(power,value) for power,value in ingredients)
    print(dict(labels=labels,coverage=coverage,directed_log2_step_subcount_lower=str(lower),
               directed_log2_step_subcount_upper=str(upper),lower_float=float(lower),interval_width=float(upper-lower),
               frozen_incoming_types=[[*key,n] for key,n in counts.items()],rational_markers=list(map(str,markers))),flush=True)
    return lower,upper


def tests():
    checks=0
    for dimension in (2,3):
        zero=(0,)*dimension;first=(1,)+(0,)*(dimension-1);last=(0,)*(dimension-1)+(1,)
        terms={zero:F(1,7),first:F(2,7),last:F(4,7)}
        exact={zero:F(1)}
        for _ in range(4):
            new={}
            for old,value in exact.items():
                for delta,probability in terms.items():
                    key=tuple(a+b for a,b in zip(old,delta));new[key]=new.get(key,F(0))+value*probability
            exact=new
        lo,hi,bits=directed_dp([(terms,4)],(5,)*dimension)
        for key in np.ndindex((5,)*dimension):
            assert F(int(lo[key]),1<<bits)<=exact.get(key,F(0))<=F(int(hi[key]),1<<bits);checks+=1
    def accepted(before,marked,after,target):
        assert not before&marked
        a=before.bit_count();e=marked.bit_count()
        if e and (a or not(target&~after)):return False
        if not a and target&after:return False
        if a and 100*(target&~after).bit_count()>=25*a+74*(before&after).bit_count():return False
        return 100*(target&~after).bit_count()>=25*e+74*(marked&after).bit_count()
    bridges=[(0,0,0,t) for t in range(16)]
    bridges.extend((u,0,v,0) for u,v in itertools.product(range(16),repeat=2))
    bridges.extend((0,e,0,15) for e in range(16));bridges.append((15,0,15,0))
    bridges.extend((u,0,15,t) for u in range(16) if u.bit_count()==2 for t in range(16))
    bridges.extend((15,0,15,t) for t in range(16))
    assert all(accepted(*case) for case in bridges)
    segregated=balanced=0;s={0,1,2,3}
    for rest in itertools.combinations(range(1,8),3):
        group={0,*rest};parts=(group,set(range(8))-group)
        segregated+=all(len(part&s) in (0,4) for part in parts)
        balanced+=all(len(part&s)==2 for part in parts)
    assert F(segregated,35)==F(math.comb(2,1),math.comb(8,4))
    assert F(balanced,35)==F(6**2,math.comb(8,4))
    print(dict(exact_fixed_point_coefficient_checks=checks,local_bridge_filter_checks=len(bridges),
               exact_bridge_partition_checks=35),flush=True)


def bridge(region_lower):
    row_den=math.comb(8192,2100);state_den=math.comb(512,256)
    prefix=F(state_den*math.comb(256,44)**2*math.comb(7168,2100),row_den)
    suffix=F(math.comb(128,64)*6**128*math.comb(7680,1844),state_den**2*row_den)
    prefix_lower=-log_upper(-1,prefix);suffix_lower=-log_upper(-1,suffix)
    lower=prefix_lower+254*region_lower+suffix_lower
    A=512+254*4096+7936;D=44+254*704;C=254*432+1844
    score=F(A,4)-C-F(37,50)*D
    assert score==F(91398,5) and score>18171
    burst_bound=2*(16*128-(2100//512)*128-min(2100%512,128))
    assert burst_bound==2968 and burst_bound<3702
    print(dict(prefix_log2_count_lower=str(prefix_lower),suffix_log2_probability_lower=str(suffix_lower),
               full_log2_expected_filtered_count_lower=str(lower),lower_float=float(lower),
               A=A,C=C,D=D,score=str(score),required_geometric_offset=18171,
               empty_initial_state=True,terminal_deaths=0,free_terminal_weight=512,
               every_regional_quota=True,all_local_filters=True,
               all_full_width_burst_scores_at_most=burst_bound,burst_good_for_counted_row_set=True,
               exceeds_remaining_field_turnover_cutoff=C+D>255463,
               probability_of_bad_geometry_not_lower_bounded=True,whole_code_distance=False),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--propose',action='store_true');args=parser.parse_args()
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests();left=step(131,proposal=args.propose);right=step(132,proposal=args.propose)
        lower=12*left[0]+4*right[0]-log_upper(1,math.comb(8192,2100))
        upper=12*left[1]+4*right[1]+log_upper(-1,math.comb(8192,2100))
        print(dict(directed_regional_log2_subcount_lower=str(lower),directed_regional_log2_subcount_upper=str(upper),
                   lower_float=float(lower),upper_float=float(upper),actual_regional_quota=True,
                   initial_terminal_not_imposed=True,whole_code_distance=False,files_written=0),flush=True)
        bridge(lower)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
