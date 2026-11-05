"""Find and exactly verify a harmless type counted by the geometry moment.

Bounded proposal only: entropy tilts suggest integer counts, then exact integer
equalities and directed logs verify their contribution. No codeword claim.
Serial,60s watchdog,32MiB numeric arrays,20GiB free, no files.
"""
from fractions import Fraction as F
import math
import random
import itertools
import os
from pathlib import Path
import shutil
import threading
import time
import argparse
from quartet_two_update import log_upper

FROZEN_COUNTS=[[0,2,6,11,7,0,2,3,4,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,2,1,0,8,2,0,0,1,0,2,1,0,0,4,4,1,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,1,0,0,6,3,1,0,0,10,5,1,0,0,0,0,0,0,0,0,2,1,0,0,0,4,1,0,0,0,0,1,0,0,2,1,0,0,0,6,2,1,0,0,7,2,0,0,0],[0,2,7,11,6,0,2,4,4,0,1,1,0,0,0,0,1,0,1,0,2,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,3,2,0,6,3,0,0,0,0,1,1,0,0,4,2,1,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,2,0,0,6,3,1,0,0,12,4,1,0,0,0,0,1,0,0,0,3,1,0,0,0,4,1,0,0,0,0,0,0,0,1,1,0,0,0,6,3,0,0,0,7,2,1,0,0]]
FROZEN_MOVES=[[[43,54,44,53],[62,115,91,97],[55,102,73,92],[15,112,41,76],[47,80,49,75],[3,8,4,7],[18,110,20,105],[50,116,77,96],[16,111,40,81],[1,11,2,10],[67,106,87,98],[26,82,45,66]],[[33,115,67,93],[2,61,7,49],[18,117,54,76],[75,116,80,111],[3,53,8,43],[55,110,81,92],[40,97,41,96],[72,106,73,105],[1,20,6,16],[26,112,66,77],[27,102,50,87]]]


def records():
    result={}
    for a in range(5):
        before=(1<<a)-1
        for after in range(16):
            v=after.bit_count();h=(before&after).bit_count();d=a-h
            for target in range(16):
                if a==0 and target&after:continue
                if a and 100*(target&~after).bit_count()>=25*a+74*h:continue
                t=target.bit_count();c=t if a else (target&after).bit_count()
                key=(a,v,t,c,d)
                result[key]=result.get(key,0)+math.comb(4,a)
    return result


def propose(types):
    """Propose two integer types without an optional solver dependency."""
    import numpy as np
    keys=list(types);size=len(keys)
    features=np.array(keys,dtype=float)
    matrix=np.vstack((np.ones(size),features.T))
    logs=np.log(np.array(list(types.values()),dtype=float))
    rng=np.random.default_rng(20261027);counts=[]
    for block in range(2):
        target=np.array([256,256,131+block,33+block,44-block],dtype=float)
        parameters=np.zeros(5)
        for _ in range(200):
            probabilities=logs+features@parameters
            probabilities=np.exp(probabilities-max(probabilities));probabilities/=sum(probabilities)
            mean=probabilities@features;centered=features-mean
            step=np.linalg.solve((centered.T*probabilities)@centered+np.eye(5)*1e-12,mean-target/128)
            if np.max(abs(mean-target/128))<1e-10:break
            parameters-=step/max(1.,np.max(abs(step))/2.)
        assert np.max(abs(mean-target/128))<1e-7
        center=128*probabilities
        desired=np.r_[128,target].astype(np.int64)
        best=None;best_error=float('inf')
        for trial in range(10000):
            base=np.floor(center).astype(np.int64)
            if trial%4:base+=(rng.random(size)<(center-base)).astype(np.int64)
            basis=rng.choice(size,6,replace=False)
            try:change=np.linalg.solve(matrix[:,basis],desired-matrix@base)
            except np.linalg.LinAlgError:continue
            if np.max(abs(change))>1000 or np.max(abs(change-np.rint(change)))>1e-7:continue
            base[basis]+=np.rint(change).astype(np.int64)
            if np.min(base)<0 or not np.array_equal(matrix@base,desired):continue
            error=float(sum((base-center)**2/(center+.1)))
            if error<best_error:best=base.copy();best_error=error
        if best is None:raise RuntimeError('bounded integer proposal found no witness')
        counts.append(best)
        print(dict(block=block,proposal_error=best_error),flush=True)
    return keys,np.array(counts)




def disjoint_moves(types,keys,ns):
    """Sum independent histogram variations with identical five statistics."""
    pairs={}
    for i in range(len(keys)):
        for j in range(i+1,len(keys)):
            key=tuple(a+b for a,b in zip(keys[i],keys[j]))
            pairs.setdefault(key,[]).append((i,j))
    candidates=[]
    for choices in pairs.values():
        for left in range(len(choices)):
            i,j=choices[left]
            for k,l in choices[left+1:]:
                if len({i,j,k,l})<4:continue
                indices=(i,j,k,l);gain=F(0)
                for shift in range(-min(ns[i],ns[j]),min(ns[k],ns[l])+1):
                    ratio=F(1)
                    for index,delta in zip(indices,(shift,shift,-shift,-shift)):
                        ratio*=F(math.factorial(ns[index]),math.factorial(ns[index]+delta))*F(types[keys[index]])**delta
                    gain+=ratio
                if gain>1:candidates.append((float(gain),indices,gain))
    used=set();product=F(1);chosen=[]
    for _,indices,gain in sorted(candidates,reverse=True):
        if used.intersection(indices):continue
        used.update(indices);product*=gain;chosen.append(indices)
    return product,chosen


def verify_moves(types,keys,ns,moves):
    used=set();product=F(1)
    for indices in moves:
        assert len(set(indices))==4 and not used.intersection(indices)
        used.update(indices);i,j,k,l=indices
        assert tuple(a+b for a,b in zip(keys[i],keys[j]))==tuple(a+b for a,b in zip(keys[k],keys[l]))
        gain=F(0)
        for shift in range(-min(ns[i],ns[j]),min(ns[k],ns[l])+1):
            ratio=F(1)
            for index,delta in zip(indices,(shift,shift,-shift,-shift)):
                assert ns[index]+delta>=0
                ratio*=F(math.factorial(ns[index]),math.factorial(ns[index]+delta))*F(types[keys[index]])**delta
            gain+=ratio
        assert gain>=1;product*=gain
    return product


def verify(types,keys,counts,all_moves):
    coefficients=[];totals=[]
    for block in range(2):
        ns=[int(x) for x in counts[block]];assert min(ns)>=0 and sum(ns)==128
        sums=[sum(n*key[j] for n,key in zip(ns,keys)) for j in range(5)]
        a,v,t,c,d=sums;assert a==v==256 and t==131+block
        coefficient=math.factorial(128)//math.prod(math.factorial(n) for n in ns)
        coefficient*=math.prod(types[key]**n for n,key in zip(ns,keys))
        moves=all_moves[block];gain=verify_moves(types,keys,ns,moves)
        extended=coefficient*gain;assert extended.denominator==1
        coefficients.append(extended.numerator);totals.append(sums)
        print(dict(block=block,disjoint_histogram_move_count=len(moves),log2_gain_upper=str(log_upper(1,gain))),flush=True)
    combined=[12*x+4*y for x,y in zip(*totals)]
    a,v,t,c,d=combined
    assert a==v==4096 and t==2100 and 50*t-100*c-74*d==0
    assert F(1,4)*a-c-F(37,50)*d==-26
    numerator=coefficients[0]**12*coefficients[1]**4
    denominator=math.comb(512,256)**16*math.comb(8192,2100)
    lower=-log_upper(-1,numerator)-log_upper(1,denominator)
    upper=log_upper(1,numerator)+log_upper(-1,denominator)
    print(dict(exact_regional_type_counts=combined,exact_reference_excess_score=0,
               exact_weighted_cut_deficit=-26,
               directed_log2_expected_count_lower=str(lower),
               directed_log2_expected_count_upper=str(upper),lower_float=float(lower),
               harmless_histories=True,global_minimum_cut_not_asserted=True,
               all_actual_regional_quotas=True,whole_code_distance=False,files_written=0),flush=True)
    prefix=F(math.comb(512,256)*math.comb(7680,2100),math.comb(8192,2100))
    global_lower=255*lower-log_upper(-1,prefix)
    assert 64-255*26==-6566
    print(dict(full_history_log2_count_lower=str(global_lower),lower_float=float(global_lower),
               first_state_empty=True,terminal_state_weight=256,full_history_score=-6566,
               full_reference_excess=0,lower_log2_markov_expression_at_logb4=str(global_lower-4*14469)),flush=True)
    return lower,upper


def structural_tests():
    """Independent local enumeration and exact component-erasure toy checks."""
    explicit={}
    for before,after,target in itertools.product(range(16),repeat=3):
        a=before.bit_count();h=(before&after).bit_count()
        if a==0 and target&after:continue
        if a and 100*(target&~after).bit_count()>=25*a+74*h:continue
        key=(a,after.bit_count(),target.bit_count(),target.bit_count() if a else (target&after).bit_count(),a-h)
        explicit[key]=explicit.get(key,0)+1
    assert explicit==records()
    rng=random.Random(20261027);checked=0;positive_components=0;birth_gate_checks=0
    width=4;steps=3;vertices=width*steps
    for _ in range(16):
        inactive=[rng.randrange(1<<width) for _ in range(steps)]
        adjacency=[0]*vertices
        def join(mask):
            for i in range(vertices):
                if mask>>i&1:adjacency[i]|=mask^(1<<i)
        for t in range(steps):
            for j in range(width):
                if t:join((1<<(t*width+j))|(1<<((t-1)*width+j)))
                if inactive[t]>>j&1:join((1<<(t*width+j))|(((1<<width)-1)<<((t-1)*width) if t else 0))
        def score(mask):
            previous=0;c=d=0
            for t in range(steps):
                current=(mask>>(t*width))&((1<<width)-1)
                c+=(inactive[t]&(((1<<width)-1) if previous else current)).bit_count()
                d+=(previous&~current).bit_count();previous=current
            return 25*mask.bit_count()-100*c-74*d
        def components(mask):
            while mask:
                part=mask&-mask;pending=part
                while pending:
                    bit=pending&-pending;pending^=bit;i=bit.bit_length()-1
                    new=adjacency[i]&mask&~part;part|=new;pending|=new
                yield part;mask&=~part
        best=0;witness=0
        for mask in range(1<<vertices):
            value=score(mask)
            assert value==sum(score(part) for part in components(mask))
            holes=((1<<vertices)-1)^mask;previous=(1<<width)-1;uncovered=births=0
            for t in range(steps):
                current=(holes>>(t*width))&((1<<width)-1)
                if previous==(1<<width)-1:uncovered+=(current&inactive[t]).bit_count()
                if t:births+=(current&~previous).bit_count()
                previous=current
            full=25*vertices-100*sum(x.bit_count() for x in inactive)
            assert value==full+100*uncovered-25*holes.bit_count()-74*births
            if value>best or (value==best and mask.bit_count()<witness.bit_count()):best=value;witness=mask
            checked+=1
        for part in components(witness):
            assert score(part)>0;positive_components+=1
        for mask in range(1<<vertices):
            if score(mask)!=best:continue
            holes=((1<<vertices)-1)^mask
            for t in range(1,steps):
                current=(holes>>(t*width))&((1<<width)-1)
                previous=(holes>>((t-1)*width))&((1<<width)-1)
                if not current&~previous:continue
                assert t+1<steps
                after=(holes>>((t+1)*width))&((1<<width)-1)
                assert current==(1<<width)-1 and after&inactive[t+1]
                birth_gate_checks+=1
        for _ in range(128):
            left=rng.randrange(1<<vertices);right=rng.randrange(1<<vertices)
            assert score(left)+score(right)<=score(left|right)+score(left&right)
    print(dict(local_mask_cases=4096,component_identities=checked,
               supermodular_checks=2048,positive_optimal_components=positive_components,
               optimal_hole_birth_gates=birth_gate_checks),flush=True)


def hole_gate_tests():
    """Exact mask DP checks optimal triples, including real internal births."""
    checks=births=0
    for inactive in ([0,0,0,0,15,15,0],[0,0,0,0,15,15,15,0],[1,3,7,2,15,15,0]):
        steps=len(inactive)
        def gain(t,previous,current):
            uncovered=(current&inactive[t]).bit_count() if previous==15 else 0
            return 100*uncovered-25*current.bit_count()-(74*(current&~previous).bit_count() if t else 0)
        forward=[[gain(0,15,v) for v in range(16)]]
        for t in range(1,steps):forward.append([max(forward[-1][u]+gain(t,u,v) for u in range(16)) for v in range(16)])
        suffix=[[0]*16 for _ in range(steps)]
        for t in range(steps-2,-1,-1):suffix[t]=[max(gain(t+1,u,v)+suffix[t+1][v] for v in range(16)) for u in range(16)]
        optimum=max(forward[-1])
        for t in range(1,steps-1):
            for before,current,after in itertools.product(range(16),repeat=3):
                if forward[t-1][before]+gain(t,before,current)+gain(t+1,current,after)+suffix[t+1][after]!=optimum:continue
                checks+=1
                if current&~before:
                    assert current==15 and after&inactive[t+1];births+=1
                    new=current&~before
                    assert 100*(after&inactive[t+1]).bit_count()>=25*new.bit_count()+74*(new&~after).bit_count()
        for before,current in itertools.product(range(16),repeat=2):
            if forward[-2][before]+gain(steps-1,before,current)==optimum:assert not current&~before
    assert births>0
    print(dict(optimal_triple_gate_checks=checks,nonvacuous_internal_births=births,
               terminal_births_excluded=True),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--propose',action='store_true');args=parser.parse_args()
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start()
    start=time.monotonic()
    try:
        types=records();keys=list(types);counts=FROZEN_COUNTS;moves=FROZEN_MOVES
        if args.propose:
            keys,counts=propose(types);moves=[disjoint_moves(types,keys,ns)[1] for ns in counts]
            print(dict(proposed_counts=[list(map(int,ns)) for ns in counts],proposed_moves=moves),flush=True)
        structural_tests();hole_gate_tests();verify(types,keys,counts,moves)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
