"""Directed ensemble bound for all zero-coverage support components.

Actual regional quotas; arbitrary fixed quartet partitions; all row sets.
This is a restricted geometry certificate, not whole-code distance.
No files, 30s watchdog, <16MiB numeric storage, >=20GiB output-volume free.
"""
from fractions import Fraction as F
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
import argparse
from quartet_two_update import log_upper


def power_upper(exponent,bits=80):
    exponent=F(exponent);assert exponent>=0
    num,den=exponent.numerator,exponent.denominator
    target=1<<(num+bits*den)
    lo=0;hi=1<<(bits+(num+den-1)//den+1)
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid**den>=target:hi=mid
        else:lo=mid
    assert (hi-1)**den<target<=hi**den
    return F(hi,1<<bits)


def bound(slope=F(13,50),B=25604,r=2100,logb=F(5,8)):
    width=512;steps=4096;rows=8192;read=17
    tau=F(37,50);alpha=3*slope-tau
    assert 0<=alpha<=slope
    active=1-F(r,rows);count=width*(steps-4);boundary=4*slope*width
    small=power_upper(read*logb*alpha);large=power_upper(read*logb*slope)
    # Positive-coefficient form avoids reversing an interval by subtraction.
    moment=1-active**13+(active**13-active**17)*small+active**17*large
    upper=-logb*(B-boundary)+F(count,read)*log_upper(1,moment)+log_upper(1,math.comb(rows,r))
    print(dict(slope=str(slope),death_weight=str(tau),inactive_rows=r,
               offset=B,log2_base=str(logb),boundary=str(boundary),read_degree=read,
               directed_log2_union_tail_upper=str(upper),upper_float=float(upper),
               actual_regional_quotas=True,all_row_sets=True,arbitrary_quartet_partitions=True,
               zero_inactive_coverage_only=True,whole_code_distance=False,files_written=0),flush=True)
    return upper


def column_moment(p,alpha_weight,continuation_weight,steps=4096,onset=3):
    """Upper moment for a complete vertical sequence; exact hidden-bit law."""
    state_bits=80;coefficient_bits=64;scale=1<<coefficient_bits
    def ceil_scaled(value,bits):return -((-value.numerator*(1<<bits))//value.denominator)
    assert onset>=1
    stride=onset+1
    state=[0]*(2*stride);state[0]=ceil_scaled(1-p,state_bits);state[stride]=ceil_scaled(p,state_bits)
    choices=((0,0,1-p),(1,0,p-p**4),(1,1,p**4))
    transitions=[]
    for a in (0,1):
        for length in range(stride):
            entries=[]
            for next_a,next_full,probability in choices:
                free=a*next_full
                out_length=min(onset,length+1) if free else 0
                weight=alpha_weight if free and length==onset-1 else continuation_weight if free and length==onset else F(1)
                entries.append((stride*next_a+out_length,ceil_scaled(probability*weight,coefficient_bits)))
            transitions.append(entries)
    for _ in range(steps-1):
        after=[0]*(2*stride)
        for source,entries in enumerate(transitions):
            for destination,coefficient in entries:after[destination]+=(state[source]*coefficient+scale-1)//scale
        state=after
    total=0
    for a in (0,1):
        for length in range(stride):
            weight=alpha_weight if a and length==onset-1 else continuation_weight if a and length==onset else F(1)
            coefficient=ceil_scaled(weight,coefficient_bits)
            total+=(state[stride*a+length]*coefficient+scale-1)//scale
    return F(total,1<<state_bits)


def column_bound(slope=F(1,4),B=18171,r=2100,logb=F(1),tau=F(37,50)):
    ratio=tau/slope;onset=-((-ratio.numerator)//ratio.denominator)
    alpha=onset*slope-tau;assert 0<=alpha<=slope
    p=1-F(r,8192);wa=power_upper(4*logb*alpha);wl=power_upper(4*logb*slope)
    moment=column_moment(p,wa,wl,onset=onset)
    common=128*log_upper(1,moment)+log_upper(1,math.comb(8192,r))
    upper=-logb*(B-512*tau)+common
    sufficient=512*tau+(common+60)/logb
    integer_offset=-((-sufficient.numerator)//sufficient.denominator)
    print(dict(column_read_degree=4,whole_vertical_sequence=True,slope=str(slope),log2_base=str(logb),
               death_weight=str(tau),run_onset=onset,
               offset=str(B),directed_log2_tail_upper=str(upper),upper_float=float(upper),
               sufficient_integer_offset_for_60_bits=integer_offset,actual_regional_quotas=True,
               all_row_sets_of_size=r,zero_coverage_only=True,whole_code_distance=False),flush=True)
    return upper,integer_offset


def column_tests():
    p=F(2,3);wa=F(3,2);wl=F(2);checks=0
    categories=((0,0,1-p),(1,0,p-p**4),(1,1,p**4))
    for steps in range(1,7):
        exact=F(0)
        for path in itertools.product(categories,repeat=steps):
            probability=math.prod(item[2] for item in path);length=0;weight=F(1)
            for t in range(steps):
                free=path[t][0]*(path[t+1][1] if t+1<steps else 1)
                if free:
                    if length==2:weight*=wa
                    elif length>=3:weight*=wl
                    length+=1
                else:length=0
            exact+=probability*weight;checks+=1
        upper=column_moment(p,wa,wl,steps)
        assert exact<=upper<exact+F(1,10**15)
        if steps==3:assert exact==1+(wa-1)*p**9
    print(dict(exact_column_paths=checks,free_terminal_column_checked=True),flush=True)


def penalty_tests():
    """Check a safe field-cutoff penalty, not a new distance certificate."""
    checked=0;p=F(2,3);wa=F(5,4);wl=F(3,2)
    categories=((0,0,1-p),(1,0,p-p**4),(1,1,p**4))
    for onset in (2,4,5):
        for steps in range(1,7):
            exact=F(0)
            for path in itertools.product(categories,repeat=steps):
                probability=math.prod(item[2] for item in path);length=0;weight=F(1)
                for t in range(steps):
                    free=path[t][0]*(path[t+1][1] if t+1<steps else 1)
                    if free:
                        if length==onset-1:weight*=wa
                        elif length>=onset:weight*=wl
                        length+=1
                    else:length=0
                exact+=probability*weight;checked+=1
            upper=column_moment(p,wa,wl,steps,onset)
            assert exact<=upper<exact+F(1,10**15)
    identities=0
    for tau in (F(37,50),F(3,4),F(19,25)):
        slope=F(1,4);ratio=tau/slope;onset=-((-ratio.numerator)//ratio.denominator)
        alpha=onset*slope-tau
        for bits in itertools.product((0,1),repeat=12):
            length=0;value=F(0)
            for bit in bits:
                if bit:length+=1
                else:value+=max(F(0),slope*length-tau);length=0
            value+=slope*length
            windows=lambda size:sum(all(bits[t:t+size]) for t in range(13-size))
            statistic=alpha*windows(onset)+(slope-alpha)*windows(onset+1)
            assert statistic<=value<=statistic+tau;identities+=1
    print(dict(generalized_exact_column_paths=checked,generalized_run_identities=identities),flush=True)
    cutoff=255463;B=18171;old_score=F(91398,5);old_turnover=290432
    for eta in (F(0),F(1,250),F(1,100),F(1,50)):
        target=B-eta*cutoff
        upper,offset=column_bound(B=target,tau=F(37,50)+eta)
        print(dict(field_penalty=str(eta),geometric_target=str(target),
                   zero_coverage_sufficient_offset=offset,positive_coverage_allowance=str(target-offset),
                   frozen_subcount_excess=str(old_score-eta*old_turnover-target),
                   field_cutoff_respected_by_actual_histories=True,
                   no_cutoff_assumption_on_geometric_maximizer=True,
                   all_component_bound_not_established=True),flush=True)


def coverage_set_screen():
    """Count paid input sets, not histories; Holder relaxation may still fail."""
    population=8;inactive=3;t=F(1,2);pi=t/(1+t);p=1-F(inactive,population)
    actual=F(0);thinned=F(0)
    def value(mask):
        return 1<<sum((mask>>j)&7==7 for j in range(population-2))
    for positions in itertools.combinations(range(population),inactive):
        imask=sum(1<<j for j in positions);base=((1<<population)-1)^imask
        for subset in range(1<<inactive):
            count=subset.bit_count();cleared=sum(1<<positions[j] for j in range(inactive) if (subset>>j)&1)
            v=value(base|cleared)
            actual+=t**count*v
            thinned+=pi**count*(1-pi)**(inactive-count)*v
    actual/=math.comb(population,inactive);thinned/=math.comb(population,inactive)
    assert actual==(1+t)**inactive*thinned
    pprime=p+(1-p)*pi;product=F(0)
    for mask in range(1<<population):
        count=mask.bit_count();product+=pprime**count*(1-pprime)**(population-count)*value(mask)
    assert thinned<=product
    print(dict(exact_thinning_label_subsets=math.comb(population,inactive)*(1<<inactive),
               exact_product_comparison_masks=1<<population,positive_polynomial_comparison=True),flush=True)
    N=1<<21;R=256*2100;p=1-F(2100,8192);tau=F(37,50);alpha=F(1,100)
    row_union=log_upper(1,math.comb(8192,2100));B=18171
    minimum=(713983-604080+3)//4;assert minimum==27476
    for count in (minimum,111572,255463):
        t=F(count,R-count);pprime=p+(1-p)*t/(1+t)
        moment=column_moment(pprime,power_upper(4*alpha),F(2))
        upper=-(B+count-512*tau)+R*log_upper(1,1+t)+log_upper(-count,t)
        upper+=128*log_upper(1,moment)+row_union
        print(dict(covered_input_count=count,paid_marker=str(t),effective_active_probability=str(pprime),
                   directed_log2_fixed_count_tail_upper=str(upper),upper_float=float(upper),
                   actual_regional_quotas=True,all_row_sets_of_size=2100,read_degree=4,
                   histories_not_enumerated=True,full_count_range_not_covered=True,whole_code_distance=False),flush=True)


def field_coverage_screen():
    """Keep field constraints before maximizing runs; compare a false idealization."""
    from quartet_rank_score import selected_prior_weights
    weights=selected_prior_weights(6092)
    a=F(39,512);death=F(381,512);constant=111127
    assert weights['activity']<=a and weights['deaths']<=-death
    assert weights['constant']+512*max(F(0),weights['terminal'])<=constant
    R=537600;p=1-F(2100,8192);rank_floor=713983;activity_floor=743225;turnover_cap=255463
    rows=log_upper(1,math.comb(8192,2100))
    # Three discriminating profiles, not a complete cover or parameter search.
    for count,mu,nu,eta,ideal_logb in ((27476,F(0),F(1,4),F(0),F(4)),
                                     (111572,F(1,5),F(0),F(0),F(3)),
                                     (255463,F(0),F(0),F(1),F(4))):
        slope=a+mu+nu;tau=death+eta;ratio=tau/slope;onset=-((-ratio.numerator)//ratio.denominator)
        alpha=onset*slope-tau
        target=(1+mu+eta)*count+mu*rank_floor+nu*activity_floor-eta*turnover_cap-constant
        markers=(F(count,R-count),F(1,5),F(1,6)) if count==111572 else (F(count,R-count),)
        for t in markers:
            pprime=p+(1-p)*t/(1+t)
            for read,logb in ((4,F(1)),(1,ideal_logb)):
                moment=column_moment(pprime,power_upper(read*logb*alpha),power_upper(read*logb*slope),onset=onset)
                upper=-logb*(target-512*tau)+R*log_upper(1,1+t)+log_upper(-count,t)
                upper+=F(512,read)*log_upper(1,moment)+rows
                print(dict(covered=count,rank_multiplier=str(mu),activity_multiplier=str(nu),turnover_multiplier=str(eta),
                           paid_marker=str(t),slope=str(slope),death_weight=str(tau),run_onset=onset,run_threshold=str(target),
                           read_degree=read,log2_base=str(logb),log2_expression=str(upper),expression_float=float(upper),
                           actual_ensemble_upper_bound=(read==4),hypothetical_independent_columns=(read==1),
                           full_count_range_not_covered=True,whole_code_distance=False),flush=True)


def tests():
    checked=0
    for slope in (F(1,4),F(13,50)):
        tau=F(37,50);alpha=3*slope-tau
        for bits in itertools.product((0,1),repeat=12):
            length=0;optimum=F(0)
            for value in bits:
                if value:length+=1
                else:optimum+=max(F(0),slope*length-tau);length=0
            optimum+=slope*length
            w3=sum(all(bits[t:t+3]) for t in range(10))
            w4=sum(all(bits[t:t+4]) for t in range(9))
            statistic=alpha*w3+(slope-alpha)*w4
            assert 0<=optimum-statistic<=tau
            truncated=sum(alpha*all(bits[t:t+3])+(slope-alpha)*all(bits[t:t+4]) for t in range(8))
            assert optimum<=truncated+4*slope
            checked+=1
    for population in (8,16):
        for inactive in range(population+1):
            p=1-F(inactive,population)
            for size in range(population+1):
                probability=F(math.comb(population-inactive,size),math.comb(population,size)) if size<=population-inactive else F(0)
                assert probability<=p**size
    print(dict(exact_free_run_identities=checked,exact_quota_monomial_tests=True),flush=True)


def birth_gate_probability(width,selected,births):
    """Exact necessary partition gate; not independent of the other kernels."""
    import numpy as np
    assert width%4==0 and 0<=births<=width-selected
    current=np.zeros((selected+1,births+1),dtype=object);current[0,0]=1
    terms=[(0,0,1)]+[(a,0,math.comb(4,a)) for a in range(1,5)]+[(0,b,math.comb(4,b)) for b in range(1,5)]
    for group in range(width//4):
        nxt=np.zeros_like(current)
        for a,b,multiplicity in terms:
            if a>selected or b>births:continue
            nxt[a:,b:]+=multiplicity*current[:selected+1-a,:births+1-b]
        current=nxt
    numerator=int(current[selected,births]);denominator=math.comb(width,selected)*math.comb(width-selected,births)
    assert 0<=numerator<=denominator
    return F(numerator,denominator)


def birth_gate_audit():
    checks=0
    for selected in range(9):
        for births in range(9-selected):
            numerator=0
            for s in itertools.combinations(range(8),selected):
                remaining=[j for j in range(8) if j not in s]
                for b in itertools.combinations(remaining,births):
                    if all(not(set(s)&set(range(4*g,4*g+4))) or not(set(b)&set(range(4*g,4*g+4))) for g in range(2)):numerator+=1
            assert birth_gate_probability(8,selected,births)==F(numerator,math.comb(8,selected)*math.comb(8-selected,births));checks+=1
    print(dict(exhaustive_colored_partition_gate_checks=checks),flush=True)
    for selected,births in ((256,43),(256,44),(256,64),(128,32),(384,32)):
        probability=birth_gate_probability(512,selected,births)
        upper=log_upper(1,probability);lower=-log_upper(-1,probability)
        print(dict(selected=selected,hole_births=births,directed_log2_partition_gate_lower=str(lower),
                   directed_log2_partition_gate_upper=str(upper),upper_float=float(upper),
                   inactive_child_requirement_not_imposed=True,
                   not_independent_of_prior_kernel=True,whole_code_distance=False),flush=True)


def marked_local(canonical=False,birth_budget=False):
    """Exact local geometry polynomial, retaining the preceding hole births."""
    result={}
    for a in range(5):
        before=(1<<a)-1
        for e in range(5-a):
            if e and a:continue
            marked=((1<<e)-1)<<a
            multiplicity=math.comb(4,a)*math.comb(4-a,e)
            for after,target in itertools.product(range(16),repeat=2):
                if e and not(target&~after):continue
                if birth_budget and 100*(target&~after).bit_count()<25*e+74*(marked&after).bit_count():continue
                if canonical:
                    if a==0 and target&after:continue
                    h=(before&after).bit_count()
                    if a and 100*(target&~after).bit_count()>=25*a+74*h:continue
                key=(a,e,after.bit_count(),(before&~after).bit_count(),target.bit_count(),
                     target.bit_count() if a else (target&after).bit_count())
                result[key]=result.get(key,0)+multiplicity
    return result


def marked_tests(canonical=False,birth_budget=False):
    explicit={};checks=0
    for colors in itertools.product(range(3),repeat=4):
        before=sum(1<<j for j,c in enumerate(colors) if c==1)
        births=sum(1<<j for j,c in enumerate(colors) if c==2)
        for after,target in itertools.product(range(16),repeat=2):
            checks+=1
            if births and (before or not(target&~after)):continue
            if canonical:
                if not before and target&after:continue
                if before and 100*(target&~after).bit_count()>=25*before.bit_count()+74*(before&after).bit_count():continue
            if birth_budget and 100*(target&~after).bit_count()<25*births.bit_count()+74*(births&after).bit_count():continue
            key=(before.bit_count(),births.bit_count(),after.bit_count(),(before&~after).bit_count(),
                 target.bit_count(),target.bit_count() if before else (target&after).bit_count())
            explicit[key]=explicit.get(key,0)+1
    assert explicit==marked_local(canonical,birth_budget)
    print(dict(marked_local_cases=checks,distinct_marked_terms=len(explicit),
               geometry_state_count=math.comb(514,2),canonical=canonical,birth_budget=birth_budget,
               full_transfer_not_run=True),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--columns',action='store_true');parser.add_argument('--penalty',action='store_true');parser.add_argument('--coverage-sets',action='store_true');parser.add_argument('--field-coverage',action='store_true');args=parser.parse_args()
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(30,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        if args.field_coverage:
            field_coverage_screen()
        elif args.coverage_sets:
            coverage_set_screen()
        elif args.penalty:
            column_tests();penalty_tests()
        elif args.columns:
            column_tests();column_bound(B=15877,logb=F(1))
        else:
            tests();bound();bound(F(1,4),18171,logb=F(9,16));birth_gate_audit();marked_tests()
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
