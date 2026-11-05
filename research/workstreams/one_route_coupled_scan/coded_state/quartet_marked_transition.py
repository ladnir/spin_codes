"""Joint marked-quartet audit for fixed partial-state regional histories.

Floating-point tilts propose parameters; exact rationals verify upper bounds.
No files, serial,60s watchdog,<16MiB numeric arrays,>=20GiB free.
"""
from fractions import Fraction as F
from functools import lru_cache
import itertools
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
from quartet_two_update import log_upper
from quartet_run_components import marked_local,power_upper,marked_tests


@lru_cache(maxsize=None)
def exp2_upper(exponent):
    exponent=F(exponent);integer=exponent.numerator//exponent.denominator
    return F(2)**integer*power_upper(exponent-integer)


def proposed_tilt(terms,target,logb,reference):
    dimension=len(target)
    keys=list(terms);features=np.array([key[:dimension] for key in keys],float)
    logs=np.array([math.log(terms[key])+math.log(2)*logb*(float(reference)*key[4]-key[5]-.74*key[3]) for key in keys])
    target=np.array(list(map(float,target)));parameters=np.zeros(dimension)
    for iteration in range(160):
        raw=logs+features@parameters;offset=max(raw);weights=np.exp(raw-offset);total=sum(weights);probabilities=weights/total
        mean=probabilities@features;gradient=mean-target
        if max(abs(gradient))<1e-10:break
        centered=features-mean;covariance=(centered.T*probabilities)@centered
        step=np.linalg.solve(covariance+np.eye(dimension)*1e-12,gradient)
        step/=max(1.,max(abs(step))/2.)
        objective=offset+math.log(total)-parameters@target
        scale=1.
        for _ in range(30):
            candidate=parameters-scale*step;values=logs+features@candidate;shift=max(values)
            updated=shift+math.log(sum(np.exp(values-shift)))-candidate@target
            if updated<=objective+1e-13:break
            scale/=2
        parameters=candidate
    assert max(abs(gradient))<1e-6
    return tuple(F(math.exp(x)).limit_denominator(1<<20) for x in parameters)


def region_bound(u,births,logb=4,canonical=False,birth_budget=False):
    terms=marked_local(canonical,birth_budget);reference=F(u,512)
    target=(F(u,128),F(births,128),F(u,128),F(births,128),F(2100,2048))
    markers=proposed_tilt(terms,target,logb,reference)
    value=F(0)
    for (a,e,v,d,l,c),multiplicity in terms.items():
        coefficient=multiplicity*exp2_upper(logb*(reference*l-c-F(37,50)*d))
        value+=coefficient*math.prod(marker**power for marker,power in zip(markers,(a,e,v,d,l)))
    upper=2048*log_upper(1,value)
    for marker,power in zip(markers,(16*u,16*births,16*u,16*births,2100)):upper+=log_upper(-power,marker)
    upper+=log_upper(-16,math.comb(512,u))+log_upper(-16,math.comb(512-u,births))+log_upper(-1,math.comb(8192,2100))
    print(dict(selected=u,inherited_hole_births=births,outgoing_deaths=births,log2_base=logb,
               directed_regional_log2_moment_upper=str(upper),upper_float=float(upper),
               exact_regional_quota=2100,constant_state_at_all_16_updates=True,
               rational_markers=list(map(str,markers)),joint_birth_child_gate=True,
               old_canonical_filter=canonical,new_hole_deletion_budget=birth_budget,
               whole_code_distance=False,files_written=0),flush=True)
    return upper


def coverage_prefix(u,births,coverage,birth_budget=False):
    terms=marked_local(True,birth_budget)
    target=(F(u,128),F(births,128),F(u,128),F(births,128),F(2100,2048),F(coverage,2048))
    markers=proposed_tilt(terms,target,0,F(0));assert markers[5]<=1
    value=sum(F(multiplicity)*math.prod(marker**power for marker,power in zip(markers,key)) for key,multiplicity in terms.items())
    upper=2048*log_upper(1,value)
    for marker,power in zip(markers,(16*u,16*births,16*u,16*births,2100,coverage)):upper+=log_upper(-power,marker)
    upper+=log_upper(-16,math.comb(512,u))+log_upper(-16,math.comb(512-u,births))+log_upper(-1,math.comb(8192,2100))
    print(dict(selected=u,births=births,regional_coverage_cap=str(coverage),
               directed_log2_count_prefix_upper=str(upper),upper_float=float(upper),
               rational_markers=list(map(str,markers)),joint_birth_gate=True,old_canonical_filter=True,
               new_hole_deletion_budget=birth_budget,constant_profile_only=True,whole_code_distance=False),flush=True)
    return upper


def tests(canonical=False,birth_budget=False):
    terms=marked_local(canonical,birth_budget);goal=(2,2,2,1,3);coefficient=F(0)
    for left,m1 in terms.items():
        for right,m2 in terms.items():
            if tuple(x+y for x,y in zip(left[:5],right[:5]))==goal:coefficient+=F(m1*m2,1<<(left[5]+right[5]))
    expected=coefficient/(math.comb(8,2)*math.comb(6,2)*math.comb(8,3))
    before={0,1};births={2,3};total=F(0);partitions=0
    for rest in itertools.combinations(range(1,8),3):
        first={0,*rest};groups=(first,set(range(8))-first);partitions+=1
        for out in itertools.combinations(range(8),2):
            after=set(out)
            if len(before-after)!=1:continue
            for labels in itertools.combinations(range(8),3):
                target=set(labels);covered=0;valid=True
                for group in groups:
                    if births&group and (before&group or not((target-after)&group)):valid=False;break
                    if canonical:
                        a=len(before&group);h=len(before&after&group)
                        if not a and target&after&group:valid=False;break
                        if a and 100*len((target-after)&group)>=25*a+74*h:valid=False;break
                    if birth_budget and 100*len((target-after)&group)<25*len(births&group)+74*len(births&after&group):valid=False;break
                    covered+=len(target&group) if before&group else len(target&after&group)
                if valid:total+=F(1,1<<covered)
    assert partitions==35 and expected==total/(35*math.comb(8,3))
    print(dict(exact_two_quartet_marked_coefficient_verified=True,partitions=35,
               weighted_joint_expected_count=str(expected),canonical=canonical,
               birth_budget=birth_budget,files_written=0),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        for canonical,birth_budget in itertools.product((False,True),repeat=2):
            marked_tests(canonical,birth_budget);tests(canonical,birth_budget)
        region_bound(256,44,4,True,True);region_bound(256,44,8,True,True)
        coverage_prefix(256,44,530,True);coverage_prefix(256,44,472,True);coverage_prefix(256,44,432,True)
        print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
