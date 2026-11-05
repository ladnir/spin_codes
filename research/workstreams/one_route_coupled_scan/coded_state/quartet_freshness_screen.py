"""FLOAT fixed-count trajectories retaining previous fully-active quartets.

Not a distance certificate.90s,one thread,<32MiB,no files; Peach only.
"""
from functools import lru_cache
from collections import Counter
import math
import os
from pathlib import Path
import shutil
import sys
import threading
import time
import numpy as np
from quartet_history_screen import occupation,column_perron,point_mass_upper,logchoose,coefficient_cauchy_upper


LOGB=4;SLOPE=707/2560;TAU=381/512
CENTER=408


def candidate():
    return column_perron(LOGB,SLOPE,TAU,CENTER)


def local_entry(a,previous_full,l,k,logJ,h):
    ds=(h[1],h[2],h[3],h[4]*2**(LOGB*(3*SLOPE-TAU)),h[4]*2**(LOGB*SLOPE))
    lo=-80.;hi=80.
    for _ in range(55):
        x=(lo+hi)/2;z=2**x
        i=max(range(1,5),key=lambda i:(1+ds[i]*z)/h[i])
        derivative=(512-a)*ds[0]*z/(1+ds[0]*z)
        derivative+=(a-previous_full)*ds[1]*z/(1+ds[1]*z)
        derivative+=previous_full*ds[i]*z/(1+ds[i]*z)-k
        if derivative<0:lo=x
        else:hi=x
    x=(lo+hi)/2;z=2**x
    value=(512-a)*math.log2(1+ds[0]*z)
    value+=(a-previous_full)*math.log2((1+ds[1]*z)/h[1])
    value+=previous_full*max(math.log2((1+ds[i]*z)/h[i]) for i in range(1,5))-k*x
    variances=[d*z/(1+d*z)**2 for d in ds]
    lower_variance=(512-a)*variances[0]+(a-previous_full)*variances[1]+previous_full*min(variances[1:])
    return logJ-logchoose(512,k)+(l-k)*math.log2(h[1])+value+math.log2(point_mass_upper(lower_variance,True))


@lru_cache(maxsize=64)
def matrix(a,l):
    h=candidate()
    old=occupation(a);new=occupation(l)
    logs=np.array([[local_entry(a,previous_full,l,k,logJ,h) for k,J,logJ in new]
                   for previous_full,Jold,logJold in old])
    scale=float(np.max(logs));weights=np.exp2(logs-scale)
    assert np.all(weights>0), 'Diagnostic matrix underflow'
    return weights,scale


def evaluate(word):
    h=candidate();first=word[0]
    initial=np.array([logJ+first*math.log2(h[1]) for k,J,logJ in occupation(first)])
    total=float(np.max(initial));vec=np.exp2(initial-total)
    for a,l in zip(word,word[1:]):
        norm=float(np.sum(vec));assert norm>0;vec/=norm;total+=math.log2(norm)
        transition,scale=matrix(a,l);vec=vec@transition;total+=scale
    last=word[-1];rewards=(SLOPE,2*SLOPE,3*SLOPE,SLOPE+TAU)
    ratios=[LOGB*r-math.log2(hi) for r,hi in zip(rewards,h[1:])]
    terminal=np.array([(last-k)*ratios[0]+k*max(ratios) for k,J,logJ in occupation(last)])
    scale=float(np.max(terminal));total+=scale+math.log2(float(vec@np.exp2(terminal-scale)))
    total+=44*(logchoose(6527,435)-logchoose(8192,2100))
    total+=212*(logchoose(6528,436)-logchoose(8192,2100))
    return total-LOGB*165556+logchoose(8192,2100)


@lru_cache(maxsize=2048)
def scalar_step(a,l):
    h=candidate()
    return coefficient_cauchy_upper(a,l,LOGB,SLOPE,TAU,h,True,point_mass='refined')[1]


def scalar_path(word):
    h=candidate()
    total=logchoose(512,word[0])+word[0]*math.log2(h[1])
    total+=sum(count*scalar_step(a,l) for (a,l),count in Counter(zip(word,word[1:])).items())
    rewards=(SLOPE,2*SLOPE,3*SLOPE,SLOPE+TAU)
    total+=word[-1]*max(LOGB*r-math.log2(hi) for r,hi in zip(rewards,h[1:]))
    total+=44*(logchoose(6527,435)-logchoose(8192,2100))
    total+=212*(logchoose(6528,436)-logchoose(8192,2100))
    return total-LOGB*165556+logchoose(8192,2100)


def maximizing_words():
    counts=list(range(392,425));n=len(counts);lo=-16;center=16
    logs=np.array([[scalar_step(a,l) for l in counts] for a in counts])
    state=np.full((1,n),-np.inf);state[0,center]=0;history=[]
    for step in range(1,17):
        after=np.full((32*step+1,n),-np.inf);previous=np.full(after.shape,-1,dtype=np.int16)
        for j,l in enumerate(counts):
            values=state+logs[:,j][None,:]
            best=np.max(values,axis=1);arg=np.argmax(values,axis=1)
            offset=l-408-lo
            after[offset:offset+len(state),j]=best
            previous[offset:offset+len(state),j]=arg
        state=after;history.append(previous)
    result=[]
    for target in (-1,0):
        total=target;current=center;word=[]
        optimum=state[total-16*lo,current]
        for step in range(16,0,-1):
            old=int(history[step-1][total-step*lo,current]);assert old>=0
            word.append(counts[current]);total-=counts[current]-408;current=old
        word.reverse();assert total==0 and current==center
        replay=sum(logs[a-392,l-392] for a,l in zip([408]+word[:-1],word))
        assert abs(replay-optimum)<1e-8
        result.append(word)
    return result


def main():
    global CENTER
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        if '--centers' in sys.argv:
            low=[424]*7+[407]+[392]*7+[408]
            high=[424]*7+[408]+[392]*7+[408]
            word=low*44+high*212
            for CENTER in (404,408,412,416,420):
                matrix.cache_clear();scalar_step.cache_clear()
                print(dict(hard_profile='regional_max_weight_path',potential_center=CENTER,
                           float_log2_freshness_path=evaluate(word),
                           constant_potential_across_counts=True,
                           directed_certificate=False,whole_code_distance=False,
                           seconds=time.monotonic()-start),flush=True)
            return
        for label,pair in (('near_flat',None),('alternating_400_416',(400,416)),('alternating_392_424',(392,424))):
            word=[]
            for region in range(256):
                local=[408]*16 if pair is None else list(pair)*8
                if region<44:local[-1]-=1
                assert sum(local)==6527+(region>=44)
                word.extend(local)
            result=evaluate(word)
            print(dict(profile=label,float_log2_fixed_path_freshness_envelope=result,
                       float_log2_fixed_path_previous_envelope=scalar_path(word),
                       previous_full_quartet_count_summed=True,regional_quotas_exact=True,
                       covered=111572,only_fixed_input_count_path=True,
                       directed_certificate=False,whole_code_distance=False,
                       seconds=time.monotonic()-start),flush=True)
        low,high=maximizing_words();word=low*44+high*212
        print(dict(profile='regional_max_weight_path',low_region_word=low,high_region_word=high,
                   float_log2_fixed_path_previous_envelope=scalar_path(word),
                   float_log2_fixed_path_freshness_envelope=evaluate(word),
                   regional_quotas_exact=True,covered=111572,
                   directed_certificate=False,whole_code_distance=False),flush=True)
        print(dict(seconds=time.monotonic()-start,files_written=0,matrix_count=matrix.cache_info().currsize),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
