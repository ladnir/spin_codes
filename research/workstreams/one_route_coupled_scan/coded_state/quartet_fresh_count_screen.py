"""FLOAT count-path screen retaining the actual number of fresh coordinates.

No distance certificate. Serial Peach only; 90s, <64MiB arrays, no files.
"""
import math
import os
from pathlib import Path
import shutil
import sys
import threading
import time
from functools import lru_cache
import numpy as np
from quartet_history_screen import occupation, column_perron, logchoose

LOGB=4; LAMBDA=707/2560; TAU=381/512
BIN=32
H=np.array(column_perron(LOGB,LAMBDA,TAU))
D=np.array((H[1],H[2],H[3],H[4]*2**(LOGB*(3*LAMBDA-TAU)),H[4]*2**(LOGB*LAMBDA)))


def buckets(a):
    return [(r,min(a,r+BIN-1)) for r in range(0,a+1,BIN)]


def active_coefficient(a,rlo,rhi,u):
    """Uniform Cauchy bound for n1 in [rlo,rhi], n2+n3+n4=a-n1.

    Maximize both the old type and fresh count BEFORE selecting z.
    The auxiliary point-mass factor uses the minimum possible variance.
    """
    lo=np.full(u.shape,-80.);hi=np.full(u.shape,80.)
    for _ in range(55):
        x=(lo+hi)/2;z=np.exp2(x)
        logs=np.log2((1+D[2:,None]*z)/H[2:,None])
        i=np.argmax(logs,axis=0)+2
        fresh=np.log2((1+D[1]*z)/H[1]);old=np.max(logs,axis=0)
        r=np.where(fresh>=old,rhi,rlo)
        derivative=r*D[1]*z/(1+D[1]*z)+(a-r)*D[i]*z/(1+D[i]*z)-u
        lo=np.where(derivative<0,x,lo);hi=np.where(derivative>=0,x,hi)
    x=(lo+hi)/2;z=np.exp2(x)
    fresh=np.log2((1+D[1]*z)/H[1]);old=np.max(np.log2((1+D[2:,None]*z)/H[2:,None]),axis=0)
    result=np.maximum(rlo*fresh+(a-rlo)*old,rhi*fresh+(a-rhi)*old)-u*x
    vf=D[1]*z/(1+D[1]*z)**2
    vo=np.min(D[2:,None]*z/(1+D[2:,None]*z)**2,axis=0)
    variance=np.minimum(rlo*vf+(a-rlo)*vo,rhi*vf+(a-rhi)*vo)
    mass=np.minimum(1.,np.sqrt(math.pi/(8*variance)))
    for theta in (.5,.75,1.,1.25):
        c=1-theta*theta/12
        mass=np.minimum(mass,1/np.sqrt(2*math.pi*variance*c)+np.exp(-variance*theta*theta*c/2))
    return result+np.log2(mass)


@lru_cache(maxsize=64)
def active_coefficients(a):
    # k-v0 has only a+1 possible values. Its bound does not depend on l,
    # so compute it once per incoming bucket rather than for every (k,v0).
    u=np.arange(a+1,dtype=float)
    return np.array([active_coefficient(a,rlo,rhi,u) for rlo,rhi in buckets(a)])


@lru_cache(maxsize=32)
def transition(a,l):
    # v0 is selected from the known inactive-input coordinates. Retaining
    # it makes the outgoing fresh count exactly l-k+v0, not merely >=l-k.
    terms=[]
    for k,J,logJ in occupation(l):
        for v0 in range(max(0,k-a),min(k,512-a)+1):
            rnext=l-k+v0
            base=logJ-logchoose(512,k)+logchoose(512-a,v0)
            base+=v0*math.log2(D[0])+(l-k)*math.log2(H[1])
            terms.append((k-v0,rnext//BIN,base))
    u=np.array([x[0] for x in terms],dtype=int)
    dest=np.array([x[1] for x in terms],dtype=int)
    base=np.array([x[2] for x in terms])
    rows=[]
    for coeff in active_coefficients(a):
        logs=base+coeff[u]
        row=np.full(len(buckets(l)),-np.inf)
        np.logaddexp2.at(row,dest,logs)
        rows.append(row)
    logs=np.array(rows);scale=float(np.max(logs));weights=np.exp2(logs-scale)
    assert np.all(weights[np.isfinite(logs)]>0), 'matrix underflow'
    return weights,scale


def evaluate(word):
    first=word[0];vec=np.zeros(len(buckets(first)));vec[first//BIN]=1
    total=logchoose(512,first)+first*math.log2(H[1])
    for a,l in zip(word,word[1:]):
        norm=float(np.sum(vec));assert norm>0;vec/=norm;total+=math.log2(norm)
        matrix,scale=transition(a,l);vec=vec@matrix;total+=scale
    last=word[-1]
    rewards=LOGB*np.array((LAMBDA,2*LAMBDA,3*LAMBDA,LAMBDA+TAU))-np.log2(H[1:])
    terminal=np.array([max(r*rewards[0]+(last-r)*max(rewards[1:]) for r in ends) for ends in buckets(last)])
    scale=float(np.max(terminal));total+=scale+math.log2(float(vec@np.exp2(terminal-scale)))
    total+=44*(logchoose(6527,435)-logchoose(8192,2100))
    total+=212*(logchoose(6528,436)-logchoose(8192,2100))
    return total-LOGB*165556+logchoose(8192,2100)


def rotation_region(middle):
    """Sum 15 distinct cyclic rotations of the first 15 counts, ending408.

    Return both the incoming408 matrix and the true zero-state initial row.
    Every local word has the required regional total; no regional reset.
    """
    core=[424]*7+[middle]+[392]*7
    entries=[]
    for shift in range(15):
        word=core[shift:]+core[:shift]+[408]
        first=word[0];mat,scale=transition(408,first)
        vec=np.zeros(len(buckets(first)));vec[first//BIN]=1
        vscale=logchoose(512,first)+first*math.log2(H[1])
        mat=mat.copy()
        for a,l in zip(word,word[1:]):
            step,ss=transition(a,l)
            mat=mat@step;vmax=float(np.max(mat));assert vmax>0
            mat/=vmax;scale+=ss+math.log2(vmax)
            vec=vec@step;vmax=float(np.max(vec));assert vmax>0
            vec/=vmax;vscale+=ss+math.log2(vmax)
        entries.append((mat,scale,vec,vscale))
    common=max(max(s,vs) for m,s,v,vs in entries)
    matrix=sum(m*2**(s-common) for m,s,v,vs in entries)
    initial=sum(v*2**(vs-common) for m,s,v,vs in entries)
    assert np.all(matrix>=0) and np.any(initial>0)
    paid=435 if middle==407 else 436
    common+=logchoose(6092+paid,paid)-logchoose(8192,2100)
    with np.errstate(divide='ignore'):
        return np.log2(matrix)+common,np.log2(initial)+common


def rotation_sum():
    low,initial_low=rotation_region(407)
    high,initial_high=rotation_region(408)
    state=np.full((45,len(buckets(408))),-np.inf)
    state[0]=initial_high;state[1]=initial_low
    def multiply(vector,matrix):
        values=vector[:,:,None]+matrix[None,:,:]
        peak=np.max(values,axis=1);safe=np.where(np.isfinite(peak),peak,0.)
        with np.errstate(divide='ignore'):
            return safe+np.log2(np.sum(np.exp2(values-safe[:,None,:]),axis=1))
    for _ in range(255):
        after=multiply(state,high)
        after[1:]=np.logaddexp2(after[1:],multiply(state[:-1],low))
        state=after
    rewards=LOGB*np.array((LAMBDA,2*LAMBDA,3*LAMBDA,LAMBDA+TAU))-np.log2(H[1:])
    terminal=np.array([max(r*rewards[0]+(408-r)*max(rewards[1:]) for r in ends) for ends in buckets(408)])
    logs=state[44]+terminal;peak=float(np.max(logs))
    return peak+math.log2(float(np.sum(np.exp2(logs-peak))))-LOGB*165556+logchoose(8192,2100)


def tests():
    import itertools
    # Exact independent check of the regional count/hidden-state summation.
    matrices=(((1,2),(3,1)),((2,1),(1,4)))
    initial=((2,3),(1,2));terminal=(3,2)
    dp=[[0,0] for _ in range(7)];dp[0]=list(initial[0]);dp[1]=list(initial[1])
    for _ in range(5):
        after=[[0,0] for _ in range(7)]
        for count,row in enumerate(dp):
            for flag,m in enumerate(matrices):
                if count+flag<=6:
                    for j in range(2):after[count+flag][j]+=sum(row[i]*m[i][j] for i in range(2))
        dp=after
    direct=[0]*7
    for word in itertools.product(range(2),repeat=6):
        row=initial[word[0]]
        for flag in word[1:]:row=tuple(sum(row[i]*matrices[flag][i][j] for i in range(2)) for j in range(2))
        direct[sum(word)]+=sum(x*y for x,y in zip(row,terminal))
    assert direct==[sum(x*y for x,y in zip(row,terminal)) for row in dp]
    for middle in (407,408):
        core=[424]*7+[middle]+[392]*7
        words={tuple(core[s:]+core[:s]+[408]) for s in range(15)}
        assert len(words)==15 and all(sum(word)==6120+middle for word in words)
    print(dict(exact_matrix_paths_checked=64,regional_words_checked=30),flush=True)


def main():
    global BIN
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests()
        if '--rotations' in sys.argv:
            BIN=8
            print(dict(profile='15_rotations_per_region_all_balanced_assignments',
                       fresh_count_bucket_width=BIN,float_log2_family_upper=rotation_sum(),
                       log2_input_count_paths=256*math.log2(15)+logchoose(256,44),
                       directed_certificate=False,whole_code_distance=False,
                       seconds=time.monotonic()-start),flush=True)
            return
        low=[424]*7+[407]+[392]*7+[408]
        high=[424]*7+[408]+[392]*7+[408]
        for BIN in (32,16,8):
            transition.cache_clear();active_coefficients.cache_clear()
            print(dict(profile='clumped_424_392',fresh_count_bucket_width=BIN,
                       float_log2_fixed_path_upper=evaluate(low*44+high*212),
                       directed_certificate=False,whole_code_distance=False,
                       seconds=time.monotonic()-start),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
