"""Exact subcount audit of the joint marked partial-state family.

Bounded numerical proposals, exact integer verification. No codeword claim.
Serial,60s watchdog,<16MiB numeric storage,>=20GiB free,no files.
"""
from fractions import Fraction as F
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
from quartet_two_update import log_upper
from quartet_run_components import marked_local
from quartet_marked_transition import proposed_tilt
from quartet_neutral_type import disjoint_moves,verify_moves


def propose(types,label_count,coverage=27):
    keys=list(types);size=len(keys);features=np.array(keys,dtype=float)
    target=(256,44,256,44,label_count,coverage)
    markers=proposed_tilt(types,tuple(F(v,128) for v in target),0,F(0))
    weights=np.array([math.log(types[key])+sum(power*math.log(float(marker)) for marker,power in zip(markers,key)) for key in keys])
    weights=np.exp(weights-max(weights));weights/=sum(weights);center=128*weights
    matrix=np.vstack((np.ones(size),features.T));desired=np.array((128,*target),dtype=np.int64)
    rng=np.random.default_rng(20261030+label_count);best=None;best_error=float('inf')
    for trial in range(20000):
        base=np.floor(center).astype(np.int64)
        if trial%4:base+=(rng.random(size)<center-base).astype(np.int64)
        basis=rng.choice(size,7,replace=False)
        try:change=np.linalg.solve(matrix[:,basis],desired-matrix@base)
        except np.linalg.LinAlgError:continue
        if max(abs(change))>1000 or max(abs(change-np.rint(change)))>1e-7:continue
        base[basis]+=np.rint(change).astype(np.int64)
        if min(base)<0 or not np.array_equal(matrix@base,desired):continue
        error=float(sum((base-center)**2/(center+.1)))
        if error<best_error:best=base.copy();best_error=error
    if best is None:raise RuntimeError('bounded proposal found no integer witness')
    ns=list(map(int,best));assert sum(ns)==128
    assert tuple(sum(n*key[j] for n,key in zip(ns,keys)) for j in range(6))==target
    print(dict(proposal_labels=label_count,coverage=coverage,error=best_error),flush=True)
    return ns


def optimize_histogram(types,keys,ns):
    """Exact improving two-type exchanges; floats choose candidates only."""
    pairs={}
    for i in range(len(keys)):
        for j in range(i,len(keys)):
            key=tuple(a+b for a,b in zip(keys[i],keys[j]));pairs.setdefault(key,[]).append((i,j))
    moves=[]
    for choices in pairs.values():
        for t,left in enumerate(choices):
            for right in choices[t+1:]:
                if set(left)&set(right):continue
                moves.append((*left,*right));moves.append((*right,*left))
    indices=np.array(moves,dtype=np.int16);assert indices.nbytes<4*(1<<20)
    counts=np.array(ns,dtype=np.int64);multiplicities=np.array([types[key] for key in keys],dtype=float)
    i,j,k,l=indices.T;ratio=multiplicities[k]*multiplicities[l]/(multiplicities[i]*multiplicities[j])
    gain=F(1);iterations=0
    for iterations in range(2000):
        values=ratio*counts[i]*(counts[j]-(i==j))/((counts[k]+1)*(counts[l]+1+(k==l)))
        chosen=int(np.argmax(values))
        if values[chosen]<=1+1e-12:break
        change={}
        for index,delta in zip(indices[chosen],(-1,-1,1,1)):change[int(index)]=change.get(int(index),0)+delta
        exact=F(1)
        for index,delta in change.items():
            assert counts[index]+delta>=0
            exact*=F(math.factorial(int(counts[index])),math.factorial(int(counts[index]+delta)))*F(types[keys[index]])**delta
        assert exact>1
        for index,delta in change.items():counts[index]+=delta
        gain*=exact
    assert np.array_equal(np.array(keys).T@counts,np.array(keys).T@np.array(ns)) and sum(counts)==128
    print(dict(histogram_exchanges=iterations,exchange_candidates=len(indices),exact_log2_gain_lower=str(-log_upper(-1,gain))),flush=True)
    return list(map(int,counts))


def audit():
    types=marked_local(True,True);keys=list(types);coefficients=[];counts=[];all_moves=[]
    for labels in (131,132):
        ns=optimize_histogram(types,keys,propose(types,labels));counts.append(ns)
        _,moves=disjoint_moves(types,keys,ns);all_moves.append(moves)
        gain=verify_moves(types,keys,ns,moves)
        coefficient=math.factorial(128)//math.prod(math.factorial(n) for n in ns)
        coefficient*=math.prod(types[key]**n for n,key in zip(ns,keys))
        coefficient*=gain;assert coefficient.denominator==1
        coefficients.append(coefficient.numerator)
        print(dict(labels=labels,moves=len(moves),gain_log2_upper=str(log_upper(1,gain))),flush=True)
    numerator=coefficients[0]**12*coefficients[1]**4
    denominator=(math.comb(512,256)*math.comb(256,44))**16*math.comb(8192,2100)
    lower=-log_upper(-1,numerator)-log_upper(1,denominator)
    upper=log_upper(1,numerator)+log_upper(-1,denominator)
    print(dict(regional_A=4096,regional_D=704,regional_C=432,regional_labels=2100,
               regional_score=str(F(4096,4)-432-F(37,50)*704),
               directed_log2_subcount_lower=str(lower),directed_log2_subcount_upper=str(upper),
               lower_float=float(lower),actual_initial_terminal_not_imposed=True,
               whole_code_distance=False,files_written=0),flush=True)
    print(dict(replay_counts=counts,replay_moves=all_moves),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(60,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:audit();print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
