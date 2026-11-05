"""Lower diagnostic for the actual geometric Markov partition, not failure.

Exact one-coordinate mean identities; FLOAT full-width Jensen expression.
Serial Peach,90s,<64MiB,no files. Does not prove or refute code distance.
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
from quartet_window_transfer import kernel
from quartet_history_screen import logchoose,logsum


def mean_reward(counts,width,alpha,lam):
    """Expected nonterminal reward, with no reward before three free bits."""
    p=[F(a,width) for a in counts]
    g=[F(math.comb(a,4),math.comb(width,4)) for a in counts]
    total=F(0)
    for t in range(3,len(counts)):
        total+=alpha*p[t-3]*g[t-2]*g[t-1]*g[t]
        if t>=4:total+=(lam-alpha)*p[t-4]*g[t-3]*g[t-2]*g[t-1]*g[t]
    return width*total


def exact_mean(counts,width):
    probabilities={(width,0,0,0,0):F(1)};reward=F(0)
    for l in counts:
        after=defaultdict(F)
        for n,p in probabilities.items():
            for (out,a,e),q in kernel(n):
                if a!=l:continue
                probability=p*q/math.comb(width,l)
                after[out]+=probability;reward+=probability*e
        assert sum(after.values())==1
        probabilities=dict(after)
    return reward


def tests():
    checked=0
    for counts in itertools.product(range(5),repeat=5):
        assert mean_reward(counts,4,F(1),F(2))==exact_mean(counts,4);checked+=1
    for counts in ((8,8,8,8,8),(6,8,7,8,6),(8,4,8,8,8),(3,5,7,8,4),(0,8,8,8,8)):
        assert mean_reward(counts,8,F(1),F(2))==exact_mean(counts,8);checked+=1
    print(dict(exact_coordinate_mean_paths=checked),flush=True)


def patterns(total):
    """A finite, explicitly deduplicated subset of regional count paths."""
    found=set()
    for full in range(13):
        remaining=total-512*full
        if not 0<=remaining<=512*(16-full):continue
        lo,high=divmod(remaining,16-full)
        path=(512,)*full+(lo+1,)*high+(lo,)*(16-full-high)
        for shift in range(16):found.add(path[shift:]+path[:shift])
    for empty in range(4):
        lo,high=divmod(total,16-empty)
        if lo+(high>0)>512:continue
        path=(0,)*empty+(lo+1,)*high+(lo,)*(16-empty-high)
        for shift in range(16):found.add(path[shift:]+path[:shift])
    assert all(len(p)==16 and sum(p)==total for p in found)
    return sorted(found)


def diagnostic():
    lam=F(707,2560);tau=F(381,512);alpha=3*lam-tau;terms={}
    for total in (6527,6528):
        paths=patterns(total);values=[]
        for path in paths:
            entropy=sum(logchoose(512,a) for a in path)
            mean=float(mean_reward(path,512,alpha,lam))
            values.append(entropy+4*mean)
        j=max(range(len(paths)),key=lambda i:values[i])
        terms[total]=logsum(values)+logchoose(total,total-6092)-logchoose(8192,2100)
        print(dict(total_allowed=total,distinct_paths=len(paths),best_path=paths[j],
                   expected_internal_reward=float(mean_reward(paths[j],512,alpha,lam)),
                   float_regional_lower=terms[total]),flush=True)
    lower=212*terms[6528]+44*terms[6527]+logchoose(256,44)
    lower+=logchoose(8192,2100)-4*165556
    print(dict(float_geometric_markov_partition_lower=lower,covered_count=111572,
               inactive_rows=2100,not_failure_lower=True,directed=False),flush=True)


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        tests();diagnostic();print(dict(seconds=time.monotonic()-start,files_written=0),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
