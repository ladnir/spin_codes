"""FLOAT bounded-band test of a second potential on the fresh-count buckets.

Serial Peach,90s,<64MiB arrays,no files. Not a distance certificate.
"""
import math
import os
from pathlib import Path
import shutil
import threading
import time
import numpy as np
import quartet_fresh_count_screen as fresh
from quartet_history_screen import logchoose
from quartet_count_band import evaluate_band


def perron(a):
    matrix,scale=fresh.transition(a,a)
    v=np.ones(len(matrix))
    for _ in range(512):
        after=matrix@v;after/=float(np.max(after))
        assert np.all(after>0)
        error=float(np.max(np.abs(np.log2(after/v))))
        v=after
        if error<1e-12:return v
    raise AssertionError('fresh-count potential proposal did not converge')


def main():
    assert shutil.disk_usage(Path(__file__).resolve().anchor).free>=20*(1<<30)
    timer=threading.Timer(90,lambda:os._exit(124));timer.daemon=True;timer.start();start=time.monotonic()
    try:
        fresh.BIN=8;fresh.transition.cache_clear();fresh.active_coefficients.cache_clear()
        counts=list(range(392,425));vs={a:perron(a) for a in counts}
        print(dict(stage='potentials',count=len(vs),seconds=time.monotonic()-start),flush=True)
        logs=np.empty((len(counts),len(counts)))
        for i,a in enumerate(counts):
            for j,l in enumerate(counts):
                matrix,scale=fresh.transition(a,l)
                ratios=(matrix@vs[l])/vs[a]
                assert np.all(ratios>0)
                logs[i,j]=scale+math.log2(float(np.max(ratios)))
        initial=[logchoose(512,l)+l*math.log2(fresh.H[1])+math.log2(vs[l][l//fresh.BIN]) for l in counts]
        rewards=fresh.LOGB*np.array((fresh.LAMBDA,2*fresh.LAMBDA,3*fresh.LAMBDA,fresh.LAMBDA+fresh.TAU))-np.log2(fresh.H[1:])
        terminal=[]
        for l in counts:
            base=np.array([max(r*rewards[0]+(l-r)*max(rewards[1:]) for r in ends) for ends in fresh.buckets(l)])
            terminal.append(float(np.max(base-np.log2(vs[l]))))
        print(dict(stage='local_coefficients',count=logs.size,seconds=time.monotonic()-start),flush=True)
        for left,right in ((400,416),(396,420),(392,424)):
            sl=slice(left-counts[0],right-counts[0]+1)
            ordered,total=evaluate_band(logs[sl,sl],counts[sl],fresh.LOGB,fresh.LAMBDA,fresh.TAU,
                                       [fresh.H]*len(counts[sl]),initial[sl],terminal[sl])
            print(dict(input_count_band=[left,right],float_log2_all_balanced_assignments=total,
                       float_log2_ordered_assignment=ordered,fresh_count_bucket_width=fresh.BIN,
                       all_paths_inside_band=True,directed_certificate=False,whole_code_distance=False,
                       seconds=time.monotonic()-start),flush=True)
    finally:timer.cancel()


if __name__=='__main__':main()
