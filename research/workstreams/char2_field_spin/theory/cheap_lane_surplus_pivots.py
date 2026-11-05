"""FLOAT surplus-pivot marker on the frozen four-packet tile operator.

Let q=256 and S=65280, the actual refresh-coefficient support size. A
retained independent pivot contracts P_q by at least1/q; a dependent
pivot survives with probability at most1/S. Starting from dimension
d0<=128r, a nonzero witness permits at most128r-1 independent pivots.
If R pivots are charged in total, the surplus factor is at most
(q/S)^max(0,R-128r+1). For q/S<=t<=1 it is at most
t^(R-128r+1), including when this expression exceeds one.

The inactive-feedback submatrix gives the deterministic pivot floor
rank(M)>=min(16-j,h) for incoming support h and outgoing support j.
Therefore multiply each RAW regular cap by t^min(16-j,h). This is an
upper bound because t<=1. The terminal remains unchanged; no additional
output pivot ranks are used. Apply the marker before raw tile4 convolution
and alpha, then add -alpha*(128r-1)log2(t) to the global log2 bound.
The lemma does not condition on the refresh b coefficients separately.

Both t=1 and.99 use the same frozen raw gamma/joint/singleton tables.
The repeated fully paired marginal is artificial; these FLOAT scores are
not support-class bounds or certificates. One48-second invocation runs
both marks serially after a20GiB free-space preflight. Estimated peak
memory<256MiB. Run python -B; stdout only, no frozen edits or data files.
"""

import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
from types import SimpleNamespace

sys.dont_write_bytecode = True
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
TILES_SHA256 = "a96bcfc95c927501264f9787692e509b03f64a8f968570bc0ac631f581d8f34e"


def load_tiles():
    path = Path(__file__).with_name("cheap_lane_mask_tiles.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=TILES_SHA256:
        raise RuntimeError("frozen mask-tiles helper hash mismatch")
    spec=importlib.util.spec_from_file_location("frozen_mask_tiles",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def marker_test():
    # Exact rational check, independent of the floating local implementation.
    from fractions import Fraction
    c=Fraction(1,255)
    checks=0
    for t in (c,Fraction(99,100),Fraction(1)):
        for dimension in (1,2,8,32):
            for pivots in range(65):
                assert c**max(0,pivots-dimension+1)<=t**(pivots-dimension+1)
                checks+=1
    return checks


def main():
    if len(sys.argv)!=1:
        raise ValueError("fixed diagnostic accepts no parameters")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on actual source volume")
    start=time.monotonic()
    tiles=load_tiles(); grouped=tiles.load_grouped(); singleton=grouped.load_singleton()
    paired=singleton.load_block().load_paired(start+48)
    tests={**tiles.toy_test(grouped,paired),**grouped.toy_test(paired),
           "exact_surplus_marker_inequalities":marker_test()}
    args=SimpleNamespace(r=256,z=.85,alpha=.523,activity_mark=1.0)
    filtered,middle=paired.imports(); np=paired.np
    m,terminal,degrees,weights=grouped.raw_tables(args,singleton,paired,filtered)
    ranks=np.array([[min(16-j,h) for j in range(17)] for h in range(17)])
    assert np.all(ranks[0]==0) and np.all(ranks[:,16]==0)
    assert np.array_equal(ranks[:,0],np.arange(17))
    results=[]
    for t in (1.0,.99):
        paired.guard(); stage=time.monotonic()
        assert 256/65280<=t<=1
        factors=np.power(t,ranks)
        marked=m*factors[None,:,:]
        if t==1:
            assert np.array_equal(marked,m)
        assert np.all(marked<=m)
        a,b=tiles.raw_polynomials(marked,terminal,degrees,weights,np,middle)
        a,b=tiles.double_pair(a,b,2,128,paired)
        a=tiles.jensen(a,128,[math.comb(16,j) for j in range(17)],args.alpha,np)
        b=tiles.jensen(b,128,[1],args.alpha,np)
        result=grouped.evaluate(args,a,b,5,paired,filtered,middle)
        chain=result["artificial_chain"]
        correction=-args.alpha*(128*args.r-1)*math.log2(t)
        results.append({"pivot_mark":t,"regional_log2":result["regional_log2"],
            "raw_nonquiet_log2":chain["nonquiet_log2_bound"],
            "global_dimension_log2_correction":correction,
            "adjusted_nonquiet_log2":chain["nonquiet_log2_bound"]+correction,
            "adjusted_quiet_log2":chain["quiet_log2_bound"]+correction,
            "adjusted_restart_floor_log2":result["artificial_restart_floor_log2"]+correction,
            "largest_first_active_column":chain["largest_first_active_columns"][0]["column_zero_based"],
            "seconds":time.monotonic()-stage})
    print(json.dumps({"scope":"FLOAT artificial repeated-marginal surplus-pivot diagnostic",
        "whole_code_certificate":False,"actual_support_class_bound":False,
        "parameters":vars(args),"tile_packets":4,"pivot_rank_floor":"min(16-outgoing,incoming)",
        "results":results,
        "adjusted_score_change":results[1]["adjusted_nonquiet_log2"]-results[0]["adjusted_nonquiet_log2"],
        "tests":tests,"estimated_peak_bytes":64*(1<<20),
        "seconds":time.monotonic()-start,"generated_bytes":0,
        "tiles_source_sha256":TILES_SHA256},indent=2))


if __name__=="__main__":
    main()
