"""Bracket the quiet-mark stationary point of one frozen FLOAT control.

The control fixes z=.85, alpha=.05, exact regional route coefficients,
consecutive256-group source profiles, and actual packet states0/full16.
Birth and floor marks remain1. It is a positive submajorant, not an upper
bound for the full source class and not a distance certificate.

Let x be the natural logarithm of the raw quiet mark. The corrected score
is convex in x, since all coefficients and chain sums are positive.
Its derivative is alpha*(E_x[q]-24575)/log(2), where E_x[q] is the
hierarchical tilted quiet count from the pinned positive-dual helper.
We bracket its zero and use at most16 safeguarded secant evaluations.
This is one-dimensional root finding, not a parameter scan. The bracket
and all numerical comparisons are FLOAT, not directed verification.

Cache geometry and raw local caps once. Each evaluation changes only the
regular0->0 coefficient, then recomputes the positive two-state dual
convolution and honest terminal. No finite differences are used.

One48-second job,20GiB free,128MiB peak, stdout only. No frozen edits,
generated files, bytecode caches, other parameter changes, or retries.
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

sys.dont_write_bytecode=True
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
PIN="84b6b757a3fc4e85f28c0c3e9f59f6063ae8daf70cc3704aca6feb9ebe27600b"
Q_TARGET=24575
MAX_EVALUATIONS=16
LN2=math.log(2)


def main():
    if len(sys.argv)!=1:
        raise ValueError("fixed one-dimensional root job accepts no scan arguments")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on the actual source volume")
    started=time.monotonic(); deadline=started+48
    import numpy as np
    from scipy.special import logsumexp
    path=Path(__file__).with_name("cheap_lane_branch_sensitivity.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:
        raise RuntimeError("frozen positive-dual helper hash mismatch")
    spec=importlib.util.spec_from_file_location("stationarity_dual",path)
    dual=importlib.util.module_from_spec(spec); spec.loader.exec_module(dual)
    dual.np=np; dual.lse=logsumexp; dual.DEADLINE=deadline
    pool,degree,pilot=dual.load()
    checked=dual.tests(pilot)
    singleton,filtered=pilot.imports()
    shapes,index,degrees,coeff,records,route_checks=pilot.geometry()
    m,t,_=pilot.raw_tables(shapes,index,records,degree.Z,singleton,filtered)
    raw,end=degree.raw_polynomials(m,t,degrees,coeff,pilot)
    raw=raw[:,:,[0,16],:][:,:,:,[0,16]].copy()
    end=end[:,:,[0,16],:].copy()
    del m,t
    old_x=math.log(degree.QUIET)
    source=degree.ALPHA*(2*math.log(256/255)+256*(160*math.log(256)-32*math.log(255)))
    output=-degree.ALPHA*pilot.LIMIT*math.log(degree.Z)
    normalizers=[math.log(math.comb(4096,f)*math.comb(4096-f,s))
                 for f,s in ((128,0),(127,2))]
    results=[]

    def evaluate(x):
        dual.guard()
        if len(results)>=MAX_EVALUATIONS:
            raise RuntimeError("16-evaluation cap reached")
        logs=raw.copy(); logs[:,:,0,0]+=x-old_x
        a,b=dual.seed(logs),dual.seed(end,True)
        for _ in range(7):
            b,a=dual.convolution(a,b),dual.convolution(a,a)
        regions=[]; tails=[]
        for (f,s),norm in zip(((128,0),(127,2)),normalizers):
            regions.append((degree.ALPHA*(a[0][s,f]-norm),degree.ALPHA*a[1][s,f]))
            tails.append((degree.ALPHA*(b[0][s,f]-norm),degree.ALPHA*b[1][s,f]))
        paths,derivatives,_=dual.chain(regions,tails)
        means=derivatives/degree.ALPHA
        identity=float(means[0]+means[1]+means[2]-means[3])
        assert abs(identity-32767)<1e-6
        score=(source+output-degree.ALPHA*Q_TARGET*x+math.log(4096)+paths)/LN2
        value={"log_quiet_mark":x,"quiet_mark":math.exp(x),"log2_score":score,
            "births_k":float(means[0]),"quiet_q":float(means[1]),
            "active_incoming_a":float(means[2]),"final_active_probability":float(means[3]),
            "quiet_residual":float(means[1]-Q_TARGET),
            "slopes_bits_per_ln_raw_mark":{
                "birth_le_512":degree.ALPHA*float(means[0]-512)/LN2,
                "quiet_le_24575":degree.ALPHA*float(means[1]-Q_TARGET)/LN2,
                "floor_4a_le_61503":degree.ALPHA*float(4*means[2]-61503)/LN2},
            "count_identity":identity}
        results.append(value)
        print(json.dumps({"evaluation":len(results),"log_quiet_mark":x,
            "quiet_residual":value["quiet_residual"],"log2_score":score,
            "seconds":time.monotonic()-started}),flush=True)
        return value

    old=evaluate(old_x)
    assert abs(old["log2_score"]-11617.188069039783)<1e-6
    assert abs(old["quiet_q"]-25579.97680245402)<1e-5
    assert old["quiet_residual"]>0
    hi=old; lo=None
    # Deterministic exponential bracketing changes only the same scalar x.
    for step in (1.,2.,4.,8.):
        candidate=evaluate(old_x-step)
        if candidate["quiet_residual"]<=0:
            lo=candidate; break
        hi=candidate
    if lo is None:
        raise RuntimeError("root not bracketed within the predeclared steps; no wider search")
    while len(results)<MAX_EVALUATIONS:
        xl,xh=lo["log_quiet_mark"],hi["log_quiet_mark"]
        if xh-xl<1e-7 or min(abs(lo["quiet_residual"]),abs(hi["quiet_residual"]))<1e-6:
            break
        fl,fh=lo["quiet_residual"],hi["quiet_residual"]
        trial=(xl*fh-xh*fl)/(fh-fl)
        # Keep both sides progressing even when the derivative is highly curved.
        trial=max(xl+.1*(xh-xl),min(xh-.1*(xh-xl),trial))
        point=evaluate(trial)
        if point["quiet_residual"]<=0:
            lo=point
        else:
            hi=point
    ordered=sorted(results,key=lambda row:row["log_quiet_mark"])
    assert all(a["quiet_q"]<=b["quiet_q"]+1e-5 for a,b in zip(ordered,ordered[1:]))
    assert lo["quiet_residual"]<=0<=hi["quiet_residual"]
    best=min(results,key=lambda row:row["log2_score"])
    dual.guard()
    peak=pool.peak_working_set()
    assert peak<=128*(1<<20)
    print(json.dumps({"scope":"FLOAT one-dimensional stationarity of an exact regional0/full16 submajorant",
        "whole_code_certificate":False,"full_source_class_bound":False,"directed":False,
        "fixed_parameters":{"z":degree.Z,"alpha":degree.ALPHA,"birth_mark":1,"floor_mark":1},
        "evaluations":len(results),"maximum_evaluations":MAX_EVALUATIONS,
        "old":old,"root_bracket":{"lower":lo,"upper":hi,
            "log_mark_width":hi["log_quiet_mark"]-lo["log_quiet_mark"]},
        "best":best,"improvement_bits":old["log2_score"]-best["log2_score"],
        "tests":{**checked,"route_normalizers":route_checks,"old_control_reproduction":"passed",
            "tilted_quiet_monotonicity":"passed","count_identity_every_evaluation":"passed"},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "dual_pin":PIN,"exact_control_pin":dual.PIN,"seconds":time.monotonic()-started,
        "peak_working_set_bytes":peak,"generated_bytes":0},indent=2))


if __name__=="__main__":
    main()
