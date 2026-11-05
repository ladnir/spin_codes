"""Directed verification of the size256 quiet-window support-class bound.

The hash-pinned screen supplies only witnesses: rational rank marks,
dyadic u,v,t obtained from exact float.as_integer_ratio(), cutoffs, and
multinomial mode counts. Optimization is not trusted. This wrapper checks
the mode exchange inequalities with integers, then recomputes every
polynomial, mixture factor, Chernoff tail, regional coefficient bound,
source factor, and profile union using Arb directed interval arithmetic.

The event is: some nonzero codeword with exactly256 active source groups
has, in every region, all regular steps with a nonzero endpoint state contained
in one non-wrapping44-packet window. The source set is arbitrary, not
necessarily consecutive. We select83 quiet packets outside each window
and omit each region's final packet. There are85 choices per region.
This is a restricted path/support class, not the whole-code distance event.

The E-defined quiet kernel has an E-scalar exact-support witness whenever
it has an extension-field exact-support witness: no union of at most256
proper GF256-linear subspaces covers a vector space. Scalar multiplication
then gives at least255 such witnesses. The quiet-only source factor is
(beta**256/255)**alpha, beta=256**160/255**32, alpha=-log(y)/log(256).

Each of16640 ordered paired profiles receives a complete bound at each
of three common marks y=3/4,2/3,3/5. The minimum is taken only between
complete valid bounds, then all profiles are summed with exact integer
source-count upper bounds. The verifier checks U<2**(-76) directly.

One48-second stdout-only run;192 or256-bit Arb;20GiB disk guard;estimated
peak below128MiB. No generated files, bytecode, encoder changes, or search.
"""

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
PIN="473428234351fd0813a4d5cf3e9bcf4db96179ef7e1cd0b884c232cd0e8a2a73"
DEADLINE=0.0


def guard():
    if time.monotonic()>DEADLINE:
        raise RuntimeError("bounded48-second verification exhausted")


def load_smoothing():
    path=Path(__file__).with_name("cheap_lane_quiet_profile_smoothing.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=PIN:
        raise RuntimeError("frozen smoothing helper hash mismatch")
    spec=importlib.util.spec_from_file_location("frozen_quiet_profile_smoothing",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ratio(value):
    if not math.isfinite(value) or value<0:
        raise ValueError("witness must be finite and nonnegative")
    return value.as_integer_ratio()


def rational(pair,arb):
    return arb(pair[0])/pair[1]


def mode_value(u_ratio,v_ratio,n,counts,arb):
    if len(counts)!=3 or any(type(k) is not int or k<0 for k in counts) or sum(counts)!=n:
        raise ValueError("invalid multinomial mode counts")
    denominator=math.lcm(u_ratio[1],v_ratio[1])
    weights=[u_ratio[0]*(denominator//u_ratio[1]),v_ratio[0]*(denominator//v_ratio[1]),denominator]
    if any(k>0 and w==0 for k,w in zip(counts,weights)):
        raise ValueError("positive count in a zero-probability category")
    if not all(counts[i]*weights[j]<=(counts[j]+1)*weights[i]
               for i in range(3) for j in range(3) if i!=j and counts[i]>0):
        raise ArithmeticError("claimed mode fails an exact exchange inequality")
    u,v=rational(u_ratio,arb),rational(v_ratio,arb)
    result=arb(math.comb(n,counts[0])*math.comb(n-counts[0],counts[1]))
    if counts[0]:
        result*=u**counts[0]
    if counts[1]:
        result*=v**counts[1]
    return result/(1+u+v)**n


def region_value(name,row,smoothing,arb,digest):
    fraction=Fraction(name); y=arb(fraction.numerator)/fraction.denominator
    ur,vr=ratio(row["u"]),ratio(row["v"])
    u,v=rational(ur,arb),rational(vr,arb)
    cutoff=row["smoothing"]["cutoff"]
    if cutoff not in smoothing.CUTOFFS:
        raise ValueError("cutoff outside fixed witness set")
    _,mode0=smoothing.exact_mode(row["u"],row["v"],1440)
    mode=row["smoothing"]["mode_counts"]
    a0=mode_value(ur,vr,1440,mode0,arb)
    am=mode_value(ur,vr,1440+cutoff,mode,arb)
    up=[u**i for i in range(17)]; vp=[v**i for i in range(17)]; yp=[y**i for i in range(9)]
    polynomial=arb(0)
    for f in range(17):
        for s in range(17-f):
            if (f and ur[0]==0) or (s and vr[0]==0):
                continue
            orientation=sum(arb(math.comb(s,j))*yp[min(4,f+j)+min(4,f+s-j)] for j in range(s+1))
            polynomial+=arb(math.comb(16,f)*math.comb(16-f,s))*orientation*up[f]*vp[s]/2**s
    b=1+y*y*u+y*v
    theta=b**16/polynomial
    p=y*y*(1+u+v)/b
    if not (theta>0 and theta<1 and p>0 and p<1):
        raise ArithmeticError("mixture parameter intervals are not strictly between zero and one")
    tr=None
    if cutoff:
        tr=ratio(row["smoothing"]["t"])
        if not 0<tr[0]<tr[1]:
            raise ValueError("Chernoff witness must satisfy0<t<1")
        t=rational(tr,arb)
        delta=(1-theta+theta*(1-p+p*t)**16)**166/t**(cutoff-1)
        if not delta<1:
            delta=arb(1)
        if not a0>am:
            raise ArithmeticError("directed mode monotonicity gap is not positive")
        cap=am+delta*(a0-am)
    else:
        cap=a0
    full,single=row["full"],row["singles"]
    if single!=256-2*full or not 0<=full<=128:
        raise ValueError("wrong source histogram")
    normalizer=math.comb(4096,full)*math.comb(4096-full,single)
    divisor=arb(normalizer)
    if full:
        if ur[0]==0:
            raise ValueError("zero full marker at positive degree")
        divisor*=u**full
    if single:
        if vr[0]==0:
            raise ValueError("zero singleton marker at positive degree")
        divisor*=v**single
    regional=polynomial**166*(1+u+v)**1440*cap/divisor
    if not regional.is_finite() or not regional>0:
        raise ArithmeticError("nonpositive or nonfinite regional interval")
    witness={"mark":[fraction.numerator,fraction.denominator],"full":full,"single":single,
             "u":ur,"v":vr,"t":tr,"cutoff":cutoff,"mode0":mode0,"mode":mode}
    digest.update(json.dumps(witness,sort_keys=True,separators=(",",":")).encode()+b"\n")
    return regional**128


def verify(precision):
    global DEADLINE
    if precision not in (192,256):
        raise ValueError("precision must be192 or256")
    free=shutil.disk_usage(Path(__file__).resolve().anchor).free
    if free<20*(1<<30):
        raise RuntimeError("requires20GiB free on the actual source volume")
    started=time.monotonic(); DEADLINE=started+48
    from flint import arb,ctx
    ctx.prec=precision
    smoothing=load_smoothing(); wrappers,marks,screen=smoothing.load_counts(DEADLINE)
    tests=smoothing.tests(); counter=wrappers.DeletionCounts(screen.load_counts())
    digest=hashlib.sha256(); data=[]
    q=arb(256); beta=q**160/arb(255)**32
    for name,y_float in marks.MARKS:
        fraction=Fraction(name); y=arb(fraction.numerator)/fraction.denominator
        rows=screen.terms(y_float); regional=[]
        for full in range(129):
            row=smoothing.strengthen(screen.optimize(rows,full),rows,y_float,screen)
            regional.append(region_value(name,row,smoothing,arb,digest))
            guard()
        alpha=-y.log()/q.log()
        source=(beta**256/255)**alpha*arb(85)**256
        data.append((name,source,regional))
    bound=arb(0); worst=arb(0); profile_count=0; usage={name:0 for name,_ in marks.MARKS}
    for first in range(129):
        for second in range(129):
            if first+second>255:
                continue
            count=counter.profile(4096,256,first,second)
            if count==0:
                continue
            best=None; selected=None
            for name,source,regional in data:
                candidate=arb(count)*source*regional[first]*regional[second]
                if best is None or candidate.upper()<best.upper():
                    best,selected=candidate,name
            if not best.is_finite() or not best>0:
                raise ArithmeticError("invalid complete-profile interval")
            bound+=best; usage[selected]+=1; profile_count+=1
            if best.upper()>worst.upper():
                worst=best
        guard()
    if profile_count!=16640:
        raise ArithmeticError("incomplete actual-profile coverage")
    passed=bool(bound<arb(2)**(-76))
    if not passed:
        raise ArithmeticError("directed bound did not prove76 bits")
    guard()
    return {"schema":"cheap-lane-quiet-window-all-profiles-v1",
        "scope":"arbitrary exact256-group supports; endpoint-live regular steps inside44-packet regional windows",
        "whole_code_certificate":False,"restricted_class_certificate":True,"directed":True,
        "precision_bits":precision,"target_margin_bits":76,"proved_below_target":passed,
        "source_size":256,"ordered_profiles":profile_count,"regional_witnesses":387,
        "regions":256,"quiet_packets_per_region":83,"window_choices_per_region":85,
        "marks":[name for name,_ in marks.MARKS],"mark_usage":usage,
        "log2_union_interval":str(bound.log()/arb(2).log()),
        "worst_profile_log2_interval":str(worst.log()/arb(2).log()),
        "selected_rational_witness_sha256":digest.hexdigest(),"exact_mode_checks":774,
        "tests":tests,"smoothing_helper_sha256":PIN,
        "count_wrapper_sha256":smoothing.PIN,"marks_wrapper_sha256":wrappers.PIN,
        "profile_screen_sha256":marks.PIN,"source_count_helper_sha256":screen.PIN,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python_executable":sys.executable,"command":"python -B cheap_lane_quiet_profile_verify.py --precision "+str(precision),
        "free_bytes_at_start":free,"seconds":time.monotonic()-started,"generated_bytes":0,
        "estimated_peak_bytes":128*(1<<20)}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--precision",type=int,choices=(192,256),default=192)
    print(json.dumps(verify(parser.parse_args().precision),indent=2))
