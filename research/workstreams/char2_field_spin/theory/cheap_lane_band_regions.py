"""Shared cyclic regional tables with separate actual-weight extraction.

one_table_many builds a positive cyclic256 table once at fixed rational
markers. Six updates T=R*T, R=R*R give H^64 and H^63*J, where J=M*T.
The local coefficients are already pair bounds: this module applies no
alpha, source factor, support count, or probability-law change.

For every requested actual weight r and valid full-pair count f, extract
the residue (f,r-2*f) and divide by that target's own exact denominator
binom(4096,f)*binom(4096-f,r-2*f)*u^f*w^(r-2*f). A normalized matrix at
another weight cannot supply this coefficient by scalar rescaling.
Positive residue aliases and upward integer rounding only enlarge the
coefficient. Regular and terminal exponents remain separate.

The caller owns the numerical deadline, 512MiB memory guard, and20GiB
free-space preflight. The function returns only normalized Arb entries
and compact metadata; it retains no large polynomial arrays. Import is
inert. Standalone execution runs small exact Fraction tests only, with
its own resource guards, stdout output, and zero generated files.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import gc
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
INTEGER_PIN="9dc583ae54caed0a78f71e409cc37aa398713b419f65538370f10c3f876f40bf"
REFERENCE_PIN="a2a2a6d784b5ae23f26b93e2c4e3caec87f7667598f8bd76d70ab3b49614e9a0"
MODULUS,SLOTS,MEMORY=256,4096,512*(1<<20)


def load_integer():
    path=Path(__file__).with_name("cheap_lane_integer_regions.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=INTEGER_PIN:
        raise RuntimeError("frozen integer arithmetic source mismatch")
    spec=importlib.util.spec_from_file_location("band_regions_integer",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def target_denominator(weight,full,marker,slots=SLOTS):
    """Exact normalizer for one actual target, never a residue normalizer."""
    single=weight-2*full
    if not 0<=full<=slots or not 0<=single<=slots-full:
        raise ValueError("invalid actual matching degree")
    u,w=marker
    if u<=0 or w<=0: raise ValueError("strictly positive rational markers required")
    return F(comb(slots,full)*comb(slots-full,single))*u**full*w**single


def one_table_many(v,degrees,item,marker,weights,integer,bits,guard,progress):
    """Return ({r:{f:(flatR4,flatV2)}}, report), with no global factors."""
    if not isinstance(weights,tuple) or not 1<=len(weights)<=17:
        raise ValueError("weights must be a tuple of1..17 distinct integers")
    if any(type(r) is not int or not 1<=r<=4096 for r in weights) or len(set(weights))!=len(weights):
        raise ValueError("distinct actual weights1..4096 required")
    if type(bits) is not int or bits not in (128,192,256):
        raise ValueError("regional mantissa bits must be128/192/256")
    center,lo,hi,kind=item
    if kind!="cyclic" or any(type(x) is not int for x in (center,lo,hi)) or not 0<=lo<=hi<=2048:
        raise ValueError("item must be(center,lo,hi,'cyclic'),0<=lo<=hi<=2048")
    if len(marker)!=2 or any(not isinstance(x,(F,int)) or isinstance(x,bool) for x in marker):
        raise ValueError("markers must be exact rational values")
    u,w=map(F,marker)
    if min(u,w)<=0: raise ValueError("strictly positive rational markers required")
    integer_path=Path(integer.__file__)
    if hashlib.sha256(integer_path.read_bytes()).hexdigest()!=INTEGER_PIN:
        raise RuntimeError("unrecognized integer arithmetic implementation")
    emit=(lambda phase,**values:None) if progress is None else progress
    guard(); up=[u**k for k in range(65)]; wp=[w**k for k in range(65)]
    scaled={}
    for (f,s),row in degrees.items():
        guard()
        if type(f) is not int or type(s) is not int or min(f,s)<0 or f+s>64 or len(row)!=6:
            raise ValueError("six pair channels with nonnegative total degree<=64 required")
        if any(not isinstance(x,(F,int)) or isinstance(x,bool) or x<0 for x in row):
            raise ValueError("exact nonnegative rational local caps required")
        scaled[f,s]=tuple(up[f]*wp[s]*x for x in row)
    if not scaled: raise ValueError("nonempty local pair coefficient table required")
    arithmetic=integer.IntegerArithmetic(255,255,cyclic_full=True,cyclic_single=True,
        bits=bits,guard=guard)
    # Reference DC is evaluated before upward quantization. Cyclic reduction
    # preserves DC; subsequent positive ceilings can only increase it.
    plain=[v.exact(sum((row[k] for row in scaled.values()),F(0))) for k in range(6)]
    plain_h=[[plain[0],plain[1]],[plain[2],plain[3]]]; plain_j=[[plain[4]],[plain[5]]]
    evaluated=(v.matpow(plain_h,64),v.matmul(v.matpow(plain_h,63),plain_j))
    rr,tt=arithmetic.pack(scaled)
    emit("band_integer_table_packed",center=center,first_full=lo,last_full=hi,
        weights=weights,matrix_exponent=rr.exponent,terminal_exponent=tt.exponent,
        arithmetic=arithmetic.report())
    del scaled,up,wp,plain,plain_h,plain_j
    for step in range(6):
        guard(); tt=arithmetic.multiply(rr,tt); rr=arithmetic.multiply(rr,rr)
        emit("band_integer_pair_doubling_complete",center=center,doubling=step+1,
            paired_factors=1<<(step+1),matrix_exponent=rr.exponent,terminal_exponent=tt.exponent)
    dc_checks=0; dc_ratios=[]
    for matrix,ordinary in zip((rr,tt),evaluated):
        for i,row in enumerate(arithmetic.mass(matrix)):
            for j,value in enumerate(row):
                guard(); dc=v.exact(value)
                if not (dc>=ordinary[i][j] or dc.overlaps(ordinary[i][j])):
                    raise ArithmeticError("cyclic rounded DC failed upper comparison")
                dc_checks+=1
                if dc>0 and ordinary[i][j]>0:
                    dc_ratios.append(str((dc/ordinary[i][j]).log()/v.arb(2).log()))
    result={}; coverage=[]
    for r in weights:
        current={}; last=min(hi,r//2)
        for f in range(lo,last+1):
            guard(); s=r-2*f
            denominator=target_denominator(r,f,(u,w))
            current[f]=tuple(tuple(v.exact(x/denominator)
                for row in arithmetic.extract(matrix,f%MODULUS,s%MODULUS) for x in row)
                for matrix in (rr,tt))
            assert tuple(map(len,current[f]))==(4,2)
        result[r]=current
        coverage.append({"weight":r,"count":len(current),
            "first_full":min(current) if current else None,"last_full":max(current) if current else None,
            "minimum_single":r-2*last if current else None,
            "maximum_single":r-2*lo if current else None})
        emit("band_weight_extraction_complete",center=center,**coverage[-1])
    report={"scope":"one shared cyclic regional table, separate actual-weight targets",
        "center":center,"first_full":lo,"last_full":hi,"method":"cyclic",
        "weights":weights,"coverage":coverage,"extracted_targets":sum(map(len,result.values())),
        "u":str(u),"v":str(w),"full_residue":MODULUS,"singleton_residue":MODULUS,
        "pair_doublings":6,"regular_object":"H^64","terminal_object":"H^63 J",
        "regular_packets":128,"terminal_updates":127,"own_target_denominators":True,
        "arithmetic":arithmetic.report(),"cyclic_DC_upper_checks":dc_checks,
        "DC_vs_unrounded_log2_ratios":dc_ratios,"alpha_applied":False,
        "source_factor_applied":False,"support_count_applied":False,
        "whole_band_certificate":False,"generated_bytes":0,
        "integer_arithmetic_pin":INTEGER_PIN,"reference_one_table_pin":REFERENCE_PIN}
    del rr,tt,matrix,ordinary,evaluated,arithmetic
    gc.collect(); guard()
    return result,report


def exact_tests(guard=lambda:None):
    """Small rational extraction controls; no FLINT or regional job."""
    checks=0; marker=(F(2,3),F(3,5)); slots=4
    # Changing r changes the selected singleton coefficient, not just its
    # denominator. Reusing the r2 numerator fails at both neighboring targets.
    scalar={(1,0):F(2),(1,1):F(7),(1,2):F(11)}
    scaled={degree:x*marker[0]**degree[0]*marker[1]**degree[1]
            for degree,x in scalar.items()}
    old_numerator=scaled[1,0]
    for r in (2,3,4):
        guard(); s=r-2; den=target_denominator(r,1,marker,slots)
        normalized=scaled[1,s]/den
        assert normalized==scalar[1,s]/(comb(slots,1)*comb(slots-1,s)); checks+=1
        if r!=2: assert old_numerator/den!=normalized; checks+=1
    # Independent noncommuting matrix-polynomial multiplication. Fold only
    # after exact products, then check every actual neighboring extraction.
    h=(( {(0,0):F(1),(1,0):F(2),(0,1):F(3)}, {(0,0):F(2),(0,1):F(1)} ),
       ( {(1,0):F(1)}, {(0,0):F(3),(0,2):F(2)} ))
    j=(( {(0,0):F(1),(0,1):F(7)}, ),( {(0,0):F(5),(1,0):F(2)}, ))
    def multiply(a,b):
        out=[]
        for row in a:
            values=[]
            for column in range(len(b[0])):
                p=defaultdict(F)
                for k,entry in enumerate(row):
                    for (f,s),x in entry.items():
                        for (ff,ss),y in b[k][column].items(): p[f+ff,s+ss]+=x*y
                values.append(dict(p))
            out.append(tuple(values))
        return tuple(out)
    hh,hj=multiply(h,h),multiply(h,j)
    assert hj[0][0] != multiply(((j[0][0],j[1][0]),),h)[0][0]; checks+=1
    for matrix in (hh,hj):
        for row in matrix:
            for coefficients in row:
                guard(); folded=defaultdict(F)
                for (f,s),x in coefficients.items():
                    folded[f%3,s%3]+=x*marker[0]**f*marker[1]**s
                assert sum(folded.values())==sum((x*marker[0]**f*marker[1]**s
                    for (f,s),x in coefficients.items()),F(0)); checks+=1
                for r in (1,2,3,4):
                    for f in range(r//2+1):
                        s=r-2*f; den=target_denominator(r,f,marker,slots)
                        normalized=folded[f%3,s%3]/den
                        exact=coefficients.get((f,s),F(0))/(comb(slots,f)*comb(slots-f,s))
                        assert normalized>=exact; checks+=1
    # Valid-domain clipping must not invent negative singleton targets.
    for r in (1,2047,2048,2049,4096):
        for lo,hi in ((0,0),(1,4),(1000,1024),(2048,2048)):
            values=tuple(range(lo,min(hi,r//2)+1))
            assert all(0<=r-2*f<=SLOTS-f for f in values); checks+=1
    return {"exact_fraction_checks":checks,"neighbor_requires_fresh_singleton_extraction":True,
        "cyclic_alias_domination":True,"noncommuting_terminal_control":True,
        "valid_weight_degree_clipping":True,"no_numerical_module_imported":True}


def main():
    parser=argparse.ArgumentParser(description="Small exact extraction tests only; no regional computation")
    parser.add_argument("--seconds",type=int,choices=(30,180),default=30)
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    integer=load_integer(); peak=integer.peak_memory_reader(); started=time.monotonic()
    def guard():
        if time.monotonic()-started>args.seconds: raise TimeoutError("band extraction test deadline")
        if peak()>MEMORY: raise MemoryError("512MiB band extraction test budget exceeded")
    tests=exact_tests(guard); guard()
    print(json.dumps({"scope":"small exact extraction tests only","tests":tests,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
