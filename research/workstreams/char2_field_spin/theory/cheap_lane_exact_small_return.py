"""Genuine one/two-input return refinement; selected-profile FLOAT replay.

In one actual AES-GF256 emission block, A has rows
(1,a,a(a+1)/6,a^2(a+1)/6), a=0..15, C=A^t, and S=I+AC.
For one or two allowed input points I, C_I is injective. Hence the
nonzero-syndrome moment is exactly q^-|I| sum_{x!=0} z^wt(S_I x).
The one-column weights are15 at points0,2,3 and16 elsewhere.

For two columns, each of the16 nonzero row forms has a projective kernel
line. If their multiplicities are n_l and m lines occur, the exact moment
is (q-1)/q^2 [(q+1-m)z^16 + sum_l z^(16-n_l)]. Compute all120 pairs,
then maximize this POSITIVE polynomial at the fixed z. This uniform bound
does not assume independent point choices in paired blocks.

Intersect these genuine moments with the frozen D1/D2 bounds. Keep every
quiet Q_l unchanged. The positive quiet/nonquiet recurrence then gives a
new raw return E=kappa*q^-16*R under the SAME actual-mask cone. Change only
E, then raise local caps to alpha before route averaging. No extra source,
gamma, probability, or projective factor is introduced.

Exact field/rank and rational moment tests are separate from the numerical
profile replay. The latter uses the frozen FLOAT residue512 engine, own-f
rational markers and the same indexed source count/source factor before
and after refinement. Unresolved residues are reported, never clipped or
replaced. Only the two named r1024 profiles are considered; no complete
weight or code certificate.120s/256MiB/20GiB; stdout; no generated files.
"""

import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations, combinations_with_replacement
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
CALIBRATION_PIN = "e9f90bdda919247df80e71e5d2555f24ef7ca12323d95b1df622b22556327c56"
Q, WEIGHT, MODULUS = 256, 1024, 512
PROFILES = ((384,500),(384,384))
SECONDS, MEMORY = 120, 256*(1 << 20)
LN2 = math.log(2)


def load_calibration():
    path=Path(__file__).with_name("cheap_lane_broad_profile_calibration.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=CALIBRATION_PIN:
        raise RuntimeError("frozen broad-profile calibration hash mismatch")
    spec=importlib.util.spec_from_file_location("small_return_calibration",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def gf_mul(a,b):
    value=0
    while b:
        if b&1: value^=a
        b>>=1; a<<=1
        if a&256: a^=0x11B
    return value


def gf_pow(a,n):
    result=1
    while n:
        if n&1: result=gf_mul(result,a)
        a=gf_mul(a,a); n>>=1
    return result


def gf_dot(x,y):
    result=0
    for a,b in zip(x,y): result^=gf_mul(a,b)
    return result


def gf_product(a,b):
    return [[gf_dot(row,col) for col in zip(*b)] for row in a]


def gf_rank(matrix,inverses):
    a=[list(row) for row in matrix]
    rank=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(rank,len(a)) if a[i][col]),None)
        if pivot is None: continue
        a[rank],a[pivot]=a[pivot],a[rank]
        scale=inverses[a[rank][col]]
        a[rank]=[gf_mul(x,scale) for x in a[rank]]
        for i in range(len(a)):
            if i!=rank and a[i][col]:
                scale=a[i][col]
                a[i]=[x^gf_mul(scale,y) for x,y in zip(a[i],a[rank])]
        rank+=1
        if rank==len(a): break
    return rank


def actual_matrices():
    """Construct and independently verify the exact implemented field basis."""
    inverses=[0]+[gf_pow(a,254) for a in range(1,Q)]
    assert all(gf_mul(a,inverses[a])==1 for a in range(1,Q))
    a=[]
    for point in range(16):
        b=gf_mul(gf_mul(point,point^1),inverses[6])
        a.append([1,point,b,gf_mul(point,b)])
    c=[list(row) for row in zip(*a)]
    assert gf_product(c,a)==[[0]*4 for _ in range(4)]
    s=gf_product(a,c)
    for i in range(16): s[i][i]^=1
    identity=[[int(i==j) for j in range(16)] for i in range(16)]
    assert gf_product(s,s)==identity and gf_rank(s,inverses)==16
    weights=[sum(row[j]!=0 for row in s) for j in range(16)]
    assert weights==[15,16,15,15]+[16]*12
    assert [(i,j) for i in range(16) for j in range(16) if s[i][j]==0]==[(0,0),(2,3),(3,2)]
    return a,c,s,inverses,weights


def pair_enumerators(c,s,inverses,guard):
    """Return exact projective weight histograms, not scalar-vector samples."""
    records=[]
    for first,second in combinations(range(16),2):
        guard()
        assert gf_rank([[row[first],row[second]] for row in c],inverses)==2
        kernels=Counter()
        for row in s:
            left,right=row[first],row[second]
            assert left or right
            # label t denotes span(1,t); Q denotes span(0,1).
            label=gf_mul(left,inverses[right]) if right else Q
            kernels[label]+=1
        histogram=[0]*17
        histogram[16]=Q+1-len(kernels)
        for multiplicity in kernels.values(): histogram[16-multiplicity]+=1
        assert sum(kernels.values())==16 and sum(histogram)==Q+1
        assert histogram[0]==0  # S_I is injective.
        records.append(((first,second),tuple(histogram)))
    assert len(records)==120
    return records


def polynomial(histogram,z):
    return sum((F(number)*z**weight for weight,number in enumerate(histogram) if number),F(0))


def exact_tests(c,s,inverses,weights,records,guard):
    """Direct projective directions independently replay selected histograms."""
    pairs=dict(records); checked=0
    for first,second in ((0,1),(0,2),(2,3),(4,15),(8,9)):
        guard(); direct=[0]*17
        for t in range(Q):
            y=[row[first]^gf_mul(t,row[second]) for row in s]
            direct[sum(value!=0 for value in y)]+=1
        direct[sum(row[second]!=0 for row in s)]+=1
        assert tuple(direct)==pairs[first,second]
        # C_I has no nonzero kernel; every direction is a nonzero syndrome.
        for x0,x1 in ((1,0),(0,1),(1,1),(1,173)):
            assert any(gf_mul(row[first],x0)^gf_mul(row[second],x1) for row in c)
        for z in (F(1,3),F(2,3),F(1)):
            actual=F(Q-1,Q**2)*polynomial(direct,z)
            assert 0<=actual<1 and actual<=((z+(1-z)/Q)**2-F(1,Q**2))
            assert actual==F(Q-1,Q**2)*polynomial(pairs[first,second],z)
            checked+=1
    for z in (F(1,3),F(2,3),F(1)):
        assert max(F(Q-1,Q)*z**weight for weight in weights)==F(Q-1,Q)*z**15
        checked+=1
    # All remaining257 projective lines represent q-1 nonzero scalar vectors.
    assert (Q+1)*(Q-1)==Q**2-1
    return {"AES_inverse_checks":255,"CA_zero":True,"S_squared_identity":True,
            "S_rank":16,"column_weights":weights,"independent_C_pairs":120,
            "direct_projective_directions":5*(Q+1),"rational_moment_checks":checked}


def refined_blocks(z,quiet,old,records):
    one=F(Q-1,Q)*z**15
    pair_values=[(pair,F(Q-1,Q**2)*polynomial(histogram,z)) for pair,histogram in records]
    two=max(value for _,value in pair_values)
    result=list(old); result[1]=min(old[1],one); result[2]=min(old[2],two)
    assert quiet[1]==F(1,Q) and quiet[2]==F(1,Q**2)
    assert all(0<=new<=before for new,before in zip(result,old))
    return result,{"old_D1":str(old[1]),"new_D1":str(result[1]),
        "old_D2":str(old[2]),"new_D2":str(result[2]),
        "D2_maximizing_pairs":[pair for pair,value in pair_values if value==two],
        "D2_distinct_histograms":len({histogram for _,histogram in records})}


def return_cap(occ,quiet,nonquiet,kappa):
    qprod,active=F(1),F(0)
    for ell in occ:
        active,qprod=active*(quiet[ell]+nonquiet[ell])+qprod*nonquiet[ell],qprod*quiet[ell]
    return kappa*active/Q**16


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z",type=F,default=F(304397809,625000000))
    parser.add_argument("--alpha",type=F,default=F(119,1280))
    args=parser.parse_args(); z,alpha=args.z,args.alpha
    if not 0<z<1 or not 0<alpha<=1:
        parser.error("require0<z<1,0<alpha<=1")
    if max(n.bit_length() for x in (z,alpha) for n in (x.numerator,x.denominator))>64:
        parser.error("rational witnesses limited to64bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+SECONDS
    calibration=load_calibration(); broad=calibration.load_broad()
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"return_neighbor")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"return_old_calibration")
    gradient=old.load_gradient(); anchor=gradient.load_anchor(); peak=anchor.peak_memory_reader()

    def guard():
        if time.monotonic()>deadline: raise TimeoutError("120-second return-refinement deadline; no automatic retry")
        if peak()>MEMORY: raise MemoryError("256MiB return-refinement budget exceeded")

    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,**values}),flush=True)

    a,c,s,inverses,weights=actual_matrices()
    records=pair_enumerators(c,s,inverses,guard)
    tests=exact_tests(c,s,inverses,weights,records,guard)
    tests["matrix_sha256"]=hashlib.sha256(bytes(value for row in s for value in row)).hexdigest()
    import numpy as np
    from flint import ctx
    ctx.prec=256
    tests["FFT_exact_toy_coefficients"]=fft.exact_toy(np,guard)
    adapter=anchor.load_adapter(); core=adapter.load_verifier()
    core.Z,core.ALPHA=z,alpha

    class LocalVerifier(core.Verifier):
        def guard(self): guard()

    verifier=LocalVerifier(2,deadline)
    _,quiet,old_d,_,_=core.rational_caps()
    new_d,block_report=refined_blocks(z,quiet,old_d,records)
    shapes=list(combinations_with_replacement(range(17),4))
    caps=verifier.local_caps(shapes); exponent=verifier.exact(alpha)
    powered={occ:tuple(verifier.positive(value)**exponent for value in row) for occ,row in caps.items()}
    for row in powered.values():
        for value in row:
            approximate=float(value.mid())
            if not math.isfinite(approximate) or approximate<0 or (value>0 and approximate==0):
                raise ArithmeticError("positive powered cap not representable in FLOAT")
    baseline=old.float_coefficients(powered,np,guard)
    changes=0; largest_ratio=F(1)
    for occ in shapes:
        guard()
        before=return_cap(occ,quiet,old_d,core.KAPPA)
        after=return_cap(occ,quiet,new_d,core.KAPPA)
        assert 0<=after<=before
        adapter.exact_point(caps[occ][2],before,verifier.fmpq)
        if after<before:
            changes+=1
            if after: largest_ratio=max(largest_ratio,before/after)
        row=list(powered[occ]); row[2]=verifier.positive(verifier.exact(after))**exponent
        if after and (not math.isfinite(float(row[2].mid())) or float(row[2].mid())<=0):
            raise ArithmeticError("refined powered return is not representable in FLOAT")
        powered[occ]=tuple(row)
    del caps
    refined=old.float_coefficients(powered,np,guard); del powered
    other=[0,1,3,4,5]
    assert np.array_equal(baseline[:,:,other],refined[:,:,other])
    assert np.all(refined[:,:,2]<=baseline[:,:,2]*(1+2e-13))
    replays={"baseline":calibration.Replay(baseline,WEIGHT,cells,fft,screen,np,guard),
             "refined":calibration.Replay(refined,WEIGHT,cells,fft,screen,np,guard)}
    counts=screen.load("cheap_lane_fast_profile_counts.py",screen.COUNT_PIN,"return_counts")
    counts.DEADLINE=deadline
    h,c_h,eta=anchor.source_factor(WEIGHT)
    log_source=float((exponent*(verifier.exact(eta).log()+WEIGHT*verifier.exact(core.BETA).log()
                               -core.LIMIT*verifier.exact(z).log())).mid())
    cf,cs=F(1562175,990029),F(1140210,880327)
    progress("exact_small_return_caps_ready",tests=tests,block_caps=block_report,
        changed_shapes=changes,max_raw_return_gain_bits=(math.log(largest_ratio.numerator)
                                                       -math.log(largest_ratio.denominator))/LN2)
    reports=[]
    for fulls in PROFILES:
        guard(); number=counts.profile_count(4096,WEIGHT,*fulls)
        assert number>0
        markers=tuple(broad.marker_curve(WEIGHT,f,cf,cs) for f in fulls)
        report={"fulls":fulls,"singles":[WEIGHT-2*f for f in fulls],
            "markers":[[str(x) for x in marker] for marker in markers],
            "log2_exact_count":math.log(number)/LN2,"comparison":{}}
        for label,replay in replays.items():
            entries=[replay.marginal(f,marker,MODULUS) for f,marker in zip(fulls,markers)]
            if any(entry["status"]!="resolved" for entry in entries):
                value={"status":"unresolved","reasons":[x.get("reason") for x in entries]}
            else:
                chain=screen.paired_chains(entries[0]["region"],entries[1]["region"])[0]
                value={"status":"resolved","log2_bound_FLOAT":(log_source+math.log(number)+chain)/LN2,
                    "log2_conditional_chain_FLOAT":chain/LN2,
                    "FFT_controls":[x["FFT_controls"] for x in entries]}
            report["comparison"][label]=value
        before,after=report["comparison"]["baseline"],report["comparison"]["refined"]
        if before["status"]==after["status"]=="resolved":
            gain=before["log2_bound_FLOAT"]-after["log2_bound_FLOAT"]
            assert gain>=-2e-5
            report["genuine_return_refinement_gain_bits_FLOAT"]=gain
        reports.append(report); progress("selected_profile_complete",**report)
    guard()
    print(json.dumps({"scope":"genuine return-cap improvement; only two selected ordered profiles replayed",
        "directed_global_bound":False,"negative_is_not_a_certificate":True,
        "whole_weight_certificate":False,"whole_code_certificate":False,
        "weight":WEIGHT,"z":str(z),"alpha":str(alpha),"modulus":MODULUS,
        "changed_channel":"raw E only, before local alpha; Q,A,L,T,U unchanged",
        "fixed_cone_and_source":True,"quiet_subtraction":"none; positive nonzero moments already exclude x0",
        "uniform_pair_policy":"max over120 actual input pairs; no independent point-average claim",
        "source_factor":{"h":h,"c_h_exact":str(c_h),"eta_exact":str(eta),
            "powered_source_log2_FLOAT":log_source/LN2},
        "tests":tests,"block_caps":block_report,"changed_shapes":changes,"results":reports,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "calibration_pin":CALIBRATION_PIN,"core_pin":anchor.CORE_PIN,"count_pin":screen.COUNT_PIN,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
