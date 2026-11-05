"""FLOAT concavity correction for the two selected actual r1024 profiles.

For raw two-packet path contributions x,y and 0<alpha<=1,
 (x+y)^alpha <= x^alpha+y^alpha-(2-2^alpha)min(x^alpha,y^alpha).
Write the powered paths as ab and cd. For each fixed positive rho,
 min(rho*a,c)*min(b/rho,d) <= min(ab,cd).
Take these minima for EACH occupancy before averaging its route weights.
Positive polynomial convolution then lower-bounds the integrated minimum.
The maximum over a finite rho bank remains a lower bound coefficientwise.
Subtract this correction from the old independently powered pair cap.
Optionally intersect the result with the frozen degree-Jensen pair cap.

This helper exports corrected_pair_logs(powered,...) for later return-cap
refinements. It changes no source factor, terminal, route law, or alpha
placement. The pair engine has64 regular pairs or63 plus M*T, not T*M.
All numerical outputs are FLOAT diagnostics, not directed certificates.
A directed implementation must subtract a LOWER bound on the correction.
Default rho bank: powers of2 from-16 through16 in steps of4. An explicit
--rho-radius may extend it to32. No optimizer or automatic widening.
180s/256MiB/20GiB; stdout only, no generated files or frozen edits.
"""

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import combinations_with_replacement, product
import json
import math
from math import comb, isqrt
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode = True
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
PAIR_PIN = "6286943ea9ccfb6acc1f0470117e2510ad2625c80c0bd51cef57482f97183c09"
SECONDS, MEMORY, WEIGHT = 180, 256*(1 << 20), 1024
PROFILES = ((384,500),(384,384))
LN2 = math.log(2)


def load_pair():
    path = Path(__file__).with_name("cheap_lane_binary_pair_pooling.py")
    if hashlib.sha256(path.read_bytes()).hexdigest() != PAIR_PIN:
        raise RuntimeError("frozen pair helper hash mismatch")
    spec = importlib.util.spec_from_file_location("concavity_pair",path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def route_geometry(shapes,np,guard):
    """Exact dyadic orientation weights; no minima have yet been averaged."""
    indices = {occ:i for i,occ in enumerate(shapes)}
    half = []
    for f in range(17):
        for s in range(17-f):
            for left in range(s+1):
                a,b = f+left,f+s-left
                if a <= b:
                    half.append((f,s,a,b,comb(16,f)*comb(16-f,s)*comb(s,left)
                                 *(1 if a == b else 2)/2**s))
    fd,sd,oi,weights = [],[],[],[]
    for i,(f,s,a,b,w) in enumerate(half):
        guard()
        for j in range(i,len(half)):
            ff,ss,c,d,ww = half[j]
            fd.append(f+ff); sd.append(s+ss)
            oi.append(indices[tuple(sorted((a,b,c,d)))])
            weights.append(w*ww*(1 if i == j else 2))
    result = (np.asarray(fd),np.asarray(sd),np.asarray(oi),np.asarray(weights))
    masses = np.zeros((33,33))
    np.add.at(masses,(result[0],result[1]),result[3])
    for f in range(33):
        for s in range(33-f):
            exact = comb(32,f)*comb(32-f,s)
            assert abs(masses[f,s]-exact) <= 2e-13*exact
    return result


def route_average(values,geometry,np,guard):
    guard(); fd,sd,oi,weights = geometry
    out = np.zeros((33,33,values.shape[1]))
    np.add.at(out,(fd,sd),weights[:,None]*values[oi])
    return out


def cross_convolution(left,right,np,guard):
    """Log-polynomial outer product of two left and three right channels."""
    n = left.shape[0]-1
    assert right.shape[:2] == left.shape[:2]
    out = np.full((2*n+1,2*n+1,6),-math.inf)
    for f in range(n+1):
        guard()
        for s in range(n+1-f):
            a = left[f,s]
            terms = np.stack((a[0]+right[:,:,0],a[0]+right[:,:,1],
                a[1]+right[:,:,0],a[1]+right[:,:,1],
                a[0]+right[:,:,2],a[1]+right[:,:,2]),axis=2)
            target = out[f:f+n+1,s:s+n+1]
            np.logaddexp(target,terms,out=target)
    return out


def subtract_correction(old,lower,alpha,np):
    """FLOAT subtraction only; its ratio is at most(2-2^alpha)/2."""
    result = old.copy(); finite = np.isfinite(old)
    if np.any(np.isfinite(lower) & ~finite):
        raise ArithmeticError("correction at a zero old coefficient")
    ratio = np.exp(lower[finite]-old[finite])
    if np.any(ratio > .5+2e-11):
        raise ArithmeticError("integrated minimum exceeds half the old pair")
    c = 2-math.exp2(alpha)
    if not 0 <= c < 1: raise ValueError("require0<alpha<=1")
    if c:
        result[finite] += np.log1p(-c*ratio)
    assert np.all(result <= old)
    return result


def corrected_pair_logs(powered,alpha,exponents,pair,old,np,guard,raw=None):
    """Reusable constructor. Mapping values are six powered local caps.

    Keys are sorted four-block occupancies. raw, if supplied, has the same
    keys and UNPOWERED caps; it enables the independent Jensen minimum.
    Returns log-polynomial tables plus source-independent diagnostics.
    """
    shapes = sorted(powered)
    if raw is not None and set(raw) != set(powered):
        raise ValueError("raw/powered occupancy mismatch")
    values = np.asarray([[float(x.mid()) if hasattr(x,"mid") else float(x)
                          for x in powered[occ]] for occ in shapes])
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ArithmeticError("invalid powered local caps")
    geometry = route_geometry(shapes,np,guard)
    single = old.float_coefficients(powered,np,guard)
    control = route_average(values,geometry,np,guard)
    assert np.allclose(single,control,rtol=3e-13,atol=0)
    old_pair = pair.pair_convolution(pair.as_logs(single,np),np,guard)
    lower = np.full(old_pair.shape,-math.inf)
    for e in exponents:
        guard(); rho = math.ldexp(1.0,e)
        mins = np.column_stack((np.minimum(rho*values[:,0],values[:,1]),
            np.minimum(rho*values[:,2],values[:,3]),
            np.minimum(values[:,0]/rho,values[:,2]),
            np.minimum(values[:,1]/rho,values[:,3]),
            np.minimum(values[:,4]/rho,values[:,5])))
        positive = np.column_stack(((values[:,0]>0)&(values[:,1]>0),
            (values[:,2]>0)&(values[:,3]>0),(values[:,0]>0)&(values[:,2]>0),
            (values[:,1]>0)&(values[:,3]>0),(values[:,4]>0)&(values[:,5]>0)))
        if not np.isfinite(mins).all() or np.any(positive & (mins == 0)):
            raise ArithmeticError("positive local rho-minimum under/overflow")
        coeff = pair.as_logs(route_average(mins,geometry,np,guard),np)
        candidate = cross_convolution(coeff[:,:,:2],coeff[:,:,2:],np,guard)
        np.maximum(lower,candidate,out=lower)
    corrected = subtract_correction(old_pair,lower,alpha,np)
    tables = {"old_pair":old_pair,"concavity_pair":corrected}
    if raw is not None:
        raw_coeff = old.float_coefficients(raw,np,guard)
        raw_pair = pair.pair_convolution(pair.as_logs(raw_coeff,np),np,guard)
        grouped = pair.degree_power(raw_pair,alpha,64,np)
        tables["old_Jensen_minimum"] = np.minimum(old_pair,grouped)
        tables["concavity_Jensen_minimum"] = np.minimum(corrected,grouped)
    finite = np.isfinite(old_pair)
    report = {"rho_exponents":list(exponents),"rho_count":len(exponents),
        "finite_pair_entries":int(np.count_nonzero(finite)),
        "improved_entries":int(np.count_nonzero(corrected < old_pair-1e-12)),
        "largest_correction_gain_bits_FLOAT":float(np.max(old_pair[finite]-corrected[finite])/LN2),
        "correction_constant_FLOAT":2-math.exp2(alpha),"route_records":len(geometry[0])}
    return tables,single,report


def sqrt_interval(x):
    scale = 1 << 80
    k = isqrt(x.numerator*scale*scale//x.denominator)
    lo = F(k,scale)
    return lo,lo if lo*lo == x else F(k+1,scale)


def exact_tests(pair,np,guard):
    """Exact scalar square checks and noncommuting weighted route toys."""
    checks = 0
    for a,b in product((F(0),F(1,7),F(2,3),F(5,2)),repeat=2):
        hi,lo = max(a,b),min(a,b)
        # At alpha1/2, the squared RHS minus a²+b² is
        # 2*lo*(hi-lo)*(sqrt(2)-1), hence nonnegative exactly.
        assert 2*lo*(hi-lo) >= 0
        for rho in (F(1,16),F(1,2),F(1),F(4),F(16)):
            for c,d in ((F(1,5),F(3,7)),(F(0),F(2)),(F(5),F(1,11))):
                witness = min(rho*a,c)*min(b/rho,d)
                assert 0 <= witness <= min(a*b,c*d)
                checks += 1
    # Values below are rational POWERED entries at alpha1/2; raw entries
    # are their squares. Each category is a distinct route, with its weight.
    atoms = (((0,0),F(2),tuple(map(F,(1,2,3,5,7,11)))),
             ((1,0),F(3),tuple(map(F,(2,3,5,7,11,13)))),
             ((0,1),F(1),tuple(map(F,(3,5,7,11,13,17)))))
    old_exact,raw_exact,route_mass = defaultdict(lambda:[F(0)]*6),defaultdict(lambda:[F(0)]*6),defaultdict(F)
    target_upper,point_min = defaultdict(lambda:[F(0)]*6),defaultdict(lambda:[F(0)]*6)
    single = np.zeros((2,2,6)); raw = np.zeros_like(single)
    paths = ((0,0,1,2),(0,1,1,3),(2,0,3,2),(2,1,3,3),(0,4,1,5),(2,4,3,5))
    for (f,s),weight,row in atoms:
        single[f,s] += float(weight)*np.asarray(list(map(float,row)))
        raw[f,s] += float(weight)*np.asarray([float(x*x) for x in row])
    for ca,cb in product(atoms,repeat=2):
        da,wa,a = ca; db,wb,b = cb; degree = (da[0]+db[0],da[1]+db[1]); w = wa*wb
        route_mass[degree] += w
        for j,(ai,bi,ci,di) in enumerate(paths):
            x,y = a[ai]*b[bi],a[ci]*b[di]
            old_exact[degree][j] += w*(x+y)
            raw_exact[degree][j] += w*(x*x+y*y)
            point_min[degree][j] += w*min(x,y)
            target_upper[degree][j] += w*sqrt_interval(x*x+y*y)[1]
    old_logs = pair.pair_convolution(pair.as_logs(single,np),np,guard)
    correction = np.full(old_logs.shape,-math.inf)
    best_exact = defaultdict(lambda:[F(0)]*6)
    for rho in (F(1,16),F(1,2),F(1),F(4),F(16)):
        left = np.zeros((2,2,2)); right = np.zeros((2,2,3))
        for (f,s),w,row in atoms:
            left[f,s] += float(w)*np.asarray(list(map(float,(min(rho*row[0],row[1]),min(rho*row[2],row[3])))))
            right[f,s] += float(w)*np.asarray(list(map(float,(min(row[0]/rho,row[2]),min(row[1]/rho,row[3]),min(row[4]/rho,row[5])))))
        actual = cross_convolution(pair.as_logs(left,np),pair.as_logs(right,np),np,guard)
        np.maximum(correction,actual,out=correction)
        enumerated = defaultdict(lambda:[F(0)]*6)
        for ca,cb in product(atoms,repeat=2):
            da,wa,a = ca; db,wb,b = cb; degree=(da[0]+db[0],da[1]+db[1])
            for j,(ai,bi,ci,di) in enumerate(paths):
                enumerated[degree][j] += wa*wb*min(rho*a[ai],a[ci])*min(b[bi]/rho,b[di])
        for degree,row in enumerated.items():
            for j,x in enumerate(row):
                assert abs(math.exp(actual[degree][j])-float(x)) <= 3e-12*max(1,float(x))
                best_exact[degree][j] = max(best_exact[degree][j],x)
                checks += 1
    corrected = subtract_correction(old_logs,correction,.5,np)
    assert np.array_equal(subtract_correction(old_logs,correction,1.0,np),old_logs)
    sqrt2lo,_ = sqrt_interval(F(2))
    for degree,row in old_exact.items():
        for j,x in enumerate(row):
            low = best_exact[degree][j]
            assert 0 <= low <= point_min[degree][j] <= x/2
            # These toy entries give strict inequalities, permitting exact
            # rational root enclosures to verify BOTH conditional caps.
            corrected_lower = x-2*low+sqrt2lo*low
            assert target_upper[degree][j] <= corrected_lower
            jensen_upper = sqrt_interval(route_mass[degree]*raw_exact[degree][j])[1]
            assert float(target_upper[degree][j]) <= float(jensen_upper)+1e-12
            assert math.exp(corrected[degree][j])+1e-9 >= float(target_upper[degree][j])
            assert math.exp(corrected[degree][j]) <= float(x)+1e-9
            checks += 5
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--z",type=F,default=F(304397809,625000000))
    parser.add_argument("--alpha",type=F,default=F(119,1280))
    parser.add_argument("--modulus",type=int,choices=(256,512),default=512)
    parser.add_argument("--rho-radius",type=int,choices=(4,8,12,16,20,24,28,32),default=16)
    args = parser.parse_args()
    if not 0 < args.z < 1 or not 0 < args.alpha <= 1:
        parser.error("require0<z<1 and0<alpha<=1")
    if max(n.bit_length() for x in (args.z,args.alpha) for n in (x.numerator,x.denominator)) > 64:
        parser.error("rational inputs limited to64 bits")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free < 20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+SECONDS
    pair=load_pair(); small=pair.load_return(); calibration=small.load_calibration(); broad=calibration.load_broad()
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"concavity_neighbor")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"concavity_old")
    anchor=old.load_gradient().load_anchor(); peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("180-second concavity deadline; no retry")
        if peak()>MEMORY: raise MemoryError("256MiB concavity budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,**values}),flush=True)
    import numpy as np
    from flint import ctx
    ctx.prec=256
    tests=exact_tests(pair,np,guard)+pair.exact_tests(np,guard)
    a,c,s,inverses,weights=small.actual_matrices()
    enumerators=small.pair_enumerators(c,s,inverses,guard)
    field_tests=small.exact_tests(c,s,inverses,weights,enumerators,guard)
    adapter=anchor.load_adapter(); core=adapter.load_verifier(); core.Z,core.ALPHA=args.z,args.alpha
    class LocalVerifier(core.Verifier):
        def guard(self): guard()
    verifier=LocalVerifier(2,deadline)
    _,quiet,old_d,_,_=core.rational_caps()
    new_d,block_report=small.refined_blocks(args.z,quiet,old_d,enumerators)
    original=verifier.local_caps(list(combinations_with_replacement(range(17),4)))
    raw,powered={},{}; exponent=verifier.exact(args.alpha)
    for occ,row in original.items():
        guard(); adapter.exact_point(row[2],small.return_cap(occ,quiet,old_d,core.KAPPA),verifier.fmpq)
        refined=list(row); refined[2]=verifier.exact(small.return_cap(occ,quiet,new_d,core.KAPPA))
        raw[occ]=tuple(refined); powered[occ]=tuple(verifier.positive(x)**exponent for x in refined)
        for x in raw[occ]+powered[occ]:
            value=float(x.mid())
            if not math.isfinite(value) or value<0 or (x>0 and value==0):
                raise ArithmeticError("positive local cap under/overflow")
    del original
    exponents=tuple(range(-args.rho_radius,args.rho_radius+1,4))
    tables,single,controls=corrected_pair_logs(powered,float(args.alpha),exponents,pair,old,np,guard,raw=raw)
    del raw,powered
    progress("correction_ready",controls=controls,tests=tests,field_tests=field_tests)
    replays={name:pair.PairReplay(logs,args.modulus,calibration,cells,fft,screen,np,guard) for name,logs in tables.items()}
    reference=calibration.Replay(single,WEIGHT,cells,fft,screen,np,guard)
    counts=screen.load("cheap_lane_fast_profile_counts.py",screen.COUNT_PIN,"concavity_counts"); counts.DEADLINE=deadline
    h,c_h,eta=anchor.source_factor(WEIGHT)
    logsource=float((exponent*(verifier.exact(eta).log()+WEIGHT*verifier.exact(core.BETA).log()-core.LIMIT*verifier.exact(args.z).log())).mid())
    cf,cs=F(1562175,990029),F(1140210,880327)
    for fulls in PROFILES:
        guard(); number=counts.profile_count(4096,WEIGHT,*fulls); assert number>0
        markers=tuple(broad.marker_curve(WEIGHT,f,cf,cs) for f in fulls)
        source_count=logsource+math.log(number)
        reports={name:replay.profile(fulls,markers,source_count) for name,replay in replays.items()}
        refs=[reference.marginal(f,m,args.modulus) for f,m in zip(fulls,markers)]
        control={"status":"unresolved"}
        if all(x["status"]=="resolved" for x in refs):
            score=(source_count+screen.paired_chains(refs[0]["region"],refs[1]["region"])[0])/LN2
            control={"status":"resolved","log2_bound_FLOAT":score}
            if reports["old_pair"]["residue"]["status"]=="resolved":
                assert abs(score-reports["old_pair"]["residue"]["log2_bound_FLOAT"])<2e-5
        for kind in ("ordinary","residue"):
            before,after=reports["old_pair"][kind],reports["concavity_pair"][kind]
            if before["status"]==after["status"]=="resolved":
                assert after["log2_bound_FLOAT"]<=before["log2_bound_FLOAT"]+2e-5
        progress("selected_profile_complete",fulls=fulls,markers=[[str(x) for x in m] for m in markers],
            log2_exact_count=math.log(number)/LN2,exact_count_sha256=hashlib.sha256(str(number).encode()).hexdigest(),
            frozen_single_packet_control=control,comparison=reports)
    guard()
    print(json.dumps({"scope":"two selected actual ordered r1024 profiles; FLOAT concavity correction",
        "directed":False,"whole_code_certificate":False,"complete_weight_bound":False,
        "z":str(args.z),"alpha":str(args.alpha),"modulus":args.modulus,"profiles":PROFILES,
        "source_factor":{"h":h,"c_h":str(c_h),"eta":str(eta),"log2_powered_source_FLOAT":logsource/LN2},
        "updates":32767,"terminal_outputs":1,"regular_pairs":64,"terminal_regular_pairs":63,
        "controls":controls,"tests":tests,"field_tests":field_tests,"refined_blocks":block_report,
        "pair_pin":PAIR_PIN,"return_pin":pair.RETURN_PIN,"count_pin":screen.COUNT_PIN,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__ == "__main__":
    main()
