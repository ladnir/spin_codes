"""Selected FLOAT control of rho-grid and first-occupancy separation losses.

For powered left-row values a,c and right-column values b,d, put
L(rho)=min(rho*a,c)min(b/rho,d). Always L(rho)<=min(ab,cd), and
rho=c/a gives equality when a,c>0. A zero a or c has exact minimum0.
Partition FIRST-packet occupancies, separately for its two incoming rows,
by a fixed rho label. Convolve each restricted left-min polynomial with
the full right-min polynomial for that label, and SUM the disjoint
buckets. This is a lower correction bound, without route Jensen or
renormalization. The terminal uses right column(T0,U), never a regular
transition in place of the output-only terminal.

Compare frozen old9 global maxima, dense33 global maxima, and dense33
bucket sums. Take maxima of finished lower correction bounds before
subtracting(2-2^alpha); keep the same degree-Jensen ceiling. Buckets use
nearest powers2^-16..2^16 with ties lower and clipped endpoints. Assignment
uses exact rational binary64 representatives of POWERED caps; any fixed
assignment is valid, so no proof depends on these being exact ratios.
Inside the grid, nearest-left-ratio rounding loses at most sqrt(2) from
the pointwise ideal correction; a spacing4 grid can lose a factor4.
Max AFTER route averaging is an additional loss not removed by density.

The fixed actual profile is r2048,(768,768), z23/100,alpha13/128.
Build genuine residual-rank caps once, fit markers on old9 ONLY, then
replay all three tables at the same rounded markers and m256/m512.
Charge identical exact run512 mass throughout: selected-profile
diagnostic, NOT a run cover or complete-size certificate.

FLOAT only. Directed successors must subtract LOWER correction bounds
from the same coherent upper old table.180s/512MiB/20GiB, stdout only;
no numerical jobs on import, no generated files or frozen edits.
"""

from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from math import comb
import os
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
CONTROL_PIN="5b5e18d2b62b8dbbf5614b5b07c258b4aa1fe6e53a554513518906c931212518"
RANK_PIN="3efa28f205a01cafe47a8058ca4d3a69e065d5e22e75ea00a8fe79ed4a5dcf7c"
Z,ALPHA=F(23,100),F(13,128)
OLD_EXPONENTS=tuple(range(-16,17,4))
DENSE_EXPONENTS=tuple(range(-16,17))
SECONDS,MEMORY,LN2=180,512*(1<<20),math.log(2)
PATHS=((0,0,1,2),(0,1,1,3),(2,0,3,2),(2,1,3,3),(0,4,1,5),(2,4,3,5))


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def bucket_exponent(a,c,low=-16,high=16):
    """Exact nearest-log dyadic for positive rational c/a; tie goes lower."""
    if a==0 or c==0: return None
    if a<0 or c<0: raise ValueError("nonnegative powered atoms required")
    ratio=c/a; left,right=ratio.numerator**2,ratio.denominator**2
    lo,hi=low,high
    while lo<hi:
        middle=(lo+hi)//2; power=2*middle+1
        below=(left<=right*(1<<power)) if power>=0 else (left*(1<<(-power))<=right)
        if below: hi=middle
        else: lo=middle+1
    return lo


def exact_tests(pair,correction,np,guard):
    """Fraction witnesses and six-channel bucket convolution controls."""
    checks=0
    vals=(F(0),F(1,7),F(2,3),F(5,2))
    for a in vals:
        for c in vals:
            for b in vals:
                for d in vals:
                    ideal=min(a*b,c*d)
                    if a and c:
                        rho=c/a
                        assert min(rho*a,c)*min(b/rho,d)==ideal
                        e=bucket_exponent(a,c,-2,2); witness=min(F(2)**e*a,c)*min(b/F(2)**e,d)
                        assert 0<=witness<=ideal
                    else: assert ideal==0
                    checks+=1
    atoms=(((0,0),F(2),tuple(map(F,(1,4,3,2,5,7)))),
           ((1,0),F(3),tuple(map(F,(3,1,0,5,2,11)))),
           ((0,1),F(1),tuple(map(F,(2,3,5,0,13,3)))))
    expected=defaultdict(lambda:[F(0)]*6)
    ideal=defaultdict(lambda:[F(0)]*6); old=defaultdict(lambda:[F(0)]*6)
    exact_ratio=defaultdict(lambda:[F(0)]*6)
    buckets={degree:tuple(bucket_exponent(row[2*i],row[2*i+1],-2,2) for i in range(2))
             for degree,w,row in atoms}
    for da,wa,a in atoms:
        for db,wb,b in atoms:
            degree=(da[0]+db[0],da[1]+db[1]); weight=wa*wb
            for j,(ai,bi,ci,di) in enumerate(PATHS):
                target=min(a[ai]*b[bi],a[ci]*b[di])
                ideal[degree][j]+=weight*target
                old[degree][j]+=weight*(a[ai]*b[bi]+a[ci]*b[di])
                e=buckets[da][ai//2]
                if e is not None:
                    rho=F(2)**e
                    expected[degree][j]+=weight*min(rho*a[ai],a[ci])*min(b[bi]/rho,b[di])
                    exactrho=a[ci]/a[ai]
                    exact_ratio[degree][j]+=weight*min(exactrho*a[ai],a[ci])*min(b[bi]/exactrho,b[di])
    computed=np.full((3,3,6),-math.inf)
    for e in range(-2,3):
        left=np.zeros((2,2,2)); right=np.zeros((2,2,3)); rho=F(2)**e
        for (f,s),weight,row in atoms:
            for i in range(2):
                if buckets[f,s][i]==e: left[f,s,i]=float(weight*min(rho*row[2*i],row[2*i+1]))
            right[f,s]=list(map(float,(weight*min(row[0]/rho,row[2]),
                weight*min(row[1]/rho,row[3]),weight*min(row[4]/rho,row[5]))))
        lower=correction.cross_convolution(pair.as_logs(left,np),pair.as_logs(right,np),np,guard)
        np.logaddexp(computed,lower,out=computed)
    for degree,row in expected.items():
        for j,value in enumerate(row):
            assert value<=ideal[degree][j]<=old[degree][j]/2
            assert exact_ratio[degree][j]==ideal[degree][j]
            got=math.exp(float(computed[degree][j])) if math.isfinite(computed[degree][j]) else 0.0
            assert abs(got-float(value))<=3e-12*max(1,float(value))
            checks+=3
    # Exact nearest-ratio comparison, including ties at rational powers2.
    for e in range(-16,17):
        assert bucket_exponent(F(1),F(2)**e)==e; checks+=1
    return {"exact_fraction_and_bucket_convolution_checks":checks,
        "direct_exact_ratio_reproduces_integrated_minimum":True,
        "two_rows_have_distinct_partitions":True,"zero_left_atoms_omit_only_zero_minima":True}


def bucket_pair_logs(powered,alpha,pair,old,correction,np,guard,raw,progress):
    """Build three matched upper tables from one coherent cap map."""
    baseline,single,baseline_report=correction.corrected_pair_logs(powered,alpha,OLD_EXPONENTS,
        pair,old,np,guard,raw=raw)
    del single
    old_pair=baseline["old_pair"]
    jensen_ceiling=baseline["old_Jensen_minimum"]
    old9=baseline["concavity_Jensen_minimum"]
    shapes=sorted(powered)
    values=np.asarray([[float(x.mid()) if hasattr(x,"mid") else float(x)
                        for x in powered[occ]] for occ in shapes])
    if not np.isfinite(values).all() or np.any(values<0):
        raise ArithmeticError("invalid coherent powered table")
    assignments=np.full((len(shapes),2),-100,dtype=int)
    counts=[defaultdict(int),defaultdict(int)]; outside=[0,0]; omitted=[0,0]
    for rowindex,row in enumerate(values):
        guard()
        for i in range(2):
            a,c=(F(*float(row[2*i+j]).as_integer_ratio()) for j in range(2))
            e=bucket_exponent(a,c)
            if e is None: omitted[i]+=1; continue
            assignments[rowindex,i]=e; counts[i][e]+=1
            outside[i]+=c<a*F(2)**(-16) or c>a*F(2)**16
    assert all(sum(counts[i].values())+omitted[i]==len(shapes) for i in range(2))
    geometry=correction.route_geometry(shapes,np,guard)
    global_lower=np.full(old_pair.shape,-math.inf)
    bucket_lower=np.full(old_pair.shape,-math.inf)
    for e in DENSE_EXPONENTS:
        guard(); rho=math.ldexp(1.0,e)
        left=np.column_stack((np.minimum(rho*values[:,0],values[:,1]),
                              np.minimum(rho*values[:,2],values[:,3])))
        right=np.column_stack((np.minimum(values[:,0]/rho,values[:,2]),
            np.minimum(values[:,1]/rho,values[:,3]),np.minimum(values[:,4]/rho,values[:,5])))
        if not np.isfinite(left).all() or not np.isfinite(right).all():
            raise ArithmeticError("rho-scaled local minimum overflow")
        positive_left=np.column_stack(((values[:,0]>0)&(values[:,1]>0),(values[:,2]>0)&(values[:,3]>0)))
        positive_right=np.column_stack(((values[:,0]>0)&(values[:,2]>0),
            (values[:,1]>0)&(values[:,3]>0),(values[:,4]>0)&(values[:,5]>0)))
        if np.any(positive_left&(left==0)) or np.any(positive_right&(right==0)):
            raise ArithmeticError("positive rho minimum underflow")
        restricted=np.where(assignments==e,left,0.0)
        coefficients=pair.as_logs(correction.route_average(np.column_stack((left,right,restricted)),geometry,np,guard),np)
        candidate=correction.cross_convolution(coefficients[:,:,:2],coefficients[:,:,2:5],np,guard)
        np.maximum(global_lower,candidate,out=global_lower)
        candidate=correction.cross_convolution(coefficients[:,:,5:7],coefficients[:,:,2:5],np,guard)
        np.logaddexp(bucket_lower,candidate,out=bucket_lower)
        if (e+16)%4==0: progress("dense_rho_complete",exponent=e)
    # Every candidate has the SAME integrated-minimum target. The sum is
    # only over the disjoint first-occupancy buckets, never over all rho.
    dense=correction.subtract_correction(old_pair,global_lower,alpha,np)
    joint_lower=np.maximum(global_lower,bucket_lower)
    bucket=correction.subtract_correction(old_pair,joint_lower,alpha,np)
    dense=np.minimum(np.minimum(dense,jensen_ceiling),old9)
    bucket=np.minimum(np.minimum(bucket,jensen_ceiling),dense)
    finite=np.isfinite(old_pair)
    assert np.all(bucket<=dense) and np.all(dense<=old9)
    corrections={"old9_global_max":old9,"dense33_global_max":dense,"dense33_bucket_sum":bucket}
    report={"old":baseline_report,"dense_exponents":list(DENSE_EXPONENTS),
        "bucket_assignment":"nearest dyadic to exact rational binary64 powered c/a; ties lower; clipped endpoints",
        "bucket_counts_by_left_row":[dict(sorted(row.items())) for row in counts],
        "positive_ratios_outside_grid_by_row":outside,"zero_minimum_omissions_by_row":omitted,
        "all_left_occupancies_accounted_for":True,"six_MM_MT_channels":True,
        "degree_entry_improvements":{name:int(np.count_nonzero(table<old9-1e-12)) for name,table in corrections.items()},
        "maximum_local_gain_bits_from_old9_FLOAT":{name:float(np.max(old9[finite]-table[finite])/LN2)
            for name,table in corrections.items()},
        "old_Jensen_ceiling_unchanged":True,"final_lower_is_max_of_valid_lower_bounds":True}
    return corrections,report


def main():
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+SECONDS
    control=load("cheap_lane_pair_profile_calibration.py",CONTROL_PIN,"rho_selected_control")
    rank=load("cheap_lane_residual_rank_return.py",RANK_PIN,"rho_rank_return")
    sizes=control.load("cheap_lane_pair_sizes_screen.py",control.SIZES_PIN,"rho_sizes")
    calibration=control.load("cheap_lane_broad_profile_calibration.py",control.CALIBRATION_PIN,"rho_moments")
    bank=sizes.load("cheap_lane_pair_bank_screen.py",sizes.BANK_PIN,"rho_bank")
    broad=sizes.load("cheap_lane_broad_peeling_screen.py",sizes.BROAD_PIN,"rho_broad")
    composition=bank.load_composition()
    correction=composition.load("cheap_lane_pair_concavity.py",composition.CORRECTION_PIN,"rho_correction")
    triple=composition.load("cheap_lane_exact_triple_return.py",composition.TRIPLE_PIN,"rho_triple")
    pair=correction.load_pair()
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"rho_neighbor")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    old=screen.load("cheap_lane_cauchy_marker_calibration.py",screen.CALIBRATION_PIN,"rho_old")
    anchor=old.load_gradient().load_anchor(); peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("180-second rho control; no automatic retry")
        if peak()>MEMORY: raise MemoryError("512MiB rho-control budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,"peak_process_bytes":peak(),**values}),flush=True)
    import numpy as np
    from flint import ctx
    ctx.prec=256
    tests=exact_tests(pair,correction,np,guard)
    tests["residual_rank"]=rank.exact_tests(guard)
    tests["positive_moments"]=calibration.exact_tests(np,screen,old)
    tests["FFT_exact_toys"]=fft.exact_toy(np,guard)
    caps=rank.augment_caps(composition.build_caps(Z,ALPHA,triple,anchor,guard,deadline),Z,ALPHA,guard)
    v,core=caps["verifier"],caps["core"]
    tables,local_report=bucket_pair_logs(caps["powered"],float(ALPHA),pair,old,correction,np,guard,caps["raw"],progress)
    del caps["raw"],caps["powered"]
    r,f,s,c=control.WEIGHT,control.FULL,control.SINGLE,control.RUNS
    missing,c_h,eta=anchor.source_factor(r)
    log_source=float((v.exact(ALPHA)*(v.exact(eta).log()+r*v.exact(core.BETA).log()-core.LIMIT*v.exact(Z).log())).mid())
    numerator=8192*comb(r-1,c-1)*comb(8192-r-1,c-1); run_mass,rem=divmod(numerator,c); assert rem==0
    source_count=log_source+math.log(run_mass)
    constant=source_count-256*math.log(comb(4096,f)*comb(4096-f,s))
    objective=control.Objective(tables["old9_global_max"],constant,calibration,np,guard)
    own=sizes.markers(r,f)
    best,history,status=control.coordinate_fit(objective,np.asarray(list(map(broad.logq,own))),np,guard)
    fitted=tuple(max(F(1,1<<30),F(math.exp(float(x))).limit_denominator(1<<20)) for x in best["x"])
    progress("old9_markers_fixed",markers=list(map(str,fitted)),proposals=len(history),status=status)
    results=[]
    for name,table in tables.items():
        replays=[]
        for modulus in control.MODULI:
            row=control.replay(table,fitted,modulus,source_count,broad,cells,fft,screen,np,guard)
            row["modulus"]=modulus; replays.append(row)
            progress("rho_matched_replay",stage=name,**row)
        results.append({"stage":name,"replays":replays})
    guard()
    print(json.dumps({"scope":"three matched selected r2048 profile768,768 FLOAT bounds, each charged run512 mass",
        "directed":False,"complete_weight_bound":False,"run_class_cover":False,"whole_code_certificate":False,
        "parameters":{"weight":r,"fulls":[f,f],"singles":[s,s],"z":str(Z),"alpha":str(ALPHA),
            "markers":list(map(str,fitted)),"moduli":control.MODULI,"seconds":SECONDS,"memory_budget_bytes":MEMORY},
        "same_caps_source_count_alpha_markers_and_terminal":True,"no_optimization_after_changing_pair_bound":True,
        "results":results,"local_rho_report":local_report,"tests":tests,
        "fit_on_old9_only":{"status":status,"evaluations":len(history),"history":history},
        "source":{"log2_factor_FLOAT":log_source/LN2,"log2_run_mass":math.log(run_mass)/LN2,
            "h":missing,"c_h":str(c_h),"eta":str(eta)},
        "residual_rank_receipt":caps["residual_rank_receipt"],"field_tests":caps["tests"],
        "control_pin":CONTROL_PIN,"rank_pin":RANK_PIN,"correction_pin":composition.CORRECTION_PIN,
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
