"""Matched selected-profile control of finer rational left-ratio buckets.

Grid: 2^e*(1,5/4,3/2,7/4), e=-16..15, plus2^16. Each powered left
row(a,c) is assigned exactly once, using (c/a)^2 versus the exact product
of adjacent grid points; ties go lower and endpoints are clipped. Only
occupied buckets are convolved against the unrestricted right side.
Sum disjoint bucket contributions, subtract(2-2^alpha) times that lower
bound, then intersect the FINISHED upper cap with the frozen coarse cap.
Reassignment alone need not dominate the coarse bucket lower bound.

For later bank use fine_bucket_pair_logs(...,coarse_builder=...) returns
all coarse tables plus fine129_bucket_sum. It never changes source/count,
alpha placement or terminal. Both tables use identical raw/powered caps.

The standalone control delegates to the pinned generic selected driver,
default full957,z23/100,alpha3/32. It fits ONLY the coarse table. At those
same rounded markers, m256/m512 additionally replay the fine table. The
ordinary/inherited/own comparisons and final selected report remain the
coarse controls, with explicit nested fine comparisons. This is not a
complete-weight bound, parameter search or certificate. Private hooks and
argv are restored.180s/512MiB/20GiB, stdout only, no generated files.
"""

from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode=True
SELECTED_PIN="97e237880824fec84527a73d400a811006ed7a492aeda774004bee0e54ade633"
GRID=tuple(F(2)**e*m for e in range(-16,16) for m in (F(1),F(5,4),F(3,2),F(7,4)))+(F(1<<16),)
LN2=math.log(2)
PATHS=((0,0,1,2),(0,1,1,3),(2,0,3,2),(2,1,3,3),(0,4,1,5),(2,4,3,5))


def load_selected():
    path=Path(__file__).with_name("cheap_lane_rho_profile_calibration.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=SELECTED_PIN:
        raise RuntimeError("frozen selected calibration mismatch")
    spec=importlib.util.spec_from_file_location("fine_rho_selected",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def nearest_bucket(a,c,grid=GRID):
    """Exact nearest-log grid index for rational powered c/a."""
    if a<0 or c<0: raise ValueError("nonnegative powered atoms required")
    if a==0 or c==0: return None
    ratio=c/a; square=ratio*ratio; low,high=0,len(grid)-1
    while low<high:
        middle=(low+high)//2
        if square<=grid[middle]*grid[middle+1]: high=middle
        else: low=middle+1
    return low


def exact_tests(pair,correction,np,guard):
    checks=0
    assert len(GRID)==129 and tuple(sorted(set(GRID)))==GRID
    assert max(b/a for a,b in zip(GRID,GRID[1:]))==F(5,4)
    for i,rho in enumerate(GRID):
        assert nearest_bucket(F(1),rho)==i; checks+=1
    for ratio,expected in ((F(1,8),0),(F(2),0),(F(8),1),(F(128),2)):
        assert nearest_bucket(F(1),ratio,(F(1),F(4),F(16)))==expected; checks+=1
    assert nearest_bucket(F(0),F(1)) is None and nearest_bucket(F(1),F(0)) is None
    atoms=(((0,0),F(2),(F(1),F(5,4),F(2),F(7,2),F(3),F(7))),
           ((1,0),F(3),(F(3),F(1),F(0),F(5),F(2),F(11))),
           ((0,1),F(1),(F(2),F(3),F(5),F(0),F(13),F(3))))
    assigned={degree:tuple(nearest_bucket(row[2*i],row[2*i+1]) for i in range(2))
              for degree,w,row in atoms}
    expected=defaultdict(lambda:[F(0)]*6)
    ideal=defaultdict(lambda:[F(0)]*6); old=defaultdict(lambda:[F(0)]*6)
    for da,wa,a in atoms:
        for db,wb,b in atoms:
            degree=da[0]+db[0],da[1]+db[1]
            for j,(ai,bi,ci,di) in enumerate(PATHS):
                ideal[degree][j]+=wa*wb*min(a[ai]*b[bi],a[ci]*b[di])
                old[degree][j]+=wa*wb*(a[ai]*b[bi]+a[ci]*b[di])
                index=assigned[da][ai//2]
                if index is not None:
                    rho=GRID[index]
                    witness=min(rho*a[ai],a[ci])*min(b[bi]/rho,b[di])
                    expected[degree][j]+=wa*wb*witness
                    if GRID[0]<=a[ci]/a[ai]<=GRID[-1]:
                        assert 5*witness*witness>=4*min(a[ai]*b[bi],a[ci]*b[di])**2
                        checks+=1
    occupied=sorted({i for row in assigned.values() for i in row if i is not None})
    computed=np.full((3,3,6),-math.inf)
    for index in occupied:
        guard(); rho=GRID[index]; left=np.zeros((2,2,2)); right=np.zeros((2,2,3))
        for (f,s),weight,row in atoms:
            for i in range(2):
                if assigned[f,s][i]==index:
                    left[f,s,i]=float(weight*min(rho*row[2*i],row[2*i+1]))
            right[f,s]=list(map(float,(weight*min(row[0]/rho,row[2]),
                weight*min(row[1]/rho,row[3]),weight*min(row[4]/rho,row[5]))))
        lower=correction.cross_convolution(pair.as_logs(left,np),pair.as_logs(right,np),np,guard)
        np.logaddexp(computed,lower,out=computed)
    for degree,row in expected.items():
        for j,value in enumerate(row):
            assert value<=ideal[degree][j]<=old[degree][j]/2
            actual=math.exp(float(computed[degree][j])) if math.isfinite(computed[degree][j]) else 0.0
            assert abs(actual-float(value))<=3e-12*max(1,float(value)); checks+=2
    return {"fine_fraction_controls":checks,"six_MM_MT_channels":True,
        "exact_adjacent_product_assignment":True,"zero_atoms_and_clipping":True,
        "rational_tie_control_uses_grid_1_4_16":True,
        "inside_grid_witness_squared_at_least_four_fifths_ideal_squared":True}


def fine_bucket_pair_logs(powered,alpha,pair,old,correction,np,guard,raw,progress,*,coarse_builder):
    """Compatible constructor; fine129_bucket_sum is an additional table."""
    tables,coarse_report=coarse_builder(powered,alpha,pair,old,correction,np,guard,raw,progress)
    coarse=tables["dense33_bucket_sum"]
    shapes=sorted(powered)
    values=np.asarray([[float(x.mid()) if hasattr(x,"mid") else float(x)
                        for x in powered[occ]] for occ in shapes])
    if not np.isfinite(values).all() or np.any(values<0):
        raise ArithmeticError("invalid coherent fine-bucket powered table")
    assignment=np.full((len(shapes),2),-1,dtype=int)
    counts=[defaultdict(int),defaultdict(int)]; zero=[0,0]; clipped=[0,0]
    for rowindex,row in enumerate(values):
        guard()
        for i in range(2):
            a,c=(F(*float(row[2*i+j]).as_integer_ratio()) for j in range(2))
            index=nearest_bucket(a,c)
            if index is None: zero[i]+=1; continue
            assignment[rowindex,i]=index; counts[i][index]+=1
            clipped[i]+=c<a*GRID[0] or c>a*GRID[-1]
    assert all(sum(counts[i].values())+zero[i]==len(shapes) for i in range(2))
    occupied=sorted(set(counts[0])|set(counts[1]))
    geometry=correction.route_geometry(shapes,np,guard)
    single=old.float_coefficients(powered,np,guard)
    old_pair=pair.pair_convolution(pair.as_logs(single,np),np,guard)
    del single
    lower=np.full(old_pair.shape,-math.inf)
    positive_left=np.column_stack(((values[:,0]>0)&(values[:,1]>0),(values[:,2]>0)&(values[:,3]>0)))
    positive_right=np.column_stack(((values[:,0]>0)&(values[:,2]>0),
        (values[:,1]>0)&(values[:,3]>0),(values[:,4]>0)&(values[:,5]>0)))
    for position,index in enumerate(occupied):
        guard(); rho=float(GRID[index])
        left=np.column_stack((np.minimum(rho*values[:,0],values[:,1]),
                              np.minimum(rho*values[:,2],values[:,3])))
        right=np.column_stack((np.minimum(values[:,0]/rho,values[:,2]),
            np.minimum(values[:,1]/rho,values[:,3]),np.minimum(values[:,4]/rho,values[:,5])))
        if not np.isfinite(left).all() or not np.isfinite(right).all() or np.any(positive_left&(left==0)) or np.any(positive_right&(right==0)):
            raise ArithmeticError("positive fine-rho minimum under/overflow")
        restricted=np.where(assignment==index,left,0.0)
        coeff=pair.as_logs(correction.route_average(np.column_stack((restricted,right)),geometry,np,guard),np)
        candidate=correction.cross_convolution(coeff[:,:,:2],coeff[:,:,2:],np,guard)
        np.logaddexp(lower,candidate,out=lower)
        if position%8==0 or position+1==len(occupied):
            progress("fine_rho_bucket_complete",occupied_done=position+1,occupied_total=len(occupied))
    fine=np.minimum(correction.subtract_correction(old_pair,lower,alpha,np),coarse)
    assert np.all(fine<=coarse)
    finite=np.isfinite(coarse)
    report={"coarse":coarse_report,"grid_size":len(GRID),"occupied_buckets":len(occupied),
        "grid":"2^e*(1,5/4,3/2,7/4),e=-16..15,plus2^16",
        "counts_by_left_row":[{str(GRID[i]):n for i,n in sorted(row.items())} for row in counts],
        "zero_minimum_omissions_by_row":zero,"ratios_outside_grid_by_row":clipped,
        "improved_entries":int(np.count_nonzero(fine<coarse-1e-12)),
        "maximum_local_gain_bits_from_coarse_FLOAT":float(np.max(coarse[finite]-fine[finite])/LN2),
        "finished_upper_intersected_with_frozen_coarse":True,
        "all_left_occupancies_accounted_for":True,"right_occupancies_unrestricted":True}
    tables["fine129_bucket_sum"]=fine
    return tables,report


def main():
    selected=load_selected(); cache={}; restorations=[]
    original_load=selected.load
    def replace(obj,name,value):
        restorations.append((obj,name,getattr(obj,name))); setattr(obj,name,value)
    def intercept_rho_load(original,name,pin,label):
        module=original(name,pin,label)
        if name=="cheap_lane_pair_profile_calibration.py":
            old_fit,old_replay=module.coordinate_fit,module.replay
            def fit(*args,**kwargs):
                result=old_fit(*args,**kwargs)
                cache["fitted"]=tuple(max(F(1,1<<30),F(math.exp(float(x))).limit_denominator(1<<20)) for x in result[0]["x"])
                return result
            def replay(logs,markers,modulus,*args,**kwargs):
                result=old_replay(logs,markers,modulus,*args,**kwargs)
                if tuple(markers)==cache.get("fitted"):
                    key=(modulus,tuple(markers))
                    if key not in cache.setdefault("fine_replays",{}):
                        extra=old_replay(cache["fine"],markers,modulus,*args,**kwargs)
                        for kind in ("ordinary","residue"):
                            before,after=result[kind],extra[kind]
                            if before["status"]==after["status"]=="resolved":
                                gain=before["log2_selected_run_mass_bound_FLOAT"]-after["log2_selected_run_mass_bound_FLOAT"]
                                assert gain>=-2e-5
                                after["gain_bits_from_coarse_same_markers_FLOAT"]=gain
                        cache["fine_replays"][key]=extra
                    result["fine_same_coarse_fitted_markers"]=cache["fine_replays"][key]
                return result
            replace(module,"coordinate_fit",fit); replace(module,"replay",replay)
        return module
    def intercepted_load(name,pin,label):
        module=original_load(name,pin,label)
        if name=="cheap_lane_rho_bucket_control.py":
            original_builder,original_tests,original_rho_load=module.bucket_pair_logs,module.exact_tests,module.load
            def builder(*args,**kwargs):
                tables,report=fine_bucket_pair_logs(*args,**kwargs,coarse_builder=original_builder)
                cache["fine"]=tables["fine129_bucket_sum"]
                return tables,report  # dense33_bucket_sum stays coarse for fitting.
            def tests(pair,correction,np,guard):
                return {"coarse":original_tests(pair,correction,np,guard),
                    "fine":exact_tests(pair,correction,np,guard)}
            replace(module,"bucket_pair_logs",builder); replace(module,"exact_tests",tests)
            replace(module,"load",lambda name,pin,label:intercept_rho_load(original_rho_load,name,pin,label))
        return module
    original_argv=sys.argv[:]
    try:
        replace(selected,"load",intercepted_load)
        arguments=sys.argv[1:]
        for option,value in (("--full","957"),("--z","23/100"),("--alpha","3/32")):
            if not any(x==option or x.startswith(option+"=") for x in arguments): arguments.extend((option,value))
        sys.argv=[sys.argv[0],*arguments]
        selected.main()
        if "fine_replays" in cache:
            print(json.dumps({"phase":"fine_rho_control_receipt","scope":"selected profile only; fine at coarse-fitted markers",
                "directed":False,"whole_code_certificate":False,"additional_parameter_fit":False,
                "fitted_markers":list(map(str,cache["fitted"])),"moduli_replayed":sorted(k[0] for k in cache["fine_replays"]),
                "selected_driver_pin":SELECTED_PIN,"generated_bytes":0,
                "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}),flush=True)
    finally:
        sys.argv=original_argv
        for obj,name,old in reversed(restorations): setattr(obj,name,old)


if __name__=="__main__":
    main()
