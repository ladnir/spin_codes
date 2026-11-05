"""Exact off-diagonal ratio-bucket correction, with a FLOAT-only control.

Use the frozen rational129 grid and its ordered nearest-log intervals.
For powered left row(a,c) and right column(b,d), put r=c/a,s=b/d.
If the left bucket i exceeds the right bucket j, min(ab,cd)=ab;
if i<j it equals cd. Only i=j uses the old rho witness. Thus sum
A_i*prefix(B)_i+C_i*suffix(D)_i+L_i*R_i over disjoint left buckets,
where L_i=min(rho_i*a,c), R_i=min(b/rho_i,d). Every term is positive.
The partitions are per incoming row and per outgoing column; the third
right column is the output-only terminal(T0,U). Zero minima are omitted.
Clipped edge buckets remain valid, though their diagonal may be loose.

ordered_bucket_pair_logs exports the finished ordered129_bucket_sum,
intersected with the frozen fine129 upper cap. Both use coherent caps.
The local route geometry is aggregated once per scalar channel, avoiding
one full geometry pass for every bucket. No occupancy-pair list is built.

This is FLOAT reconnaissance. Bucket comparisons are exact for the binary64
atoms, not certified comparisons of the underlying Arb powers. A directed
successor must certify ordering of coherent powers; midpoint comparisons
must not justify exact off-diagonal terms in a certificate.

Standalone: default f915,z23/100,alpha3/32. Fit the frozen N4 finished cap
only, then compare the ordered cap at those same rational markers using
m256/m512. Source/count/normalizer/terminal stay in the pinned selected
driver. No complete-weight claim, jobs on import, files, or frozen edits.
Parent execution only; inherited180s/512MiB/20GiB, stdout only.
"""

from collections import defaultdict
from fractions import Fraction as F
import builtins
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode=True
FINE_PIN="ca8a73442ff87ca46ce64192dd999967716cf145c65a9712553c0e81a9a69701"
PATHS=((0,0,1,2),(0,1,1,3),(2,0,3,2),(2,1,3,3),(0,4,1,5),(2,4,3,5))
LN2=math.log(2)


def load_fine():
    path=Path(__file__).with_name("cheap_lane_rho_fine_buckets.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=FINE_PIN:
        raise RuntimeError("frozen fine-bucket helper mismatch")
    spec=importlib.util.spec_from_file_location("ordered_rho_fine",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def bucket_polynomials(values,geometry,grid,nearest,np,guard):
    """Route the2 left and3 right partitions into positive degree polynomials."""
    guard(); n=len(grid)
    left_assignment=np.full((len(values),2),-1,dtype=int)
    right_assignment=np.full((len(values),3),-1,dtype=int)
    left_atom=np.zeros((len(values),6)); right_atom=np.zeros((len(values),9))
    counts=[[0]*n for _ in range(5)]; clipped=[0]*5; omitted=[0]*5
    for k,row in enumerate(values):
        guard()
        for i in range(2):
            a,c=(F(*float(row[2*i+j]).as_integer_ratio()) for j in range(2))
            index=nearest(a,c,grid)
            if index is None: omitted[i]+=1; continue
            left_assignment[k,i]=index; counts[i][index]+=1
            clipped[i]+=c<a*grid[0] or c>a*grid[-1]
            rho=float(grid[index])
            left_atom[k,i],left_atom[k,2+i]=float(a),float(c)
            left_atom[k,4+i]=min(rho*float(a),float(c))
        for j,(bi,di) in enumerate(((0,2),(1,3),(4,5))):
            b,d=(F(*float(row[i]).as_integer_ratio()) for i in (bi,di))
            index=nearest(d,b,grid)  # Right ratio is b/d, not d/b.
            if index is None: omitted[2+j]+=1; continue
            right_assignment[k,j]=index; counts[2+j][index]+=1
            clipped[2+j]+=b<d*grid[0] or b>d*grid[-1]
            rho=float(grid[index])
            right_atom[k,j],right_atom[k,3+j]=float(b),float(d)
            right_atom[k,6+j]=min(float(b)/rho,float(d))
    for i in range(5): assert sum(counts[i])+omitted[i]==len(values)
    if not np.isfinite(left_atom).all() or not np.isfinite(right_atom).all():
        raise ArithmeticError("ordered-ratio atom overflow")
    if np.any((left_assignment>=0)&(left_atom[:,4:]==0)) or np.any((right_assignment>=0)&(right_atom[:,6:]==0)):
        raise ArithmeticError("positive ordered-ratio witness underflow")
    fd,sd,oi,weights=geometry
    # cross_convolution iterates a TOTAL-degree triangle, not a rectangle.
    width=int(np.max(fd+sd))+1
    assert np.all(fd>=0) and np.all(sd>=0) and np.all(fd+sd<width)
    left=np.zeros((n,width,width,6)); right=np.zeros((n,width,width,9))
    for assignment,atoms,out,channels in ((left_assignment,left_atom,left,2),
                                        (right_assignment,right_atom,right,3)):
        for channel in range(3*channels):
            guard(); ids=assignment[oi,channel%channels]; valid=ids>=0
            mass=weights[valid]*atoms[oi[valid],channel]
            if np.any(mass<=0) or not np.isfinite(mass).all():
                raise ArithmeticError("positive bucket route atom under/overflow")
            np.add.at(out,(ids[valid],fd[valid],sd[valid],channel),mass)
    return left,right,{"counts_by_left_row":counts[:2],"counts_by_right_column":counts[2:],
        "clipped_by_row_or_column":clipped,"zero_minimum_omissions":omitted,
        "all_occupancies_accounted_for":True,"geometry_passes":15}


def ordered_lower(left,right,pair,correction,np,guard,progress):
    """Three positive convolutions per occupied left bucket, no subtraction."""
    a=pair.as_logs(left,np); b=pair.as_logs(right,np)
    del left,right
    buckets,width=a.shape[:2]
    prefix=np.full((buckets,width,width,3),-math.inf)
    suffix=np.full_like(prefix,-math.inf)
    running=np.full((width,width,3),-math.inf)
    for i in range(buckets):
        guard(); prefix[i]=running
        np.logaddexp(running,b[i,:,:,:3],out=running)
    running.fill(-math.inf)
    for i in range(buckets-1,-1,-1):
        guard(); suffix[i]=running
        np.logaddexp(running,b[i,:,:,3:6],out=running)
    lower=np.full((2*width-1,2*width-1,6),-math.inf)
    occupied=[i for i in range(buckets) if np.isfinite(a[i,:,:,:2]).any()]
    products=0
    for position,i in enumerate(occupied):
        guard()
        for lhs,rhs in ((a[i,:,:,:2],prefix[i]),(a[i,:,:,2:4],suffix[i]),
                        (a[i,:,:,4:6],b[i,:,:,6:9])):
            if not np.isfinite(rhs).any(): continue
            contribution=correction.cross_convolution(lhs,rhs,np,guard)
            np.logaddexp(lower,contribution,out=lower); products+=1
        if position%8==0 or position+1==len(occupied):
            progress("ordered_ratio_bucket_complete",left_buckets_done=position+1,
                     left_buckets_total=len(occupied),cross_convolutions=products)
    return lower,{"occupied_left_buckets":len(occupied),"cross_convolutions":products}


def exact_tests(fine,pair,correction,np,guard):
    """Independent rational all-route six-channel reference, including ties."""
    grid=(F(1),F(4),F(16)); checks=0
    atoms=(((0,0),F(2),(F(1),F(2),F(3),F(6),F(4),F(8))),
           ((1,0),F(3),(F(4),F(1),F(0),F(5),F(2),F(0))),
           ((0,1),F(1),(F(2),F(10),F(6),F(2),F(0),F(7))),
           ((1,1),F(5),(F(1),F(128),F(1),F(1,128),F(3),F(7))))
    expected=defaultdict(lambda:[F(0)]*6)
    one_sided=defaultdict(lambda:[F(0)]*6)
    ideal=defaultdict(lambda:[F(0)]*6)
    exact_offdiag=ties=zeros=clipped=0
    for da,wa,a in atoms:
        for db,wb,b in atoms:
            degree=da[0]+db[0],da[1]+db[1]
            for j,(ai,bi,ci,di) in enumerate(PATHS):
                aa,cc,bb,dd=a[ai],a[ci],b[bi],b[di]
                wanted=min(aa*bb,cc*dd)
                i=fine.nearest_bucket(aa,cc,grid); k=fine.nearest_bucket(dd,bb,grid)
                if i is None or k is None:
                    witness=old=F(0); zeros+=1; assert wanted==0
                else:
                    rho=grid[i]; old=min(rho*aa,cc)*min(bb/rho,dd)
                    witness=aa*bb if k<i else cc*dd if k>i else old
                    if k!=i: assert witness==wanted; exact_offdiag+=1
                    clipped+=cc<aa*grid[0] or cc>aa*grid[-1] or bb<dd*grid[0] or bb>dd*grid[-1]
                    ties+=any((cc/aa)**2==x*y for x,y in zip(grid,grid[1:]))
                assert 0<=old<=witness<=wanted; checks+=1
                expected[degree][j]+=wa*wb*witness
                one_sided[degree][j]+=wa*wb*old; ideal[degree][j]+=wa*wb*wanted
    values=np.asarray([list(map(float,row)) for _,_,row in atoms])
    geometry=(np.asarray([d[0] for d,_,_ in atoms]),np.asarray([d[1] for d,_,_ in atoms]),
              np.arange(len(atoms)),np.asarray([float(w) for _,w,_ in atoms]))
    left,right,report=bucket_polynomials(values,geometry,grid,fine.nearest_bucket,np,guard)
    computed,work=ordered_lower(left,right,pair,correction,np,guard,lambda *a,**k:None)
    for degree,row in expected.items():
        for j,value in enumerate(row):
            assert one_sided[degree][j]<=value<=ideal[degree][j]
            actual=math.exp(float(computed[degree][j])) if math.isfinite(computed[degree][j]) else 0.0
            assert abs(actual-float(value))<=4e-12*max(1,float(value)); checks+=2
    assert exact_offdiag and ties and zeros and clipped
    return {"exact_fraction_checks":checks,"exact_offdiagonal_cases":exact_offdiag,
        "boundary_ties":ties,"zero_cases":zeros,"clipped_cases":clipped,
        "all_six_MM_MT_channels":True,"one_sided_pointwise_dominated":True,
        "routing":report,"convolutions":work}


def ordered_bucket_pair_logs(powered,alpha,pair,old,correction,np,guard,raw,progress,*,coarse_builder):
    """Compatible bank constructor, adding ordered129_bucket_sum."""
    fine=load_fine()
    tables,old_report=fine.fine_bucket_pair_logs(powered,alpha,pair,old,correction,np,
        guard,raw,progress,coarse_builder=coarse_builder)
    previous=tables["fine129_bucket_sum"]
    shapes=sorted(powered); rows=[]
    for occ in shapes:
        row=[]
        for value in powered[occ]:
            number=float(value.mid()) if hasattr(value,"mid") else float(value)
            if not math.isfinite(number) or number<0:
                raise ArithmeticError("invalid powered ordered-bucket atom")
            if hasattr(value,"is_zero") and not value.is_zero() and number==0:
                raise ArithmeticError("positive powered atom underflow before ratio classification")
            row.append(number)
        rows.append(row)
    values=np.asarray(rows); del rows
    geometry=correction.route_geometry(shapes,np,guard)
    left,right,routing=bucket_polynomials(values,geometry,fine.GRID,fine.nearest_bucket,np,guard)
    lower,work=ordered_lower(left,right,pair,correction,np,guard,progress)
    single=old.float_coefficients(powered,np,guard)
    old_pair=pair.pair_convolution(pair.as_logs(single,np),np,guard); del single
    ordered=np.minimum(correction.subtract_correction(old_pair,lower,alpha,np),previous)
    assert np.all(ordered<=previous)
    finite=np.isfinite(previous)
    tables["ordered129_bucket_sum"]=ordered
    return tables,{"frozen_N4":old_report,"grid_size":len(fine.GRID),
        "routing":routing,"work":work,"same_grid_as_N4":True,
        "exact_offdiagonal_ratio_cases":True,"diagonal_uses_old_rho_witness":True,
        "finished_upper_intersected_with_N4":True,"directed":False,
        "improved_entries":int(np.count_nonzero(ordered<previous-1e-12)),
        "maximum_local_gain_bits_from_N4_FLOAT":float(np.max(previous[finite]-ordered[finite])/LN2)}


def main():
    fine=load_fine(); selected=fine.load_selected(); cache={}; restore=[]
    def replace(obj,name,value):
        exists=hasattr(obj,name); restore.append((obj,name,exists,getattr(obj,name,None))); setattr(obj,name,value)
    original_load=selected.load
    def load_control(original,name,pin,label):
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
                    key=modulus,tuple(markers)
                    if key not in cache.setdefault("replays",{}):
                        extra=old_replay(cache["ordered"],markers,modulus,*args,**kwargs)
                        for kind in ("ordinary","residue"):
                            before,after=result[kind],extra[kind]
                            if before["status"]==after["status"]=="resolved":
                                gain=before["log2_selected_run_mass_bound_FLOAT"]-after["log2_selected_run_mass_bound_FLOAT"]
                                assert gain>=-2e-5; after["gain_bits_from_N4_same_markers_FLOAT"]=gain
                        cache["replays"][key]=extra
                    result["ordered_same_N4_fitted_markers"]=cache["replays"][key]
                return result
            replace(module,"coordinate_fit",fit); replace(module,"replay",replay)
        return module
    def intercept(name,pin,label):
        module=original_load(name,pin,label)
        if name=="cheap_lane_rho_bucket_control.py":
            old_builder,old_tests,old_load=module.bucket_pair_logs,module.exact_tests,module.load
            def builder(*args,**kwargs):
                tables,report=ordered_bucket_pair_logs(*args,**kwargs,coarse_builder=old_builder)
                cache["ordered"]=tables["ordered129_bucket_sum"]
                tables["dense33_bucket_sum"]=tables["fine129_bucket_sum"]
                report["selected_fit_key_is_N4"]=True
                return tables,report
            def tests(pair,correction,np,guard):
                return {"coarse":old_tests(pair,correction,np,guard),
                    "N4":fine.exact_tests(pair,correction,np,guard),
                    "ordered":exact_tests(fine,pair,correction,np,guard)}
            replace(module,"bucket_pair_logs",builder); replace(module,"exact_tests",tests)
            replace(module,"load",lambda name,pin,label:load_control(old_load,name,pin,label))
        return module
    def emit(*values,**kwargs):
        if len(values)==1 and isinstance(values[0],str):
            row=json.loads(values[0])
            if isinstance(row,dict) and "pair_bound" in row:
                row["pair_bound"]="N4 fine129 control; nested ordered129 at SAME N4-fitted markers"
                row["ordered_wrapper_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
                row["fine_pin"]=FINE_PIN; row["directed"]=False
                values=(json.dumps(row,indent=2),)
        builtins.print(*values,**kwargs)
    argv=sys.argv[:]
    try:
        replace(selected,"load",intercept); replace(selected,"print",emit)
        arguments=sys.argv[1:]
        for option,value in (("--full","915"),("--z","23/100"),("--alpha","3/32")):
            if not any(x==option or x.startswith(option+"=") for x in arguments): arguments.extend((option,value))
        sys.argv=[sys.argv[0],*arguments]; selected.main()
    finally:
        sys.argv=argv
        for obj,name,exists,value in reversed(restore):
            if exists: setattr(obj,name,value)
            else: delattr(obj,name)


if __name__=="__main__":
    main()
