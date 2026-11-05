"""A fixed W=5/4 deficit improvement for the ordered FLOAT pair cap.

For powered path products P,Q, put w=max(P,Q)/min(P,Q). Their exact
pooled value is P+Q-h_alpha(w)*min(P,Q), where
h_alpha(w)=w+1-(w^(1/alpha)+1)^alpha. This function increases on w>=1.
The frozen ordered cap charges kappa=2-2^alpha times a lower minimum.
Here Delta=h_alpha(5/4)-kappa is charged additionally only on certified
off-diagonal bucket rectangles. Exact rational extrema of the SAME
binary64 powered atoms certify each rectangle. No uncertain rectangle,
zero atom, or same-bucket pair receives the extra charge.

The two certified cases are max(s)<=min(r)/W and min(s)>=W*max(r),
where r=c/a for a left row and s=b/d for a right column. They give the
exact minimum ab and cd respectively. Positive prefix/suffix polynomials
sum these cases without an occupancy-pair list. The finished weighted
upper cap is intersected with the frozen ordered upper cap.

This is FLOAT reconnaissance, not directed arithmetic. In particular,
the extrema certify binary64 representatives, and h/Delta are FLOAT.
A directed successor must freeze coherent powered caps or certify their
ordering, and subtract a certified LOWER deficit. No new alpha, source
factor, route normalization, or terminal update is introduced.

weighted_bucket_pair_logs is a bank-compatible constructor. Standalone
uses the pinned selected driver: fit ordered129 once, then replay the
weighted cap at exactly those rational markers. Default f920,z23/100,
alpha3/32; inherited180s/512MiB/20GiB guards, stdout only. No jobs on
import, generated files, or edits to frozen sources.
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
ORDERED_PIN="f816addd50969a83ad951567944d905938ece966c9260175e3b4a08b1be53255"
THRESHOLD=F(5,4)
CHANNELS=((0,0),(0,1),(1,0),(1,1),(0,2),(1,2))
PATHS=((0,0,1,2),(0,1,1,3),(2,0,3,2),(2,1,3,3),(0,4,1,5),(2,4,3,5))
LN2=math.log(2)


def load_ordered():
    path=Path(__file__).with_name("cheap_lane_rho_ordered_buckets.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=ORDERED_PIN:
        raise RuntimeError("frozen ordered helper mismatch")
    spec=importlib.util.spec_from_file_location("weighted_ordered_frozen",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def extrema(values,grid,nearest,guard):
    """Exact extrema of positive binary64 ratios, separately per row/column."""
    left=[[None]*len(grid) for _ in range(2)]
    right=[[None]*len(grid) for _ in range(3)]
    def add(table,channel,index,value):
        old=table[channel][index]
        table[channel][index]=(value,value) if old is None else (min(old[0],value),max(old[1],value))
    for row in values:
        guard(); exact=tuple(F(*float(x).as_integer_ratio()) for x in row)
        if any(x<0 for x in exact): raise ValueError("nonnegative atoms required")
        for i in range(2):
            a,c=exact[2*i:2*i+2]; index=nearest(a,c,grid)
            if index is not None: add(left,i,index,c/a)
        for j,(bi,di) in enumerate(((0,2),(1,3),(4,5))):
            b,d=exact[bi],exact[di]; index=nearest(d,b,grid)
            if index is not None: add(right,j,index,b/d)
    return left,right


def certified_bands(left,right,guard):
    """Return inclusive prefix/suffix cutoffs; -1/B denote empty bands."""
    buckets=len(left[0]); prefix=[[[-1]*3 for _ in range(2)] for _ in range(buckets)]
    suffix=[[[buckets]*3 for _ in range(2)] for _ in range(buckets)]
    rectangles=0
    for i in range(buckets):
        guard()
        for row in range(2):
            interval=left[row][i]
            if interval is None: continue
            lo,hi=interval
            for col in range(3):
                before=[j for j in range(i) if right[col][j] is not None
                        and THRESHOLD*right[col][j][1]<=lo]
                after=[j for j in range(i+1,buckets) if right[col][j] is not None
                       and right[col][j][0]>=THRESHOLD*hi]
                p=max(before,default=-1); s=min(after,default=buckets)
                # Shared ordered buckets make each certified set a prefix or
                # suffix. Check that property rather than silently assuming it.
                assert all(right[col][j] is None or THRESHOLD*right[col][j][1]<=lo for j in range(p+1))
                assert all(right[col][j] is None or right[col][j][0]>=THRESHOLD*hi for j in range(s,buckets))
                assert p<i<s
                prefix[i][row][col]=p; suffix[i][row][col]=s
                rectangles+=len(before)+len(after)
    return prefix,suffix,rectangles


def channel_convolution(left,right,np,guard):
    """Six positive convolutions with a separate right polynomial per channel."""
    n=left.shape[0]-1
    assert left.shape==(n+1,n+1,2) and right.shape==(n+1,n+1,6)
    out=np.full((2*n+1,2*n+1,6),-math.inf)
    for f in range(n+1):
        guard()
        for s in range(n+1-f):
            a=left[f,s]
            terms=right+np.asarray((a[0],a[0],a[1],a[1],a[0],a[1]))
            target=out[f:f+n+1,s:s+n+1]
            np.logaddexp(target,terms,out=target)
    return out


def extra_lower(left,right,bands,pair,np,guard,progress):
    """Positive ab/cd terms from certified, disjoint off-diagonal rectangles."""
    a=pair.as_logs(left,np); b=pair.as_logs(right,np)
    buckets,width=a.shape[:2]; pc,sc,rectangles=bands
    prefix=np.full((buckets+1,width,width,3),-math.inf)
    suffix=np.full_like(prefix,-math.inf)
    for i in range(buckets):
        guard(); np.logaddexp(prefix[i],b[i,:,:,:3],out=prefix[i+1])
    for i in range(buckets-1,-1,-1):
        guard(); np.logaddexp(suffix[i+1],b[i,:,:,3:6],out=suffix[i])
    result=np.full((2*width-1,2*width-1,6),-math.inf); products=0
    occupied=[i for i in range(buckets) if np.isfinite(a[i,:,:,:2]).any()]
    for position,i in enumerate(occupied):
        guard()
        for cuts,bank,lhs,offset in ((pc,prefix,a[i,:,:,:2],1),(sc,suffix,a[i,:,:,2:4],0)):
            rhs=np.empty((width,width,6))
            for channel,(row,col) in enumerate(CHANNELS):
                rhs[:,:,channel]=bank[cuts[i][row][col]+offset,:,:,col]
            if not np.isfinite(rhs).any(): continue
            contribution=channel_convolution(lhs,rhs,np,guard)
            np.logaddexp(result,contribution,out=result); products+=1
        if position%8==0 or position+1==len(occupied):
            progress("weighted_ratio_bucket_complete",left_buckets_done=position+1,
                     left_buckets_total=len(occupied),extra_convolutions=products)
    return result,{"certified_row_column_bucket_rectangles":rectangles,"extra_convolutions":products}


def deficit_constants(alpha):
    if not math.isfinite(alpha) or not 0<alpha<=1: raise ValueError("require0<alpha<=1")
    if alpha==1: return 0.0,0.0,0.0
    w=float(THRESHOLD); kappa=2-math.exp2(alpha)
    # Stable form, including very small alpha: no exp(log(w)/alpha).
    h=1-w*math.expm1(alpha*math.log1p(math.exp(-math.log(w)/alpha)))
    delta=h-kappa
    if delta < -2e-14 or not -2e-14<=h<=1+2e-14:
        raise ArithmeticError("invalid FLOAT deficit constants")
    return kappa,max(0.0,delta),h


def subtract_weighted(old,base,extra,alpha,np):
    kappa,delta,h=deficit_constants(alpha)
    if np.any(np.isfinite(extra)&~np.isfinite(base)):
        raise ArithmeticError("extra deficit outside the ordered minimum")
    common=np.isfinite(extra)&np.isfinite(base)
    if np.any(extra[common]>base[common]+2e-11):
        raise ArithmeticError("extra minimum exceeds the ordered minimum")
    deficit=np.full(old.shape,-math.inf)
    if kappa: deficit=base+math.log(kappa)
    if delta: deficit=np.logaddexp(deficit,extra+math.log(delta))
    finite=np.isfinite(old)
    if np.any(np.isfinite(deficit)&~finite): raise ArithmeticError("deficit at zero old coefficient")
    ratio=np.exp(deficit[finite]-old[finite])
    if np.any(ratio>.5+2e-11): raise ArithmeticError("weighted deficit exceeds half the old pair")
    out=old.copy(); out[finite]+=np.log1p(-ratio)
    return out,{"threshold_exact":str(THRESHOLD),"kappa_FLOAT":kappa,
                "h_threshold_FLOAT":h,"additional_delta_FLOAT":delta}


def exact_tests(ordered,fine,pair,correction,np,guard):
    """Rational reference enumerates toy pairs, never production occupancies."""
    grid=(F(1),F(4),F(16)); checks=0
    atoms=(((0,0),F(2),(F(1),F(2),F(3),F(6),F(4),F(8))),
           ((1,0),F(3),(F(4),F(1),F(0),F(5),F(2),F(0))),
           ((0,1),F(1),(F(2),F(10),F(6),F(2),F(0),F(7))),
           ((1,1),F(5),(F(1),F(128),F(1),F(1,128),F(3),F(7))))
    values=np.asarray([list(map(float,row)) for _,_,row in atoms])
    ex=extrema(values,grid,fine.nearest_bucket,guard); bands=certified_bands(*ex,guard)
    expected=defaultdict(lambda:[F(0)]*6); ideal=defaultdict(lambda:[F(0)]*6)
    selected=zero=ties=0
    for da,wa,a in atoms:
        for db,wb,b in atoms:
            degree=da[0]+db[0],da[1]+db[1]
            for channel,(ai,bi,ci,di) in enumerate(PATHS):
                aa,cc,bb,dd=a[ai],a[ci],b[bi],b[di]
                i=fine.nearest_bucket(aa,cc,grid); j=fine.nearest_bucket(dd,bb,grid)
                target=min(aa*bb,cc*dd); ideal[degree][channel]+=wa*wb*target
                if i is None or j is None:
                    assert target==0; zero+=1; continue
                row,col=CHANNELS[channel]
                low=j<=bands[0][i][row][col]; high=j>=bands[1][i][row][col]
                assert not(low and high)
                if low or high:
                    witness=aa*bb if low else cc*dd
                    ratio=max((cc/aa)/(bb/dd),(bb/dd)/(cc/aa))
                    assert ratio>=THRESHOLD and witness==target and i!=j
                    expected[degree][channel]+=wa*wb*witness; selected+=1
                ties+=any((cc/aa)**2==x*y for x,y in zip(grid,grid[1:]))
                checks+=2
    geometry=(np.asarray([d[0] for d,_,_ in atoms]),np.asarray([d[1] for d,_,_ in atoms]),
              np.arange(len(atoms)),np.asarray([float(w) for _,w,_ in atoms]))
    left,right,routing=ordered.bucket_polynomials(values,geometry,grid,fine.nearest_bucket,np,guard)
    extra,work=extra_lower(left,right,bands,pair,np,guard,lambda *a,**k:None)
    base,_=ordered.ordered_lower(left,right,pair,correction,np,guard,lambda *a,**k:None)
    for degree,row in ideal.items():
        for channel,total in enumerate(row):
            want=expected[degree][channel]
            got=math.exp(float(extra[degree][channel])) if math.isfinite(extra[degree][channel]) else 0.0
            assert 0<=want<=total and abs(got-float(want))<=5e-12*max(1,float(want)); checks+=2
            if math.isfinite(extra[degree][channel]): assert extra[degree][channel]<=base[degree][channel]+2e-12
    assert selected and zero and ties
    assert deficit_constants(1.0)[1]<=2e-14
    for alpha in (1/2,3/32,13/128):
        k,d,h=deficit_constants(alpha); assert 0<=k<=h<=1 and d>=0
    return {"exact_fraction_checks":checks,"certified_toy_pairs":selected,
        "zero_pairs":zero,"boundary_ties":ties,"no_double_count":True,
        "all_six_MM_MT_channels":True,"routing":routing,"extra_work":work}


def weighted_bucket_pair_logs(powered,alpha,pair,old,correction,np,guard,raw,progress,*,coarse_builder):
    """Add weighted_ordered129 to the unmodified ordered/fine table bank."""
    ordered=load_ordered(); fine=ordered.load_fine(); cache={}
    original=ordered.ordered_lower
    def capture(left,right,*args,**kwargs):
        lower,work=original(left,right,*args,**kwargs)
        cache.update(left=left,right=right,base=lower)
        return lower,work
    ordered.ordered_lower=capture
    try:
        tables,report=ordered.ordered_bucket_pair_logs(powered,alpha,pair,old,correction,np,
            guard,raw,progress,coarse_builder=coarse_builder)
    finally:
        ordered.ordered_lower=original
    values=[]
    for shape in sorted(powered):
        row=[float(x.mid()) if hasattr(x,"mid") else float(x) for x in powered[shape]]
        if any(not math.isfinite(x) or x<0 for x in row): raise ArithmeticError("invalid powered atom")
        values.append(row)
    bounds=extrema(values,fine.GRID,fine.nearest_bucket,guard)
    bands=certified_bands(*bounds,guard)
    extra,work=extra_lower(cache["left"],cache["right"],bands,pair,np,guard,progress)
    single=old.float_coefficients(powered,np,guard)
    old_pair=pair.pair_convolution(pair.as_logs(single,np),np,guard)
    candidate,constants=subtract_weighted(old_pair,cache["base"],extra,alpha,np)
    previous=tables["ordered129_bucket_sum"]; weighted=np.minimum(candidate,previous)
    assert np.all(weighted<=previous)
    finite=np.isfinite(previous)
    tables["weighted_ordered129"]=weighted
    return tables,{"ordered":report,"constants_FLOAT":constants,"work":work,
        "exact_extrema_of_binary64_atoms":True,"same_bucket_extra":False,
        "zero_atoms_excluded_from_ratios":True,"finished_ordered_upper_retained":True,
        "improved_entries":int(np.count_nonzero(weighted<previous-1e-12)),
        "maximum_local_gain_bits_from_ordered_FLOAT":float(np.max(previous[finite]-weighted[finite])/LN2),
        "directed":False}


def main():
    ordered=load_ordered(); fine=ordered.load_fine(); selected=fine.load_selected()
    cache={}; restore=[]
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
                        extra=old_replay(cache["weighted"],markers,modulus,*args,**kwargs)
                        for kind in ("ordinary","residue"):
                            before,after=result[kind],extra[kind]
                            if before["status"]==after["status"]=="resolved":
                                gain=before["log2_selected_run_mass_bound_FLOAT"]-after["log2_selected_run_mass_bound_FLOAT"]
                                assert gain>=-2e-5; after["gain_bits_from_ordered_same_markers_FLOAT"]=gain
                        cache["replays"][key]=extra
                    result["weighted_same_ordered_fitted_markers"]=cache["replays"][key]
                return result
            replace(module,"coordinate_fit",fit); replace(module,"replay",replay)
        return module
    def intercept(name,pin,label):
        module=original_load(name,pin,label)
        if name=="cheap_lane_rho_bucket_control.py":
            old_builder,old_tests,old_load=module.bucket_pair_logs,module.exact_tests,module.load
            def builder(*args,**kwargs):
                tables,report=weighted_bucket_pair_logs(*args,**kwargs,coarse_builder=old_builder)
                cache["weighted"]=tables["weighted_ordered129"]
                tables["dense33_bucket_sum"]=tables["ordered129_bucket_sum"]
                report["selected_fit_key_is_unweighted_ordered"]=True
                return tables,report
            def tests(pair,correction,np,guard):
                return {"coarse":old_tests(pair,correction,np,guard),
                    "N4":fine.exact_tests(pair,correction,np,guard),
                    "ordered":ordered.exact_tests(fine,pair,correction,np,guard),
                    "weighted":exact_tests(ordered,fine,pair,correction,np,guard)}
            replace(module,"bucket_pair_logs",builder); replace(module,"exact_tests",tests)
            replace(module,"load",lambda name,pin,label:load_control(old_load,name,pin,label))
        return module
    def emit(*values,**kwargs):
        if len(values)==1 and isinstance(values[0],str):
            row=json.loads(values[0])
            if isinstance(row,dict) and "pair_bound" in row:
                row["pair_bound"]="ordered129 control; fixed W=5/4 weighted correction at SAME ordered-fitted markers"
                row["weighted_wrapper_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
                row["ordered_pin"]=ORDERED_PIN; row["directed"]=False
                values=(json.dumps(row,indent=2),)
        builtins.print(*values,**kwargs)
    argv=sys.argv[:]
    try:
        replace(selected,"load",intercept); replace(selected,"print",emit)
        arguments=sys.argv[1:]
        for option,value in (("--full","920"),("--z","23/100"),("--alpha","3/32")):
            if not any(x==option or x.startswith(option+"=") for x in arguments): arguments.extend((option,value))
        sys.argv=[sys.argv[0],*arguments]; selected.main()
    finally:
        sys.argv=argv
        for obj,name,exists,value in reversed(restore):
            if exists: setattr(obj,name,value)
            else: delattr(obj,name)


if __name__=="__main__":
    main()
