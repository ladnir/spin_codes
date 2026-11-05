"""Directed count-scaled envelopes for bands of ACTUAL source sizes.

Inputs are own-normalized regional R_(r,f), V_(r,f), at common z/alpha,
and source scales s_r=Omega_r^(alpha/256). For a fixed d=a+b let
B_d={r:N_r(r-d)>0}. Multiply each r's regions by
 s_r * N_r(r-d)^(1/256),
then take entrywise maxima over B_d. The finite chain has255 regular
factors and one terminal factor. Thus |B_d| times its maximum over
ordered a+b=d bounds the sum of the counted, source-scaled r-chains.
The count exponent is1/256, NOT alpha/256. The caller applies the common
z^(-alpha*L) once; no further source/count/alpha factor is needed.

This uniform-in-r choice need not improve the frozen count-sum envelope.
Only minima of COMPLETED valid degree bounds may be used to compare them.
The count scaling is fixed across every profile pair at its degree d.
Actual f indices only: hole-indexed regions require different regrouping.

selected_envelope accepts sparse original mappings for selected controls.
counted_join rejects incomplete original coverage and covers every degree.
Both return dictionaries; selected output contains envelope, active_weights,
counts_by_weight, countfactor and report. Complete output contains bound,
degree_bounds and report. The caller owns time/memory/disk guards. At most
17 actual weights1..4096; no computation, file output or tests on import.
"""

from fractions import Fraction as F
import hashlib
import importlib.util
from math import comb
from pathlib import Path
import sys
import time

sys.dont_write_bytecode=True
BASE_PIN="127a6a82fcb785bcb339115f0e0ea8d951841340ce32540655b3e86b06a20503"
POINTS,FACTORS=8192,256


def load_base():
    path=Path(__file__).with_name("cheap_lane_band_run_mass.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=BASE_PIN:
        raise RuntimeError("frozen band helper mismatch")
    spec=importlib.util.spec_from_file_location("counted_frozen_band",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def _scale(region,value):
    r,v=region
    return ((value*r[0],value*r[1],value*r[2],value*r[3]),
            (value*v[0],value*v[1]))


def _maximum(old,new):
    a,x=old; b,y=new
    return ((a[0].max(b[0]),a[1].max(b[1]),a[2].max(b[2]),a[3].max(b[3])),
            (x[0].max(y[0]),x[1].max(y[1])))


def _prepare(regions,scales,join,guard,complete):
    from collections.abc import Mapping
    from flint import arb
    base=load_base()
    if not isinstance(regions,Mapping) or not isinstance(scales,Mapping):
        raise TypeError("regions and scales must be weight-indexed mappings")
    weights=base._weights(regions)
    if set(scales)!=set(weights):
        raise ValueError("source scales must match the actual weights exactly")
    scaled={}; counts={}; intersections=0
    for r in weights:
        guard(); rows=regions[r]; scale=scales[r]
        if not isinstance(rows,Mapping):
            raise TypeError("each weight needs a sparse actual-f mapping")
        if not isinstance(scale,arb) or not scale.is_finite() or not scale>0:
            raise ArithmeticError("source scale must be a finite positive Arb enclosure")
        if any(type(f) is not int or not 0<=f<=r//2 for f in rows):
            raise ValueError("invalid actual full-pair index")
        if complete and set(rows)!=set(range(r//2+1)):
            raise ValueError("complete join requires every original own-r marginal")
        current={}
        for f,region in rows.items():
            guard(); clean,number=base._region(region,arb); intersections+=number
            value,number=base._region(_scale(clean,scale),arb); intersections+=number
            current[f]=value
        scaled[r]=current
        numbers=tuple(join.run_counts(points=POINTS,weight=r))
        if (len(numbers)!=r+1 or numbers[0]!=0
                or any(type(n) is not int or n<0 for n in numbers)
                or sum(numbers)!=comb(POINTS,r)):
            raise ArithmeticError("invalid exact indexed run counts")
        counts[r]=numbers
    return weights,scaled,counts,intersections


def _make(prepared,d,guard):
    """Stream one d-specific envelope; no all-degree matrix cache."""
    from flint import arb
    weights,scaled,counts,_=prepared
    active=tuple(r for r in weights if 1<=r-d<=r and counts[r][r-d]>0)
    numbers={r:counts[r][r-d] for r in active}; result={}
    if not active:
        return result,active,numbers
    full=max(active)//2; lo=max(0,d-full); hi=min(full,d)
    for r in active:
        guard(); number=numbers[r]
        # Directed positive root of an EXACT count, never a FLOAT exponent.
        root=arb(1) if number==1 else (arb(number).log()/FACTORS).exp()
        if not root.is_finite() or not root>0:
            raise ArithmeticError("invalid directed count root")
        rows=scaled[r]
        # Both indices must be valid for THIS r before taking a cross-r max.
        for f in range(max(lo,d-r//2),min(hi,r//2)+1):
            if f%64==0: guard()
            region=rows.get(f)
            if region is None: continue
            value=_scale(region,root)
            result[f]=value if f not in result else _maximum(result[f],value)
    return result,active,numbers


def selected_envelope(regions,scales,d,join,guard):
    """Sparse selected-d control. No claim that omitted marginals vanish.

    envelope contains only supplied indices capable of participating at d.
    An empty B_d returns an empty envelope and countfactor0. Otherwise a
    selected profile's counted bound is countfactor*two_orientations(...),
    followed by the caller's common output correction. A profile control
    is not a maximum over its degree class or a complete band certificate.
    """
    if type(d) is not int or d<0:
        raise ValueError("degree sum must be a nonnegative integer")
    prepared=_prepare(regions,scales,join,guard,False)
    envelope,active,numbers=_make(prepared,d,guard)
    guard()
    return {"envelope":envelope,"active_weights":active,"counts_by_weight":numbers,
        "countfactor":len(active),"report":{
            "scope":"selected count-scaled degree envelope only",
            "actual_weights":prepared[0],"degree_sum":d,"active_weights":active,
            "supplied_envelope_indices":tuple(sorted(envelope)),
            "count_root_denominator":FACTORS,"count_root_includes_alpha":False,
            "source_scale_already_in_each_regional_factor":True,
            "additional_degree_factor":len(active),"input_intersections":prepared[3],
            "output_weight_correction_applied":False,"complete_band_bound":False,
            "whole_code_certificate":False}}


def counted_join(regions,scales,join,guard,progress=None):
    """Complete actual-weight band; return bound before common z correction.

    The source/count-scaled matrices are rebuilt for each d, but exact run
    counts and source-scaled original regions are cached once. Both terminal
    orders share each unordered matrix power. No interpolation is used.
    """
    from flint import arb
    started=time.monotonic()
    prepared=_prepare(regions,scales,join,guard,True)
    weights,_,counts,intersections=prepared
    full=max(weights)//2; total=arb(0); degree_bounds=[]
    ordered=unordered=0; observed_count=0; active_sizes=[]
    for d in range(max(weights)):
        guard(); envelope,active,numbers=_make(prepared,d,guard)
        maximum=arb(0); active_sizes.append(len(active)); observed_count+=sum(numbers.values())
        lo,hi=max(0,d-full),min(full,d)
        if active:
            if set(envelope)!=set(range(lo,hi+1)):
                raise ArithmeticError("missing complete d-specific envelope coverage")
            for a in range(lo,min(hi,d//2)+1):
                b=d-a
                if unordered%128==0: guard()
                left,right=join.two_orientations(envelope[a],envelope[b],diagonal=a==b)
                maximum=maximum.max(left)
                if a!=b: maximum=maximum.max(right)
                ordered+=1 if a==b else 2; unordered+=1
        term=len(active)*maximum
        if not term.is_finite() or term.upper()<0:
            raise ArithmeticError("invalid counted band contribution")
        term=term.nonnegative_part(); total+=term; degree_bounds.append(term)
        if progress is not None and (d==0 or (d+1)%64==0 or d+1==max(weights)):
            progress("directed_counted_band",degree_sum=d,active_weights=len(active),
                ordered_envelope_pairs=ordered,unordered_powers=unordered,
                seconds=time.monotonic()-started)
    expected=(full+1)**2-(1 if max(weights)%2==0 else 0)
    if ordered!=expected or observed_count!=sum(comb(POINTS,r) for r in weights):
        raise ArithmeticError("complete pair/count checksum failed")
    guard()
    return {"bound":total,"degree_bounds":tuple(degree_bounds),"report":{
        "scope":"complete directed count-scaled run envelope for the listed actual sizes",
        "actual_weights":weights,"points":POINTS,"actual_full_indices_not_holes":True,
        "all_original_weight_marginals_present":True,"ordered_envelope_pairs":ordered,
        "unordered_powers":unordered,"degree_entries":len(degree_bounds),
        "maximum_active_weights":max(active_sizes),"all_band_counts_sum_to_binomials":True,
        "band_count_sum_sha256":hashlib.sha256(str(observed_count).encode()).hexdigest(),
        "source_scale_already_in_each_regional_factor":True,"new_alpha_applied":False,
        "count_root_denominator":FACTORS,"count_root_includes_alpha":False,
        "degree_factor":"number of actual sizes with positive run mass at d",
        "regional_order":"e0(Rbar_d,a Rbar_d,b)^127 Rbar_d,a Vbar_d,b; reverse keeps Vbar_d,a",
        "regional_factors":FACTORS,"packet_updates":32767,
        "input_intersections":intersections,"output_weight_correction_applied":False,
        "caller_final_correction":"multiply by common z^(-alpha*L) exactly once",
        "guaranteed_improvement_over_old_envelope":False,
        "before_output_correction_log2_interval":str(total.log()/arb(2).log()) if total>0 else "zero",
        "seconds":time.monotonic()-started,"whole_code_certificate":False}}


def exact_tests(join,guard):
    """Exact tiny chains with perfect-power counts; parent execution only."""
    from flint import arb
    checks=0
    def toy(r,f):
        return ((F(r+f+1,19),F(f,23),F(2,29),F(r+1,31)),
                (F(r+f,37),F(0) if f==0 else F(f+2,41)))
    def maximum(regions,scales,d):
        keys=set().union(*(rows.keys() for rows in regions.values()))
        return {f:tuple(tuple(max((scales[r]*rows[f][i][j]
                for r,rows in regions.items() if f in rows),default=F(0))
                for j in range(width)) for i,width in enumerate((4,2))) for f in keys}
    def chain(left,right,bits):
        a,_=left; b,v=right
        for k in range((1<<(bits+1))-1):
            m=a if k%2==0 else b
            v=(m[0]*v[0]+m[1]*v[1],m[2]*v[0]+m[3]*v[1])
        return v[0]
    def max_chain(rows,d,bits):
        return max((chain(rows[a],rows[d-a],bits) for a in rows if d-a in rows),default=F(0))
    original={2:{f:toy(2,f) for f in (0,1)},3:{f:toy(3,f) for f in (0,1)}}
    source_roots={2:F(2,3),3:F(3,5)}
    for bits in (1,2,3):
        guard(); factors=1<<(bits+1)
        for d in (0,1,2):
            # r=2 has no positive run count at d=2; its matrices cannot enter.
            active=tuple(r for r in original if r>d)
            rows={r:original[r] for r in active}
            count_roots={r:F(r+d+1,1) for r in active}
            numbers={r:count_roots[r]**factors for r in active}
            scales={r:source_roots[r]*count_roots[r] for r in active}
            envelope=maximum(rows,scales,d)
            separate=sum((numbers[r]*source_roots[r]**factors*max_chain(rows[r],d,bits)
                          for r in active),F(0))
            bound=len(active)*max_chain(envelope,d,bits)
            assert separate<=bound; checks+=1
            for r in active:
                one=maximum({r:rows[r]},{r:scales[r]},d)
                assert max_chain(one,d,bits)==numbers[r]*source_roots[r]**factors*max_chain(rows[r],d,bits)
                checks+=1
            for a in envelope:
                b=d-a
                if b in envelope:
                    forward,reverse=join.two_orientations(envelope[a],envelope[b],bits)
                    assert forward==chain(envelope[a],envelope[b],bits)
                    assert reverse==chain(envelope[b],envelope[a],bits); checks+=2
    # Exact positive count checksum is separate from the synthetic root toys.
    base=load_base()
    for points in (6,8,10):
        weights=(2,3); values=base._band_counts(weights,points,join,guard)
        assert sum(values)==sum(comb(points,r) for r in weights); checks+=1
    toarb=lambda x:arb(x.numerator)/arb(x.denominator)
    directed={r:{f:tuple(tuple(toarb(x) for x in row) for row in value)
                 for f,value in rows.items()} for r,rows in original.items()}
    scales={r:toarb(source_roots[r]) for r in original}
    selected=selected_envelope(directed,scales,2,join,guard)
    assert selected["active_weights"]==(3,) and selected["countfactor"]==1
    assert set(selected["envelope"])=={1}; checks+=2
    n=join.run_counts(weight=3)[1]
    direct=join.two_orientations(directed[3][1],directed[3][1],diagonal=True)[0]
    value=join.two_orientations(selected["envelope"][1],selected["envelope"][1],diagonal=True)[0]
    assert value.overlaps(n*scales[3]**FACTORS*direct); checks+=1
    for r in (2,3):
        answer=counted_join({r:directed[r]},{r:scales[r]},join,guard)
        old=join.directed_join(tuple(directed[r][f] for f in range(r//2+1)),guard,weight=r)
        assert answer["bound"].overlaps(scales[r]**FACTORS*old["bound"]); checks+=1
    both=counted_join(directed,scales,join,guard)
    separate=sum((scales[r]**FACTORS*join.directed_join(
        tuple(directed[r][f] for f in range(r//2+1)),guard,weight=r)["bound"] for r in (2,3)),arb(0))
    assert both["bound"]>=separate or both["bound"].overlaps(separate); checks+=1
    sparse={2:{0:directed[2][0]},3:{1:directed[3][1]}}
    assert selected_envelope(sparse,scales,3,join,guard)["countfactor"]==0; checks+=1
    # At d=2, r3 cannot contribute f0: its partner f2 is invalid for r3.
    huge=((arb(1000),)*4,(arb(1000),)*2); small=directed[2][0]
    clipped=selected_envelope({3:{0:huge},4:{0:small,2:small}},
                             {3:arb(1),4:arb(1)},2,join,guard)
    root4=(arb(join.run_counts(weight=4)[2]).log()/FACTORS).exp()
    assert set(clipped["envelope"])=={0,2}
    assert clipped["envelope"][0][0][0].overlaps(root4*small[0][0]); checks+=2
    try: counted_join(sparse,scales,join,guard)
    except ValueError: checks+=1
    else: raise AssertionError("incomplete original coverage accepted")
    try: selected_envelope({2:{2:directed[2][0]}},{2:scales[2]},0,join,guard)
    except ValueError: checks+=1
    else: raise AssertionError("invalid actual full index accepted")
    return {"exact_fraction_and_public_interface_checks":checks,
        "noncommuting_terminal_controls":True,"perfect_power_count_toys":True,
        "singleton_reproduction":True,"sum_of_individual_bounds_control":True,
        "partial_envelope_rejected_by_complete_join":True,
        "actual_large_matrix_joins":0,"tests_require_parent_execution":True}
