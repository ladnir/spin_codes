"""Directed ordinary-marker selection control for all actual4096 supports.

Keep the frozen local ordered-pair table, z1/16,alpha119/1280 and31 fixed
positive rational markers. Every marker gives an ordinary H64/H63J upper
bound for EVERY target, after division by that target's own denominator.
Select each of the six entries independently by the smallest certified
upper endpoint, retaining the selected enclosure. This is selection among
valid majorants, not intersection of enclosures of an unknown probability.

The target denominator starts at v^4096 and advances by
 (u/v^2)*(4096-2f)*(4095-2f)/(f+1)^2.
No advance is made after f2048. Exactf0 is then restored. The unchanged
positive directed join covers all4,198,400 ordered profiles and charges
the actual source factor once. No extra alpha or marker union factor.

This measures selection loss for one source size, not whole-code closure
or the success of a fixed-cell size-band envelope. Numerical receipts are
indexed separately; source alone is not evidence. Parent execution only:
300s/512MiB/20GiB free, -B, stdout, no generated files or work on import.
"""

import argparse
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
ORDINARY_PIN="cb2cc15b9209a95cdfd57b51d98baed99aa205a82099595075c2e72563b826ed"
MEMORY=512*(1<<20)
CHECK_FULL=(0,1,14,36,512,1024,1980,2013,2048)
CHANNELS=("R00","R01","R10","R11","V0","V1")


def load_ordinary():
    path=Path(__file__).with_name("cheap_lane_ordered_ordinary_directed.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=ORDINARY_PIN:
        raise RuntimeError("frozen ordinary directed source mismatch")
    spec=importlib.util.spec_from_file_location("marker_envelope_ordinary",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def candidate_wins(current,candidate):
    """Upper-endpoint choice only; never intersect unrelated majorants."""
    return current is None or candidate.upper()<current.upper()


def exact_tests(guard):
    """Rational recurrence and positive products, independent of Arb helpers."""
    checks=0
    for slots in (2,4,8,12):
        for u,w in ((F(1),F(1)),(F(2,3),F(5,7)),(F(11,4),F(3,2))):
            denominator=w**slots
            for f in range(slots//2+1):
                guard(); s=slots-2*f
                direct=F(comb(slots,f)*comb(slots-f,s))*u**f*w**s
                assert denominator==direct and denominator>0; checks+=1
                if f<slots//2:
                    denominator*=u/(w*w)*F(s*(s-1),(f+1)**2)
            assert denominator==F(comb(slots,slots//2))*u**(slots//2); checks+=1

    def mm(a,b):
        return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(len(b))),F(0))
                           for j in range(len(b[0]))) for i in range(len(a)))
    def le(a,b):
        return all(x<=y for ar,br in zip(a,b) for x,y in zip(ar,br))
    base=((F(1,7),F(2,9)),(F(3,11),F(5,13)))
    terminal=((F(7,17),),(F(11,19),))
    candidates=[]; terminals=[]
    for marker in range(3):
        candidates.append(tuple(tuple(base[i][j]+F((marker+i+2*j)%3+1,17)
                                      for j in range(2)) for i in range(2)))
        terminals.append(tuple((terminal[i][0]+F((2*marker+i)%3+1,23),) for i in range(2)))
    selected=tuple(tuple(min(x[i][j] for x in candidates) for j in range(2)) for i in range(2))
    selected_t=tuple((min(x[i][0] for x in terminals),) for i in range(2))
    assert le(base,selected) and le(terminal,selected_t); checks+=2
    br,bt=base,terminal; sr,st=selected,selected_t
    cr=list(candidates); ct=list(terminals)
    for _ in range(3):
        guard(); bt=mm(br,bt); br=mm(br,br)
        st=mm(sr,st); sr=mm(sr,sr)
        assert le(br,sr) and le(bt,st); checks+=2
        for marker in range(3):
            ct[marker]=mm(cr[marker],ct[marker]); cr[marker]=mm(cr[marker],cr[marker])
            assert le(sr,cr[marker]) and le(st,ct[marker]); checks+=2
    assert mm(candidates[0],candidates[1])!=mm(candidates[1],candidates[0]); checks+=1

    class Interval:
        def __init__(self,lo,hi): self.lo,self.hi=F(lo),F(hi)
        def upper(self): return self.hi
    wide=Interval(0,3); narrow=Interval(2, F(5,2)); tie=Interval(1,F(5,2))
    assert candidate_wins(None,wide) and candidate_wins(wide,narrow)
    assert not candidate_wins(narrow,wide) and not candidate_wins(narrow,tie)
    checks+=4
    return {"exact_fraction_checks":checks,"recurrence_includes_both_endpoints":True,
        "no_zero_final_advance":True,"independent_entry_min_positive_products":True,
        "terminal_minimum_independent":True,"upper_endpoint_not_midpoint_selection":True}


def build_envelope(v,coefficients,ordinary,large,extraction,guard,progress):
    """Stream marker candidates; retain only best and nearest-cell controls."""
    nearest=[None]*2049
    for center,lo,hi,_ in large.schedule():
        for f in range(lo,hi+1):
            assert nearest[f] is None
            nearest[f]=center
    assert all(x is not None for x in nearest)
    best=[[None]*6 for _ in range(2049)]
    selected=[[None]*6 for _ in range(2049)]
    assigned=[None]*2049
    recurrence_checks=0; candidate_entries=0
    check_full=set(CHECK_FULL)
    for center,marker in large.WITNESSES.items():
        guard(); h,j=ordinary.evaluate_pairs(v,coefficients,marker,guard)
        regular=v.matpow(h,64); terminal=v.matmul(v.matpow(h,63),j)
        numerator=tuple(v.positive(x) for matrix in (regular,terminal) for row in matrix for x in row)
        assert len(numerator)==6
        u,w=map(F,marker)
        denominator=v.exact(w)**4096
        ratio_marker=v.exact(u/(w*w))
        updates=0
        for f in range(2049):
            guard()
            if not denominator.is_finite() or not denominator>0:
                raise ArithmeticError("denominator recurrence lost strict positivity")
            if f in check_full:
                direct=v.exact(extraction.target_denominator(4096,f,marker))
                assert denominator.overlaps(direct)
                recurrence_checks+=1
            own=tuple(v.positive(x/denominator) for x in numerator)
            if nearest[f]==center:
                assert assigned[f] is None
                assigned[f]=own
            for k,candidate in enumerate(own):
                if not candidate.is_finite():
                    raise ArithmeticError("nonfinite normalized ordinary candidate")
                if candidate_wins(best[f][k],candidate):
                    best[f][k]=candidate; selected[f][k]=center; updates+=1
            candidate_entries+=6
            if f<2048:
                s=4096-2*f
                denominator*=ratio_marker*v.exact(F(s*(s-1),(f+1)**2))
        progress("ordinary_marker_all_targets",center=center,u=str(u),v=str(w),
                 targets=2049,entry_selections=updates)
        del h,j,regular,terminal,numerator,denominator,ratio_marker,own
        gc.collect()
    assert recurrence_checks==len(large.WITNESSES)*len(CHECK_FULL)
    assert candidate_entries==len(large.WITNESSES)*2049*6
    for f in range(2049):
        assert assigned[f] is not None and all(x is not None for x in best[f])
        assert all(x.upper()<=y.upper() for x,y in zip(best[f],assigned[f]))
    histograms={name:{str(center):0 for center in large.WITNESSES} for name in CHANNELS}
    changed=0; strictly_smaller=0
    for f in range(1,2049):
        for k,name in enumerate(CHANNELS):
            histograms[name][str(selected[f][k])]+=1
            changed+=selected[f][k]!=nearest[f]
            strictly_smaller+=best[f][k].upper()<assigned[f][k].upper()
    assert all(sum(row.values())==2048 for row in histograms.values())
    regions=[(tuple(row[:4]),tuple(row[4:])) for row in best]
    return regions,{"markers":len(large.WITNESSES),"targets_per_marker":2049,
        "candidate_entries":candidate_entries,"exact_denominator_overlap_checks":recurrence_checks,
        "crosschecked_full_targets":list(CHECK_FULL),"histograms_exclude_exact_f0":True,
        "selected_marker_histograms":histograms,"entries_changed_from_nearest_marker":changed,
        "strictly_smaller_upper_endpoints_than_nearest":strictly_smaller,
        "entry_selection":"smallest certified upper endpoint; retain chosen enclosure",
        "all_selected_uppers_at_most_nearest_candidate_uppers":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    parser.add_argument("--precision",type=int,choices=(256,384),default=384)
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    ordinary=load_ordinary(); large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"envelope_anchor")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"envelope_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"envelope_normalizers")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"envelope_join")
    screen=large.load("cheap_lane_ordered_large_screen.py",large.SCREEN_PIN,"envelope_boundary_tests")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("marker envelope deadline; no incomplete cover")
        if peak()>MEMORY: raise MemoryError("marker envelope512MiB budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    tests=exact_tests(guard)
    tests["ordinary_fraction_controls"]=ordinary.exact_tests(guard)
    tests["frozen_boundary_extraction_join"]=large.exact_tests(screen,extraction,join,guard)
    built=local.build_local_pairs(large.Z,large.ALPHA,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]; logtwo=v.arb(2).log()
    assert large.WEIGHT==4096 and large.POINTS==8192 and len(large.WITNESSES)==31
    regions,selection=build_envelope(v,built["coefficients"],ordinary,large,extraction,guard,progress)
    boundary,boundary_report=large.singleton_boundary(v,built["rational_coefficients"],join,guard)
    for row,exactrow in zip(regions[0],boundary):
        for upper,exact in zip(row,exactrow):
            assert upper>=exact or upper.overlaps(exact)
    regions[0]=boundary
    tests["selected_ordinary_dominates_exact_f0_entries"]=6
    counts=join.run_counts(points=8192,weight=4096)
    assert sum(counts)==comb(8192,4096) and counts[4096]==2
    h,c_h,eta=anchor.source_factor(4096)
    assert h==16 and c_h==F(1,17)
    source=(v.exact(eta)*v.exact(core.BETA)**4096*v.exact(large.Z)**(-core.LIMIT))**v.exact(large.ALPHA)
    selected=join.two_orientations(regions[2013],regions[36])
    controls=[{"full":full,"runs":2047,"not_a_run_or_size_cover":True,
        "log2_selected_with_full_run_mass_interval":str((source*counts[2047]*value).log()/logtwo)}
        for full,value in zip(((2013,36),(36,2013)),selected)]
    joined=join.directed_join(regions,guard,progress,weight=4096,points=8192)
    assert joined["report"]["ordered_profiles"]==4198400
    total=v.positive(source*joined["bound"])
    worst=max(range(1,4097),key=lambda c:joined["run_bounds"][c-1].upper())
    guard()
    print(json.dumps({"scope":"complete actual4096 ordinary entrywise marker envelope",
        "directed_arithmetic":True,"all_profiles_at_weight_covered":True,
        "whole_code_certificate":False,"other_weights_covered":False,
        "fixed_cell_size_band_claim":False,"parameters":{"weight":4096,
            "z":str(large.Z),"alpha":str(large.ALPHA),"precision":built["report"]["working_precision"],
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "regional_method":"ordinary H64/H63J; all31 own-target candidates; independent six-entry selection",
        "selection":selection,"exact_f0_replacement":True,"residue_projection":False,
        "source_once_counts_outside_alpha":True,"new_local_alpha":False,"marker_union_factor":False,
        "log2_union_bound_interval":str(total.log()/logtwo),
        "actual4096_below_2_pow_minus_40":bool(total<v.arb(2)**(-40)),
        "actual4096_below_2_pow_minus_60":bool(total<v.arb(2)**(-60)),
        "worst_run":{"runs":worst,"log2_interval":str((source*joined["run_bounds"][worst-1]).log()/logtwo)},
        "selected_controls":controls,"boundary":boundary_report,"join":joined["report"],
        "source_log2_interval":str(source.log()/logtwo),"joint_missing":{"h":h,"c_h":str(c_h),"eta":str(eta)},
        "tests":tests,"local_coefficients":built["report"],"generated_bytes":0,
        "elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),
        "nonnegative_intersections":v.intersections,"pins":{"ordinary":ORDINARY_PIN,
            "large":ordinary.LARGE_PIN,"local":large.LOCAL_PIN,"normalizers":large.REGION_PIN,
            "join":large.JOIN_PIN,"anchor":large.ANCHOR_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
