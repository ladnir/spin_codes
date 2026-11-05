"""Complete FLOAT profile diagnosis of the fixed ordinary r4096 bound.

Rebuild the hash-pinned directed local ordered-pair coefficients and the
same31 rational marker evaluations as ordered_ordinary_directed. Each
ordinary H^64/H^63J numerator is divided by its own exact actual target
denominator. Replace only f0 by the independent top-degree projection.
Then export the positive Arb regional entries as outward normalized
FLOAT tuples and run the frozen complete ordered-profile join.

The export avoids direct conversion of tiny or huge regional entries.
Use an integer log scale and outward FLOAT upper endpoints for exp of
the centered log enclosures. Exact zeros remain zero; ambiguous signs,
nonfinite values, or positive underflow abort. Subsequent FLOAT matrix
arithmetic is not directed, so even a negative total is NOT a certificate.

Optional --edge-full1..8 replaces f0..edge_full by exact FULL/EMPTY
projections from a separately pinned helper. It changes no other region,
source factor, count, marker witness, or local cap. Its dependency pin and
exact tests are checked before local construction. Edge0 reproduces
the ordinary diagnostic with only the existing exact f0 replacement.

All2049 marginals and4,198,400 ordered profiles are retained. The source
factor uses actual4096 once; exact run masses sum to C(8192,4096).
Leading FLOAT channels receive separate directed selected-chain checks,
but those checks do not cover their whole run classes or the source size.
No numerical work on import. Parent-only300s,512MiB,20GiB free,-B,
stdout only, zero generated files, and no frozen edits.
"""

import argparse
from fractions import Fraction as F
import gc
import hashlib
import importlib.util
import json
import math
from math import comb
from pathlib import Path
import shutil
import sys
import time

sys.dont_write_bytecode=True
ORDINARY_PIN="cb2cc15b9209a95cdfd57b51d98baed99aa205a82099595075c2e72563b826ed"
EDGE_PIN="2fbe2e1220f083914aa4419939e84927f7c61573163f6d56250e250c7a5a5c42"
MEMORY=512*(1<<20)
LN2=math.log(2)


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if not pin or hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def upper_float(value,v):
    """Finite binary64 upper endpoint; conversion rounding is checked exactly."""
    if not value.is_finite(): raise ArithmeticError("nonfinite export enclosure")
    endpoint=value.upper().fmpq(); result=float(endpoint)
    if not math.isfinite(result): raise ArithmeticError("FLOAT endpoint overflow")
    if endpoint!=0 and result==0: raise ArithmeticError("positive or signed FLOAT endpoint underflow")
    while v.fmpq(*result.as_integer_ratio())<endpoint:
        result=math.nextafter(result,math.inf)
        if not math.isfinite(result): raise ArithmeticError("outward FLOAT endpoint overflow")
    return result


def normalized_upper(items,v):
    """Export nonnegative R4 or V2 without converting its absolute scale."""
    if len(items) not in (2,4): raise ValueError("expected flat R4 or V2")
    logs=[]
    for value in items:
        if not value.is_finite() or value.upper()<0:
            raise ArithmeticError("invalid regional entry for FLOAT export")
        if value==0: logs.append(None)
        elif value>0: logs.append(value.log())
        else: raise ArithmeticError("regional sign is not resolved for FLOAT export")
    finite=[value for value in logs if value is not None]
    if not finite: raise ArithmeticError("zero regional matrix object is unsupported")
    scale=math.ceil(max(upper_float(value,v) for value in finite))
    if abs(scale)>1<<53: raise ArithmeticError("integer log scale is not exactly representable")
    exact_scale=v.exact(F(scale)); entries=[]
    for value in logs:
        if value is None: entries.append(0.0); continue
        entry=upper_float((value-exact_scale).exp(),v)
        if not 0<entry<=1+1e-12:
            raise ArithmeticError("positive normalized entry underflow or scale failure")
        entries.append(entry)
    return tuple(entries),float(scale)


def export_region(region,v):
    if len(region)!=2 or tuple(map(len,region))!=(4,2):
        raise ValueError("expected a normalized flat Arb R4/V2 pair")
    # The exact target denominator was already divided out. The frozen
    # FLOAT join must therefore subtract0, not a second marker denominator.
    return normalized_upper(region[0],v),normalized_upper(region[1],v),0.0,F(1),F(1)


def export_tests(v,guard):
    checks=0
    for power in (-2048,0,2048):
        guard(); base=F(1<<power) if power>=0 else F(1,1<<(-power))
        values=(v.exact(base),v.arb(0),v.exact(3*base),v.exact(base/7))
        entries,scale=normalized_upper(values,v)
        assert entries[1]==0.0 and scale==int(scale)
        factor=v.exact(F(int(scale))).exp()
        for actual,coefficient in zip(values,entries):
            upper=v.exact(F(coefficient))*factor
            assert upper>=actual or upper.overlaps(actual)
            checks+=1
        region=export_region((values,values[:1]+values[2:3]),v)
        assert region[2:]==(0.0,F(1),F(1)); checks+=1
    return {"Arb_export_domination_checks":checks,"exact_zeros_preserved":True,
        "tiny_and_large_absolute_scales":True,"already_removed_denominator_not_reapplied":True,
        "subsequent_FLOAT_join_not_certified":True}


def build_ordinary_regions(v,coefficients,large,ordinary,extraction,guard,progress):
    """All ordinary normalized marginals, before any endpoint replacement."""
    regions=[None]*2049; reports=[]; assignments={}
    for center,lo,hi,_ in large.schedule():
        guard(); marker=large.WITNESSES[center]
        h,j=ordinary.evaluate_pairs(v,coefficients,marker,guard)
        regular=v.matpow(h,64); terminal=v.matmul(v.matpow(h,63),j)
        numerator=tuple(tuple(v.positive(x) for row in matrix for x in row)
                        for matrix in (regular,terminal))
        for f in range(lo,hi+1):
            guard(); denominator=v.exact(extraction.target_denominator(4096,f,marker))
            assert regions[f] is None and denominator>0
            regions[f]=tuple(tuple(v.positive(x/denominator) for x in row) for row in numerator)
            assignments[f]=center
        report={"center":center,"first_full":lo,"last_full":hi,
            "markers":list(map(str,marker)),"marginals":hi-lo+1,
            "method":"ordinary H64/H63J, own exact denominator already removed"}
        reports.append(report); progress("ordinary_profile_cell_complete",**report)
        del h,j,regular,terminal,numerator,denominator
        gc.collect()
    assert all(entry is not None for entry in regions) and len(assignments)==2049
    return regions,reports,assignments


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edge-full",type=int,choices=range(9),default=0)
    parser.add_argument("--precision",type=int,choices=(256,384),default=384)
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    args=parser.parse_args()
    if args.edge_full and not EDGE_PIN:
        parser.error("edge projection helper is not yet pinned; use --edge-full0")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    ordinary=load("cheap_lane_ordered_ordinary_directed.py",ORDINARY_PIN,"profile_ordinary")
    edge=load("cheap_lane_edge_regions_directed.py",EDGE_PIN,"profile_exact_edges") if args.edge_full else None
    large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"profile_anchor")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"profile_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"profile_normalizers")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"profile_selected_checks")
    old_screen=large.load("cheap_lane_ordered_large_screen.py",large.SCREEN_PIN,"profile_large_screen")
    sizes=old_screen.load("cheap_lane_pair_sizes_screen.py",old_screen.SIZES_PIN,"profile_sizes")
    broad=sizes.load("cheap_lane_broad_peeling_screen.py",sizes.BROAD_PIN,"profile_batched")
    neighbor=broad.load("cheap_lane_neighbor_screen.py",broad.NEIGHBOR_PIN,"profile_neighbor")
    cells=neighbor.load_cells(); fft=cells.load_fft(); screen=fft.load_screen()
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("ordinary profile deadline; no incomplete cover")
        if peak()>MEMORY: raise MemoryError("512MiB ordinary profile budget exceeded")
    def progress(phase,**values):
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**values}),flush=True)
    import numpy as np
    tests={"ordinary_rational":ordinary.exact_tests(guard),
        "boundary_extraction_join":large.exact_tests(old_screen,extraction,join,guard)}
    if edge is not None: tests["exact_FULL_EMPTY_projection"]=edge.exact_tests(guard)
    built=local.build_local_pairs(large.Z,large.ALPHA,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]; logtwo=v.arb(2).log()
    tests["FLOAT_export"]=export_tests(v,guard)
    regions,tables,assignments=build_ordinary_regions(v,built["coefficients"],large,ordinary,
                                                     extraction,guard,progress)
    boundary,boundary_report=large.singleton_boundary(v,built["rational_coefficients"],join,guard)
    regions[0]=boundary
    edge_report=None
    if args.edge_full:
        refined,edge_report=edge.edge_regions(v,built["rational_coefficients"],args.edge_full,
                                             guard,progress=progress)
        if tuple(sorted(refined))!=tuple(range(args.edge_full+1)):
            raise ArithmeticError("edge helper returned an incomplete or different projection range")
        for f,entry in refined.items():
            if len(entry)!=2 or tuple(map(len,entry))!=(4,2): raise ValueError("invalid edge region shape")
            for upper_row,projected_row in zip(regions[f],entry):
                for upper,projected in zip(upper_row,projected_row):
                    if not (upper>=projected or upper.overlaps(projected)):
                        raise ArithmeticError("exact edge exceeds its ordinary own-target bound")
            if f==0:
                if not all(a.overlaps(b) for ar,br in zip(entry,boundary) for a,b in zip(ar,br)):
                    raise ArithmeticError("edge f0 does not overlap frozen exact endpoint")
            regions[f]=entry
    converted=[export_region(region,v) for region in regions]
    assert len(converted)==2049
    h,c_h,eta=anchor.source_factor(4096)
    assert h==16 and c_h==F(1,17)
    source=(v.exact(eta)*v.exact(core.BETA)**4096*v.exact(large.Z)**(-core.LIMIT))**v.exact(large.ALPHA)
    log_source=upper_float(source.log(),v)
    progress("ordinary_profile_marginals_ready",marginals=2049,edge_full=args.edge_full,
        local_working_precision=built["report"]["working_precision"])
    result=sizes.complete_join(4096,converted,[],broad,screen,np,log_source,guard,progress)
    assert result["ordered_profiles"]==4198400 and result["actual_weight"]==4096
    # Rename only the frozen join's method label: these are primary ordinary
    # bounds, not resolved FFT residues. Counts and numerical values are unchanged.
    categories=result["method_categories"]
    result["method_categories"]={"ordinary_or_exact_projection":categories["resolved_residue"],
                                 "uses_ordinary_fallback":categories["uses_ordinary_fallback"]}
    counts=join.run_counts(points=8192,weight=4096)
    assert sum(counts)==comb(8192,4096)
    leaders=[]
    for row in result["largest_run_bounds"]:
        guard(); a,b=row["actual_full"]; c=4096-a-b
        assert c==row["runs"]
        forward,reverse=join.two_orientations(regions[a],regions[b],diagonal=a==b)
        charged=v.positive(source*counts[c]*forward)
        marker_rows=[]
        for f in (a,b):
            method="exact FULL/EMPTY projection" if args.edge_full and f<=args.edge_full else \
                   "exact top-degree f0" if f==0 else "ordinary Cauchy"
            marker_rows.append({"full":f,"single":4096-2*f,"method":method,
                "ordinary_center":assignments[f],
                "ordinary_markers":list(map(str,large.WITNESSES[assignments[f]])),
                "marker_used_for_final_region":method=="ordinary Cauchy"})
        leaders.append({"full":[a,b],"runs":c,"marginals":marker_rows,
            "log2_selected_chain_with_full_run_mass_interval":str(charged.log()/logtwo),
            "reverse_log2_selected_chain_with_full_run_mass_interval":str((source*counts[c]*reverse).log()/logtwo),
            "FLOAT_log2_run_maximum":row["log2_run_envelope_FLOAT"],
            "FLOAT_minus_directed_selected_midpoint_bits":row["log2_run_envelope_FLOAT"]-float((charged.log()/logtwo).mid()),
            "not_a_directed_run_or_size_cover":True})
    guard()
    print(json.dumps({"scope":"complete actual4096 FLOAT ordinary-profile diagnostic; selected leaders checked with Arb",
        "directed_complete_join":False,"negative_is_not_a_certificate":True,
        "whole_code_certificate":False,"other_weights_covered":False,
        "parameters":{"weight":4096,"z":str(large.Z),"alpha":str(large.ALPHA),
            "edge_full":args.edge_full,"working_precision":built["report"]["working_precision"],
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "result":result,"leader_details":leaders,"tables":tables,
        "boundary":boundary_report,"edge_projection":edge_report,
        "ordinary_fallbacks":0,"ordinary_is_primary_method_not_failure_fallback":True,
        "source_once_counts_outside_alpha":True,"normalizers_removed_before_FLOAT_export":True,
        "FLOAT_export":"upper log enclosures, integer scale, outward positive binary64 entries; exact zeros preserved",
        "source_log2_interval":str(source.log()/logtwo),
        "joint_missing":{"h":h,"c_h":str(c_h),"eta":str(eta)},
        "tests":tests,"local_coefficients":built["report"],"nonnegative_intersections":v.intersections,
        "pins":{"ordinary":ORDINARY_PIN,"edge":EDGE_PIN,"large":ordinary.LARGE_PIN,
            "local":large.LOCAL_PIN,"normalizers":large.REGION_PIN,"selected_join":large.JOIN_PIN,
            "sizes":old_screen.SIZES_PIN,"batched":sizes.BROAD_PIN,"anchor":large.ANCHOR_PIN},
        "seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
