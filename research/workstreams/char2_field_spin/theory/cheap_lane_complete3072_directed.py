"""Complete directed actual3072 replay of the frozen selected-region witness.

Fix z1/10,alpha119/1280 and the current ordered-pair local constructor.
Rebuild the same31-marker bank as selected_region_ab and REQUIRE its
recorded table/bank hashes. Fitting proposes exact rational witnesses
only; it cannot change the immutable local bound. Arb precision is384,
the minimum enforced by that local constructor, not a256/384 control.

Unique normalized-cell assignment covers every integer full count0..1536.
Build31 positive cyclic256 tables sequentially, with six H64/H63J
doublings each. Each actual target(f,3072-2f) removes its OWN multinomial
and marker denominator. Keep1537 normalized R4/V2 tuples; discard each
large polynomial table before the next. No r4096 endpoint substitution,
ordinary fallback, new alpha, source factor, or count enters this stage.

The frozen complete join covers2,362,368 ordered matching pairs, retains
both terminal orders and charges exact indexed run masses summing to
C(8192,3072). Multiply ONCE by(eta3072*beta^3072*z^(-L))^alpha afterward.
This is all actual3072 supports/histories only, not whole-code closure.
A positive final bound remains a valid but unsuccessful upper bound.

Parent execution only:1500/1800s,512MiB,20GiB free,-B,stdout,no generated
files. No work or proof-provider imports occur on module import.
"""

import argparse
from fractions import Fraction as F
import gc
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import time
from types import MappingProxyType

sys.dont_write_bytecode=True
SELECTED_PIN="cf4d025277a281aecfd00b105bd7f791ca00c9864838b4849f97d1da7a2de26b"
COUNTS_PIN="2468da0c91e6e44c7b707a1158252bb695818980de33eb7f270d620748b4f028"
TABLE_SHA="636b295acf8ba93a20f76a983ab4723c290f9967a252a04b9562148a3a7c539d"
BANK_SHA="3e4843f4abd6be8730f72a5f77c02f4b4e554bf32a834464bda3df8643f551e7"
WEIGHT,POINTS,PRECISION,MEMORY=3072,8192,384,512*(1<<20)
Z,ALPHA=F(1,10),F(119,1280)
ORDERED_PROFILES,UNORDERED_POWERS=2362368,1181952


def load_selected():
    path=Path(__file__).with_name("cheap_lane_selected_region_ab.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=SELECTED_PIN:
        raise RuntimeError("frozen selected-region baseline mismatch")
    spec=importlib.util.spec_from_file_location("complete3072_selected",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def schedule(cells,count_module,guard):
    """Exact frozen tie-to-lower assignment; no guessed rounded intervals."""
    clean=count_module._cells(cells)
    ranges=count_module.marginal_ranges(clean,WEIGHT,guard)
    out=[]; covered=[]; digest=hashlib.sha256()
    for cell,interval in zip(clean,ranges):
        guard()
        if interval is None:
            raise ArithmeticError("fixed3072 bank unexpectedly has an empty cell")
        center,low,high,u,w=cell; lo,hi=interval
        for f in range(lo,hi+1):
            x=F(2048*f,WEIGHT)
            if not low<=x<=high or (out and not low<x):
                raise ArithmeticError("assigned target leaves its unique normalized cell")
        out.append(((center,lo,hi,"cyclic"),(u,w)))
        covered.extend(range(lo,hi+1))
        digest.update((str(center)+":"+str(lo)+":"+str(hi)+"\n").encode("ascii"))
    if len(out)!=31 or covered!=list(range(WEIGHT//2+1)):
        raise ArithmeticError("incomplete or duplicate full-count assignment")
    return tuple(out),{"tables":len(out),"marginals":len(covered),
        "first_full":covered[0],"last_full":covered[-1],"tie_policy":"lower cell",
        "assignment_sha256":digest.hexdigest(),"all_actual_singletons_nonnegative":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits",type=int,choices=(192,256),default=192)
    parser.add_argument("--precision",type=int,choices=(384,),default=384)
    parser.add_argument("--seconds",type=int,choices=(1500,1800),default=1500)
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    selected=load_selected()
    if (selected.WEIGHT,selected.POINTS,selected.Z,selected.ALPHA)!=(WEIGHT,POINTS,Z,ALPHA):
        raise ArithmeticError("selected baseline parameters changed")
    run=selected.load_run(); actual=run.load_actual(); support=actual.load_support()
    normalized=support.load_normalized(); band=normalized.load_band()
    count_module=band.load("cheap_lane_cell_run_counts.py",COUNTS_PIN,"complete3072_assignment")
    fitter=band.load(support.FITTER_FILE,support.FITTER_PIN,"complete3072_fitter")
    reference=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"complete3072_reference")
    ordinary=reference.load("cheap_lane_ordered_ordinary_directed.py",reference.ORDINARY_PIN,
                            "complete3072_ordinary_controls")
    large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"complete3072_source")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"complete3072_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"complete3072_regions")
    integer=large.load("cheap_lane_integer_regions.py",large.INTEGER_PIN,"complete3072_integer")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"complete3072_join")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("complete3072 deadline; no incomplete coverage claim")
        if peak()>MEMORY: raise MemoryError("complete3072 exceeds512MiB budget")
    def progress(phase,**values):
        if phase=="actual_marker_fit" and values.get("center") not in (0,512,1024):
            return
        keys=("center","status","evaluations","first_full","last_full","tables_completed",
              "run","ordered_profiles","unordered_powers","working_precision")
        kept={key:values[key] for key in keys if key in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**kept}),flush=True)
    counts,own=selected.own_tests(extraction,join,anchor,guard)
    tests={"selected_controls":own,"extraction":extraction.exact_tests(guard),
        "integer":integer.exact_tests(guard),"ordered_terminal":join.exact_tests(guard),
        "ordinary":ordinary.exact_tests(guard),"bank":support.wrapper_tests(reference.BANK,guard)}
    built=local.build_local_pairs(Z,ALPHA,anchor,guard,deadline,PRECISION,progress=progress)
    v,core=built["verifier"],built["core"]; logtwo=v.arb(2).log()
    if built["report"]["working_precision"]!=PRECISION:
        raise ArithmeticError("frozen working precision no longer384")
    degrees=built["rational_coefficients"]; arb_identity=dict(built["coefficients"])
    table_hash=support.coefficient_fingerprint(degrees,guard)
    if table_hash!=TABLE_SHA or built["report"]["rational_coefficient_sha256"]!=TABLE_SHA:
        raise ArithmeticError("local table does not match selected3072 witness")
    proposed,fit_report=actual.fit_actual_markers(MappingProxyType(degrees),WEIGHT,fitter,guard,progress)
    fitted,public_bank,bank_hash=support.validate_bank(proposed,reference.BANK)
    if bank_hash!=BANK_SHA:
        raise ArithmeticError("regenerated bank does not match selected3072 witness")
    if support.coefficient_fingerprint(degrees,guard)!=TABLE_SHA:
        raise ArithmeticError("fitter changed the exact local table")
    if len(arb_identity)!=len(built["coefficients"]) or any(
            row is not built["coefficients"][key] for key,row in arb_identity.items()):
        raise ArithmeticError("fitter changed original Arb coefficients")
    rows,coverage=schedule(normalized.normalized_cells(fitted),count_module,guard)
    tests["complete_assignment"]={**coverage,"expected_ordered_profiles":ORDERED_PROFILES,
        "expected_shared_powers":UNORDERED_POWERS}
    if (WEIGHT//2+1)**2-1!=ORDERED_PROFILES:
        raise ArithmeticError("ordered-profile control constant is incorrect")
    regions=[None]*(WEIGHT//2+1); reports=[]
    for index,(item,marker) in enumerate(rows):
        guard(); center,lo,hi,_=item
        progress("complete3072_table_start",center=center,first_full=lo,last_full=hi)
        partial,report=extraction.one_table_many(v,degrees,item,marker,(WEIGHT,),
            integer,args.bits,guard,None)
        own_regions=partial[WEIGHT]
        if tuple(own_regions)!=tuple(range(lo,hi+1)):
            raise ArithmeticError("regional extraction lost an assigned integer target")
        for f,entry in own_regions.items():
            guard()
            if regions[f] is not None or tuple(map(len,entry))!=(4,2):
                raise ArithmeticError("duplicate target or invalid regional shape")
            regions[f]=entry
        reports.append({"center":center,"first_full":lo,"last_full":hi,
            "marginals":len(own_regions),"u":str(marker[0]),"v":str(marker[1]),
            "pair_doublings":report["pair_doublings"],
            "cyclic_DC_upper_checks":report["cyclic_DC_upper_checks"],
            "DC_vs_unrounded_log2_ratios":report["DC_vs_unrounded_log2_ratios"],
            "arithmetic":report["arithmetic"]})
        del partial,own_regions; gc.collect(); guard()
        progress("complete3072_table_complete",center=center,first_full=lo,last_full=hi,
                 tables_completed=index+1)
    if len(reports)!=31 or any(region is None for region in regions):
        raise ArithmeticError("incomplete3072 marginal cover; join forbidden")
    h,c_h,eta=anchor.source_factor(WEIGHT)
    if h!=12 or c_h!=F(1,13): raise ArithmeticError("actual3072 source interface changed")
    source=(v.exact(eta)*v.exact(core.BETA)**WEIGHT*v.exact(Z)**(-core.LIMIT))**v.exact(ALPHA)
    # Reuse already-extracted entries for the matched selected witness.
    # This is a broad inequality control, not equality to rounded decimals.
    for full,expected_center in ((1186,813),(1452,980)):
        assigned=next(item[0] for item,_ in rows if item[1]<=full<=item[2])
        if assigned!=expected_center:
            raise ArithmeticError("selected control no longer uses its recorded marker")
    selected_values=join.two_orientations(regions[1186],regions[1452])
    selected_controls=[]
    for fulls,value in zip(((1186,1452),(1452,1186)),selected_values):
        charged=v.positive(source*v.exact(F(counts[434]))*value)
        if not charged<v.arb(2)**(-3379):
            raise ArithmeticError("matched selected434-run diagnostic failed its recorded threshold")
        selected_controls.append({"fulls":list(fulls),"runs":434,
            "log2_selected_chain_with_full_run_mass_interval":str(charged.log()/logtwo),
            "below_2_pow_minus_3379":True,"not_added_separately_to_join":True})
    tests["matched_selected_control_from_complete_tables"]=True
    joined=join.directed_join(regions,guard,progress,weight=WEIGHT,points=POINTS)
    if (joined["report"]["ordered_profiles"]!=ORDERED_PROFILES
            or joined["report"]["unordered_powers"]!=UNORDERED_POWERS
            or not joined["report"]["all_run_counts_sum_to_binomial"]):
        raise ArithmeticError("complete ordered run-count coverage control failed")
    total=v.positive(source*joined["bound"])
    if not total>0: raise ArithmeticError("positive complete upper bound required")
    leaders=sorted(range(1,WEIGHT+1),key=lambda c:joined["run_bounds"][c-1].upper(),reverse=True)[:8]
    run_reports=[{"runs":c,"log2_source_charged_interval":
        str((source*joined["run_bounds"][c-1]).log()/logtwo),
        "indexed_run_count_sha256":hashlib.sha256(str(counts[c]).encode("ascii")).hexdigest()}
        for c in leaders]
    guard()
    print(json.dumps({"scope":"all actual3072-group indexed source supports and all histories at35%",
        "directed_arithmetic":True,"all_profiles_at_weight_covered":True,
        "whole_code_certificate":False,"log2_union_bound_interval":str(total.log()/logtwo),
        "actual3072_below_2_pow_minus_40":bool(total<v.arb(2)**(-40)),
        "actual3072_below_2_pow_minus_60":bool(total<v.arb(2)**(-60)),
        "parameters":{"weight":WEIGHT,"points":POINTS,"z":str(Z),"alpha":str(ALPHA),
            "regional_moduli":[256,256],"regional_integer_bits":args.bits,
            "working_Arb_precision":PRECISION,"seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "coverage":coverage,"join":joined["report"],"largest_runs":run_reports,
        "selected434_run_controls":selected_controls,
        "source_log2_interval":str(source.log()/logtwo),
        "joint_missing":{"h":h,"c_h":str(c_h),"eta":str(eta)},
        "source_once_after_join":True,"counts_outside_alpha":True,
        "region_normalizers_already_removed":True,"both_terminal_orders_retained":True,
        "regular_packets":128,"terminal_updates":127,"packet_updates":32767,
        "ordinary_fallbacks":0,"boundary_replacement":False,"degree_reversals":0,
        "coefficient_projection_exact":False,"positive_cyclic_aliases_remain":True,
        "shared_bits_not_relative_coefficient_precision":True,
        "local_table_builds":1,"local_coefficient_sha256":table_hash,
        "positive_rational_bank_sha256":bank_hash,"baseline_hashes_enforced":True,
        "positive_rational_bank":public_bank,"tables":reports,
        "fit":{"scope":"same frozen FLOAT proposal procedure; final witnesses are exact rationals",
            "successful_proposal_evaluations_recorded":fit_report["successful_proposal_evaluations_recorded"],
            "fit_failures":fit_report["fit_failures"]},
        "tests":tests,"local_coefficients":built["report"],
        "source_count_checksum":own["run_count_sum_sha256"],
        "nonnegative_intersections":v.intersections,
        "elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "pins":{"selected":SELECTED_PIN,"counts":COUNTS_PIN,"run":selected.RUN_PIN,
            "actual":run.ACTUAL_PIN,"support":actual.SUPPORT_PIN,"fitter":support.FITTER_PIN,
            "normalized":support.NORMALIZED_PIN,"band":normalized.BAND_PIN,
            "reference":band.HOLE_PIN,"ordinary":reference.ORDINARY_PIN,
            "local":large.LOCAL_PIN,"anchor":large.ANCHOR_PIN,"regions":large.REGION_PIN,
            "integer":large.INTEGER_PIN,"join":large.JOIN_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
