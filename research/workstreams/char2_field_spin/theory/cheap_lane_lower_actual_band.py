"""Ordinary directed actual-size bands for225<=r<=4096.

The old641 lower guard came from a hole-fitting interface, not from the
actual-support bound. This adapter supplies its own lower-domain fitting
targets and bounded proposal loop. It leaves every frozen helper unchanged.
The complete actual-degree pair table is built once at the requested z and
alpha. No degree reversal, scalar-message restriction, or new alpha occurs.

The fitted bank contains31 positive rational witnesses. FLOAT fitting is
only witness selection; all proof numerators use the preserved directed
coefficients. Each actual target has its own multinomial/marker denominator.
The256 denominator tilts supply beta^(alpha*r). The actual eta_r^alpha and
indexed run-interval masses enter once, and the output factor enters once.
Cell count intervals can overlap: they are upper counts, not exact joint
profile counts. Both off-diagonal terminal orders remain in the union.

Importable validate, fit_actual_markers, load_dependencies and exact_tests
also support a later residue/replay adapter. evaluate_band consumes already
regional UNNORMALIZED R4/V2 numerators; it does not power them again.
The CLI currently builds ordinary H^64/H^63J numerators only. --checks-only
runs exact interface tests, with no local table, fitting, or Arb evaluation.
The exact old3072 witness additionally checks frozen fitter, ordinary
numerators and frozen band assembly, not merely hashes. Other anchors do
not claim that replay.

Source/receipt scope is distinct: source alone certifies no size. Parent
execution only, -B, stdout,300/600seconds,512MiB and20GiB free; no files.
There is no provider import or numerical work at module import.
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
from types import MappingProxyType, SimpleNamespace


GENERATOR_PIN="ed6745fb570200c27300a0da5e4c179ae4c998483c026592163f1428272fd619"
RUN_PIN="19914e8bc698c410971f8edb0ab7f884cb0b83955bee2a6933e9791aa9dfcfc2"
ACTUAL_PIN="ca72174f0fc19848e02de159093f476eb61369629d8dc6f1271332171e7ae635"
COUNTS_PIN="2468da0c91e6e44c7b707a1158252bb695818980de33eb7f270d620748b4f028"
OLD_TABLE="636b295acf8ba93a20f76a983ab4723c290f9967a252a04b9562148a3a7c539d"
OLD_BANK="3e4843f4abd6be8730f72a5f77c02f4b4e554bf32a834464bda3df8643f551e7"
MINIMUM,MAXIMUM,POINTS,REFERENCE=225,4096,8192,2048
MEMORY,PRECISION=512*(1<<20),384
TEST_WEIGHTS=(225,254,255,256,257,258,507,508,511,512,515,516,
              640,641,767,768,1023,3072)


def _load(filename,pin,label):
    path=Path(__file__).with_name(filename)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen lower-domain dependency mismatch: "+filename)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def load_dependencies():
    """Load pinned inert interfaces only; invoke no builder, fitter or test."""
    generator=_load("cheap_lane_recalibrated_residue_max_band.py",GENERATOR_PIN,"lower_generator")
    run=_load("cheap_lane_calibrated_run_band.py",RUN_PIN,"lower_run_band")
    if (generator.RUN_PIN,generator.COUNTS_PIN,run.ACTUAL_PIN,run.COUNTS_PIN)!=(
            RUN_PIN,COUNTS_PIN,ACTUAL_PIN,COUNTS_PIN):
        raise ArithmeticError("frozen domain/count provenance changed")
    actual=run.load_actual(); support=actual.load_support()
    normalized=support.load_normalized(); band=normalized.load_band()
    counts=band.load("cheap_lane_cell_run_counts.py",COUNTS_PIN,"lower_cell_counts")
    fitter=band.load(support.FITTER_FILE,support.FITTER_PIN,"lower_fit_utilities")
    geometry=band.load(normalized.GEOMETRY_FILE,normalized.GEOMETRY_PIN,"lower_geometry")
    reference=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"lower_reference_cells")
    ordinary=reference.load("cheap_lane_ordered_ordinary_directed.py",reference.ORDINARY_PIN,"lower_ordinary")
    large=ordinary.load_large()
    source=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"lower_source")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"lower_normalizer")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"lower_terminal")
    pins={"generator":GENERATOR_PIN,"run":RUN_PIN,"actual":ACTUAL_PIN,"counts":COUNTS_PIN,
        "support":actual.SUPPORT_PIN,"fitter":support.FITTER_PIN,"normalized":support.NORMALIZED_PIN,
        "band":normalized.BAND_PIN,"geometry":normalized.GEOMETRY_PIN,"reference":band.HOLE_PIN,
        "ordinary":reference.ORDINARY_PIN,"source":large.ANCHOR_PIN,"normalizer":large.REGION_PIN,
        "join":large.JOIN_PIN,"local":large.LOCAL_PIN}
    return SimpleNamespace(generator=generator,run=run,actual=actual,support=support,
        normalized=normalized,band=band,counts=counts,fitter=fitter,geometry=geometry,
        reference=reference,ordinary=ordinary,large=large,source=source,
        extraction=extraction,join=join,pins=pins)


def rational64(text):
    if not isinstance(text,str) or len(text)>128 or "e" in text.lower():
        raise argparse.ArgumentTypeError("use a bounded exact fraction or decimal")
    try: value=F(text)
    except (ValueError,ZeroDivisionError) as error:
        raise argparse.ArgumentTypeError("expected an exact rational") from error
    if max(abs(value.numerator).bit_length(),value.denominator.bit_length())>64:
        raise argparse.ArgumentTypeError("rational witnesses limited to64 bits")
    return value


def validate(anchor,first,last,z,alpha):
    if any(type(x) is not int for x in (anchor,first,last)):
        raise TypeError("anchor and band endpoints must be integer actual sizes")
    if not MINIMUM<=anchor<=MAXIMUM or not MINIMUM<=first<=last<=MAXIMUM:
        raise ValueError("actual anchor and band must lie in225..4096")
    if any(isinstance(x,bool) or not isinstance(x,(F,int)) for x in (z,alpha)):
        raise TypeError("z and alpha must be exact rationals")
    if not 0<z<1 or not 0<alpha<=1:
        raise ValueError("require0<z<1 and0<alpha<=1")


def _interior_target(actual,center):
    """Proposal algebra extends to r2; this is not the public proof domain."""
    if type(actual) is not int or not 2<=actual<=4096:
        raise ValueError("interior proxy algebra requires2..4096")
    if type(center) is not int or not 0<=center<=1024:
        raise ValueError("normalized center must be0..1024")
    full=F(actual*center,REFERENCE)
    if center==0: full=F(1,2)
    elif center==1024: full=F(actual-1,2)
    single=actual-2*full
    if min(full,single,4096-actual+full)<=0:
        raise ArithmeticError("proposal targets must have positive full/single/empty counts")
    return full,single


def fit_target(actual,center):
    if type(actual) is not int or not MINIMUM<=actual<=MAXIMUM:
        raise ValueError("actual fitting anchor must be225..4096")
    return _interior_target(actual,center)


def fit_actual_markers(rational_degrees,actual,fitter,guard,progress=None):
    """Frozen proposal arithmetic with own targets; no frozen global changes."""
    fit_target(actual,0)
    template,control,moments,sizes,pins=fitter._dependencies(guard)
    import numpy as np
    table,fingerprint,positives,zero_count=fitter._prepare(rational_degrees,guard)
    logs=np.full((65,65,6),-np.inf,dtype=float)
    for (f,s),row in table.items():
        guard(); logs[f,s]=tuple(fitter._log_fraction(value) for value in row)
        if row!=rational_degrees[f,s]:
            raise ArithmeticError("actual degree labels changed during proposal conversion")
    if np.isnan(logs).any() or np.isposinf(logs).any():
        raise ArithmeticError("invalid actual-degree proposal log table")
    del table
    objective=control.Objective(logs,0.0,moments,np,guard)
    bank=[]; reports=[]; failures=[]; evaluations=0
    for center,lo,hi,_,_ in template:
        guard(); full,single=fit_target(actual,center)
        objective.targets=256*np.asarray([float(full),float(single)])
        seed=tuple(fitter._positive_rational(value) for value in sizes.markers(actual,full))
        marker=seed; history=[]; rounded=None
        status="positive rational seed retained"
        try:
            best,history,status=control.coordinate_fit(objective,
                np.asarray([fitter._log_fraction(value) for value in seed]),np,guard)
            marker=tuple(fitter._positive_rational(F.from_float(math.exp(float(x)))) for x in best["x"])
            rounded=objective.evaluate(np.asarray([fitter._log_fraction(value) for value in marker]))
        except (ArithmeticError,ValueError) as error:
            marker=seed; rounded=None
            status="FLOAT fit failed; positive rational seed retained"
            failures.append({"center":center,"reason":type(error).__name__+": "+str(error)})
        if len(history)>fitter.MAX_PROPOSALS:
            raise ArithmeticError("frozen coordinate fit exceeded12 proposals")
        if not all(value>0 and value.denominator<=fitter.DENOMINATOR_LIMIT for value in marker):
            raise ArithmeticError("invalid exported rational witness")
        bank.append((center,lo,hi,*marker)); evaluations+=len(history)
        report={"center":center,"status":status,"evaluations":len(history)}
        if rounded is not None:
            report["rounded_region_means_FLOAT"]=[float(x/256) for x in rounded["means"]]
        reports.append(report)
        if progress: progress("actual_marker_fit",center=center,status=status,evaluations=len(history))
    if tuple(row[:3] for row in bank)!=tuple(row[:3] for row in template):
        raise ArithmeticError("fitting changed reference coverage")
    return tuple(bank),{"scope":"FLOAT proposals only; exact positive rational markers exported",
        "certificate":False,"actual_anchor":actual,"coordinate_weight":actual,"degree_reversals":0,
        "local_cap_builds":0,"local_degree_rows":2145,"input_utility_sha256":fingerprint,
        "positive_coefficients_per_channel":positives,"zero_coefficients":zero_count,
        "successful_proposal_evaluations_recorded":evaluations,"fits":reports,"fit_failures":failures,
        "frozen_coordinate_evaluation_limit":fitter.MAX_PROPOSALS,
        "marker_denominator_limit":fitter.DENOMINATOR_LIMIT,"positive_floor":"1/1048576",
        "objective_numerator":"e0 H^16383 J","objective_targets":"256*(actual_full,actual_single)",
        "marker_independent_objective_constant":0,"endpoint_proxies_only_for_fitting":True,
        "proxy_monotonicity_not_assumed":True,"source_or_count_factor_applied":False,"pins":pins}


def exact_tests(deps,guard):
    """Exact small/interface tests only: no fitting, local build or Arb calls."""
    d=deps; checks=0; empty_total=0; source_rows=[]; digest=hashlib.sha256()
    template=d.reference.BANK; cells=d.normalized.normalized_cells(template)
    for r in TEST_WEIGHTS:
        guard(); validate(r,r,r,F(1,10),F(119,1280)); checks+=1
        for center,_,_,_,_ in template:
            full,single=fit_target(r,center)
            assert 2*full+single==r and min(full,single,4096-r+full)>0; checks+=1
            if r>=641:
                assert (full,single)==d.actual.fit_target(r,center); checks+=1
        ranges=d.counts.marginal_ranges(cells,r,guard)
        assigned=[f for interval in ranges if interval is not None for f in range(interval[0],interval[1]+1)]
        assert assigned==list(range(r//2+1)); checks+=1
        empty_total+=sum(interval is None for interval in ranges)
        for index,interval in enumerate(ranges):
            if interval is None: continue
            for f in range(interval[0],interval[1]+1):
                guard(); x=F(2048*f,r)
                expected=next(i for i,cell in enumerate(cells) if x<=cell[2])
                assert expected==index and cells[index][1]<=x<=cells[index][2]; checks+=1
        counts=d.join.run_counts(points=POINTS,weight=r)
        assert len(counts)==r+1 and counts[0]==0 and sum(counts)==comb(POINTS,r); checks+=1
        for c in sorted({1,r//2,r}):
            want=F(POINTS,c)*comb(r-1,c-1)*comb(POINTS-r-1,c-1)
            assert want.denominator==1 and counts[c]==want; checks+=1
        digest.update((str(r)+":"+str(sum(counts))+"\n").encode("ascii"))
        h,c,eta=d.source.source_factor(r)
        want_c=1-F(r,256*(r//256+1))
        want_eta=(1+F(sum(comb(r,j) for j in range(1,r//256+1)),255))/(255*want_c)
        assert (h,c,eta)==(r//256,want_c,want_eta) and 0<=h<r and c>0; checks+=1
        if r<256:
            assert eta==F(256,255*(256-r)); checks+=1
        source_rows.append({"actual":r,"h":h,"c":str(c),"eta":str(eta)})
        marker=(F(2,3),F(3,5))
        for f in sorted({0,1,r//4,r//2}):
            single=r-2*f; empty=4096-r+f
            assert min(f,single,empty)>=0 and f+single+empty==4096; checks+=1
            expected=F(comb(4096,single)*comb(4096-single,f))*marker[0]**f*marker[1]**single
            assert d.extraction.target_denominator(r,f,marker)==expected; checks+=1
        a,b=r//4,r//3; x,y=F(2048*a,r),F(2048*b,r)
        assert d.geometry.vertices(r,r,x,x,y,y,guard)==((F(r),F(a),F(b)),); checks+=1
        del counts,assigned,ranges
    assert empty_total>0 and d.source.source_factor(254)[2]==F(128,255)<1; checks+=1
    for center,_,_,_,_ in template:
        f,s=_interior_target(2,center)
        assert 2*f+s==2 and min(f,s,4094+f)>0; checks+=1
    assert _interior_target(2,0)==_interior_target(2,1024)==(F(1,2),F(1)); checks+=1
    assert fit_target(225,0)[0]>fit_target(225,3)[0]; checks+=1
    for bad in (224,4097,True,F(225)):
        try: fit_target(bad,512)
        except ValueError: checks+=1
        else: raise AssertionError("unsupported fitting anchor accepted")
    for args in ((224,225,225,F(1,2),F(1,2)),(225,226,225,F(1,2),F(1,2)),
            (4097,225,225,F(1,2),F(1,2)),(225,225,4097,F(1,2),F(1,2)),
            (True,225,225,F(1,2),F(1,2)),(225,225,225,0.5,F(1,2)),
            (225,225,225,F(1),F(1,2)),(225,225,225,F(1,2),F(0))):
        try: validate(*args)
        except (ValueError,TypeError): checks+=1
        else: raise AssertionError("invalid lower-domain parameters accepted")
    validate(4096,225,226,F(1,2),F(1)); checks+=1
    for r in TEST_WEIGHTS:
        assert 256*F(119,1280)*F(r,256)==F(119,1280)*r; checks+=1
    return {"lower_domain_exact_checks":checks,"weights":list(TEST_WEIGHTS),
        "empty_assigned_cells_retained":empty_total,"odd_sizes_and_source_jumps_checked":True,
        "unclipped_h0_eta_checked":True,"r2_proxy_algebra_only":True,
        "crossing_endpoint_proxy_targets_allowed":True,"actual_sizes_not_anchor_used_in_source":True,
        "indexed_support_totals_sha256":digest.hexdigest(),"source_controls":source_rows,
        "run_count_small_exact_controls":d.counts.exact_tests(d.join,guard),
        "ordinary_fraction_controls":d.ordinary.exact_tests(guard),
        "no_local_table_or_fit_or_Arb_evaluation":True}


def evaluate_band(deps,v,core,cells,numerators,marker_logs,first,last,z,alpha,guard,progress=None):
    """Return native log bound plus receipt; numerator inputs are regional."""
    d=deps; validate(first,first,last,z,alpha)
    den=d.actual.actual_denominators(v,core.BETA,alpha,d.band,guard)
    cover=d.actual.actual_coverage(d.normalized,cells,first,last,guard)
    bounds,geometry=d.run.cell_bounds(v,d.band,d.geometry,den,numerators,marker_logs,cells,
        first,last,"actual",d.join,guard,progress)
    counts=d.counts.source_counts(v,d.source,d.join,cells,alpha,first,last,"actual",guard)
    old_weight,_=d.band.outer_source_sum(v,d.source,alpha,first,last,guard)
    if not counts["full_weight"].overlaps(old_weight):
        raise ArithmeticError("actual eta/count totals disagree")
    joined=d.counts.combine(v,counts,bounds,guard)
    maximum=None
    for row in bounds.values():
        for _,value in row: maximum=value if maximum is None else maximum.max(value)
    if not joined["old_log_bound"].overlaps(old_weight.log()+maximum):
        raise ArithmeticError("complete old max envelope changed")
    logtwo=v.arb(2).log(); output=-v.exact(alpha)*core.LIMIT*v.exact(z).log()
    final=output+joined["minimum_log_bound"]
    count_report=dict(counts["report"])
    count_report["indexed_support_count_decimal_digits"]=len(count_report.pop("indexed_support_count"))
    geometry["largest_conditional_cells"]=geometry["largest_conditional_cells"][:4]
    report={"actual_first":first,"actual_last":last,"sizes":last-first+1,
        "all_integer_sizes_and_indexed_supports_covered":True,
        "log2_old_global_max_interval":str((output+joined["old_log_bound"])/logtwo),
        "log2_cell_run_sum_interval":str((output+joined["run_log_bound"])/logtwo),
        "log2_union_bound_interval":str(final/logtwo),
        "below_2_pow_minus_40":bool(final<-40*logtwo),
        "below_2_pow_minus_60":bool(final<-60*logtwo),
        "log2_eta_weighted_support_count_interval":str(old_weight.log()/logtwo),
        "log2_output_correction_interval":str(output/logtwo),
        "coverage":{k:x for k,x in cover.items() if k!="marginals_by_center"},
        "geometry":geometry,"counts":count_report,"join":joined["report"],
        "frozen_source_sum_and_complete_max_controls":True}
    del den,bounds,counts,joined
    gc.collect(); guard()
    return final,report


def old3072_control(deps,built,bank_hash,new_numerators,new_log_bound,guard,progress=None):
    """Independent frozen fitter/ordinary bank replay at the old witness."""
    d=deps; v,core=built["verifier"],built["core"]
    proposed,_=d.actual.fit_actual_markers(MappingProxyType(built["rational_coefficients"]),
        3072,d.fitter,guard,progress)
    bank,_,reference_hash=d.support.validate_bank(proposed,d.reference.BANK)
    if reference_hash!=bank_hash or reference_hash!=OLD_BANK:
        raise ArithmeticError("frozen3072 fitted markers differ")
    cells=d.normalized.normalized_cells(bank)
    old_num,old_logs=d.normalized.numerator_tables(v,built["coefficients"],cells,d.ordinary,guard,progress)
    checks=0
    for old,new in zip(old_num,new_numerators):
        for row1,row2 in zip(old,new):
            for x,y in zip(row1,row2):
                guard()
                if not x.overlaps(y): raise ArithmeticError("frozen3072 ordinary numerator mismatch")
                checks+=1
    if checks!=186:
        raise ArithmeticError("incomplete frozen3072 six-channel comparison")
    old_value,receipt=evaluate_band(d,v,core,cells,old_num,old_logs,3072,3072,
        F(1,10),F(119,1280),guard,progress)
    if not old_value.overlaps(new_log_bound):
        raise ArithmeticError("frozen3072 ordinary singleton bound mismatch")
    # The frozen assembler returns a report, not its native Arb result.
    # Parse its interval ONLY for this regression assertion. It never enters
    # a source factor, proof numerator, bound selection, or returned bound.
    frozen=d.generator.evaluate_band(v,core,d.run,d.actual,d.normalized,d.band,
        d.geometry,d.counts,d.source,d.join,cells,old_num,old_logs,
        3072,3072,F(1,10),F(119,1280),guard,progress)
    frozen_log2=v.arb(frozen["log2_union_bound_interval"])
    if not frozen_log2.overlaps(new_log_bound/v.arb(2).log()):
        raise ArithmeticError("frozen3072 assembler report disagrees with native bound")
    return {"performed":True,"regional_entry_overlap_checks":checks,
        "frozen_fitter_replayed":True,"native_directed_log_bounds_overlap":True,
        "frozen_generator_assembly_report_overlap":True,
        "parsed_reference_interval_used_only_for_regression_assertion":True,
        "ordinary_log2_union_bound_interval":receipt["log2_union_bound_interval"],
        "not_a_comparison_to_the_stronger_residue_certificate":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor",type=int,default=512)
    parser.add_argument("--first",type=int)
    parser.add_argument("--last",type=int)
    parser.add_argument("--z",type=rational64,default=F(17447,25000))
    parser.add_argument("--alpha",type=rational64,default=F(119,1280))
    parser.add_argument("--seconds",type=int,choices=(300,600),default=300)
    parser.add_argument("--ordinary-only",action="store_true",help="accepted explicitly; ordinary is the only current mode")
    parser.add_argument("--checks-only",action="store_true")
    args=parser.parse_args()
    if args.first is None: args.first=args.anchor
    if args.last is None: args.last=args.anchor
    validate(args.anchor,args.first,args.last,args.z,args.alpha)
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    sys.dont_write_bytecode=True
    started=time.monotonic(); deadline=started+args.seconds
    d=load_dependencies(); peak=d.source.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("lower actual-band deadline; no incomplete claim")
        if peak()>MEMORY: raise MemoryError("lower actual-band512MiB budget exceeded")
    def progress(phase,**values):
        if phase=="actual_marker_fit" and values.get("center") not in (0,512,1024): return
        if phase=="run_band_geometry" and values.get("row",0)%8: return
        kept={k:values[k] for k in ("center","status","evaluations","row","vertices") if k in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**kept}),flush=True)
    tests=exact_tests(d,guard)
    if args.checks_only:
        guard()
        print(json.dumps({"scope":"exact lower-domain interface checks only","checks_only":True,
            "numerical_bound_computed":False,"whole_code_certificate":False,"local_table_builds":0,
            "marker_fits":0,"tests":tests,"pins":d.pins,"generated_bytes":0,
            "elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),
            "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},separators=(",",":")),flush=True)
        return
    local=d.large.load("cheap_lane_ordered_local_directed.py",d.large.LOCAL_PIN,"lower_coherent_local")
    built=local.build_local_pairs(args.z,args.alpha,d.source,guard,deadline,PRECISION,progress=progress)
    v,core=built["verifier"],built["core"]; degrees=built["rational_coefficients"]
    before=d.support.coefficient_fingerprint(degrees,guard); identities=dict(built["coefficients"])
    if before!=built["report"]["rational_coefficient_sha256"]:
        raise ArithmeticError("local coefficient fingerprint mismatch")
    proposed,fit=fit_actual_markers(MappingProxyType(degrees),args.anchor,d.fitter,guard,progress)
    fitted,public_bank,bank_hash=d.support.validate_bank(proposed,d.reference.BANK)
    if d.support.coefficient_fingerprint(degrees,guard)!=before or len(identities)!=len(built["coefficients"]):
        raise ArithmeticError("fitting changed exact actual caps")
    if any(row is not built["coefficients"][key] for key,row in identities.items()):
        raise ArithmeticError("fitting changed original directed coefficient rows")
    old=d.generator.old_reference(args.anchor,args.z,args.alpha)
    if old and (before,bank_hash)!=(OLD_TABLE,OLD_BANK):
        raise ArithmeticError("old3072 table/bank hashes changed")
    cells=d.normalized.normalized_cells(fitted)
    den=d.actual.actual_denominators(v,core.BETA,args.alpha,d.band,guard)
    tests["frozen_directed_actual"]=d.actual.directed_tests(v,den,d.extraction,fitted,
        args.alpha,core.BETA,d.join,d.band,guard)
    own_checks=0
    for r in TEST_WEIGHTS:
        for f in (0,r//2):
            guard(); marker=(fitted[15][3],fitted[15][4])
            exact=d.extraction.target_denominator(r,f,marker)
            got=den.value(r,f,v.exact(marker[0]).log(),v.exact(marker[1]).log())
            want=v.exact(exact).log()-v.exact(args.alpha)*v.exact(F(r,256))*v.exact(core.BETA).log()
            if not got.overlaps(want): raise ArithmeticError("lower own-r directed normalizer mismatch")
            own_checks+=1
    tests["lower_directed_own_normalizers"]=own_checks
    del den
    numerators,marker_logs=d.normalized.numerator_tables(v,built["coefficients"],cells,d.ordinary,guard,progress)
    anchor_value,singleton=evaluate_band(d,v,core,cells,numerators,marker_logs,
        args.anchor,args.anchor,args.z,args.alpha,guard,progress)
    if (args.first,args.last)==(args.anchor,args.anchor): requested=singleton
    else:
        _,requested=evaluate_band(d,v,core,cells,numerators,marker_logs,
            args.first,args.last,args.z,args.alpha,guard,progress)
    regression=old3072_control(d,built,bank_hash,numerators,anchor_value,guard,progress) if old else {
        "performed":False,"reason":"requires anchor3072,z1/10,alpha119/1280"}
    guard()
    print(json.dumps({"scope":"complete ordinary actual-size band, all indexed supports and histories at35%",
        "directed_arithmetic":True,"whole_code_certificate":False,"mode":"ordinary_only",
        "parameters":{"anchor_actual":args.anchor,"actual_first":args.first,"actual_last":args.last,
            "z":str(args.z),"alpha":str(args.alpha),"working_precision":built["report"]["working_precision"],
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "anchor_singleton_control":singleton,"requested_band":requested,
        "separate_reported_bounds_not_added":True,"ordinary_only_has_no_residue_claim":True,
        "residue_tables_built":0,"positive_rational_bank":public_bank,
        "positive_rational_bank_sha256":bank_hash,"local_coefficient_sha256":before,
        "old3072_ordinary_control":regression,"local_table_builds":1,
        "local_caps_unchanged_by_fitting":True,"degree_reversals":0,
        "own_target_denominators_in_geometry":True,"source_beta_only_in_denominator_tilt":True,
        "eta_powered_once_counts_outside_alpha":True,"output_once_per_bound":True,
        "alpha_not_reapplied_after_local_caps":True,"both_terminal_orders_retained":True,
        "cell_counts_are_overlapping_exact_run_interval_upper_counts":True,
        "fit":{"scope":"FLOAT proposals only","evaluations":fit["successful_proposal_evaluations_recorded"],
            "fit_failures":fit["fit_failures"]},"tests":tests,
        "local_checks":{"tests":built["report"]["tests"],"field_tests":built["report"]["field_tests"]},
        "pins":d.pins,"generated_bytes":0,"elapsed_seconds":time.monotonic()-started,
        "peak_process_bytes":peak(),"nonnegative_intersections":v.intersections,
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},separators=(",",":")),flush=True)


if __name__=="__main__":
    main()
