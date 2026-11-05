"""Recalibrated actual-size bands with one coherent regional witness bank.

Build actual-occupancy local pair caps at the requested exact z,alpha.
The frozen bounded FLOAT fitter proposes31 positive rational markers at
the requested anchor. It does not change the cap table or provide any
proof numerator. Reuse this bank for the anchor singleton and one complete
inclusive actual-size band, including odd sizes and source-factor jumps.

Default mode builds cyclic256 H^64/H^63 J tables and exports their six
exact entrywise residue maxima per marker. Each marked target coefficient
is below its positive alias and its maximum. Mode --ordinary-only instead
uses directed H(u,v)^64/H(u,v)^63 J(u,v), with no cyclic tables or residue
claim. It is a separate, potentially looser upper bound. Neither mode
applies another fractional power after the local table.

The finite256-region chain retains both terminal orders. Every actual
target uses its own multinomial and marker denominator through the frozen
normalized-cell geometry. beta^(alpha*r) is restored by256 denominator
tilts; eta_r^alpha weights exact run-interval count upper bounds. Counts
stay outside alpha and the output factor z^(-alpha*L) occurs once. The
two reported bands are separate bounds, not added as a disjoint union.

No provider imports or work on import. Parent execution only, stdout and
no generated files. A canonical dyadic export is identity/provenance data;
its upper-bound status comes from the completed verified generator run.
Frozen helpers remain unchanged. Source alone is not a numerical receipt.
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
BASE_PIN="b2761383e6a929978764baa1af9f777e22907c662d505f352015e16f3d8a8938"
MAXIMUM_PIN="d03a5dff0c607cac71f7b6da45c355c07654ca91bea879d726e58ce05c098036"
RUN_PIN="19914e8bc698c410971f8edb0ab7f884cb0b83955bee2a6933e9791aa9dfcfc2"
COUNTS_PIN="2468da0c91e6e44c7b707a1158252bb695818980de33eb7f270d620748b4f028"
PRECISION,MEMORY=384,512*(1<<20)


def load_base():
    path=Path(__file__).with_name("cheap_lane_residue_max_band.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=BASE_PIN:
        raise RuntimeError("frozen residue-max band source mismatch")
    spec=importlib.util.spec_from_file_location("recalibrated_maximum_base",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    if (module.MAXIMUM_PIN,module.RUN_PIN,module.COUNTS_PIN)!=(MAXIMUM_PIN,RUN_PIN,COUNTS_PIN):
        raise ArithmeticError("frozen arithmetic/count provenance changed")
    return module


def rational64(text):
    if not isinstance(text,str) or len(text)>128 or "e" in text.lower():
        raise argparse.ArgumentTypeError("use a bounded exact fraction or decimal")
    try: value=F(text)
    except (ValueError,ZeroDivisionError) as exc:
        raise argparse.ArgumentTypeError("expected an exact rational") from exc
    if max(abs(value.numerator).bit_length(),value.denominator.bit_length())>64:
        raise argparse.ArgumentTypeError("rational witnesses limited to64 bits")
    return value


def validate(anchor,first,last,z,alpha):
    if any(type(x) is not int for x in (anchor,first,last)) or not 641<=anchor<=4096 or not 641<=first<=last<=4096:
        raise ValueError("actual anchor and inclusive band must lie in641..4096")
    if any(isinstance(x,bool) or not isinstance(x,(int,F)) for x in (z,alpha)):
        raise TypeError("exact rational z and alpha required")
    if not 0<z<1 or not 0<alpha<=1:
        raise ValueError("require0<z<1 and0<alpha<=1")


def old_reference(anchor,z,alpha):
    return (anchor,z,alpha)==(3072,F(1,10),F(119,1280))


def exact_tests(guard):
    checks=0
    for anchor,first,last in ((2944,2944,2944),(3264,3200,3328),(3072,3008,3136),(641,4095,4096)):
        guard(); validate(anchor,first,last,F(1,10),F(119,1280)); checks+=1
    for anchor,first,last,z,alpha in ((640,641,641,F(1,10),F(1,2)),
            (4097,4096,4096,F(1,10),F(1,2)),(2944,2945,2944,F(1,10),F(1,2)),
            (2944,2944,4097,F(1,10),F(1,2)),(True,2944,2944,F(1,10),F(1,2)),
            (2944,2944,2944,F(1),F(1,2)),(2944,2944,2944,F(1,10),F(0))):
        try: validate(anchor,first,last,z,alpha)
        except ValueError: checks+=1
        else: raise AssertionError("invalid parameter domain accepted")
    for text,want in (("1/10",F(1,10)),("0.125",F(1,8)),("119/1280",F(119,1280))):
        assert rational64(text)==want; checks+=1
    for text in ("1/0","1e10",str(1<<65),"x","1"*129):
        try: rational64(text)
        except argparse.ArgumentTypeError: checks+=1
        else: raise AssertionError("unbounded rational accepted")
    assert old_reference(3072,F(1,10),F(119,1280)); checks+=1
    for anchor,z,alpha in ((2944,F(1,10),F(119,1280)),(3072,F(1,9),F(119,1280)),
                         (3072,F(1,10),F(3,32))):
        assert not old_reference(anchor,z,alpha); checks+=1
    return {"parameter_and_reference_checks":checks,"anchor_can_lie_outside_requested_band":True}


def evaluate_band(v,core,run,actual,normalized,band,geometry,counts_module,source,join,
                  cells,numerators,marker_logs,first,last,z,alpha,guard,progress):
    """Explicit dynamic parameters; no modification of frozen module globals."""
    den=actual.actual_denominators(v,core.BETA,alpha,band,guard)
    cover=actual.actual_coverage(normalized,cells,first,last,guard)
    bounds,geometry_report=run.cell_bounds(v,band,geometry,den,numerators,marker_logs,cells,
        first,last,"actual",join,guard,progress)
    counts=counts_module.source_counts(v,source,join,cells,alpha,first,last,"actual",guard)
    old_weight,_=band.outer_source_sum(v,source,alpha,first,last,guard)
    if not counts["full_weight"].overlaps(old_weight):
        raise ArithmeticError("actual eta/count sum disagrees with frozen source interface")
    joined=counts_module.combine(v,counts,bounds,guard)
    old_max=None
    for row in bounds.values():
        for _,value in row: old_max=value if old_max is None else old_max.max(value)
    if not joined["old_log_bound"].overlaps(old_weight.log()+old_max):
        raise ArithmeticError("complete global-max reference changed")
    output=-v.exact(alpha)*core.LIMIT*v.exact(z).log(); logtwo=v.arb(2).log()
    final=output+joined["minimum_log_bound"]
    count_report=dict(counts["report"])
    count_report["indexed_support_count_decimal_digits"]=len(count_report.pop("indexed_support_count"))
    join_report=dict(joined["report"])
    # Printed logs are re-enclosed ONLY for display. No proof input uses them.
    join_report["largest_run_weighted_cells"]=[
        {**row,"log2_weighted_cell_with_output_interval":str(
            v.arb(row["log2_weighted_cell_interval"])+output/logtwo)}
        for row in join_report["largest_run_weighted_cells"][:4]]
    geometry_report["largest_conditional_cells"]=geometry_report["largest_conditional_cells"][:4]
    result={"actual_first":first,"actual_last":last,"sizes":last-first+1,
        "all_integer_sizes_and_indexed_supports_covered":True,
        "log2_old_global_max_interval":str((output+joined["old_log_bound"])/logtwo),
        "log2_cell_run_sum_interval":str((output+joined["run_log_bound"])/logtwo),
        "log2_union_bound_interval":str(final/logtwo),
        "below_2_pow_minus_40":bool(final<-40*logtwo),
        "below_2_pow_minus_60":bool(final<-60*logtwo),
        "log2_eta_weighted_support_count_interval":str(old_weight.log()/logtwo),
        "log2_output_correction_interval":str(output/logtwo),
        "coverage":{key:value for key,value in cover.items() if key!="marginals_by_center"},
        "geometry":geometry_report,"counts":count_report,"join":join_report,
        "source_sum_and_global_max_controls_passed":True}
    del den,bounds,counts,joined
    gc.collect(); guard()
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor",type=int,default=2944)
    parser.add_argument("--first",type=int)
    parser.add_argument("--last",type=int)
    parser.add_argument("--z",type=rational64,default=F(1,10))
    parser.add_argument("--alpha",type=rational64,default=F(119,1280))
    parser.add_argument("--bits",type=int,choices=(192,256),default=192)
    parser.add_argument("--seconds",type=int,choices=(300,600,1500,1800))
    parser.add_argument("--ordinary-only",action="store_true")
    args=parser.parse_args()
    if args.first is None: args.first=args.anchor
    if args.last is None: args.last=args.anchor
    if args.seconds is None: args.seconds=300 if args.ordinary_only else 1500
    if not args.ordinary_only and args.seconds not in (1500,1800):
        parser.error("residue-max generation requires --seconds1500 or1800")
    validate(args.anchor,args.first,args.last,args.z,args.alpha)
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    base=load_base(); maximum=base.load_maximum()
    complete=maximum.load_complete(); selected=complete.load_selected(); run=selected.load_run()
    if (selected.RUN_PIN,run.COUNTS_PIN)!=(RUN_PIN,COUNTS_PIN):
        raise ArithmeticError("run/count provenance changed")
    actual=run.load_actual(); support=actual.load_support()
    normalized=support.load_normalized(); band=normalized.load_band()
    counts_module=band.load("cheap_lane_cell_run_counts.py",COUNTS_PIN,"recalibrated_max_counts")
    fitter=band.load(support.FITTER_FILE,support.FITTER_PIN,"recalibrated_max_fitter")
    geometry=band.load(normalized.GEOMETRY_FILE,normalized.GEOMETRY_PIN,"recalibrated_max_geometry")
    reference=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"recalibrated_max_reference")
    ordinary=reference.load("cheap_lane_ordered_ordinary_directed.py",reference.ORDINARY_PIN,"recalibrated_max_ordinary")
    large=ordinary.load_large()
    if (large.REGION_PIN,large.INTEGER_PIN)!=(maximum.REGIONS_PIN,maximum.INTEGER_PIN):
        raise ArithmeticError("regional arithmetic provenance changed")
    source=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"recalibrated_max_source")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"recalibrated_max_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"recalibrated_max_extraction")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"recalibrated_max_terminal")
    peak=source.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("recalibrated band deadline; no incomplete claim")
        if peak()>MEMORY: raise MemoryError("recalibrated band512MiB budget exceeded")
    def progress(phase,**values):
        if phase=="actual_marker_fit" and values.get("center") not in (0,512,1024): return
        if phase=="run_band_geometry" and values.get("row",0)%8: return
        kept={key:values[key] for key in ("center","status","evaluations","row","vertices",
            "tables_completed","actual_first","actual_last") if key in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**kept}),flush=True)
    tests={"parameters":exact_tests(guard),"base_export":base.own_tests(guard),
        "bank":support.wrapper_tests(reference.BANK,guard),"ordinary":ordinary.exact_tests(guard),
        "run_wrapper":run.wrapper_tests(guard),"run_counts":counts_module.exact_tests(join,guard),
        "actual":actual.exact_tests(reference.BANK,fitter,normalized,geometry,guard)}
    integer=None
    if not args.ordinary_only:
        integer=large.load("cheap_lane_integer_regions.py",maximum.INTEGER_PIN,"recalibrated_max_integer")
        tests["maxima"]=maximum.exact_tests(integer,guard)
        tests["integer"]=integer.exact_tests(guard)
        tests["extraction"]=extraction.exact_tests(guard)
    built=local.build_local_pairs(args.z,args.alpha,source,guard,deadline,PRECISION,progress=progress)
    v,core=built["verifier"],built["core"]; degrees=built["rational_coefficients"]
    original_arb=dict(built["coefficients"])
    table_hash=support.coefficient_fingerprint(degrees,guard)
    if table_hash!=built["report"]["rational_coefficient_sha256"]:
        raise ArithmeticError("local table fingerprint mismatch")
    proposed,fit=actual.fit_actual_markers(MappingProxyType(degrees),args.anchor,fitter,guard,progress)
    fitted,public_bank,bank_hash=support.validate_bank(proposed,reference.BANK)
    if support.coefficient_fingerprint(degrees,guard)!=table_hash or len(original_arb)!=len(built["coefficients"]):
        raise ArithmeticError("fitter changed local caps")
    if any(row is not built["coefficients"][key] for key,row in original_arb.items()):
        raise ArithmeticError("fitter changed original directed rows")
    regression=old_reference(args.anchor,args.z,args.alpha)
    if regression and (table_hash,bank_hash)!=(base.TABLE_SHA,base.BANK_SHA):
        raise ArithmeticError("old3072 table/bank regression mismatch")
    cells=normalized.normalized_cells(fitted)
    den=actual.actual_denominators(v,core.BETA,args.alpha,band,guard)
    tests["directed_actual"]=actual.directed_tests(v,den,extraction,fitted,args.alpha,core.BETA,join,band,guard)
    del den
    raw_witness=witness_sha=None; diagnostics=[]
    totals={key:0 for key in ("DC_upper_checks","legal_coefficients","zero_padding_checked",
        "matrix_products","polynomial_products","total_quantizations","rounded_coefficients")}
    totals["maximum_shift"]=0
    if args.ordinary_only:
        numerators,marker_logs=normalized.numerator_tables(v,built["coefficients"],cells,ordinary,guard,progress)
    else:
        numerators=[]; marker_logs=[]; exported=[]
        for index,(center,_,_,u,w) in enumerate(cells):
            guard(); progress("recalibrated_maximum_start",center=center,tables_completed=index)
            maxima,report=maximum.build_maximum(v,degrees,(u,w),integer,args.bits,guard,selected=None)
            dyadics=base.dyadic_export(maxima,report)
            if any(value<=0 for row in maxima for value in row):
                raise ArithmeticError("nonpositive maximum cannot enter logarithmic cell bound")
            numerators.append(tuple(tuple(v.exact(value) for value in row) for row in maxima))
            marker_logs.append((v.exact(u).log(),v.exact(w).log()))
            exported.append([center,str(u),str(w),dyadics[:4],dyadics[4:]])
            diagnostics.append({"center":center,
                "six_dyadics_sha256":hashlib.sha256(json.dumps(dyadics,separators=(",",":")).encode("ascii")).hexdigest(),
                "log2_ordinary_over_maximum_range":base.gain_range(v,report["entries"])})
            totals["DC_upper_checks"]+=report["DC_upper_checks"]
            for entry in report["entries"]:
                for key in ("legal_coefficients","zero_padding_checked"): totals[key]+=entry[key]
            for key in ("matrix_products","polynomial_products","total_quantizations","rounded_coefficients"):
                totals[key]+=report["arithmetic"][key]
            totals["maximum_shift"]=max(totals["maximum_shift"],report["arithmetic"]["maximum_shift"])
            del maxima,report; gc.collect(); guard()
            progress("recalibrated_maximum_complete",center=center,tables_completed=index+1)
        if len(numerators)!=31 or totals["DC_upper_checks"]!=186:
            raise ArithmeticError("incomplete residue maximum bank")
        raw_witness={"coordinate":"actual","moduli":[256,256],"integer_bits":args.bits,
            "z":str(args.z),"alpha":str(args.alpha),"local_coefficient_sha256":table_hash,
            "positive_rational_bank_sha256":bank_hash,"rows":exported}
        witness_sha=hashlib.sha256(json.dumps(raw_witness,sort_keys=True,separators=(",",":")).encode("ascii")).hexdigest()
        numerators,marker_logs=tuple(numerators),tuple(marker_logs)
    if len(numerators)!=31 or tuple(map(lambda x:tuple(map(len,x)),numerators))!=((4,2),)*31:
        raise ArithmeticError("incomplete six-channel numerator bank")
    progress("recalibrated_anchor_control",actual_first=args.anchor,actual_last=args.anchor)
    singleton=evaluate_band(v,core,run,actual,normalized,band,geometry,counts_module,source,join,
        cells,numerators,marker_logs,args.anchor,args.anchor,args.z,args.alpha,guard,progress)
    if (args.first,args.last)==(args.anchor,args.anchor): requested=singleton
    else:
        progress("recalibrated_requested_band",actual_first=args.first,actual_last=args.last)
        requested=evaluate_band(v,core,run,actual,normalized,band,geometry,counts_module,source,join,
            cells,numerators,marker_logs,args.first,args.last,args.z,args.alpha,guard,progress)
    guard()
    mode="ordinary_only" if args.ordinary_only else "universal_residue_maximum"
    print(json.dumps({"scope":"complete actual-size band using "+mode+" regional numerators",
        "directed_arithmetic":True,"whole_code_certificate":False,"mode":mode,
        "parameters":{"anchor_actual":args.anchor,"actual_first":args.first,"actual_last":args.last,
            "z":str(args.z),"alpha":str(args.alpha),"integer_bits":None if args.ordinary_only else args.bits,
            "working_precision":built["report"]["working_precision"],"seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "anchor_singleton_control":singleton,"requested_band":requested,
        "separate_reported_bounds_not_added":True,"ordinary_only_has_no_residue_claim":args.ordinary_only,
        "residue_tables_built":0 if args.ordinary_only else 31,"same_numerators_reused_for_both_bounds":True,
        "exact_maximum_witness":raw_witness,"maximum_witness_sha256":witness_sha,
        "witness_hash_format":"canonical sorted-key compact JSON of exact_maximum_witness",
        "witness_row_format":["center","u","v","R4_dyadics","V2_dyadics"],
        "dyadic_convention":"value=mantissa*2^exponent; channels R00,R01,R10,R11,V0,V1",
        "witness_diagnostics":diagnostics,"positive_rational_bank":public_bank,
        "local_coefficient_sha256":table_hash,"positive_rational_bank_sha256":bank_hash,
        "old3072_regression_hashes_enforced":regression,"local_table_builds":1,
        "actual_caps_rebuilt_at_requested_z_alpha":True,"local_caps_unchanged_by_fitting":True,
        "degree_reversals":0,"own_target_denominators_in_geometry":True,
        "source_beta_only_in_denominator_tilt":True,"eta_powered_once_counts_outside_alpha":True,
        "output_once_per_bound":True,"alpha_not_reapplied_after_local_caps":True,
        "finite256_region_chain_and_both_terminal_orders_retained":True,
        "completed_global_and_cell_run_bounds_minimized":True,"ordinary_maximum_numerator_minimum_used":False,
        "cell_counts_are_overlapping_exact_run_interval_upper_counts":True,
        "fit":{"scope":"FLOAT proposals only; exact positive rational witnesses replayed",
            "successful_proposal_evaluations_recorded":fit["successful_proposal_evaluations_recorded"],
            "fit_failures":fit["fit_failures"]},
        "aggregate_regional_arithmetic_checks":totals,"tests":tests,
        "local_checks":{"tests":built["report"]["tests"],"field_tests":built["report"]["field_tests"]},
        "generated_bytes":0,"elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),
        "nonnegative_intersections":v.intersections,
        "pins":{"base":BASE_PIN,"maximum":MAXIMUM_PIN,"run":RUN_PIN,"counts":COUNTS_PIN,
            "actual":run.ACTUAL_PIN,"support":actual.SUPPORT_PIN,"fitter":support.FITTER_PIN,
            "normalized":support.NORMALIZED_PIN,"band":normalized.BAND_PIN,"geometry":normalized.GEOMETRY_PIN,
            "reference":band.HOLE_PIN,"ordinary":reference.ORDINARY_PIN,"local":large.LOCAL_PIN,
            "source":large.ANCHOR_PIN,"join":large.JOIN_PIN,"extraction":large.REGION_PIN,"integer":large.INTEGER_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},separators=(",",":")),flush=True)


if __name__=="__main__":
    main()
