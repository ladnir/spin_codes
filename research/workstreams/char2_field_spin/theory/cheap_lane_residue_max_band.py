"""Directed actual-size bands from universal regional residue maxima.

Rebuild the frozen actual3072 table and31 positive rational markers at
z1/10, alpha119/1280. For each marker, build H^64 and H^63 J in the
positive cyclic256 ring with upward dyadic rounding. Replace each regional
entry by its maximum over all65536 residues. Every marked target
coefficient is below its alias and therefore below this entrywise maximum.
The maximizing residues may differ; only positive entrywise domination
is used. No additional regional power or fractional power is applied.

The frozen cell evaluator still composes the finite256-region chain.
Each actual target keeps D(r,f)=C(4096,f)C(4096-f,r-2f)u^f v^(r-2f).
Concave log-denominator minima use exact rational normalized-cell vertices.
The beta tilt contributes beta^(alpha*r) through256 region factors.
Exact run-interval masses, weighted by eta_r^alpha, upper-bound cell
profile counts; overlapping cell intervals are charged separately. Both
terminal orders are charged off diagonal. The output factor occurs once.

The SAME31 maxima supply the singleton3072 control and the requested
inclusive band. This control need not reproduce the stronger own-residue
certificate. Export six exact dyadics per marker for later read-only reuse;
no proof value is reconstructed from printed logs. Parent execution only:
1500/1800s,512MiB,20GiB free,-B,stdout only. No work/provider imports on
import. Source alone is not evidence of a numerical result.
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
MAXIMUM_PIN="d03a5dff0c607cac71f7b6da45c355c07654ca91bea879d726e58ce05c098036"
RUN_PIN="19914e8bc698c410971f8edb0ab7f884cb0b83955bee2a6933e9791aa9dfcfc2"
COUNTS_PIN="2468da0c91e6e44c7b707a1158252bb695818980de33eb7f270d620748b4f028"
TABLE_SHA="636b295acf8ba93a20f76a983ab4723c290f9967a252a04b9562148a3a7c539d"
BANK_SHA="3e4843f4abd6be8730f72a5f77c02f4b4e554bf32a834464bda3df8643f551e7"
Z,ALPHA=F(1,10),F(119,1280)
ANCHOR,PRECISION,MEMORY=3072,384,512*(1<<20)
CHANNELS=("R00","R01","R10","R11","V0","V1")


def load_maximum():
    path=Path(__file__).with_name("cheap_lane_residue_maximum.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=MAXIMUM_PIN:
        raise RuntimeError("frozen residue maximum source mismatch")
    spec=importlib.util.spec_from_file_location("residue_max_band_engine",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    if (module.Z,module.ALPHA,module.WEIGHT,module.TABLE_SHA,module.BANK_SHA)!=(Z,ALPHA,ANCHOR,TABLE_SHA,BANK_SHA):
        raise ArithmeticError("frozen maximum witness changed")
    return module


def validate_band(first,last):
    if any(type(x) is not int for x in (first,last)) or not 641<=first<=last<=4096:
        raise ValueError("inclusive actual band must lie in641..4096")


def dyadic_export(maxima,report):
    """Validate all six unnormalized exports against their separate scales."""
    if tuple(map(len,maxima))!=(4,2) or len(report["entries"])!=6:
        raise ArithmeticError("maximum engine changed its six-channel layout")
    order=(("H64",0,0),("H64",0,1),("H64",1,0),("H64",1,1),("H63J",0,0),("H63J",1,0))
    result=[]
    for value,item,expected in zip((*maxima[0],*maxima[1]),report["entries"],order):
        if (item["object"],item["row"],item["column"])!=expected:
            raise ArithmeticError("maximum entry order changed")
        mantissa,exponent=int(item["mantissa"]),item["exponent"]
        if type(exponent) is not int or mantissa<0:
            raise ArithmeticError("invalid exact maximum dyadic")
        represented=F(mantissa<<exponent) if exponent>=0 else F(mantissa,1<<(-exponent))
        if represented!=value:
            raise ArithmeticError("compact dyadic does not reproduce maximum")
        result.append([str(mantissa),exponent])
    return result


def own_tests(guard):
    checks=0
    for first,last in ((641,641),(3008,3136),(3072,3072),(4095,4096)):
        guard(); validate_band(first,last); checks+=1
    for first,last in ((640,641),(4096,4097),(3073,3072),(True,3072)):
        try: validate_band(first,last)
        except ValueError: checks+=1
        else: raise AssertionError("invalid actual band accepted")
    entries=[]; values=[]
    for index,(obj,i,j) in enumerate((("H64",0,0),("H64",0,1),("H64",1,0),("H64",1,1),("H63J",0,0),("H63J",1,0))):
        exponent=-5 if obj=="H64" else 7
        entries.append({"object":obj,"row":i,"column":j,"mantissa":str(index),"exponent":exponent})
        values.append(F(index,32) if exponent<0 else F(index*128))
    maxima=(tuple(values[:4]),tuple(values[4:])); report={"entries":entries}
    assert dyadic_export(maxima,report)==[[str(i),-5 if i<4 else 7] for i in range(6)]; checks+=1
    bad={"entries":[dict(x) for x in entries]}; bad["entries"][4]["exponent"]=-5
    try: dyadic_export(maxima,bad)
    except ArithmeticError: checks+=1
    else: raise AssertionError("terminal scale was silently replaced")
    return {"wrapper_exact_checks":checks,"separate_terminal_export_scale":True}


def gain_range(v,entries):
    """Compact DISPLAY only; these parsed logs never enter proof arithmetic."""
    values=[v.arb(item["log2_ordinary_over_maximum_interval"]) for item in entries
            if item["log2_ordinary_over_maximum_interval"] is not None]
    if len(values)!=6: return {"finite_gain_entries":len(values)}
    low=high=values[0]
    for value in values[1:]: low=low.min(value); high=high.max(value)
    return {"minimum_interval":low.str(16),"maximum_interval":high.str(16)}


def evaluate_band(v,core,run,actual,normalized,band,geometry,counts_module,anchor,join,
                  cells,numerators,marker_logs,first,last,guard,progress):
    """Return one completed bound; maxima have no target/source normalization."""
    validate_band(first,last)
    den=actual.actual_denominators(v,core.BETA,ALPHA,band,guard)
    cover=actual.actual_coverage(normalized,cells,first,last,guard)
    bounds,geometry_report=run.cell_bounds(v,band,geometry,den,numerators,marker_logs,cells,
        first,last,"actual",join,guard,progress)
    counts=counts_module.source_counts(v,anchor,join,cells,ALPHA,first,last,"actual",guard)
    old_weight,_=band.outer_source_sum(v,anchor,ALPHA,first,last,guard)
    if not counts["full_weight"].overlaps(old_weight):
        raise ArithmeticError("actual eta/count sum changed")
    joined=counts_module.combine(v,counts,bounds,guard)
    old_max=None
    for row in bounds.values():
        for _,value in row: old_max=value if old_max is None else old_max.max(value)
    if not joined["old_log_bound"].overlaps(old_weight.log()+old_max):
        raise ArithmeticError("old global-max control changed")
    output=-v.exact(ALPHA)*core.LIMIT*v.exact(Z).log(); logtwo=v.arb(2).log()
    final=output+joined["minimum_log_bound"]
    count_report=dict(counts["report"])
    count_report["indexed_support_count_decimal_digits"]=len(count_report.pop("indexed_support_count"))
    # Exact checksums remain; neither the enormous count nor all cell arrays
    # are needed in the stdout receipt. Leader logs below include output.
    join_report=dict(joined["report"])
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
        "source_sum_and_old_max_controls_passed":True}
    del bounds,counts,den,joined
    gc.collect(); guard()
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first",type=int,default=3008)
    parser.add_argument("--last",type=int,default=3136)
    parser.add_argument("--bits",type=int,choices=(192,256),default=192)
    parser.add_argument("--seconds",type=int,choices=(1500,1800),default=1500)
    args=parser.parse_args(); validate_band(args.first,args.last)
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    maximum=load_maximum(); complete=maximum.load_complete(); selected=complete.load_selected()
    if selected.RUN_PIN!=RUN_PIN: raise ArithmeticError("run-band provenance changed")
    run=selected.load_run()
    if run.COUNTS_PIN!=COUNTS_PIN: raise ArithmeticError("run-count provenance changed")
    actual=run.load_actual(); support=actual.load_support()
    normalized=support.load_normalized(); band=normalized.load_band()
    counts_module=band.load("cheap_lane_cell_run_counts.py",COUNTS_PIN,"maximum_band_counts")
    fitter=band.load(support.FITTER_FILE,support.FITTER_PIN,"maximum_band_fitter")
    geometry=band.load(normalized.GEOMETRY_FILE,normalized.GEOMETRY_PIN,"maximum_band_geometry")
    reference=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"maximum_band_reference")
    ordinary=reference.load("cheap_lane_ordered_ordinary_directed.py",reference.ORDINARY_PIN,"maximum_band_ordinary")
    large=ordinary.load_large()
    if (large.REGION_PIN,large.INTEGER_PIN)!=(maximum.REGIONS_PIN,maximum.INTEGER_PIN):
        raise ArithmeticError("regional arithmetic provenance changed")
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"maximum_band_source")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"maximum_band_local")
    extraction=large.load("cheap_lane_band_regions.py",maximum.REGIONS_PIN,"maximum_band_extraction")
    integer=large.load("cheap_lane_integer_regions.py",maximum.INTEGER_PIN,"maximum_band_integer")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"maximum_band_terminal")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("residue-max band deadline; no incomplete claim")
        if peak()>MEMORY: raise MemoryError("residue-max band512MiB limit exceeded")
    def progress(phase,**values):
        if phase=="actual_marker_fit" and values.get("center") not in (0,512,1024): return
        if phase=="run_band_geometry" and values.get("row",0)%8: return
        selected_values={key:values[key] for key in
            ("center","status","evaluations","row","vertices","tables_completed","actual_first","actual_last")
            if key in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**selected_values}),flush=True)
    tests={"wrapper":own_tests(guard),"maxima":maximum.exact_tests(integer,guard),
        "integer":integer.exact_tests(guard),"extraction":extraction.exact_tests(guard),
        "ordinary":ordinary.exact_tests(guard),"run_wrapper":run.wrapper_tests(guard),
        "bank":support.wrapper_tests(reference.BANK,guard),
        "actual":actual.exact_tests(reference.BANK,fitter,normalized,geometry,guard),
        "run_counts":counts_module.exact_tests(join,guard)}
    built=local.build_local_pairs(Z,ALPHA,anchor,guard,deadline,PRECISION,progress=progress)
    v,core=built["verifier"],built["core"]; degrees=built["rational_coefficients"]
    original_arb=dict(built["coefficients"])
    table_hash=support.coefficient_fingerprint(degrees,guard)
    if table_hash!=TABLE_SHA or built["report"]["rational_coefficient_sha256"]!=TABLE_SHA:
        raise ArithmeticError("local table differs from frozen actual3072 witness")
    proposed,fit=actual.fit_actual_markers(MappingProxyType(degrees),ANCHOR,fitter,guard,progress)
    fitted,public_bank,bank_hash=support.validate_bank(proposed,reference.BANK)
    if bank_hash!=BANK_SHA: raise ArithmeticError("rational bank differs from frozen actual3072 witness")
    if support.coefficient_fingerprint(degrees,guard)!=TABLE_SHA or len(original_arb)!=len(built["coefficients"]):
        raise ArithmeticError("fitter changed the local table")
    if any(row is not built["coefficients"][key] for key,row in original_arb.items()):
        raise ArithmeticError("fitter changed original Arb rows")
    cells=normalized.normalized_cells(fitted)
    den=actual.actual_denominators(v,core.BETA,ALPHA,band,guard)
    tests["directed_actual"]=actual.directed_tests(v,den,extraction,fitted,ALPHA,core.BETA,join,band,guard)
    del den
    numerators=[]; marker_logs=[]; exported=[]; witness_diagnostics=[]
    totals={key:0 for key in ("DC_upper_checks","legal_coefficients","zero_padding_checked",
        "matrix_products","polynomial_products","total_quantizations","rounded_coefficients")}
    totals["maximum_shift"]=0
    for index,(center,_,_,u,w) in enumerate(cells):
        guard(); progress("maximum_band_table_start",center=center,tables_completed=index)
        maxima,report=maximum.build_maximum(v,degrees,(u,w),integer,args.bits,guard)
        dyadics=dyadic_export(maxima,report)
        if any(value<=0 for row in maxima for value in row):
            raise ArithmeticError("nonpositive regional maximum; log-cell replay unavailable")
        numerators.append(tuple(tuple(v.exact(value) for value in row) for row in maxima))
        marker_logs.append((v.exact(u).log(),v.exact(w).log()))
        exported.append([center,str(u),str(w),dyadics[:4],dyadics[4:]])
        witness_diagnostics.append({"center":center,
            "six_dyadics_sha256":hashlib.sha256(json.dumps(dyadics,separators=(",",":")).encode("ascii")).hexdigest(),
            "log2_ordinary_over_maximum_range":gain_range(v,report["entries"])})
        totals["DC_upper_checks"]+=report["DC_upper_checks"]
        for entry in report["entries"]:
            for key in ("legal_coefficients","zero_padding_checked"): totals[key]+=entry[key]
        for key in ("matrix_products","polynomial_products","total_quantizations","rounded_coefficients"):
            totals[key]+=report["arithmetic"][key]
        totals["maximum_shift"]=max(totals["maximum_shift"],report["arithmetic"]["maximum_shift"])
        del maxima,report; gc.collect(); guard()
        progress("maximum_band_table_complete",center=center,tables_completed=index+1)
    if len(numerators)!=31 or totals["DC_upper_checks"]!=186:
        raise ArithmeticError("incomplete universal numerator bank")
    raw_witness={"coordinate":"actual","moduli":[256,256],"integer_bits":args.bits,
        "z":str(Z),"alpha":str(ALPHA),"local_coefficient_sha256":table_hash,
        "positive_rational_bank_sha256":bank_hash,"rows":exported}
    witness_sha=hashlib.sha256(json.dumps(raw_witness,sort_keys=True,separators=(",",":")).encode("ascii")).hexdigest()
    numerators,marker_logs=tuple(numerators),tuple(marker_logs)
    progress("maximum_band_singleton_control",actual_first=ANCHOR,actual_last=ANCHOR)
    singleton=evaluate_band(v,core,run,actual,normalized,band,geometry,counts_module,anchor,join,
        cells,numerators,marker_logs,ANCHOR,ANCHOR,guard,progress)
    if (args.first,args.last)==(ANCHOR,ANCHOR): requested=singleton
    else:
        progress("maximum_band_requested_envelope",actual_first=args.first,actual_last=args.last)
        requested=evaluate_band(v,core,run,actual,normalized,band,geometry,counts_module,anchor,join,
            cells,numerators,marker_logs,args.first,args.last,guard,progress)
    guard()
    print(json.dumps({"scope":"complete actual-size band using universal residue-max regional numerators",
        "directed_arithmetic":True,"whole_code_certificate":False,
        "parameters":{"actual_first":args.first,"actual_last":args.last,"fit_anchor_actual":ANCHOR,
            "z":str(Z),"alpha":str(ALPHA),"integer_bits":args.bits,
            "working_precision":built["report"]["working_precision"],"modulus":256,
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "singleton3072_control":singleton,"requested_band":requested,
        "singleton_control_is_not_the_stronger_own_residue_certificate":True,
        "regional_tables_built":31,"same_maxima_reused_for_both_bounds":True,
        "exact_maximum_witness":raw_witness,"maximum_witness_sha256":witness_sha,
        "witness_hash_format":"SHA256 of canonical sorted-key compact JSON of exact_maximum_witness",
        "witness_row_format":["center","u","v","R4_dyadics","V2_dyadics"],
        "witness_diagnostics":witness_diagnostics,
        "dyadic_convention":"value = mantissa * 2^exponent; channels R00,R01,R10,R11,V0,V1",
        "local_coefficient_sha256":table_hash,"positive_rational_bank_sha256":bank_hash,
        "baseline_hashes_enforced":True,"local_table_builds":1,"degree_reversals":0,
        "regional_maxima_unnormalized":True,"own_target_denominators_in_geometry":True,
        "separate_regular_terminal_scales_preserved":True,"maxima_need_not_share_a_residue":True,
        "regional_powers_applied_once":True,"finite_256_region_chain_retained":True,
        "alpha_not_reapplied":True,"beta_only_in_denominator_tilt":True,
        "eta_powered_once_counts_outside_alpha":True,"output_factor_once_per_bound":True,
        "both_offdiagonal_terminal_orders_charged":True,
        "cell_count_scope":"exact run-interval masses upper-bound cell counts; overlaps are charged",
        "completed_old_and_run_upper_bounds_minimized":True,"ordinary_numerator_minimum_used":False,
        "fit":{"scope":"FLOAT-only proposals; proof uses exact rational witnesses",
            "successful_proposal_evaluations_recorded":fit["successful_proposal_evaluations_recorded"],
            "fit_failures":fit["fit_failures"]},
        "aggregate_arithmetic_and_scan_checks":totals,"tests":tests,
        "local_checks":{"exact_tests":built["report"]["tests"],"field_tests":built["report"]["field_tests"]},
        "generated_bytes":0,"elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),
        "nonnegative_intersections":v.intersections,
        "pins":{"maximum":MAXIMUM_PIN,"complete3072":maximum.COMPLETE_PIN,"selected":maximum.SELECTED_PIN,
            "run":RUN_PIN,"counts":COUNTS_PIN,"actual":run.ACTUAL_PIN,"support":actual.SUPPORT_PIN,
            "fitter":support.FITTER_PIN,"normalized":support.NORMALIZED_PIN,"band":normalized.BAND_PIN,
            "geometry":normalized.GEOMETRY_PIN,"reference":band.HOLE_PIN,"ordinary":reference.ORDINARY_PIN,
            "local":large.LOCAL_PIN,"source":large.ANCHOR_PIN,"join":large.JOIN_PIN,
            "extraction":maximum.REGIONS_PIN,"integer":maximum.INTEGER_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},separators=(",",":")),flush=True)


if __name__=="__main__":
    main()
