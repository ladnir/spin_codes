"""Complete directed actual3064..3080 band from the fixed3072 witness.

Preserve the recorded current local table, z1/10,alpha119/1280 and31
positive rational markers.  Keep the3072 marginal-to-marker assignment;
extend only the last full-count interval through1540.  For EACH actual
weight, clip its ranges at floor(r/2), freshly extract(f,r-2*f) from the
shared cyclic256 table and divide by that target's OWN exact denominator.
The same31 polynomial tables serve all17 sizes.  No interpolation of
normalized regions, hole reversal, endpoint shortcut or fallback occurs.

The frozen count-aware band join takes all own-r R4/V2 arrays and source
scales s_r=(eta_r*beta^r)^(alpha/256).  At each d it absorbs the EXACT run
count through a1/256 root, takes its valid cross-r envelope, and retains
both terminal orders.  The caller then applies only z^(-alpha*L), once.
Actual source h changes at3072; no common eta or source size is substituted.

Independently replay singleton3072 from these same arrays, without adding
that diagnostic to the band.  At192 regional bits it must overlap the
recorded full3072 interval.  At256 bits it instead checks a conservative
threshold; that is a regional quantization control, not a new Arb precision.

Parent execution only:1500/1800s,512MiB,20GiB free,-B,stdout,no data files.
Source and numerical receipts are separate.  No provider imports or work
on module import.  A complete band result is not whole-code closure.
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
from types import MappingProxyType

sys.dont_write_bytecode=True
COMPLETE_PIN="e599c0cd7ed7fd4f594f7b045aeb3702e9ded8fe0d2128ff0fd7dd6f47417b53"
COUNTED_PIN="2c481652ef8a37a14c5bd88ce6af31e0517a3f51706b7a07e3aa6447906ab74b"
FIRST,LAST,ANCHOR,POINTS,PRECISION=3064,3080,3072,8192,384
MEMORY=512*(1<<20)
EXPECTED_ORDERED,EXPECTED_POWERS=2374680,1188110
CONTROL_3072=(
    "-3356.7491176288611133963989569096663623099291764354924315215656796089412766978847349530488304140983050546941498112",
    "8.53e-110")


def load_complete():
    path=Path(__file__).with_name("cheap_lane_complete3072_directed.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=COMPLETE_PIN:
        raise RuntimeError("frozen complete3072 source mismatch")
    spec=importlib.util.spec_from_file_location("band3072_complete",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def schedule(complete,cells,count_module,weights,guard):
    """Reuse all anchor ranges, extending only the final positive witness."""
    original,anchor_report=complete.schedule(cells,count_module,guard)
    rows=list(original)
    item,marker=rows[-1]
    rows[-1]=((item[0],item[1],max(weights)//2,item[3]),marker)
    if any(rows[i]!=original[i] for i in range(len(rows)-1)):
        raise ArithmeticError("nonterminal marker interval changed")
    coverage=[]; digest=hashlib.sha256()
    for r in weights:
        guard(); observed=[]
        for item,_ in rows:
            center,lo,hi,_=item
            values=tuple(range(lo,min(hi,r//2)+1))
            if any(not 0<=r-2*f<=4096-f for f in values):
                raise ArithmeticError("invalid own-r matching target")
            observed.extend(values)
            digest.update((str(r)+":"+str(center)+":"+str(values)+"\n").encode("ascii"))
        if observed!=list(range(r//2+1)):
            raise ArithmeticError("own-r marginal assignment is incomplete or duplicated")
        coverage.append({"weight":r,"marginals":len(observed),"first_full":0,"last_full":r//2})
    if len(rows)!=31 or rows[-1][0][2]!=1540:
        raise ArithmeticError("fixed seventeen-size schedule changed")
    return tuple(rows),{"anchor_assignment":anchor_report,"own_weight_coverage":coverage,
        "marker_assignment":"fixed3072 ranges; last interval extended; each r clipped at floor(r/2)",
        "last_interval_full_limit":1540,"assignment_sha256":digest.hexdigest(),
        "all_integer_sizes":len(weights),"total_own_marginals":sum(x["marginals"] for x in coverage)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bits",type=int,choices=(192,256),default=192)
    parser.add_argument("--precision",type=int,choices=(384,),default=384)
    parser.add_argument("--seconds",type=int,choices=(1500,1800),default=1500)
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    weights=tuple(range(FIRST,LAST+1)); complete=load_complete()
    if (complete.WEIGHT,complete.POINTS,complete.PRECISION)!=(ANCHOR,POINTS,PRECISION):
        raise ArithmeticError("frozen singleton setup changed")
    z,alpha=complete.Z,complete.ALPHA
    selected=complete.load_selected(); run=selected.load_run(); actual=run.load_actual()
    support=actual.load_support(); normalized=support.load_normalized(); band=normalized.load_band()
    count_module=band.load("cheap_lane_cell_run_counts.py",complete.COUNTS_PIN,"band3072_assignment")
    fitter=band.load(support.FITTER_FILE,support.FITTER_PIN,"band3072_fitter")
    reference=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"band3072_reference")
    ordinary=reference.load("cheap_lane_ordered_ordinary_directed.py",reference.ORDINARY_PIN,
                            "band3072_ordinary_controls")
    large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"band3072_source")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"band3072_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"band3072_extraction")
    integer=large.load("cheap_lane_integer_regions.py",large.INTEGER_PIN,"band3072_integer")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"band3072_single_join")
    counted=large.load("cheap_lane_counted_band_run_mass.py",COUNTED_PIN,"band3072_counted_join")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("complete3072-band deadline; no incomplete claim")
        if peak()>MEMORY: raise MemoryError("complete3072-band exceeds512MiB")
    def progress(phase,**values):
        if phase=="actual_marker_fit" and values.get("center") not in (0,512,1024): return
        keys=("center","status","evaluations","first_full","last_full","tables_completed",
              "run","degree_sum","active_weights","ordered_profiles","ordered_envelope_pairs",
              "unordered_powers","working_precision")
        kept={key:values[key] for key in keys if key in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**kept}),flush=True)
    _,own_tests=selected.own_tests(extraction,join,anchor,guard)
    tests={"selected_controls":own_tests,"extraction":extraction.exact_tests(guard),
        "integer":integer.exact_tests(guard),"ordered_terminal":join.exact_tests(guard),
        "ordinary":ordinary.exact_tests(guard),"bank":support.wrapper_tests(reference.BANK,guard)}
    built=local.build_local_pairs(z,alpha,anchor,guard,deadline,PRECISION,progress=progress)
    v,core=built["verifier"],built["core"]; logtwo=v.arb(2).log()
    if built["report"]["working_precision"]!=PRECISION:
        raise ArithmeticError("working precision changed")
    tests["frozen_counted_public_interface"]=counted.exact_tests(join,guard)
    degrees=built["rational_coefficients"]; arb_identity=dict(built["coefficients"])
    table_hash=support.coefficient_fingerprint(degrees,guard)
    if table_hash!=complete.TABLE_SHA or built["report"]["rational_coefficient_sha256"]!=table_hash:
        raise ArithmeticError("local table differs from the recorded3072 witness")
    proposed,fit_report=actual.fit_actual_markers(MappingProxyType(degrees),ANCHOR,fitter,guard,progress)
    fitted,public_bank,bank_hash=support.validate_bank(proposed,reference.BANK)
    if bank_hash!=complete.BANK_SHA:
        raise ArithmeticError("marker bank differs from the recorded3072 witness")
    if support.coefficient_fingerprint(degrees,guard)!=table_hash:
        raise ArithmeticError("fitter changed the exact local table")
    if len(arb_identity)!=len(built["coefficients"]) or any(
            row is not built["coefficients"][key] for key,row in arb_identity.items()):
        raise ArithmeticError("fitter changed original Arb rows")
    rows,coverage=schedule(complete,normalized.normalized_cells(fitted),count_module,weights,guard)
    regions={r:{} for r in weights}; reports=[]
    for index,(item,marker) in enumerate(rows):
        guard(); center,lo,hi,_=item
        progress("complete3072_band_table_start",center=center,first_full=lo,last_full=hi)
        partial,report=extraction.one_table_many(v,degrees,item,marker,weights,
            integer,args.bits,guard,None)
        for r in weights:
            own=partial[r]
            if tuple(own)!=tuple(range(lo,min(hi,r//2)+1)) or set(own).intersection(regions[r]):
                raise ArithmeticError("shared extraction lost or duplicated own-r targets")
            regions[r].update(own)
        reports.append({"center":center,"first_full":lo,"last_full":hi,
            "own_targets":report["extracted_targets"],"u":str(marker[0]),"v":str(marker[1]),
            "pair_doublings":report["pair_doublings"],
            "cyclic_DC_upper_checks":report["cyclic_DC_upper_checks"],
            "DC_vs_unrounded_log2_ratios":report["DC_vs_unrounded_log2_ratios"],
            "arithmetic":report["arithmetic"]})
        del partial,own; gc.collect(); guard()
        progress("complete3072_band_table_complete",center=center,first_full=lo,last_full=hi,
                 tables_completed=index+1)
    for r in weights:
        if set(regions[r])!=set(range(r//2+1)):
            raise ArithmeticError("incomplete original own-r coverage; join forbidden")
    scales={}; sources={}; expected_count=0
    output=v.exact(z)**v.exact(-alpha*core.LIMIT)
    for r in weights:
        guard(); h,c_h,eta=anchor.source_factor(r)
        omega=v.exact(eta)*v.exact(core.BETA)**r
        scales[r]=omega**v.exact(alpha/F(256))
        direct=(omega*v.exact(z)**(-core.LIMIT))**v.exact(alpha)
        if not (scales[r]**256*output).overlaps(direct):
            raise ArithmeticError("own-r source scale/output factorization changed")
        number=comb(POINTS,r); expected_count+=number
        sources[str(r)]={"h":h,"c_h":str(c_h),"eta":str(eta),
            "indexed_support_count_sha256":hashlib.sha256(str(number).encode("ascii")).hexdigest()}
    if (sources["3071"]["h"],sources["3071"]["c_h"],sources["3072"]["h"],
            sources["3072"]["c_h"])!=(11,"1/3072",12,"1/13"):
        raise ArithmeticError("actual source h boundary control failed")
    tests["own_source_factorization_checks"]=len(weights)
    tests["source_h_boundary3071_3072"]=True
    singleton=join.directed_join(tuple(regions[ANCHOR][f] for f in range(ANCHOR//2+1)),
        guard,progress,weight=ANCHOR,points=POINTS)
    _,_,eta=anchor.source_factor(ANCHOR)
    single_source=(v.exact(eta)*v.exact(core.BETA)**ANCHOR*v.exact(z)**(-core.LIMIT))**v.exact(alpha)
    single_total=v.positive(single_source*singleton["bound"]); single_log=single_total.log()/logtwo
    if (singleton["report"]["ordered_profiles"]!=complete.ORDERED_PROFILES
            or singleton["report"]["unordered_powers"]!=complete.UNORDERED_POWERS):
        raise ArithmeticError("singleton own-array profile coverage changed")
    overlap=single_log.overlaps(v.arb(*CONTROL_3072))
    if args.bits==192 and not overlap:
        raise ArithmeticError("singleton3072 fails recorded192-bit overlap")
    if not single_total<v.arb(2)**(-3356):
        raise ArithmeticError("singleton3072 fails its conservative recorded threshold")
    singleton_report={"weight":ANCHOR,"log2_union_bound_interval":str(single_log),
        "recorded192_bit_interval_overlap":bool(overlap),"overlap_required":args.bits==192,
        "below_2_pow_minus_3356":True,"join":singleton["report"],
        "not_added_separately_to_band":True}
    del singleton; gc.collect(); guard()
    joined=counted.counted_join(regions,scales,join,guard,progress=progress)
    if (joined["report"]["ordered_envelope_pairs"]!=EXPECTED_ORDERED
            or joined["report"]["unordered_powers"]!=EXPECTED_POWERS
            or not joined["report"]["all_band_counts_sum_to_binomials"]
            or joined["report"]["band_count_sum_sha256"]!=
                hashlib.sha256(str(expected_count).encode()).hexdigest()):
        raise ArithmeticError("complete counted-band coverage/checksum mismatch")
    total=v.positive(output*joined["bound"])
    largest=sorted(range(len(joined["degree_bounds"])),
        key=lambda d:joined["degree_bounds"][d].upper(),reverse=True)[:8]
    top=[{"degree_sum":d,"log2_output_charged_interval":
        str((output*joined["degree_bounds"][d]).log()/logtwo)} for d in largest]
    guard()
    print(json.dumps({"scope":"complete actual3064..3080 band, all indexed supports and histories at35%",
        "directed_arithmetic":True,"all_integer_sizes_in_band":True,"whole_code_certificate":False,
        "actual_weights":weights,"log2_union_bound_interval":str(total.log()/logtwo),
        "band_below_2_pow_minus_40":bool(total<v.arb(2)**(-40)),
        "band_below_2_pow_minus_60":bool(total<v.arb(2)**(-60)),
        "parameters":{"anchor":ANCHOR,"z":str(z),"alpha":str(alpha),"regional_moduli":[256,256],
            "regional_integer_bits":args.bits,"working_Arb_precision":PRECISION,
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "coverage":coverage,"tables":reports,"singleton3072_control":singleton_report,
        "join":joined["report"],"largest_degree_terms":top,"actual_source_factors":sources,
        "log2_common_output_factor_interval":str(output.log()/logtwo),
        "source_scale":"(eta_r*beta^r)^(alpha/256)","count_root":"N_r(r-d)^(1/256)",
        "source_already_in_join":True,"output_factor_once_after_join":True,
        "no_additional_source_or_count_factor":True,"all_own_target_denominators":True,
        "fresh_actual_degree_extraction":True,"interpolated_normalized_coefficients":False,
        "boundary_replacement":False,"ordinary_fallbacks":0,"degree_reversals":0,
        "local_table_builds":1,"local_coefficient_sha256":table_hash,
        "positive_rational_bank_sha256":bank_hash,"baseline_hashes_enforced":True,
        "positive_rational_bank":public_bank,
        "fit":{"scope":"frozen3072 proposal-only fit, reused unchanged for every actual weight",
            "successful_proposal_evaluations_recorded":fit_report["successful_proposal_evaluations_recorded"],
            "fit_failures":fit_report["fit_failures"]},
        "tests":tests,"local_coefficients":built["report"],
        "nonnegative_intersections":v.intersections,"elapsed_seconds":time.monotonic()-started,
        "peak_process_bytes":peak(),"generated_bytes":0,
        "pins":{"complete3072":COMPLETE_PIN,"selected":complete.SELECTED_PIN,
            "counts":complete.COUNTS_PIN,"counted_join":COUNTED_PIN,"counted_base":counted.BASE_PIN,
            "run":selected.RUN_PIN,"actual":run.ACTUAL_PIN,"support":actual.SUPPORT_PIN,
            "fitter":support.FITTER_PIN,"normalized":support.NORMALIZED_PIN,"band":normalized.BAND_PIN,
            "reference":band.HOLE_PIN,"ordinary":reference.ORDINARY_PIN,"local":large.LOCAL_PIN,
            "source":large.ANCHOR_PIN,"extraction":large.REGION_PIN,"integer":large.INTEGER_PIN,
            "single_join":large.JOIN_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
