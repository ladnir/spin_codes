"""Choose saved marker pairs after completing each conditional cell bound.

Read the verified compact residue-max witness through the frozen validator.
Rebuild its local table and rational bank for identity checks, but build no
regional polynomial. The upper-bound validity of saved dyadics still relies
on the trusted generator receipt; a matching hash alone does not prove it.

Keep every assigned support cell and its exact rational joint vertices.
For candidate markers p,q, minimize log D_p(r,a)+log D_q(r,b) over those
vertices. Subtract128 times that certified lower bound from the appropriate
finite-chain numerator log. Choose the completed candidate with the smallest
certified upper endpoint, separately for the two terminal orders. No interval
intersection is used. The assigned marker pair remains an unchanged candidate.

Candidate neighborhoods use indices in the31-marker bank. A diagonal SUPPORT
cell is counted once, even when its chosen MARKERS differ. Both ordered marker
chains remain distinct. The original cell count bounds, beta tilt, eta powers,
and single output factor are unchanged. No marker union factor is introduced.

Numerical receipts are indexed separately; source alone is not evidence.
Parent execution only:300/600s,512MiB,20GiB free,-B,stdout,no generated files.
Import performs no numerical work or proof-provider imports.
"""

import argparse
from collections import Counter
from fractions import Fraction as F
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
REPLAY_PIN="77bd0acc081e68691810391c2cd4e9a5dad1f731f14296ec78d4de4831018867"
Z,ALPHA=F(1,10),F(119,1280)
ANCHOR,PRECISION,MEMORY=3072,384,512*(1<<20)


def load_replay():
    path=Path(__file__).with_name("cheap_lane_residue_max_replay.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=REPLAY_PIN:
        raise RuntimeError("frozen compact receipt replayer source mismatch")
    spec=importlib.util.spec_from_file_location("marker_envelope_receipt",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    if (module.Z,module.ALPHA,module.ANCHOR,module.PRECISION)!=(Z,ALPHA,ANCHOR,PRECISION):
        raise ArithmeticError("frozen receipt parameters changed")
    return module


def candidate_indices(index,radius,size=31):
    if any(type(x) is not int for x in (index,radius,size)) or not 1<=size<=31:
        raise ValueError("bounded integer marker-bank layout required")
    if not 0<=index<size or not 0<=radius<=31:
        raise ValueError("invalid marker index or neighborhood radius")
    return tuple(range(max(0,index-radius),min(size,index+radius+1)))


def minimum_endpoint_sum(left,right):
    """Exact minimum of sums of cached LOWER endpoints on matching vertices."""
    if not left or len(left)!=len(right): raise ValueError("incompatible vertex caches")
    low=left[0]+right[0]; where=0
    for index in range(1,len(left)):
        value=left[index]+right[index]
        if value<low: low,where=value,index
    return low,where


def smaller_upper(current,candidate):
    if not candidate.is_finite(): raise ArithmeticError("nonfinite completed candidate")
    return candidate.upper()<current.upper()


def numerator_logs(numerators,join,guard,progress=None):
    """Store all ORDERED finite-chain logs; only p=q permits diagonal reuse."""
    result={}; powers=0
    for p,left in enumerate(numerators):
        for q in range(p,len(numerators)):
            guard()
            forward,reverse=join.two_orientations(left,numerators[q],diagonal=p==q)
            if any(not value.is_finite() or not value>0 for value in (forward,reverse)):
                raise ArithmeticError("saved regional caps give a nonpositive or nonfinite chain")
            result[p,q]=forward.log(); result[q,p]=reverse.log(); powers+=1
        if progress: progress("marker_pair_numerators",row=p)
    if len(result)!=len(numerators)**2: raise ArithmeticError("missing ordered marker numerator")
    return result,{"unordered_chain_powers":powers,"ordered_numerator_logs":len(result),
        "distinct_terminal_orders_preserved":True}


def cell_envelope(v,band,geometry,den,numerators,logs,cells,first,last,radius,
                  join,guard,progress=None):
    """Return unchanged assigned bounds, candidate bounds, and compact provenance.

    Cached endpoint sums can round slightly further down than lower(A+B).
    The frozen assigned enclosure is therefore retained verbatim, including
    at radius0. All other candidates use valid lower bounds on the SAME sum.
    """
    h,numerator_report=numerator_logs(numerators,join,guard,progress)
    old={}; refined={}; chosen={}; empty=vertices_total=ordered=trials=changed=0
    digest=hashlib.sha256(); histogram=Counter(); leaders=[]; logtwo=v.arb(2).log()
    for i,left in enumerate(cells):
        for j in range(i,len(cells)):
            guard(); right=cells[j]
            points=geometry.vertices(first,last,left[1],left[2],right[1],right[2],guard)
            digest.update((str(i)+","+str(j)+":"+repr(points)+"\n").encode("ascii"))
            if not points: empty+=1; continue
            vertices_total+=len(points)
            low,vertex,original=band.evaluate_cell(v,den,points,logs[i],logs[j],
                numerators[i],numerators[j],join,diagonal=i==j)
            old[i,j]=original; ordered+=len(original)
            best={name:value for name,value in original}
            winners={name:(i,j,vertex,low,"frozen_assigned") for name,_ in original}
            # The precomputed ordered numerator must agree with the original
            # evaluator before any marker choice is made.
            for name,value in original:
                key=(i,j) if name=="forward" else (j,i)
                if not value.overlaps(h[key]-128*v.exact(low)):
                    raise ArithmeticError("cached finite terminal chain differs from frozen cell evaluator")
            left_indices=candidate_indices(i,radius,len(cells))
            right_indices=candidate_indices(j,radius,len(cells))
            if radius:
                # Each column keeps the same vertex order and common actual r.
                # Taking these lower endpoints separately is conservative.
                left_cache={p:tuple(band.lower_fraction(den.value(r,a,*logs[p]))
                    for r,a,b in points) for p in left_indices}
                right_cache={q:tuple(band.lower_fraction(den.value(r,b,*logs[q]))
                    for r,a,b in points) for q in right_indices}
                for p in left_indices:
                    guard()
                    for q in right_indices:
                        if (p,q)==(i,j): continue
                        trial_low,where=minimum_endpoint_sum(left_cache[p],right_cache[q])
                        correction=128*v.exact(trial_low); trials+=1
                        for name in best:
                            key=(p,q) if name=="forward" else (q,p)
                            candidate=h[key]-correction
                            if smaller_upper(best[name],candidate):
                                best[name]=candidate
                                winners[name]=(p,q,points[where],trial_low,"candidate_pair")
            row=tuple((name,best[name]) for name,_ in original)
            refined[i,j]=original if radius==0 else row
            for name,oldvalue in original:
                value=best[name]; p,q,point,lower,method=winners[name]
                if value.upper()>oldvalue.upper():
                    raise ArithmeticError("assigned candidate was not retained")
                if (p,q)!=(i,j): changed+=1
                histogram[name,p,q]+=1
                detail={"cell_indices":[i,j],"assigned_centers":[left[0],right[0]],
                    "terminal_order":name,"selected_marker_indices":[p,q],
                    "selected_marker_centers":[cells[p][0],cells[q][0]],"selection":method,
                    "minimum_vertex_actual_full_full":list(map(str,point)),
                    "certified_lower_log_denominator_sum":str(lower),
                    "log2_conditional_interval":str(value/logtwo),
                    "log2_gain_over_assigned_interval":str((oldvalue-value)/logtwo)}
                chosen[i,j,name]=detail; leaders.append((value.upper(),detail))
            if radius==0 and refined[i,j] is not original:
                raise ArithmeticError("radius0 changed a frozen assigned enclosure")
        if progress: progress("marker_envelope_geometry",row=i,center=left[0],vertices=vertices_total)
    if not old or set(old)!=set(refined): raise ArithmeticError("incomplete candidate cell envelope")
    histogram_rows=[[name,p,q,count] for (name,p,q),count in sorted(histogram.items())]
    axis_histogram=Counter()
    for (name,p,q),count in histogram.items():
        axis_histogram[name,"left",p]+=count; axis_histogram[name,"right",q]+=count
    report={"nonempty_unordered_cells":len(old),"empty_unordered_cells":empty,
        "ordered_support_cells":ordered,"rational_vertices":vertices_total,
        "vertices_sha256":digest.hexdigest(),"alternative_marker_pair_trials":trials,
        "changed_ordered_support_cells":changed,"radius":radius,
        "radius_uses_bank_indices":True,"radius0_reuses_frozen_enclosures":radius==0,
        "assigned_candidate_retained_for_every_order":True,
        "no_interval_intersections":True,"no_marker_union_factor":True,
        "diagonal_support_cells_counted_once":True,
        "common_actual_size_and_joint_vertices_preserved":True,
        "loggamma_cache_entries":len(den.gamma),"log_denominator_evaluations":den.calls,
        "numerator_cache":numerator_report,
        "selected_marker_pair_histogram_sha256":hashlib.sha256(
            json.dumps(histogram_rows,separators=(",",":")).encode("ascii")).hexdigest(),
        "most_frequent_selected_marker_pairs":[{"order":name,"indices":[p,q],
            "centers":[cells[p][0],cells[q][0]],"support_cells":count}
            for (name,p,q),count in sorted(histogram.items(),key=lambda item:(-item[1],item[0]))[:16]],
        "selected_marker_axis_histogram":[[cells[p][0],
            *(axis_histogram[name,side,p] for name in ("forward","reverse") for side in ("left","right"))]
            for p in range(len(cells))],
        "axis_histogram_columns":["center","forward_left","forward_right","reverse_left","reverse_right"],
        "largest_conditional_cells":[row for _,row in sorted(leaders,key=lambda x:x[0],reverse=True)[:6]]}
    return old,refined,chosen,report


def exact_tests(join,guard):
    """Exact toy selection, own denominators, and noncommuting terminal checks."""
    checks=0
    assert candidate_indices(0,1,3)==(0,1); checks+=1
    assert candidate_indices(1,0,3)==(1,); checks+=1
    assert candidate_indices(2,31,3)==(0,1,2); checks+=1
    for args in ((True,1,3),(0,-1,3),(3,1,3),(0,32,3),(0,1,32)):
        try: candidate_indices(*args)
        except ValueError: checks+=1
        else: raise AssertionError("invalid marker neighborhood accepted")
    aa=(F(-3,7),F(5,3),F(-1,2)); bb=(F(1,7),F(-7,3),F(3,2))
    low,index=minimum_endpoint_sum(aa,bb)
    assert low==min(a+b for a,b in zip(aa,bb)) and low==aa[index]+bb[index]; checks+=1
    regions=tuple(((F(i+2,11),F(1,7),F(3,13),F(i+5,17)),
                   (F(i+7,19),F(2*i+3,23))) for i in range(3))
    markers=((F(2,3),F(5,7)),(F(7,5),F(3,4)),(F(11,7),F(9,5)))
    def direct(left,right,bits):
        row=(F(1),F(0)); a,_=left; b,terminal=right
        for matrix in (a,b)*((1<<bits)-1)+(a,):
            row=(row[0]*matrix[0]+row[1]*matrix[2],row[0]*matrix[1]+row[1]*matrix[3])
        return row[0]*terminal[0]+row[1]*terminal[1]
    h={}; bits=2; exponent=1<<bits
    for p in range(3):
        for q in range(p,3):
            guard(); forward,reverse=join.two_orientations(regions[p],regions[q],bits=bits,diagonal=p==q)
            assert forward==direct(regions[p],regions[q],bits); checks+=1
            assert reverse==direct(regions[q],regions[p],bits); checks+=1
            h[p,q]=forward; h[q,p]=reverse
    assert h[0,1]!=h[1,0]; checks+=1
    def denominator(full,marker):
        single=6-2*full; u,w=marker
        return F(comb(8,full)*comb(8-full,single))*u**full*w**single
    ranges=((0,1),(2,),(3,)); weighted={radius:F(0) for radius in (0,1,3)}
    for i in range(3):
        for j in range(i,3):
            points=tuple((a,b) for a in ranges[i] for b in ranges[j] if a+b<=5)
            if not points: continue
            for radius in (0,1,3):
                guard(); candidates=tuple((p,q) for p in candidate_indices(i,radius,3)
                    for q in candidate_indices(j,radius,3))
                for reverse in range(1 if i==j else 2):
                    upper={}
                    for p,q in candidates:
                        own=tuple(denominator(a,markers[p])*denominator(b,markers[q]) for a,b in points)
                        numerator=h[q,p] if reverse else h[p,q]
                        upper[p,q]=numerator/min(own)**exponent
                        assert upper[p,q]==max(numerator/d**exponent for d in own); checks+=1
                    chosen=min(upper.values()); baseline=upper[i,j]
                    assert chosen<=baseline; checks+=1
                    # Each completed candidate is uniform over this support
                    # cell. Minimum AFTER completion still dominates the
                    # pointwise minimum, without a marker-count multiplier.
                    pointwise=max(min((h[q,p] if reverse else h[p,q])/
                        (denominator(a,markers[p])*denominator(b,markers[q]))**exponent
                        for p,q in candidates) for a,b in points)
                    assert chosen>=pointwise; checks+=1
                    if radius==0: assert chosen==baseline; checks+=1
                    weighted[radius]+=(i+j+1)*chosen
    assert weighted[3]<=weighted[1]<=weighted[0]; checks+=1
    return {"exact_checks":checks,"noncommuting_terminal_orders_distinguished":True,
        "diagonal_support_cell_counted_once_with_unequal_markers":True,
        "minimum_of_completed_candidates_matches_direct_toy":True,
        "radius0_matches_assigned_toy":True,"no_marker_union_factor":True}


def evaluate_band(v,core,run,actual,normalized,band,geometry,counts_module,anchor,join,
                  cells,numerators,marker_logs,first,last,radius,guard,progress=None):
    den=actual.actual_denominators(v,core.BETA,ALPHA,band,guard)
    cover=actual.actual_coverage(normalized,cells,first,last,guard)
    old,new,chosen,geometry_report=cell_envelope(v,band,geometry,den,numerators,marker_logs,
        cells,first,last,radius,join,guard,progress)
    # An independent frozen radius0 reconstruction protects the source/count
    # seam and exact geometry. It does not rebuild any regional polynomial.
    frozen,frozen_report=run.cell_bounds(v,band,geometry,den,numerators,marker_logs,cells,
        first,last,"actual",join,guard)
    if set(old)!=set(frozen) or geometry_report["vertices_sha256"]!=frozen_report["vertices_sha256"]:
        raise ArithmeticError("assigned geometry differs from frozen replay")
    for key,row in old.items():
        if tuple(name for name,_ in row)!=tuple(name for name,_ in frozen[key]):
            raise ArithmeticError("assigned terminal orders differ")
        for (_,a),(_,b) in zip(row,frozen[key]):
            if not a.overlaps(b): raise ArithmeticError("radius0 differs from frozen conditional bound")
    counts=counts_module.source_counts(v,anchor,join,cells,ALPHA,first,last,"actual",guard)
    old_weight,_=band.outer_source_sum(v,anchor,ALPHA,first,last,guard)
    if not counts["full_weight"].overlaps(old_weight): raise ArithmeticError("source/count factor changed")
    old_join=counts_module.combine(v,counts,old,guard)
    new_join=counts_module.combine(v,counts,new,guard)
    oldmax=None
    for row in frozen.values():
        for _,value in row: oldmax=value if oldmax is None else oldmax.max(value)
    if not old_join["old_log_bound"].overlaps(old_weight.log()+oldmax):
        raise ArithmeticError("frozen global-max source control changed")
    if radius==0:
        for name in ("old_log_bound","run_log_bound","minimum_log_bound"):
            if not old_join[name].overlaps(new_join[name]): raise ArithmeticError("radius0 join changed")
    baseline=old_join["minimum_log_bound"]; selected=new_join["minimum_log_bound"]
    if not smaller_upper(baseline,selected): selected=baseline
    if selected.upper()>baseline.upper(): raise ArithmeticError("completed baseline was not retained")
    logtwo=v.arb(2).log(); output=-v.exact(ALPHA)*core.LIMIT*v.exact(Z).log()
    final=output+selected
    count_report=dict(counts["report"])
    count_report["indexed_support_count_decimal_digits"]=len(count_report.pop("indexed_support_count"))
    join_report=dict(new_join["report"])
    enriched=[]
    for row in join_report["largest_run_weighted_cells"][:6]:
        i,j=row["cell_indices"]
        enriched.append({**row,"log2_weighted_cell_with_output_interval":str(
            v.arb(row["log2_weighted_cell_interval"])+output/logtwo),
            "chosen_conditional_bounds":[chosen[i,j,name] for name,_ in new[i,j]]})
    join_report["largest_run_weighted_cells"]=enriched
    return {"actual_first":first,"actual_last":last,"sizes":last-first+1,
        "all_integer_sizes_and_indexed_supports_covered":True,
        "log2_frozen_assigned_old_global_max_interval":str((output+old_join["old_log_bound"])/logtwo),
        "log2_frozen_assigned_run_sum_interval":str((output+old_join["run_log_bound"])/logtwo),
        "log2_frozen_assigned_minimum_interval":str((output+baseline)/logtwo),
        "log2_marker_choice_old_global_max_interval":str((output+new_join["old_log_bound"])/logtwo),
        "log2_marker_choice_run_sum_interval":str((output+new_join["run_log_bound"])/logtwo),
        "log2_union_bound_interval":str(final/logtwo),
        "log2_gain_over_frozen_assigned_interval":str((baseline-selected)/logtwo),
        "below_2_pow_minus_40":bool(final<-40*logtwo),
        "below_2_pow_minus_60":bool(final<-60*logtwo),
        "log2_eta_weighted_support_count_interval":str(old_weight.log()/logtwo),
        "log2_output_correction_interval":str(output/logtwo),
        "coverage":{key:value for key,value in cover.items() if key!="marginals_by_center"},
        "geometry":geometry_report,"counts":count_report,"join":join_report,
        "frozen_radius0_geometry_and_conditional_controls_passed":True,
        "source_sum_and_old_max_controls_passed":True,
        "completed_frozen_baseline_retained":True,"certified_upper_endpoint_selection":True}


def main():
    replay=load_replay()
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--witness",type=Path,required=True)
    parser.add_argument("--expect-witness-sha",type=replay.sha256_argument,required=True)
    parser.add_argument("--first",type=int,default=3008)
    parser.add_argument("--last",type=int,default=3136)
    parser.add_argument("--radius",type=int,choices=(0,1,2,31),default=1)
    parser.add_argument("--seconds",type=int,choices=(300,600),default=300)
    args=parser.parse_args()
    if not 641<=args.first<=args.last<=4096: raise ValueError("actual band must lie in641..4096")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    witness,input_report=replay.read_witness(args.witness,args.expect_witness_sha)
    generator=replay.load_generator(); maximum=generator.load_maximum()
    complete=maximum.load_complete(); selected=complete.load_selected(); run=selected.load_run()
    if (selected.RUN_PIN,run.COUNTS_PIN)!=(generator.RUN_PIN,generator.COUNTS_PIN):
        raise ArithmeticError("frozen run/count provenance changed")
    actual=run.load_actual(); support=actual.load_support()
    normalized=support.load_normalized(); band=normalized.load_band()
    counts_module=band.load("cheap_lane_cell_run_counts.py",generator.COUNTS_PIN,"choice_counts")
    fitter=band.load(support.FITTER_FILE,support.FITTER_PIN,"choice_fitter")
    geometry=band.load(normalized.GEOMETRY_FILE,normalized.GEOMETRY_PIN,"choice_geometry")
    reference=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"choice_reference")
    ordinary=reference.load("cheap_lane_ordered_ordinary_directed.py",reference.ORDINARY_PIN,"choice_ordinary")
    large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"choice_source")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"choice_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"choice_normalizers")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"choice_terminals")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("marker-envelope deadline; no incomplete claim")
        if peak()>MEMORY: raise MemoryError("marker-envelope512MiB limit exceeded")
    def progress(phase,**values):
        if phase=="actual_marker_fit" and values.get("center") not in (0,512,1024): return
        if phase in ("marker_pair_numerators","marker_envelope_geometry") and values.get("row",0)%8: return
        kept={key:values[key] for key in ("center","status","evaluations","row","vertices") if key in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**kept}),flush=True)
    tests={"marker_envelope":exact_tests(join,guard),
        "receipt_validation":replay.exact_tests(reference.BANK,guard),
        "run_wrapper":run.wrapper_tests(guard),"bank":support.wrapper_tests(reference.BANK,guard),
        "actual":actual.exact_tests(reference.BANK,fitter,normalized,geometry,guard),
        "run_counts":counts_module.exact_tests(join,guard)}
    built=local.build_local_pairs(Z,ALPHA,anchor,guard,deadline,PRECISION,progress=progress)
    v,core=built["verifier"],built["core"]; degrees=built["rational_coefficients"]
    original_arb=dict(built["coefficients"])
    if (support.coefficient_fingerprint(degrees,guard),built["report"]["rational_coefficient_sha256"])!=\
            (replay.TABLE_SHA,replay.TABLE_SHA): raise ArithmeticError("local table differs from receipt")
    proposed,fit=actual.fit_actual_markers(MappingProxyType(degrees),ANCHOR,fitter,guard,progress)
    fitted,_,bank_hash=support.validate_bank(proposed,reference.BANK)
    if bank_hash!=replay.BANK_SHA: raise ArithmeticError("regenerated bank differs from receipt")
    if support.coefficient_fingerprint(degrees,guard)!=replay.TABLE_SHA or len(original_arb)!=len(built["coefficients"]):
        raise ArithmeticError("fitter changed local coefficients")
    if any(row is not built["coefficients"][key] for key,row in original_arb.items()):
        raise ArithmeticError("fitter changed original Arb rows")
    exact,validation=replay.validate_witness(witness,args.expect_witness_sha,fitted,guard)
    numerators=tuple(tuple(tuple(v.exact(x) for x in row) for row in region) for region in exact)
    cells=normalized.normalized_cells(fitted)
    marker_logs=tuple((v.exact(cell[3]).log(),v.exact(cell[4]).log()) for cell in cells)
    den=actual.actual_denominators(v,core.BETA,ALPHA,band,guard)
    tests["directed_actual"]=actual.directed_tests(v,den,extraction,fitted,ALPHA,core.BETA,join,band,guard)
    del den
    result=evaluate_band(v,core,run,actual,normalized,band,geometry,counts_module,anchor,join,
        cells,numerators,marker_logs,args.first,args.last,args.radius,guard,progress)
    guard()
    print(json.dumps({"scope":"complete actual-size band with completed marker-pair choices from a verified maximum receipt",
        "directed_arithmetic":True,"whole_code_certificate":False,
        "trusted_verified_generator_receipt_required":True,
        "regional_upper_bounds_independently_recomputed":False,"regional_tables_built":0,
        "regional_polynomial_products":0,"local_table_rebuilt_for_provenance":True,
        "parameters":{"actual_first":args.first,"actual_last":args.last,"anchor_actual":ANCHOR,
            "z":str(Z),"alpha":str(ALPHA),"radius":args.radius,
            "recorded_integer_bits":witness["integer_bits"],
            "working_precision":built["report"]["working_precision"],"seconds":args.seconds,
            "memory_budget_bytes":MEMORY},
        "requested_band":result,"input":input_report,"validation":validation,
        "local_coefficient_sha256":replay.TABLE_SHA,"positive_rational_bank_sha256":bank_hash,
        "own_joint_target_denominators_preserved":True,"assigned_support_counts_unchanged":True,
        "alpha_not_reapplied_to_numerators":True,"degree_reversals":0,
        "beta_only_in_denominator_tilt":True,"eta_powered_once_counts_outside_alpha":True,
        "output_factor_once":True,"both_offdiagonal_terminal_orders_charged":True,
        "fit":{"scope":"frozen proposal procedure verifies exact witness markers only",
            "successful_proposal_evaluations_recorded":fit["successful_proposal_evaluations_recorded"],
            "fit_failures":fit["fit_failures"]},
        "tests":tests,"local_checks":{"tests":built["report"]["tests"],"field_tests":built["report"]["field_tests"]},
        "elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "nonnegative_intersections":v.intersections,"pins":{"replay":REPLAY_PIN,
            "generator":replay.GENERATOR_PIN,"maximum":replay.MAXIMUM_PIN,
            "complete3072":maximum.COMPLETE_PIN,"selected":maximum.SELECTED_PIN,
            "run":generator.RUN_PIN,"counts":generator.COUNTS_PIN,"actual":run.ACTUAL_PIN,
            "support":actual.SUPPORT_PIN,"fitter":support.FITTER_PIN,"normalized":support.NORMALIZED_PIN,
            "band":normalized.BAND_PIN,"geometry":normalized.GEOMETRY_PIN,"reference":band.HOLE_PIN,
            "ordinary":reference.ORDINARY_PIN,"local":large.LOCAL_PIN,"source":large.ANCHOR_PIN,
            "join":large.JOIN_PIN,"normalizers":large.REGION_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},separators=(",",":")),flush=True)


if __name__=="__main__":
    main()
