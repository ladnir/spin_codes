"""Bounded directed refinement of calibrated run-aware cell contributions.

Hold one z,alpha, coherent local table and fitted rational bank fixed.
Compute the COMPLETE ordinary/global-max and cell-run-mass controls.
Rank unordered cells by their complete weighted contribution and refine
only the largest selected cells into disjoint INTEGER run slabs.

For each actual size, each slab receives the exact indexed run mass in
the intersection of its run range with that cell's assigned integer
range. These masses sum exactly to the original cell mass before source
weighting. eta_r^alpha uses the actual size; beta remains only in the256
denominator tilts. The new geometry adds c_lo<=k-a-b<=c_hi to the frozen
normalized cell. Integer targets are covered even though adjacent closed
slabs leave gaps containing only noninteger run values.

For an off-diagonal cell, sum BOTH ordered terminal contributions within
every slab. Minimize the completed slab sum and completed original CELL
contribution, retain every unselected cell unchanged, and finally retain
the old complete global envelope too. The global output factor is applied
once. No local-cap change, refit, new alpha, or event independence is used.

Actual641..4096 and hole4097..7551 modes retain the frozen same-side
domain contract. Default singleton actual3072 uses z1/10,alpha119/1280.
Parent execution only:300s default/600s explicit,512MiB,20GiB free,-B,
stdout only. No imports of proof providers or work occur on import.
Numerical receipts are separate; source alone is not numerical evidence.
"""

import argparse
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
RUN_PIN="19914e8bc698c410971f8edb0ab7f884cb0b83955bee2a6933e9791aa9dfcfc2"
GEOMETRY_PIN="9428a131913d6d6ac590393c4b73b988889d5e327bce54362ae1f6fecdeeb7b9"
MEMORY=512*(1<<20)


def load_run():
    path=Path(__file__).with_name("cheap_lane_calibrated_run_band.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=RUN_PIN:
        raise RuntimeError("frozen run-aware band driver mismatch")
    spec=importlib.util.spec_from_file_location("slabs_frozen_run",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def integer_slabs(lo,hi,width,maximum):
    """Yield(global slab lo,hi, intersection lo,hi), all exact integers."""
    if any(type(x) is not int for x in (lo,hi,width,maximum)):
        raise TypeError("integer slab arguments required")
    if not 1<=lo<=hi<=maximum or width<1:
        raise ValueError("invalid inclusive integer slab interval")
    first=1+((lo-1)//width)*width
    for lower in range(first,hi+1,width):
        upper=min(maximum,lower+width-1)
        yield lower,upper,max(lo,lower),min(hi,upper)


def exact_tests(guard):
    checks=cases=0
    for maximum in (1,2,7,8,17):
        prefix=[0]
        for c in range(1,maximum+1): prefix.append(prefix[-1]+c*c+2*c+3)
        for width in (1,2,3,4,16,32):
            for lo in range(1,maximum+1):
                for hi in range(lo,maximum+1):
                    guard(); pieces=tuple(integer_slabs(lo,hi,width,maximum))
                    covered=[c for _,_,a,b in pieces for c in range(a,b+1)]
                    assert covered==list(range(lo,hi+1)); checks+=1
                    assert len({(a,b) for a,b,_,_ in pieces})==len(pieces); checks+=1
                    assert all(1<=a<=c<=d<=b<=maximum for a,b,c,d in pieces); checks+=1
                    mass=sum(prefix[d]-prefix[c-1] for _,_,c,d in pieces)
                    assert mass==prefix[hi]-prefix[lo-1]; checks+=1
                    # Both terminal orders are charged inside each off-diagonal slab.
                    forward=F(2,7); reverse=F(5,11); eta=F(3,13)
                    old=eta*mass*(forward+reverse)
                    split=sum(eta*(prefix[d]-prefix[c-1])*(forward+reverse)
                              for _,_,c,d in pieces)
                    assert split==old; checks+=1; cases+=1
    # Per-cell minima retain unrefined cells and both complete old controls.
    old_cells=(F(11,5),F(7,4),F(9,8),F(3,2))
    true_cells=tuple(x/3 for x in old_cells)
    candidates={0:old_cells[0]*F(2,3),2:old_cells[2]*F(4,3)}
    result=sum(min(old,candidates[i]) if i in candidates else old for i,old in enumerate(old_cells))
    assert sum(true_cells)<=result<=sum(old_cells); checks+=1
    other_complete=sum(true_cells)*2
    assert sum(true_cells)<=min(result,other_complete); checks+=1
    for bad in ((0,1,1,2),(2,1,1,2),(1,3,1,2),(1,2,0,2)):
        try: tuple(integer_slabs(*bad))
        except ValueError: checks+=1
        else: raise AssertionError("invalid run interval accepted")
    try: tuple(integer_slabs(1,2,True,2))
    except TypeError: checks+=1
    else: raise AssertionError("boolean width accepted")
    return {"exact_checks":checks,"integer_interval_cases":cases,
        "partition_mass_identity":True,"both_terminal_orders_inside_each_slab":True,
        "minimum_only_after_complete_cell_sum":True}


def original_terms(v,counts,bounds,guard):
    """Return all positive complete unordered-cell contributions, no output."""
    terms={}; total=v.exact(F(0))
    for pair,weight in zip(counts["pairs"],counts["weights"]):
        guard()
        if weight==0: continue
        row=bounds.get(pair)
        if row is None: raise ArithmeticError("positive cell count has no conditional bound")
        expected=("forward",) if pair[0]==pair[1] else ("forward","reverse")
        if tuple(name for name,_ in row)!=expected:
            raise ArithmeticError("missing ordered terminal contribution")
        value=weight*sum((logbound.exp() for _,logbound in row),v.exact(F(0)))
        if not value.is_finite() or not value>0:
            raise ArithmeticError("nonpositive original cell contribution")
        terms[pair]=value; total+=value
    if not total>0: raise ArithmeticError("empty complete cell-run sum")
    return terms,total


def selected_slab_counts(v,anchor,join,count_module,cells,selected,counts,alpha,
                         first,last,mode,width,guard):
    """Stream one exact run prefix; retain only selected positive slab masses."""
    selected=tuple(selected)
    if len(set(selected))!=len(selected) or any(pair not in counts["pairs"] for pair in selected):
        raise ValueError("invalid selected unordered cells")
    if not selected:
        return {},{"selected_cells":0,"slab_mass_checks":0,"source_sizes_streamed":0,
            "no_cells_selected":True}
    clean=count_module._cells(cells)
    maximum=last if mode=="actual" else 8192-first
    result={pair:{} for pair in selected}; original_integer={pair:0 for pair in selected}
    digest=hashlib.sha256(); source_digest=hashlib.sha256(); exact_checks=0
    exponent=v.exact(alpha)
    for r in range(first,last+1):
        guard(); k=r if mode=="actual" else 8192-r
        run=join.run_counts(points=8192,weight=k)
        prefix=count_module._run_prefix(run,comb(8192,r),guard)
        ranges=count_module.marginal_ranges(clean,k,guard)
        h,c_h,eta=anchor.source_factor(r)
        if c_h<=0 or eta<=0: raise ArithmeticError("invalid actual source factor")
        factor=v.exact(eta)**exponent
        source_digest.update((str(r)+":"+str(h)+":"+str(c_h)+":"+str(eta)+"\n").encode("ascii"))
        for pair in selected:
            guard(); left,right=ranges[pair[0]],ranges[pair[1]]
            if left is None or right is None: continue
            lo=max(1,k-left[1]-right[1]); hi=min(k,k-left[0]-right[0])
            if lo>hi: continue
            old=prefix[hi]-prefix[lo-1]; split=0
            for lower,upper,a,b in integer_slabs(lo,hi,width,maximum):
                mass=prefix[b]-prefix[a-1]
                if mass<=0: raise ArithmeticError("nonpositive exact slab mass")
                split+=mass; key=(lower,upper)
                target=result[pair]
                target[key]=target.get(key,v.exact(F(0)))+mass*factor
            if split!=old: raise ArithmeticError("integer slab masses lost or duplicated a run")
            exact_checks+=1; original_integer[pair]+=old
            digest.update((str(r)+":"+str(pair)+":"+str(lo)+":"+str(hi)+":"+str(old)+"\n").encode("ascii"))
        del run,prefix,ranges
    old_weights=dict(zip(counts["pairs"],counts["weights"]))
    weights={}; mass_hash=hashlib.sha256(); all_slabs=0
    for pair in selected:
        guard(); row=tuple((lo,hi,value) for (lo,hi),value in sorted(result[pair].items()))
        if not row: raise ArithmeticError("positive selected cell has no run slab")
        total=sum((value for _,_,value in row),v.exact(F(0)))
        if not total.overlaps(old_weights[pair]):
            raise ArithmeticError("weighted slab masses disagree with frozen cell weight")
        if any(not value.is_finite() or not value>0 for _,_,value in row):
            raise ArithmeticError("invalid positive slab weight")
        weights[pair]=row; all_slabs+=len(row)
        mass_hash.update((str(pair)+":"+str(original_integer[pair])+"\n").encode("ascii"))
    if source_digest.hexdigest()!=counts["report"]["actual_source_sha256"]:
        raise ArithmeticError("slab source factors differ from the complete run control")
    return weights,{"selected_cells":len(selected),"positive_cell_slabs":all_slabs,
        "width":width,"maximum_coordinate_run":maximum,"source_sizes_streamed":last-first+1,
        "slab_mass_checks":exact_checks,"exact_partition_before_eta":True,
        "weighted_cell_sums_overlap_frozen":True,"source_factors_match_frozen":True,
        "integer_mass_partition_sha256":digest.hexdigest(),"unweighted_selected_mass_sha256":mass_hash.hexdigest(),
        "one_actual_size_prefix_retained":True,"slabs_are_disjoint_integer_intervals":True}


def refine_cells(v,band,slab_geometry,den,numerators,marker_logs,cells,klo,khi,mode,
                 join,terms,selected,weights,guard,progress=None):
    """Retain all cells; replace selected cells only by a completed minimum."""
    chosen={}; reports=[]; slab_leaders=[]; geometry_hash=hashlib.sha256()
    geometry_count=vertex_count=orders=0; logtwo=v.arb(2).log()
    for index,pair in enumerate(selected):
        guard(); i,j=pair; left,right=cells[i],cells[j]
        refined=v.exact(F(0)); local_slabs=local_vertices=0
        for clo,chi,weight in weights[pair]:
            guard()
            points=slab_geometry.vertices(klo,khi,left[1],left[2],right[1],right[2],clo,chi,guard)
            geometry_hash.update((str(pair)+":"+str((clo,chi))+":"+repr(points)+"\n").encode("ascii"))
            if not points: raise ArithmeticError("positive run-slab mass has empty exact geometry")
            low,point,bounds=band.evaluate_cell(v,den,points,marker_logs[i],marker_logs[j],
                numerators[i],numerators[j],join,diagonal=i==j)
            expected=("forward",) if i==j else ("forward","reverse")
            if tuple(name for name,_ in bounds)!=expected:
                raise ArithmeticError("incomplete slab terminal orders")
            term=weight*sum((value.exp() for _,value in bounds),v.exact(F(0)))
            if not term.is_finite() or not term>0:
                raise ArithmeticError("nonpositive directed slab contribution")
            refined+=term; geometry_count+=1; local_slabs+=1
            vertex_count+=len(points); local_vertices+=len(points); orders+=len(bounds)
            k,a,b=point
            slab_leaders.append((term.log().upper(),{"cell_indices":[i,j],"centers":[left[0],right[0]],
                "run_slab":[clo,chi],"minimum_vertex_coordinate_full_full":list(map(str,point)),
                "actual_size_at_minimum_vertex":str(k if mode=="actual" else 8192-k),
                "run_at_minimum_vertex":str(k-a-b),"terminal_orders":len(bounds),
                "log2_weighted_slab_interval":str(term.log()/logtwo)}))
        if not refined>0: raise ArithmeticError("empty completed refined cell")
        old=terms[pair]; chosen[pair]=old.min(refined)
        reports.append({"cell_indices":[i,j],"centers":[left[0],right[0]],
            "positive_slabs":local_slabs,"vertices":local_vertices,
            "log2_original_cell_interval":str(old.log()/logtwo),
            "log2_refined_cell_interval":str(refined.log()/logtwo),
            "log2_chosen_cell_interval":str(chosen[pair].log()/logtwo)})
        if progress: progress("completed_run_slab_cell",row=index,center=left[0],vertices=vertex_count)
    total=sum((chosen.get(pair,value) for pair,value in terms.items()),v.exact(F(0)))
    if not total.is_finite() or not total>0:
        raise ArithmeticError("incomplete or nonpositive refined total")
    complete_leaders=sorted(((chosen.get(pair,value).log().upper(),pair,chosen.get(pair,value))
        for pair,value in terms.items()),key=lambda row:row[0],reverse=True)[:8]
    unrefined=[(value.log().upper(),pair,value) for pair,value in terms.items() if pair not in chosen]
    unrefined.sort(key=lambda row:row[0],reverse=True)
    largest_unrefined=None
    if unrefined:
        _,pair,value=unrefined[0]
        largest_unrefined={"cell_indices":list(pair),"centers":[cells[pair[0]][0],cells[pair[1]][0]],
            "log2_cell_interval":str(value.log()/logtwo)}
    unrefined_total=sum((value for _,_,value in unrefined),v.exact(F(0)))
    return total,{"selected_cells":len(selected),"unchanged_positive_cells":len(terms)-len(selected),
        "slab_geometry_calls":geometry_count,"slab_rational_vertices":vertex_count,
        "slab_terminal_orders":orders,"slab_vertices_sha256":geometry_hash.hexdigest(),
        "completed_cell_minima_only":True,"unselected_cells_retained":True,
        "cell_controls":reports,
        "contribution_logs_exclude_common_output_factor":True,
        "largest_remaining_unrefined_cell":largest_unrefined,
        "log2_remaining_unrefined_sum_interval":str(unrefined_total.log()/logtwo) if unrefined else None,
        "largest_weighted_slabs":[row for _,row in sorted(slab_leaders,key=lambda x:x[0],reverse=True)[:8]],
        "largest_final_cells":[{"cell_indices":list(pair),"centers":[cells[pair[0]][0],cells[pair[1]][0]],
            "refined":pair in chosen,"log2_cell_interval":str(value.log()/logtwo)} for _,pair,value in complete_leaders]}


def main():
    run=load_run(); actual=run.load_actual()
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor",type=int,default=3072)
    parser.add_argument("--first",type=int,default=3072)
    parser.add_argument("--last",type=int,default=3072)
    parser.add_argument("--coordinate",choices=("actual","hole"))
    parser.add_argument("--z",type=actual.rational64,default=F(1,10))
    parser.add_argument("--alpha",type=actual.rational64,default=F(119,1280))
    parser.add_argument("--top-cells",type=int,choices=(0,1,4,8,16,32,64,128,496),default=16)
    parser.add_argument("--slab-width",type=int,choices=(4,8,16,32,64,128),default=32)
    parser.add_argument("--seconds",type=int,choices=(300,600),default=300)
    parser.add_argument("--precision",type=int,choices=(256,384),default=384)
    args=parser.parse_args()
    mode=run.coordinate_mode(args.first,args.last,args.anchor,args.coordinate)
    if not 0<args.z<1 or not 0<args.alpha<=1:
        raise ValueError("require0<z<1 and0<alpha<=1")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    support=actual.load_support(); normalized=support.load_normalized(); band=normalized.load_band()
    counts_module=band.load("cheap_lane_cell_run_counts.py",run.COUNTS_PIN,"slabs_frozen_counts")
    fitter=band.load(support.FITTER_FILE,support.FITTER_PIN,"slabs_fitter")
    geometry=band.load(normalized.GEOMETRY_FILE,normalized.GEOMETRY_PIN,"slabs_base_geometry")
    slab_geometry=band.load("cheap_lane_run_slab_geometry.py",GEOMETRY_PIN,"slabs_geometry")
    reference=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"slabs_reference")
    ordinary=reference.load("cheap_lane_ordered_ordinary_directed.py",reference.ORDINARY_PIN,"slabs_ordinary")
    large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"slabs_source")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"slabs_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"slabs_normalizer")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"slabs_terminal")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("run-slab deadline; no incomplete claim")
        if peak()>MEMORY: raise MemoryError("run-slab512MiB budget exceeded")
    def progress(phase,**values):
        selected_values={key:values[key] for key in ("center","status","evaluations","row","vertices") if key in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**selected_values}),flush=True)
    tests={"wrapper":run.wrapper_tests(guard),"slab_partition":exact_tests(guard),
        "slab_geometry":slab_geometry.exact_tests(geometry,guard),
        "bank":support.wrapper_tests(reference.BANK,guard),"fitting_utilities":fitter.exact_tests(guard),
        "ordinary":ordinary.exact_tests(guard),"run_counts":counts_module.exact_tests(join,guard)}
    built=local.build_local_pairs(args.z,args.alpha,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]; logtwo=v.arb(2).log()
    degrees=built["rational_coefficients"]; arb_identity=dict(built["coefficients"])
    before=support.coefficient_fingerprint(degrees,guard)
    assert before==built["report"]["rational_coefficient_sha256"]
    if mode=="actual":
        proposed,fit_report=actual.fit_actual_markers(MappingProxyType(degrees),args.anchor,fitter,guard,progress)
    else:
        proposed,fit_report=fitter.fit_markers(MappingProxyType(degrees),args.anchor,guard,progress)
    fitted,public_bank,bank_hash=support.validate_bank(proposed,reference.BANK)
    after=support.coefficient_fingerprint(degrees,guard)
    if before!=after or len(arb_identity)!=len(built["coefficients"]):
        raise ArithmeticError("fitter changed the coherent local table")
    assert all(row is built["coefficients"][key] for key,row in arb_identity.items())
    cells=normalized.normalized_cells(fitted)
    if mode=="actual":
        klo,khi=args.first,args.last; coefficients=built["coefficients"]
        den=actual.actual_denominators(v,core.BETA,args.alpha,band,guard)
        cover=actual.actual_coverage(normalized,cells,klo,khi,guard)
        tests["actual_coordinates"]=actual.exact_tests(fitted,fitter,normalized,geometry,guard)
        tests["directed_bank"]=actual.directed_tests(v,den,extraction,fitted,args.alpha,core.BETA,join,band,guard)
    else:
        klo,khi=8192-args.last,8192-args.first
        coefficients={(64-f-s,s):row for (f,s),row in built["coefficients"].items()}
        assert len(coefficients)==2145 and all(row is built["coefficients"][64-f-s,s]
            for (f,s),row in coefficients.items())
        den=band.LogDenominators(v,core.BETA,args.alpha,guard)
        cover=normalized.coverage(cells,klo,khi,guard)
        tests["normalized_coordinates"]=normalized.normalized_tests(cells,fitted,geometry,guard)
        tests["directed_bank"]=band.directed_tests(v,den,extraction,fitted,args.alpha,core.BETA,join,guard)
    numerators,marker_logs=normalized.numerator_tables(v,coefficients,cells,ordinary,guard,progress)
    bounds,geometry_report=run.cell_bounds(v,band,geometry,den,numerators,marker_logs,cells,
        klo,khi,mode,join,guard,progress)
    counts=counts_module.source_counts(v,anchor,join,cells,args.alpha,args.first,args.last,mode,guard)
    old_weight,old_source_report=band.outer_source_sum(v,anchor,args.alpha,args.first,args.last,guard)
    if not counts["full_weight"].overlaps(old_weight):
        raise ArithmeticError("complete source sum differs from frozen control")
    joined=counts_module.combine(v,counts,bounds,guard)
    terms,original=original_terms(v,counts,bounds,guard)
    if not original.log().overlaps(joined["run_log_bound"]):
        raise ArithmeticError("complete original cell contributions differ from frozen join")
    old_max=None
    for row in bounds.values():
        for _,value in row: old_max=value if old_max is None else old_max.max(value)
    if not joined["old_log_bound"].overlaps(old_weight.log()+old_max):
        raise ArithmeticError("old global-max control changed")
    ordered=sorted(terms,key=lambda pair:terms[pair].log().upper(),reverse=True)
    selected=tuple(ordered[:args.top_cells])
    weights,mass_report=selected_slab_counts(v,anchor,join,counts_module,cells,selected,counts,args.alpha,
        args.first,args.last,mode,args.slab_width,guard)
    refined,refinement_report=refine_cells(v,band,slab_geometry,den,numerators,marker_logs,cells,
        klo,khi,mode,join,terms,selected,weights,guard,progress)
    if not refined<=original and not refined.overlaps(original):
        raise ArithmeticError("completed cell minima increased the run-aware sum")
    log_output=-v.exact(args.alpha)*core.LIMIT*v.exact(args.z).log()
    old_log=joined["old_log_bound"]+log_output
    run_log=joined["run_log_bound"]+log_output
    refined_log=refined.log()+log_output
    final_log=joined["minimum_log_bound"].min(refined.log())+log_output
    tests["frozen_old_and_complete_run_controls_match"]=True
    guard()
    print(json.dumps({"scope":"complete inclusive actual-size band, all indexed supports and histories at35%",
        "directed_arithmetic":True,"whole_code_certificate":False,"all_integer_sizes_in_band":True,
        "parameters":{"actual_first":args.first,"actual_last":args.last,"fit_anchor_actual":args.anchor,
            "coordinate":mode,"coordinate_first":klo,"coordinate_last":khi,
            "z":str(args.z),"alpha":str(args.alpha),"top_cells_requested":args.top_cells,
            "top_cells_selected":len(selected),"slab_width":args.slab_width,"seconds":args.seconds,
            "working_precision":built["report"]["working_precision"],"memory_budget_bytes":MEMORY},
        "log2_old_global_max_interval":str(old_log/logtwo),
        "log2_complete_run_control_interval":str(run_log/logtwo),
        "log2_complete_refined_interval":str(refined_log/logtwo),
        "log2_union_bound_interval":str(final_log/logtwo),
        "band_below_2_pow_minus_40":bool(final_log<-40*logtwo),
        "log2_output_correction_interval":str(log_output/logtwo),
        "selection":"top complete unordered-cell contributions by directed log upper endpoint",
        "selected_cell_indices":[list(pair) for pair in selected],
        "original_run_join":joined["report"],"run_counts":counts["report"],
        "old_source_control":old_source_report,"base_geometry":geometry_report,"coverage":cover,
        "slab_mass":mass_report,"refinement":refinement_report,
        "positive_rational_bank":public_bank,"positive_rational_bank_sha256":bank_hash,
        "fit":{"status":"FLOAT proposals only; same rational witnesses used by every complete control",
            "not_a_certificate":True,"report":fit_report},
        "local_table_builds":1,"alpha_explicitly_rebuilt":True,"refits_during_refinement":0,
        "local_coefficient_sha256_before_fit":before,"local_coefficient_sha256_after_fit":after,
        "arb_coefficients_preserved_by_identity":True,"directed_degree_reversals":int(mode=="hole"),
        "source_beta_only_in_denominator_tilt":True,"source_eta_powered_once":True,"counts_outside_alpha":True,
        "output_factor_once":True,"both_offdiagonal_terminal_orders_charged":True,
        "no_extra_alpha_or_event_independence":True,"unrefined_cells_retained":True,
        "tests":tests,"local_coefficients":built["report"],"elapsed_seconds":time.monotonic()-started,
        "peak_process_bytes":peak(),"generated_bytes":0,"nonnegative_intersections":v.intersections,
        "pins":{"run_band":RUN_PIN,"slab_geometry":GEOMETRY_PIN,"counts":run.COUNTS_PIN,
            "actual":run.ACTUAL_PIN,"support":actual.SUPPORT_PIN,"normalized":support.NORMALIZED_PIN,
            "fitter":support.FITTER_PIN,"band":normalized.BAND_PIN,"base_geometry":normalized.GEOMETRY_PIN,
            "local":large.LOCAL_PIN,"anchor":large.ANCHOR_PIN,"join":large.JOIN_PIN,
            "ordinary":reference.ORDINARY_PIN,"normalizer":large.REGION_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
