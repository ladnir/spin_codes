"""Matched selected-integer ordinary/residue256 diagnostic at actual r3072.

Rebuild the frozen current local pair table at z1/10,alpha119/1280 and
the same31-marker actual3072 bank as the calibrated run/slab drivers.
The FLOAT fitter only proposes rational markers.  Both directed replays
use the SAME immutable rational local coefficients and chosen markers.
Optional expected hashes reject a changed table or bank before regional
products.  No marker calibration is performed specifically for this pair.

For each selected integer full count f, s=3072-2*f, normalize by its own
D=C(4096,f)C(4096-f,s)u^f v^s.  The ordinary numerator is H(u,v)^64 or
H(u,v)^63 J(u,v).  The second numerator is the positive cyclic256
coefficient upper bound at(f mod256,s mod256), with upward integer
quantization.  This is a RESIDUE projection, not exact coefficient
extraction.  Six doublings preserve J=M*T and127 terminal updates.

Each ordered finite chain receives the identical actual source/output
factor and the exact FULL RUN-CLASS count N3072(3072-f0-f1), outside
alpha.  That count upper-bounds a selected profile's count; it is not an
exact profile count.  Neither a selected positive nor negative result
covers the whole run class, source size, or code.  Both terminal orders
are reported separately, not added into a claimed coverage bound.

At most two regional tables. Parent execution only:300/600s,512MiB,
20GiB free,-B,stdout,no generated files. No work/provider import on import.
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
RUN_PIN="19914e8bc698c410971f8edb0ab7f884cb0b83955bee2a6933e9791aa9dfcfc2"
WEIGHT,POINTS,MEMORY=3072,8192,512*(1<<20)
Z,ALPHA=F(1,10),F(119,1280)


def load_run():
    path=Path(__file__).with_name("cheap_lane_calibrated_run_band.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=RUN_PIN:
        raise RuntimeError("frozen calibrated run-band mismatch")
    spec=importlib.util.spec_from_file_location("selected_region_run",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def sha256_text(text):
    if len(text)!=64 or any(c not in "0123456789abcdefABCDEF" for c in text):
        raise argparse.ArgumentTypeError("expected64 hexadecimal SHA256 characters")
    return text.lower()


def validate_pair(left,right):
    if any(type(f) is not int or not 0<=f<=WEIGHT//2 for f in (left,right)):
        raise ValueError("full counts must be integers0..1536")
    runs=WEIGHT-left-right
    if not 1<=runs<=WEIGHT:
        raise ValueError("selected profile must have positive cyclic run count")
    return runs


def ordinary_region(v,coefficients,marker,full,ordinary,extraction,guard):
    h,j=ordinary.evaluate_pairs(v,coefficients,marker,guard)
    rr=v.matpow(h,64); tt=v.matmul(v.matpow(h,63),j)
    denominator=v.exact(extraction.target_denominator(WEIGHT,full,marker))
    result=tuple(tuple(v.positive(x/denominator) for row in matrix for x in row)
                 for matrix in (rr,tt))
    if tuple(map(len,result))!=(4,2):
        raise ArithmeticError("ordinary regional shape changed")
    return result


def own_tests(extraction,join,anchor,guard):
    """Exact domain/normalizer/count checks, with no selected-profile claim."""
    checks=0; marker=F(2,3),F(5,7)
    for f in (0,1,1282,1510,1536):
        guard(); s=WEIGHT-2*f; e=4096-WEIGHT+f
        assert min(f,s,e)>=0 and f+s+e==4096; checks+=1
        denominator=extraction.target_denominator(WEIGHT,f,marker)
        mass=F(comb(4096,f)*comb(4096-f,s))
        assert denominator==mass*marker[0]**f*marker[1]**s; checks+=1
    assert validate_pair(1282,1282)==508; checks+=1
    assert validate_pair(1536,0)==1536; checks+=1
    for bad in ((1536,1536),(-1,1),(1537,0),(True,1)):
        try: validate_pair(*bad)
        except ValueError: checks+=1
        else: raise AssertionError("invalid selected matching pair accepted")
    h,c_h,eta=anchor.source_factor(WEIGHT)
    assert h==12 and c_h==F(1,13); checks+=1
    assert eta==(1+F(sum(comb(WEIGHT,j) for j in range(1,13)),255))/(255*c_h); checks+=1
    counts=join.run_counts(points=POINTS,weight=WEIGHT)
    assert sum(counts)==comb(POINTS,WEIGHT); checks+=1
    for c in (1,52,508,1536,3072):
        guard(); numerator=POINTS*comb(WEIGHT-1,c-1)*comb(POINTS-WEIGHT-1,c-1)
        assert numerator%c==0 and counts[c]==numerator//c; checks+=1
    # Exactly128 factors from each matching, even with the output terminal.
    assert 127*2+1+1==256 and 64*256-1==16383; checks+=1
    assert 127+1==128; checks+=1
    return counts,{"checks":checks,"actual_source_h":h,"actual_source_c_h":str(c_h),
        "run_count_sum_sha256":hashlib.sha256(str(sum(counts)).encode("ascii")).hexdigest(),
        "indexed_run_count_not_profile_count":True,"own_normalizers":True,
        "finite_chain_has256_regional_factors":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left-center",type=int,default=813)
    parser.add_argument("--right-center",type=int,default=980)
    parser.add_argument("--left-full",type=int,default=1186)
    parser.add_argument("--right-full",type=int,default=1452)
    parser.add_argument("--bits",type=int,choices=(192,256),default=192)
    parser.add_argument("--precision",type=int,choices=(256,384),default=384)
    parser.add_argument("--seconds",type=int,choices=(300,600),default=300)
    parser.add_argument("--expect-table-sha",type=sha256_text)
    parser.add_argument("--expect-bank-sha",type=sha256_text)
    args=parser.parse_args()
    runs=validate_pair(args.left_full,args.right_full)
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    run=load_run(); actual=run.load_actual(); support=actual.load_support()
    normalized=support.load_normalized(); band=normalized.load_band()
    fitter=band.load(support.FITTER_FILE,support.FITTER_PIN,"selected_region_fitter")
    reference=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"selected_region_reference")
    ordinary=reference.load("cheap_lane_ordered_ordinary_directed.py",reference.ORDINARY_PIN,
                            "selected_region_ordinary")
    large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"selected_region_source")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"selected_region_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"selected_region_extraction")
    integer=large.load("cheap_lane_integer_regions.py",large.INTEGER_PIN,"selected_region_integer")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"selected_region_join")
    labels={row[0] for row in reference.BANK}
    if args.left_center not in labels or args.right_center not in labels:
        raise ValueError("centers must be labels in the frozen normalized31-marker bank")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("selected-region deadline; no incomplete result")
        if peak()>MEMORY: raise MemoryError("selected-region512MiB budget exceeded")
    def progress(phase,**values):
        keys=("center","status","evaluations","doubling","paired_factors","first_full","last_full")
        selected={key:values[key] for key in keys if key in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**selected}),flush=True)
    counts,own=own_tests(extraction,join,anchor,guard)
    tests={"own":own,"extraction":extraction.exact_tests(guard),
        "integer":integer.exact_tests(guard),"ordinary":ordinary.exact_tests(guard),
        "ordered_terminal":join.exact_tests(guard),
        "bank":support.wrapper_tests(reference.BANK,guard)}
    built=local.build_local_pairs(Z,ALPHA,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]; logtwo=v.arb(2).log()
    degrees=built["rational_coefficients"]; arb_identity=dict(built["coefficients"])
    table_hash=support.coefficient_fingerprint(degrees,guard)
    if table_hash!=built["report"]["rational_coefficient_sha256"]:
        raise ArithmeticError("local-table canonical hash mismatch")
    if args.expect_table_sha and table_hash!=args.expect_table_sha:
        raise ArithmeticError("local table differs from the requested baseline receipt")
    proposed,fit_report=actual.fit_actual_markers(MappingProxyType(degrees),WEIGHT,fitter,guard,progress)
    fitted,public_bank,bank_hash=support.validate_bank(proposed,reference.BANK)
    if args.expect_bank_sha and bank_hash!=args.expect_bank_sha:
        raise ArithmeticError("regenerated marker bank differs from the requested baseline receipt")
    if support.coefficient_fingerprint(degrees,guard)!=table_hash:
        raise ArithmeticError("fitter changed the immutable local rational table")
    if len(arb_identity)!=len(built["coefficients"]) or any(
            row is not built["coefficients"][key] for key,row in arb_identity.items()):
        raise ArithmeticError("fitter changed the original Arb table")
    # Reconstruct ordinary inputs from the exact table used by the integer
    # engine, avoiding any dependence on FLOAT or a differently rounded cap.
    coefficients={degree:tuple(v.exact(x) for x in row) for degree,row in degrees.items()}
    cells={row[0]:row for row in normalized.normalized_cells(fitted)}
    selected=((args.left_center,args.left_full),(args.right_center,args.right_full))
    ordinary_regions={}; residue_regions={}; tables=[]; target_reports=[]
    for center,full in dict.fromkeys(selected):
        guard(); cell=cells[center]; marker=cell[3:5]
        ordinary_regions[center,full]=ordinary_region(v,coefficients,marker,full,ordinary,extraction,guard)
        x=F(2048*full,WEIGHT)
        assigned=next(row[0] for index,row in enumerate(cells.values())
            if (row[1]<=x if index==0 else row[1]<x) and x<=row[2])
        target_reports.append({"center":center,"full":full,"single":WEIGHT-2*full,
            "empty":4096-WEIGHT+full,"normalized_full":str(x),"assigned_center":assigned,
            "uses_its_assigned_cell_marker":center==assigned,"u":str(marker[0]),"v":str(marker[1]),
            "full_residue":full%256,"single_residue":(WEIGHT-2*full)%256,
            "normalizer_scope":"own actual multinomial times own marker powers"})
    for center in dict.fromkeys(center for center,_ in selected):
        guard(); fulls=tuple(full for label,full in selected if label==center)
        item=(center,min(fulls),max(fulls),"cyclic")
        partial,report=extraction.one_table_many(v,degrees,item,cells[center][3:5],(WEIGHT,),
            integer,args.bits,guard,progress)
        for full in fulls: residue_regions[center,full]=partial[WEIGHT][full]
        tables.append(report); del partial; gc.collect(); guard()
    if not 1<=len(tables)<=2: raise ArithmeticError("selected adapter exceeded two regional tables")
    left,right=selected; diagonal=left==right
    ordinary_chains=join.two_orientations(ordinary_regions[left],ordinary_regions[right],diagonal=diagonal)
    residue_chains=join.two_orientations(residue_regions[left],residue_regions[right],diagonal=diagonal)
    h,c_h,eta=anchor.source_factor(WEIGHT)
    source=(v.exact(eta)*v.exact(core.BETA)**WEIGHT*v.exact(Z)**(-core.LIMIT))**v.exact(ALPHA)
    charge=source*v.exact(F(counts[runs]))
    orientations=[]
    orders=(("forward",left,right),) if diagonal else (("forward",left,right),("reverse",right,left))
    for index,(name,a,b) in enumerate(orders):
        ordinary_value=ordinary_chains[index]; residue_value=residue_chains[index]
        if not ordinary_value>0 or not residue_value>0:
            raise ArithmeticError("selected positive chain required for finite log diagnostics")
        orientations.append({"order":name,"centers":[a[0],b[0]],"actual_fulls":[a[1],b[1]],
            "actual_singles":[WEIGHT-2*a[1],WEIGHT-2*b[1]],"runs":runs,
            "log2_ordinary_conditional_interval":str(ordinary_value.log()/logtwo),
            "log2_residue_conditional_interval":str(residue_value.log()/logtwo),
            "log2_ordinary_with_full_run_mass_interval":str((charge*ordinary_value).log()/logtwo),
            "log2_residue_with_full_run_mass_interval":str((charge*residue_value).log()/logtwo),
            "ordinary_minus_residue_bits_interval":str((ordinary_value/residue_value).log()/logtwo),
            "not_a_run_or_size_cover":True})
    entry_controls=[]
    for target in dict.fromkeys(selected):
        pairs=zip((x for row in ordinary_regions[target] for x in row),
                  (x for row in residue_regions[target] for x in row))
        ratios=[str((a/b).log()/logtwo) if a>0 and b>0 else None for a,b in pairs]
        entry_controls.append({"center":target[0],"full":target[1],
            "ordinary_minus_residue_entry_bits":ratios})
    guard()
    print(json.dumps({"scope":"selected actual3072 integer profiles only; matched directed ordinary/residue256",
        "whole_code_certificate":False,"all_profiles_at_weight_covered":False,
        "run_class_covered":False,"directed_arithmetic":True,
        "parameters":{"weight":WEIGHT,"z":str(Z),"alpha":str(ALPHA),"integer_bits":args.bits,
            "requested_precision":args.precision,"working_precision":built["report"]["working_precision"],
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "selected_targets":target_reports,"ordered_comparison":orientations,
        "regional_entry_controls":entry_controls,"regional_tables":tables,
        "region_table_count":len(tables),"coefficient_projection_exact":False,
        "positive_aliases_and_upward_quantization_remain":True,
        "shared_bits_not_relative_coefficient_precision":True,
        "source_log2_interval":str(source.log()/logtwo),
        "run_count_log2_interval":str(v.exact(F(counts[runs])).log()/logtwo),
        "run_count_sha256":hashlib.sha256(str(counts[runs]).encode("ascii")).hexdigest(),
        "source":{"h":h,"c_h":str(c_h),"eta":str(eta)},
        "source_once_counts_outside_alpha":True,"same_run_mass_in_both_replays":True,
        "ordered_diagnostics_not_summed_as_coverage":True,"degree_reversals":0,
        "local_table_builds":1,"alpha_reapplied_after_local_table":False,
        "local_coefficient_sha256":table_hash,"positive_rational_bank_sha256":bank_hash,
        "expected_table_hash_checked":args.expect_table_sha is not None,
        "expected_bank_hash_checked":args.expect_bank_sha is not None,
        "positive_rational_bank":public_bank,
        "fit":{"scope":"same frozen31-marker proposal procedure; not proof arithmetic",
            "fits":fit_report["fits"],"successful_proposal_evaluations_recorded":
                fit_report["successful_proposal_evaluations_recorded"],
            "fit_failures":fit_report["fit_failures"]},
        "tests":tests,"local_coefficients":built["report"],
        "nonnegative_intersections":v.intersections,"elapsed_seconds":time.monotonic()-started,
        "peak_process_bytes":peak(),"generated_bytes":0,
        "pins":{"run":RUN_PIN,"actual":run.ACTUAL_PIN,"support":actual.SUPPORT_PIN,
            "fitter":support.FITTER_PIN,"normalized":support.NORMALIZED_PIN,"band":normalized.BAND_PIN,
            "reference":band.HOLE_PIN,"ordinary":reference.ORDINARY_PIN,"local":large.LOCAL_PIN,
            "source":large.ANCHOR_PIN,"extraction":large.REGION_PIN,"integer":large.INTEGER_PIN,
            "join":large.JOIN_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
