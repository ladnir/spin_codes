"""Directed normalized hole-band bound with a fitted rational marker bank.

Build the current directed local pair table ONCE at the selected exact
z/alpha. A read-only view of its actual-degree rational coefficients is
passed to the bounded FLOAT fitter. Fitting only proposes31 positive
rational witnesses at one actual-size anchor. No FLOAT value enters the
final bound: the unchanged Arb coefficients are reversed once, evaluated
at those rational markers, and joined through the frozen normalized-cell
geometry, directed log denominators and actual-r source-count envelope.

The anchor may differ from the requested band; it selects proof parameters,
not a setup distribution. Every integer actual size in the inclusive band,
including odd sizes, receives its own source factor. Counts stay outside
alpha. No endpoint replacement, marker union factor, new coefficient power
after fitting, or uniform full-field law is introduced. A poor fit may
produce a positive bound; it cannot invalidate a completed directed replay.

Parent execution only:300s/512MiB/20GiB free,-B,stdout,no generated files.
No work runs on import. Numerical receipts are indexed separately.
"""

import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import time
from types import MappingProxyType

sys.dont_write_bytecode=True
NORMALIZED_PIN="8686ec42dd814974efeb3ec9230d9f4da898b1df83712fec38ed1acce2b7983e"
FITTER_PIN="d926a1283b4fb016687e6d5fdc66ce59d031db79a465cc015c3a81888ca8ec0d"
FITTER_FILE="cheap_lane_hole_marker_fit.py"
MEMORY=512*(1<<20)


def load_normalized():
    path=Path(__file__).with_name("cheap_lane_normalized_hole_band.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=NORMALIZED_PIN:
        raise RuntimeError("frozen normalized hole-band mismatch")
    spec=importlib.util.spec_from_file_location("calibrated_normalized_band",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def rational64(text):
    """Bound input size before parsing and preserve the exact CLI witness."""
    if not isinstance(text,str) or len(text)>128:
        raise argparse.ArgumentTypeError("rational witness text is too long")
    if "e" in text.lower():
        raise argparse.ArgumentTypeError("use a bounded fraction or decimal, not exponent notation")
    try: value=F(text)
    except (ValueError,ZeroDivisionError) as exc:
        raise argparse.ArgumentTypeError("expected an exact rational") from exc
    if max(abs(value.numerator).bit_length(),value.denominator.bit_length())>64:
        raise argparse.ArgumentTypeError("rational witnesses limited to64 bits")
    return value


def coefficient_fingerprint(table,guard):
    """Same canonical ordering/serialization as the frozen local constructor."""
    digest=hashlib.sha256(); degrees=0; entries=0
    for f in range(65):
        guard()
        for s in range(65-f):
            row=table[f,s]
            if not isinstance(row,tuple) or len(row)!=6:
                raise TypeError("local coefficient row must be an immutable six-tuple")
            if any(not isinstance(x,F) or x<0 for x in row):
                raise TypeError("local coefficients must be nonnegative exact Fractions")
            digest.update((str((f,s))+":"+":".join(map(str,row))+"\n").encode())
            degrees+=1; entries+=6
    if len(table)!=degrees or degrees!=2145 or entries!=12870:
        raise ArithmeticError("incomplete actual-degree pair table")
    return digest.hexdigest()


def validate_bank(bank,reference):
    """Only witness values may change; normalized labels/cuts stay frozen."""
    if not isinstance(bank,(list,tuple)) or len(bank)!=31 or len(reference)!=31:
        raise ValueError("expected31 fitted marker rows")
    rows=[]; public=[]
    for row,old in zip(bank,reference):
        if not isinstance(row,(list,tuple)) or len(row)!=5:
            raise TypeError("marker row must contain center,lo,hi,u,v")
        center,lo,hi,u,w=row
        if (center,lo,hi)!=old[:3] or any(type(x) is not int for x in (center,lo,hi)):
            raise ValueError("fitter changed the normalized cell labels")
        if any(isinstance(x,bool) or not isinstance(x,(int,F)) for x in (u,w)):
            raise TypeError("fitted markers must be exact rational values")
        u,w=F(u),F(w)
        if min(u,w)<=0: raise ValueError("fitted markers must be strictly positive")
        rows.append((center,lo,hi,u,w))
        public.append({"center":center,"reference_lo":lo,"reference_hi":hi,"u":str(u),"v":str(w)})
    if len(rows)!=31: raise ValueError("incomplete fitted bank")
    serialized=json.dumps(public,sort_keys=True,separators=(",",":"))
    return tuple(rows),public,hashlib.sha256(serialized.encode("ascii")).hexdigest()


def wrapper_tests(reference,guard):
    checks=0
    for text,want in (("1/32",F(1,32)),("119/1280",F(119,1280)),("0.125",F(1,8))):
        assert rational64(text)==want; checks+=1
    for text in ("1/0","no","1e999999999",str(1<<65),"1/"+str(1<<65),"1"*129):
        try: rational64(text)
        except argparse.ArgumentTypeError: checks+=1
        else: raise AssertionError("unbounded or invalid rational accepted")
    valid,public,digest=validate_bank(reference,reference)
    assert valid==reference and len(public)==31 and len(digest)==64; checks+=1
    for bad in (reference[:-1],reference+reference[:1]):
        try: validate_bank(bad,reference)
        except ValueError: checks+=1
        else: raise AssertionError("incomplete/extra marker rows accepted")
    row=list(reference[0])
    for value in (F(0),F(-1),1.0,True):
        guard(); row[3]=value; bad=(tuple(row),)+reference[1:]
        try: validate_bank(bad,reference)
        except (ValueError,TypeError): checks+=1
        else: raise AssertionError("invalid witness accepted")
    row=list(reference[0]); row[0]+=1
    try: validate_bank((tuple(row),)+reference[1:],reference)
    except ValueError: checks+=1
    else: raise AssertionError("changed proof cell accepted")
    toy={(0,0):(F(1),)*6}; readonly=MappingProxyType(toy)
    try: readonly[0,0]=(F(2),)*6
    except TypeError: checks+=1
    else: raise AssertionError("read-only coefficient view was mutable")
    return {"wrapper_exact_checks":checks,"positive_rational_bank_required":True,
        "fixed_normalized_cell_layout":True,"fit_input_mapping_readonly":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchor",type=int,default=5120)
    parser.add_argument("--first",type=int,default=4864)
    parser.add_argument("--last",type=int,default=5376)
    parser.add_argument("--z",type=rational64,default=F(1,32))
    parser.add_argument("--alpha",type=rational64,default=F(119,1280))
    parser.add_argument("--seconds",type=int,choices=(180,300),default=300)
    parser.add_argument("--precision",type=int,choices=(256,384),default=384)
    args=parser.parse_args()
    if not 4097<=args.first<=args.last<=7551 or not 4097<=args.anchor<=7551:
        raise ValueError("actual anchor and inclusive band must lie in4097..7551")
    if not 0<args.z<1 or not 0<args.alpha<=1:
        raise ValueError("require0<z<1 and0<alpha<=1")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    normalized=load_normalized(); band=normalized.load_band()
    fitter=band.load(FITTER_FILE,FITTER_PIN,"calibrated_hole_fitter")
    geometry=band.load(normalized.GEOMETRY_FILE,normalized.GEOMETRY_PIN,"calibrated_normalized_geometry")
    old_geometry=band.load(band.GEOMETRY_FILE,band.GEOMETRY_PIN,"calibrated_rectangle_tests")
    hole=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"calibrated_reference_bank")
    ordinary=hole.load("cheap_lane_ordered_ordinary_directed.py",hole.ORDINARY_PIN,"calibrated_ordinary")
    large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"calibrated_source_interface")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"calibrated_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"calibrated_denominator")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"calibrated_terminal")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("calibrated band deadline; no incomplete claim")
        if peak()>MEMORY: raise MemoryError("calibrated band512MiB budget exceeded")
    def progress(phase,**values):
        selected={key:values[key] for key in ("center","status","evaluations","row","vertices") if key in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**selected}),flush=True)
    tests={"wrapper":wrapper_tests(hole.BANK,guard),"fitter":fitter.exact_tests(guard),
        "frozen_band":band.exact_tests(hole.BANK,old_geometry,join,guard),
        "ordinary":ordinary.exact_tests(guard)}
    built=local.build_local_pairs(args.z,args.alpha,anchor,guard,deadline,args.precision,progress=progress)
    v,core=built["verifier"],built["core"]; logtwo=v.arb(2).log()
    degrees=built["rational_coefficients"]
    before=coefficient_fingerprint(degrees,guard)
    assert before==built["report"]["rational_coefficient_sha256"]
    arb_identity={key:row for key,row in built["coefficients"].items()}
    proposed,fit_report=fitter.fit_markers(MappingProxyType(degrees),args.anchor,guard,progress=progress)
    fitted,public_bank,bank_hash=validate_bank(proposed,hole.BANK)
    after=coefficient_fingerprint(degrees,guard)
    if before!=after or len(arb_identity)!=len(built["coefficients"]):
        raise ArithmeticError("fitter changed the shared coherent local table")
    assert all(row is built["coefficients"][key] for key,row in arb_identity.items())
    cells=normalized.normalized_cells(fitted)
    klo,khi=8192-args.last,8192-args.first
    cover=normalized.coverage(cells,klo,khi,guard)
    tests["normalized_new_bank"]=normalized.normalized_tests(cells,fitted,geometry,guard)
    coefficients={(64-f-s,s):row for (f,s),row in built["coefficients"].items()}
    assert len(coefficients)==2145
    assert all(row is built["coefficients"][64-f-s,s] for (f,s),row in coefficients.items())
    den=band.LogDenominators(v,core.BETA,args.alpha,guard)
    tests["directed_new_bank"]=band.directed_tests(v,den,extraction,fitted,args.alpha,core.BETA,join,guard)
    numerators,marker_logs=normalized.numerator_tables(v,coefficients,cells,ordinary,guard,progress)
    maximum,worst,geometry_report=normalized.cell_envelope(v,band,geometry,den,numerators,marker_logs,
        cells,klo,khi,join,guard,progress)
    weighted,source_report=band.outer_source_sum(v,anchor,args.alpha,args.first,args.last,guard)
    if not weighted>0: raise ArithmeticError("nonpositive weighted source sum")
    log_count=weighted.log(); log_output=-v.exact(args.alpha)*core.LIMIT*v.exact(args.z).log()
    log_total=log_count+log_output+maximum
    guard()
    print(json.dumps({"scope":"all indexed supports in the inclusive actual-size band and all histories at35%",
        "directed_arithmetic":True,"whole_code_certificate":False,"all_integer_sizes_in_band":True,
        "parameters":{"actual_first":args.first,"actual_last":args.last,"fit_anchor_actual":args.anchor,
            "fit_anchor_holes":8192-args.anchor,"anchor_inside_band":args.first<=args.anchor<=args.last,
            "hole_first":klo,"hole_last":khi,"z":str(args.z),"alpha":str(args.alpha),
            "working_precision":built["report"]["working_precision"],"seconds":args.seconds,
            "memory_budget_bytes":MEMORY},
        "log2_union_bound_interval":str(log_total/logtwo),
        "band_below_2_pow_minus_40":bool(log_total<-40*logtwo),
        "band_below_2_pow_minus_60":bool(log_total<-60*logtwo),
        "log2_eta_weighted_support_count_interval":str(log_count/logtwo),
        "log2_output_correction_interval":str(log_output/logtwo),
        "worst_cell":worst,"geometry":geometry_report,"coverage":cover,"source_counts":source_report,
        "cell_coordinate":"x=2048*hole_full/hole_weight","continuous_cut_rule":"reference hi+1/2; endpoints0,1024",
        "positive_rational_bank":public_bank,"positive_rational_bank_sha256":bank_hash,
        "fit":{"status":"FLOAT proposals only; exact rational witnesses independently used in directed proof",
            "not_a_certificate":True,"report":fit_report},
        "local_table_builds":1,"alpha_explicitly_rebuilt":True,"no_coefficient_power_after_fitting":True,
        "local_coefficient_sha256_before_fit":before,"local_coefficient_sha256_after_fit":after,
        "arb_coefficients_preserved_by_identity":True,"actual_degree_fit_input":True,
        "directed_degree_reversals":1,"source_beta_only_in_denominator_tilt":True,
        "output_factor_once":True,"boundary_replacement":False,"extra_profile_or_marker_factor":False,
        "tests":tests,"local_coefficients":built["report"],"elapsed_seconds":time.monotonic()-started,
        "peak_process_bytes":peak(),"generated_bytes":0,"nonnegative_intersections":v.intersections,
        "pins":{"normalized":NORMALIZED_PIN,"fitter":FITTER_PIN,"band":normalized.BAND_PIN,
            "geometry":normalized.GEOMETRY_PIN,"old_geometry":band.GEOMETRY_PIN,"hole":band.HOLE_PIN,
            "ordinary":hole.ORDINARY_PIN,"local":large.LOCAL_PIN,"anchor":large.ANCHOR_PIN,
            "join":large.JOIN_PIN,"normalizer":large.REGION_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
