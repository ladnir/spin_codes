"""One-marker directed diagnostic for a universal positive residue maximum.

Rebuild the frozen actual3072 local table and31-marker bank at z1/10,
alpha119/1280, enforcing their recorded hashes. The FLOAT fitter proposes
the same positive rational bank only; all regional inputs are exact caps.
Mark every local degree before positive cyclic256 reduction. Six integer
updates T=R*T, R=R*R give H^64 and H^63 J with separate dyadic exponents.

For each of the six entries, scan all65536 legal residue positions and
take its exact nonnegative integer maximum. Every marked target coefficient
is bounded by its positive residue alias, and therefore by this maximum.
The six maxima need not occur at the same residue. Positivity permits
entrywise majorization; this does not assert independence or a probability.

Return compact exact dyadics (mantissa, exponent), plus directed maximum
versus DC diagnostics. DC controls use the SAME rational local table before
regional rounding. No own-target denominator, source factor, support count,
or further alpha is applied to the maximum numerator. An optional selected
target check uses its own denominator separately and reports six ratios;
it does not change the exported universal numerator or cover any source size.

Parent execution only:300/600s,512MiB,20GiB free,-B,stdout only. No work or
proof-provider imports occur on import. Numerical receipts are separate;
source alone is not evidence of a numerical bound.
"""

import argparse
from collections import defaultdict
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
COMPLETE_PIN="e599c0cd7ed7fd4f594f7b045aeb3702e9ded8fe0d2128ff0fd7dd6f47417b53"
SELECTED_PIN="cf4d025277a281aecfd00b105bd7f791ca00c9864838b4849f97d1da7a2de26b"
REGIONS_PIN="ab74014cbd1e2e6334ab15acaf2fcbc9dc6d9bc97b19e5e6c495f0d2eb8a5852"
INTEGER_PIN="9dc583ae54caed0a78f71e409cc37aa398713b419f65538370f10c3f876f40bf"
TABLE_SHA="636b295acf8ba93a20f76a983ab4723c290f9967a252a04b9562148a3a7c539d"
BANK_SHA="3e4843f4abd6be8730f72a5f77c02f4b4e554bf32a834464bda3df8643f551e7"
Z,ALPHA=F(1,10),F(119,1280)
WEIGHT,MODULUS,PRECISION,MEMORY=3072,256,384,512*(1<<20)


def load_complete():
    path=Path(__file__).with_name("cheap_lane_complete3072_directed.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=COMPLETE_PIN:
        raise RuntimeError("frozen complete3072 source mismatch")
    spec=importlib.util.spec_from_file_location("maximum_complete3072",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    if (module.SELECTED_PIN,module.TABLE_SHA,module.BANK_SHA)!=(SELECTED_PIN,TABLE_SHA,BANK_SHA):
        raise ArithmeticError("frozen witness provenance changed")
    return module


def scan_maxima(arithmetic,matrix,guard):
    """Exact legal-position maxima, with each matrix's own shared exponent.

    Return (matrix-shaped Fraction maxima, flat compact entry reports).
    Padding positions are checked but never included in a maximum or DC.
    """
    arithmetic._check(matrix)
    full,single=arithmetic.full,arithmetic.single
    if not arithmetic.cyclic_full or not arithmetic.cyclic_single:
        raise ValueError("universal residue scan requires two cyclic axes")
    stride=arithmetic.stride
    scale=F(1<<matrix.exponent) if matrix.exponent>=0 else F(1,1<<(-matrix.exponent))
    maxima=[]; reports=[]
    for i,row in enumerate(matrix.entries):
        result=[]
        for j,poly in enumerate(row):
            guard()
            if len(poly)>arithmetic.length:
                raise ArithmeticError("regional polynomial exceeds retained layout")
            best=mass=ties=0; position=(0,0); padding=0
            for s in range(single+1):
                if s%8==0: guard()
                base=stride*s
                for f in range(full+1):
                    value=int(poly[base+f])
                    if value<0: raise ArithmeticError("negative exact residue coefficient")
                    mass+=value
                    if value>best: best=value; position=(f,s); ties=1
                    elif value==best: ties+=1
                for f in range(full+1,stride):
                    if int(poly[base+f])!=0:
                        raise ArithmeticError("nonzero retained padding coefficient")
                    padding+=1
            if not 0<=best<=mass:
                raise ArithmeticError("residue maximum exceeds exact represented DC")
            if mass>(full+1)*(single+1)*best:
                raise ArithmeticError("residue average exceeds maximum")
            result.append(best*scale)
            reports.append({"row":i,"column":j,"mantissa":str(best),
                "exponent":matrix.exponent,"argmax_full_single":list(position),
                "maximizing_residues":ties,"DC_mantissa":str(mass),
                "legal_coefficients":(full+1)*(single+1),"zero_padding_checked":padding})
        maxima.append(tuple(result))
    return tuple(maxima),reports


def build_maximum(v,degrees,marker,integer,bits,guard,progress=None,selected=None,extraction=None):
    """Build one table; return unnormalized H64/H63J maxima and metadata."""
    if bits not in (192,256) or type(bits) is not int:
        raise ValueError("mantissa bits must be192 or256")
    if len(marker)!=2 or any(isinstance(x,bool) or not isinstance(x,(int,F)) or x<=0 for x in marker):
        raise ValueError("two positive exact rational markers required")
    if hashlib.sha256(Path(integer.__file__).read_bytes()).hexdigest()!=INTEGER_PIN:
        raise RuntimeError("integer arithmetic pin mismatch")
    u,w=map(F,marker); up=[u**f for f in range(65)]; wp=[w**s for s in range(65)]
    scaled={}
    for (f,s),row in degrees.items():
        guard()
        if type(f) is not int or type(s) is not int or min(f,s)<0 or f+s>64 or len(row)!=6:
            raise ValueError("six-channel pair degree triangle required")
        if any(isinstance(x,bool) or not isinstance(x,(int,F)) or x<0 for x in row):
            raise ValueError("exact nonnegative rational caps required")
        scaled[f,s]=tuple(up[f]*wp[s]*x for x in row)
    if len(scaled)!=2145: raise ValueError("complete degree64 local triangle required")
    plain=[v.exact(sum((row[k] for row in scaled.values()),F(0))) for k in range(6)]
    h=[[plain[0],plain[1]],[plain[2],plain[3]]]; j=[[plain[4]],[plain[5]]]
    evaluated=(v.matpow(h,64),v.matmul(v.matpow(h,63),j))
    arithmetic=integer.IntegerArithmetic(255,255,cyclic_full=True,cyclic_single=True,
        bits=bits,guard=guard)
    rr,tt=arithmetic.pack(scaled)
    del scaled,up,wp,plain,h,j
    emit=(lambda phase,**kw:None) if progress is None else progress
    emit("maximum_table_packed",matrix_exponent=rr.exponent,terminal_exponent=tt.exponent)
    for step in range(6):
        guard(); tt=arithmetic.multiply(rr,tt); rr=arithmetic.multiply(rr,rr)
        emit("maximum_pair_doubling_complete",doubling=step+1,paired_factors=1<<(step+1),
            matrix_exponent=rr.exponent,terminal_exponent=tt.exponent)
    maxima=[]; entries=[]; dc_checks=0; logtwo=v.arb(2).log()
    for name,matrix,ordinary in zip(("H64","H63J"),(rr,tt),evaluated):
        values,report=scan_maxima(arithmetic,matrix,guard); maxima.append(values)
        masses=arithmetic.mass(matrix)
        for item in report:
            guard(); i,j=item["row"],item["column"]
            raw=values[i][j]; dc=masses[i][j]; old=ordinary[i][j]
            exact_scale=integer.dyadic(matrix.exponent)
            if F(int(item["DC_mantissa"]))*exact_scale!=dc:
                raise ArithmeticError("scan DC differs from frozen arithmetic.mass")
            if F(int(item["mantissa"]))*exact_scale!=raw:
                raise ArithmeticError("maximum dyadic representation changed")
            directed_dc=v.exact(dc); directed_max=v.exact(raw)
            if not (directed_dc>=old or directed_dc.overlaps(old)):
                raise ArithmeticError("rounded cyclic DC is below unrounded ordinary mass")
            dc_checks+=1
            entries.append({"object":name,**item,
                "log2_maximum_interval":str(directed_max.log()/logtwo) if raw else None,
                "log2_rounded_DC_over_maximum":str(v.exact(dc/raw).log()/logtwo) if raw else None,
                "log2_ordinary_over_maximum_interval":str((old/directed_max).log()/logtwo)
                    if raw and old>0 else None,
                "log2_rounded_DC_over_ordinary_interval":str((directed_dc/old).log()/logtwo)
                    if dc and old>0 else None})
    selected_report=None
    if selected is not None:
        if extraction is None or type(selected) is not int or not 0<=selected<=WEIGHT//2:
            raise ValueError("selected own-target check requires a valid actual3072 full count")
        f=selected; s=WEIGHT-2*f
        denominator=extraction.target_denominator(WEIGHT,f,(u,w))
        controls=[]
        for name,matrix,maximum in zip(("H64","H63J"),(rr,tt),maxima):
            own=arithmetic.extract(matrix,f%MODULUS,s%MODULUS)
            for i,row in enumerate(own):
                for j,value in enumerate(row):
                    guard()
                    exact=int(matrix.entries[i][j][f%MODULUS+arithmetic.stride*(s%MODULUS)])*integer.dyadic(matrix.exponent)
                    if value!=exact or not 0<=value<=maximum[i][j]:
                        raise ArithmeticError("selected residue disagrees with exact scan")
                    normalized=v.exact(value/denominator)
                    controls.append({"object":name,"row":i,"column":j,
                        "log2_selected_own_normalized_interval":str(normalized.log()/logtwo) if value else None,
                        "log2_maximum_over_selected_interval":str(v.exact(maximum[i][j]/value).log()/logtwo)
                            if value else None})
        selected_report={"actual_weight":WEIGHT,"full":f,"single":s,
            "full_residue":f%MODULUS,"single_residue":s%MODULUS,
            "own_denominator_used_only_in_this_diagnostic":True,"same_live_table_no_recomputation":True,
            "entries":controls}
    result=(tuple(x for row in maxima[0] for x in row),tuple(x for row in maxima[1] for x in row))
    report={"entries":entries,"selected_control":selected_report,"DC_upper_checks":dc_checks,
        "pair_doublings":6,"regular_object":"H^64","terminal_object":"H^63 J",
        "regular_exponent":rr.exponent,"terminal_exponent":tt.exponent,
        "separate_matrix_exponents":True,"full_residue":256,"single_residue":256,
        "arithmetic":arithmetic.report(),"maxima_need_not_share_a_residue":True,
        "no_target_normalizer_on_export":True,"no_alpha_source_or_count_on_export":True}
    del rr,tt,evaluated,arithmetic
    gc.collect(); guard()
    return result,report


def exact_tests(integer,guard):
    """Exact small scans, alias domination, padding and scale separation."""
    checks=0
    arith=integer.IntegerArithmetic(2,2,True,True,bits=8,guard=guard)
    maps=({(0,0):F(1,3),(2,1):F(2,5)}, {(1,0):F(7,11)},
          {(0,2):F(3,13)}, {(0,0):F(5,17),(1,1):F(1,19)},
          {(0,1):F(1,1<<80)}, {(2,0):F(3,1<<80)})
    rr=arith.pack_matrix(((maps[0],maps[1]),(maps[2],maps[3])))
    tt=arith.pack_matrix(((maps[4],),(maps[5],)))
    assert rr.exponent!=tt.exponent; checks+=1
    def exact_product(left,right):
        out=[]
        for row in left:
            target=[]
            for j in range(len(right[0])):
                poly=defaultdict(F)
                for k,p in enumerate(row):
                    for (f,s),x in p.items():
                        for (ff,ss),y in right[k][j].items(): poly[f+ff,s+ss]+=x*y
                target.append(dict(poly))
            out.append(tuple(target))
        return tuple(out)
    exact_h=exact_product(((maps[0],maps[1]),(maps[2],maps[3])),
                          ((maps[0],maps[1]),(maps[2],maps[3])))
    exact_t=exact_product(((maps[0],maps[1]),(maps[2],maps[3])),((maps[4],),(maps[5],)))
    old_r,old_t=rr,tt; tt=arith.multiply(old_r,old_t); rr=arith.multiply(old_r,old_r)
    assert rr.exponent!=tt.exponent; checks+=1
    for matrix,truth in ((rr,exact_h),(tt,exact_t)):
        maxima,reports=scan_maxima(arith,matrix,guard)
        mass=arith.mass(matrix)
        for i,row in enumerate(matrix.entries):
            for j,_ in enumerate(row):
                extracted=[arith.extract(matrix,f,s)[i][j] for s in range(3) for f in range(3)]
                assert maxima[i][j]==max(extracted); checks+=1
                assert mass[i][j]==sum(extracted,F(0)); checks+=1
                for (f,s),value in truth[i][j].items():
                    alias=arith.extract(matrix,f%3,s%3)[i][j]
                    assert maxima[i][j]>=alias>=value; checks+=1
                item=reports[i*matrix.shape[1]+j]
                assert maxima[i][j]==int(item["mantissa"])*integer.dyadic(matrix.exponent); checks+=1
                assert item["legal_coefficients"]==9 and item["zero_padding_checked"]==6; checks+=1
    zero=arith.pack_matrix((({},),))
    values,report=scan_maxima(arith,zero,guard)
    assert values==((F(0),),) and report[0]["mantissa"]=="0"; checks+=1
    assert report[0]["maximizing_residues"]==9; checks+=1
    # A padding location must never silently compete with legal maxima.
    for location,value in ((3,100), (0,-1)):
        raw=[0]*arith.length; raw[location]=value
        invalid=integer.DyadicMatrix(((arith.poly(raw),),),0,arith.layout)
        try: scan_maxima(arith,invalid,guard)
        except ArithmeticError: checks+=1
        else: raise AssertionError("invalid retained coefficient accepted")
    # Flattening cannot replace two object scales by one shared exponent.
    a=integer.DyadicMatrix(((arith.poly([3]),),),-5,arith.layout)
    b=integer.DyadicMatrix(((arith.poly([3]),),),7,arith.layout)
    aa,_=scan_maxima(arith,a,guard); bb,_=scan_maxima(arith,b,guard)
    assert aa[0][0]==F(3,32) and bb[0][0]==384; checks+=1
    assert bb[0][0]/aa[0][0]==4096; checks+=1
    return {"exact_checks":checks,"scan_equals_all_frozen_extractions":True,
        "positive_alias_domination":True,"distinct_regular_terminal_exponents":True,
        "zero_and_invalid_padding_controls":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--center",type=int,choices=(813,980),default=813)
    parser.add_argument("--bits",type=int,choices=(192,256),default=192)
    parser.add_argument("--seconds",type=int,choices=(300,600),default=300)
    parser.add_argument("--selected-control",action="store_true")
    args=parser.parse_args()
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    complete=load_complete(); selected=complete.load_selected()
    if (selected.WEIGHT,selected.Z,selected.ALPHA)!=(WEIGHT,Z,ALPHA):
        raise ArithmeticError("frozen selected witness parameters changed")
    run=selected.load_run(); actual=run.load_actual(); support=actual.load_support()
    normalized=support.load_normalized(); band=normalized.load_band()
    fitter=band.load(support.FITTER_FILE,support.FITTER_PIN,"maximum_fitter")
    reference=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"maximum_reference")
    ordinary=reference.load("cheap_lane_ordered_ordinary_directed.py",reference.ORDINARY_PIN,"maximum_ordinary")
    large=ordinary.load_large()
    if (large.REGION_PIN,large.INTEGER_PIN)!=(REGIONS_PIN,INTEGER_PIN):
        raise ArithmeticError("regional arithmetic provenance changed")
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"maximum_resources")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"maximum_local")
    extraction=large.load("cheap_lane_band_regions.py",REGIONS_PIN,"maximum_extraction")
    integer=large.load("cheap_lane_integer_regions.py",INTEGER_PIN,"maximum_integer")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("residue maximum deadline; no incomplete result")
        if peak()>MEMORY: raise MemoryError("residue maximum512MiB limit exceeded")
    def progress(phase,**values):
        if phase=="actual_marker_fit" and values.get("center") not in (0,512,1024): return
        keys=("center","status","evaluations","doubling","paired_factors","matrix_exponent","terminal_exponent")
        kept={key:values[key] for key in keys if key in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**kept}),flush=True)
    tests={"maxima":exact_tests(integer,guard),"integer":integer.exact_tests(guard),
        "extraction":extraction.exact_tests(guard),"ordinary":ordinary.exact_tests(guard),
        "bank":support.wrapper_tests(reference.BANK,guard)}
    built=local.build_local_pairs(Z,ALPHA,anchor,guard,deadline,PRECISION,progress=progress)
    v=built["verifier"]; degrees=built["rational_coefficients"]
    arb_identity=dict(built["coefficients"])
    table_hash=support.coefficient_fingerprint(degrees,guard)
    if table_hash!=TABLE_SHA or built["report"]["rational_coefficient_sha256"]!=TABLE_SHA:
        raise ArithmeticError("coherent local table differs from recorded witness")
    proposed,fit_report=actual.fit_actual_markers(MappingProxyType(degrees),WEIGHT,fitter,guard,progress)
    fitted,public_bank,bank_hash=support.validate_bank(proposed,reference.BANK)
    if bank_hash!=BANK_SHA: raise ArithmeticError("rational bank differs from recorded witness")
    if support.coefficient_fingerprint(degrees,guard)!=TABLE_SHA or len(arb_identity)!=len(built["coefficients"]):
        raise ArithmeticError("fitter changed local coefficients")
    if any(row is not built["coefficients"][key] for key,row in arb_identity.items()):
        raise ArithmeticError("fitter changed original Arb rows")
    cells=normalized.normalized_cells(fitted)
    cell=next(row for row in cells if row[0]==args.center); marker=cell[3:5]
    target={813:1186,980:1452}[args.center] if args.selected_control else None
    if target is not None:
        x=F(2048*target,WEIGHT)
        assigned=next(row[0] for index,row in enumerate(cells)
            if (row[1]<=x if index==0 else row[1]<x) and x<=row[2])
        if assigned!=args.center: raise ArithmeticError("selected control no longer uses recorded marker")
    maxima,report=build_maximum(v,degrees,marker,integer,args.bits,guard,progress,target,extraction)
    if tuple(map(len,maxima))!=(4,2): raise ArithmeticError("six-channel maximum interface changed")
    guard()
    print(json.dumps({"scope":"one-marker universal residue-max numerator diagnostic only",
        "directed_arithmetic":True,"whole_code_certificate":False,"source_size_covered":False,
        "parameters":{"center":args.center,"u":str(marker[0]),"v":str(marker[1]),
            "z":str(Z),"alpha":str(ALPHA),"local_witness_anchor_actual":WEIGHT,
            "integer_bits":args.bits,"working_precision":built["report"]["working_precision"],
            "seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "universal_numerator":report,"exact_dyadic_convention":"value = mantissa * 2^exponent",
        "source_applied":False,"support_count_applied":False,"target_normalization_on_maxima":False,
        "alpha_reapplied_after_local_table":False,"positive_aliases_and_upward_rounding_remain":True,
        "regional_tables_built":1,"degree_reversals":0,"local_table_builds":1,
        "local_coefficient_sha256":table_hash,"positive_rational_bank_sha256":bank_hash,
        "baseline_hashes_enforced":True,"positive_rational_bank":public_bank,
        "fit":{"scope":"frozen FLOAT proposal only; exact rational witness replay",
            "successful_proposal_evaluations_recorded":fit_report["successful_proposal_evaluations_recorded"],
            "fit_failures":fit_report["fit_failures"]},
        "tests":tests,"local_coefficients":built["report"],"generated_bytes":0,
        "nonnegative_intersections":v.intersections,"elapsed_seconds":time.monotonic()-started,
        "peak_process_bytes":peak(),"pins":{"complete3072":COMPLETE_PIN,"selected":SELECTED_PIN,
            "run":selected.RUN_PIN,"actual":run.ACTUAL_PIN,"support":actual.SUPPORT_PIN,
            "fitter":support.FITTER_PIN,"normalized":support.NORMALIZED_PIN,"band":normalized.BAND_PIN,
            "reference":band.HOLE_PIN,"ordinary":reference.ORDINARY_PIN,"local":large.LOCAL_PIN,
            "resource_anchor":large.ANCHOR_PIN,"regions":REGIONS_PIN,"integer":INTEGER_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2),flush=True)


if __name__=="__main__":
    main()
