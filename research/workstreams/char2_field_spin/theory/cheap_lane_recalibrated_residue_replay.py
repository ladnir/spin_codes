"""Replay a trusted dynamically calibrated residue-maximum witness.

Require the raw <=64KiB exact_maximum_witness and its canonical SHA256.
Read z and alpha as bounded canonical rational strings from that header.
The fitting anchor is a required CLI argument because it is not stored in
the witness. Rebuild actual-occupancy caps at those z,alpha and fit the
same31 markers at that anchor; check both fingerprints and every marker.
Restore186 exact dyadics with separate regular/terminal object scales.

Upper-bound validity of these saved values is inherited from the verified
generator receipt. A matching hash proves identity, not that arbitrary
values are coefficient bounds. No regional polynomial is reconstructed.
The requested complete actual band uses the dynamic generator's unchanged
denominator geometry, actual-r eta/counts, both terminal orders and output
factor. No new alpha, target denominator or source is put into the saved
unnormalized numerators. Source alone is not a numerical receipt.

Parent execution only:300/600s,512MiB,20GiB free,-B,stdout,no generated
files. No provider imports or work on import; no frozen globals changed.
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
GENERATOR_PIN="ed6745fb570200c27300a0da5e4c179ae4c998483c026592163f1428272fd619"
PARSER_PIN="77bd0acc081e68691810391c2cd4e9a5dad1f731f14296ec78d4de4831018867"
PRECISION,MEMORY=384,512*(1<<20)


def load(name,pin,label):
    path=Path(__file__).with_name(name)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pin:
        raise RuntimeError("frozen replay dependency mismatch: "+name)
    spec=importlib.util.spec_from_file_location(label,path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def canonical_fraction(text):
    if type(text) is not str or not 1<=len(text)<=128 or not text.isascii():
        raise ValueError("bounded canonical rational string required")
    # Reject exponent/decimal syntax BEFORE Fraction can allocate a large
    # power of ten. A nonnegative64-bit component needs at most20 digits.
    parts=text.split("/")
    if len(parts) not in (1,2) or any(not 1<=len(part)<=20 or not part.isdecimal()
            or (len(part)>1 and part[0]=="0") for part in parts):
        raise ValueError("canonical bounded integer ratio required")
    try: value=F(text)
    except (ValueError,ZeroDivisionError) as exc: raise ValueError("invalid rational header") from exc
    if max(abs(value.numerator).bit_length(),value.denominator.bit_length())>64 or str(value)!=text:
        raise ValueError("header rational must be canonical and at most64 bits")
    return value


def hash_text(value):
    return type(value) is str and len(value)==64 and all(c in "0123456789abcdef" for c in value)


def parameters(witness,parser):
    if type(witness) is not dict or set(witness)!=parser.HEADER_KEYS:
        raise ValueError("raw exact_maximum_witness object required")
    if witness["coordinate"]!="actual" or type(witness["moduli"]) is not list or witness["moduli"]!=[256,256]:
        raise ValueError("actual coordinates and cyclic256 axes required")
    if any(type(x) is not int for x in witness["moduli"]):
        raise ValueError("exact integer moduli required")
    bits=witness["integer_bits"]
    if type(bits) is not int or bits not in (192,256): raise ValueError("recorded bits must be192 or256")
    z,alpha=canonical_fraction(witness["z"]),canonical_fraction(witness["alpha"])
    if not 0<z<1 or not 0<alpha<=1: raise ValueError("require0<z<1 and0<alpha<=1")
    if not all(hash_text(witness[key]) for key in ("local_coefficient_sha256","positive_rational_bank_sha256")):
        raise ValueError("canonical coefficient/bank SHA256 strings required")
    if type(witness["rows"]) is not list or len(witness["rows"])!=31:
        raise ValueError("exactly31 ordered marker rows required")
    return z,alpha,bits


def validate_witness(witness,expected_sha,bank,table_hash,bank_hash,parser,guard):
    """Identity and structure only; cap validity needs the generator receipt."""
    z,alpha,bits=parameters(witness,parser)
    if len(parser.canonical_bytes(witness))>parser.MAX_BYTES or parser.witness_hash(witness)!=expected_sha:
        raise ValueError("canonical witness size/hash mismatch")
    if (witness["local_coefficient_sha256"],witness["positive_rational_bank_sha256"])!=(table_hash,bank_hash):
        raise ValueError("rebuilt local table or anchor-fitted bank differs from witness")
    if len(bank)!=31: raise ValueError("incomplete rebuilt bank")
    numerators=[]; exponent_pairs=[]; zeros=maximum_bits=0
    for row,reference in zip(witness["rows"],bank):
        guard()
        if type(row) is not list or len(row)!=5: raise ValueError("row must be[center,u,v,R4,V2]")
        center,u,w,regular,terminal=row
        if type(center) is not int or (center,u,w)!=(reference[0],str(reference[3]),str(reference[4])):
            raise ValueError("marker order/value differs from rebuilt positive rational bank")
        if type(u) is not str or type(w) is not str or max(len(u),len(w))>80:
            raise ValueError("bounded exact marker strings required")
        if type(regular) is not list or type(terminal) is not list or (len(regular),len(terminal))!=(4,2):
            raise ValueError("six-entry R4/V2 layout required")
        values=[]; exponents=[]
        for group in (regular,terminal):
            decoded=[parser.decode_dyadic(pair,bits) for pair in group]
            if len({exponent for _,exponent in decoded})!=1:
                raise ValueError("each matrix object must retain its own shared exponent")
            values.append(tuple(value for value,_ in decoded)); exponents.append(decoded[0][1])
            zeros+=sum(value==0 for value,_ in decoded)
            maximum_bits=max(maximum_bits,max(int(pair[0]).bit_length() for pair in group))
        numerators.append(tuple(values)); exponent_pairs.append([center,*exponents])
    return tuple(numerators),{"rows":31,"exact_dyadics":186,"integer_bits":bits,
        "z":str(z),"alpha":str(alpha),"channel_order":["R00","R01","R10","R11","V0","V1"],
        "maximum_mantissa_bits":maximum_bits,"zero_entries":zeros,
        "separate_object_exponents":exponent_pairs,"regenerated_marker_rows_checked":31,
        "trusted_generator_receipt_required":True,"upper_bound_values_independently_recomputed":False}


def exact_tests(bank,parser,guard):
    """Tiny stand-in values test metadata, not their upper-bound validity."""
    checks=0
    def clone(value): return json.loads(json.dumps(value))
    witness={"coordinate":"actual","moduli":[256,256],"integer_bits":192,
        "z":"1/12","alpha":"119/1280","local_coefficient_sha256":"1"*64,
        "positive_rational_bank_sha256":"2"*64,"rows":[
            [row[0],str(row[3]),str(row[4]),[[str(6*i+j),-11] for j in range(4)],
             [[str(6*i+j),7] for j in (4,5)]] for i,row in enumerate(bank)]}
    def validate(value,digest=None):
        return validate_witness(value,parser.witness_hash(value) if digest is None else digest,
            bank,"1"*64,"2"*64,parser,guard)
    digest=parser.witness_hash(witness); decoded,report=validate(witness,digest)
    assert report["exact_dyadics"]==186 and len(decoded)==31; checks+=1
    for i,(regular,terminal) in enumerate(decoded):
        guard(); assert regular==tuple(F(6*i+j,2048) for j in range(4)); checks+=1
        assert terminal==tuple(F((6*i+j)*128) for j in (4,5)); checks+=1
    assert parser.witness_hash(parser.decode_witness(json.dumps(witness,indent=2).encode("ascii")))==digest; checks+=1
    # Identity checks reject changed values even if their schema remains valid.
    for key,value in (("z","1/10"),("alpha","3/32")):
        bad=clone(witness); bad[key]=value
        try: validate(bad,digest)
        except ValueError: checks+=1
        else: raise AssertionError("changed canonical witness accepted under original hash")
    invalid=[]
    for key,value in (("coordinate","hole"),("moduli",[128,256]),("integer_bits",True),
            ("z","0.1"),("z","1/1"),("alpha","0"),("alpha","1e-2"),
            ("local_coefficient_sha256","3"*64),("positive_rational_bank_sha256","3"*64)):
        bad=clone(witness); bad[key]=value; invalid.append(bad)
    bad=clone(witness); bad["rows"].pop(); invalid.append(bad)
    bad=clone(witness); bad["rows"][0],bad["rows"][1]=bad["rows"][1],bad["rows"][0]; invalid.append(bad)
    bad=clone(witness); bad["rows"][0][1]="2"; invalid.append(bad)
    bad=clone(witness); bad["rows"][0][3].pop(); invalid.append(bad)
    bad=clone(witness); bad["rows"][0][4][1][1]=8; invalid.append(bad)
    for pair in (["-1",0],["01",0],[str(1<<192),0],["0",parser.MAX_EXPONENT+1],["0",False]):
        bad=clone(witness); bad["rows"][0][3][0]=pair; invalid.append(bad)
    bad=clone(witness); bad["extra"]=0; invalid.append(bad)
    for bad in invalid:
        try: validate(bad)
        except ValueError: checks+=1
        else: raise AssertionError("malformed metadata accepted under fresh identity hash")
    for data in (b'{"a":1,"a":2}',b'{"a":1.0}',b'{"a":NaN}',b'{"a":123456789}',b'[]',b' '*(parser.MAX_BYTES+1)):
        try: parser.decode_witness(data)
        except ValueError: checks+=1
        else: raise AssertionError("invalid bounded JSON accepted")
    for text in ("1"*129,str(1<<65),"1/0","01/12","0.25","1e99999999999999999999"):
        try: canonical_fraction(text)
        except ValueError: checks+=1
        else: raise AssertionError("noncanonical/unbounded header rational accepted")
    return {"metadata_exact_checks":checks,"standin_values_not_claimed_coefficient_caps":True,
        "dynamic_parameters_and_separate_scales_tested":True,"no_regional_products":True}


def main():
    generator=load("cheap_lane_recalibrated_residue_max_band.py",GENERATOR_PIN,"dynamic_replay_generator")
    parser=load("cheap_lane_residue_max_replay.py",PARSER_PIN,"dynamic_replay_parser")
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--witness",type=Path,required=True)
    cli.add_argument("--expect-witness-sha",type=parser.sha256_argument,required=True)
    cli.add_argument("--anchor",type=int,required=True)
    cli.add_argument("--first",type=int)
    cli.add_argument("--last",type=int)
    cli.add_argument("--z",type=generator.rational64,help="optional expected z; header is authoritative")
    cli.add_argument("--alpha",type=generator.rational64,help="optional expected alpha; header is authoritative")
    cli.add_argument("--seconds",type=int,choices=(300,600),default=300)
    args=cli.parse_args()
    if args.first is None: args.first=args.anchor
    if args.last is None: args.last=args.anchor
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    witness,input_report=parser.read_witness(args.witness,args.expect_witness_sha)
    z,alpha,bits=parameters(witness,parser)
    generator.validate(args.anchor,args.first,args.last,z,alpha)
    if args.z is not None and args.z!=z: raise ValueError("requested z differs from witness")
    if args.alpha is not None and args.alpha!=alpha: raise ValueError("requested alpha differs from witness")
    base=generator.load_base(); maximum=base.load_maximum()
    complete=maximum.load_complete(); selected=complete.load_selected(); run=selected.load_run()
    if (selected.RUN_PIN,run.COUNTS_PIN)!=(generator.RUN_PIN,generator.COUNTS_PIN):
        raise ArithmeticError("frozen run/count provenance changed")
    actual=run.load_actual(); support=actual.load_support()
    normalized=support.load_normalized(); band=normalized.load_band()
    counts_module=band.load("cheap_lane_cell_run_counts.py",generator.COUNTS_PIN,"dynamic_replay_counts")
    fitter=band.load(support.FITTER_FILE,support.FITTER_PIN,"dynamic_replay_fitter")
    geometry=band.load(normalized.GEOMETRY_FILE,normalized.GEOMETRY_PIN,"dynamic_replay_geometry")
    reference=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"dynamic_replay_reference")
    ordinary=reference.load("cheap_lane_ordered_ordinary_directed.py",reference.ORDINARY_PIN,"dynamic_replay_ordinary")
    large=ordinary.load_large()
    source=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"dynamic_replay_source")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"dynamic_replay_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"dynamic_replay_normalizers")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"dynamic_replay_terminal")
    peak=source.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("dynamic residue replay deadline; no incomplete claim")
        if peak()>MEMORY: raise MemoryError("dynamic residue replay512MiB budget exceeded")
    def progress(phase,**values):
        if phase=="actual_marker_fit" and values.get("center") not in (0,512,1024): return
        if phase=="run_band_geometry" and values.get("row",0)%8: return
        kept={key:values[key] for key in ("center","status","evaluations","row","vertices") if key in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,"peak_process_bytes":peak(),**kept}),flush=True)
    tests={"receipt":exact_tests(reference.BANK,parser,guard),"parameters":generator.exact_tests(guard),
        "bank":support.wrapper_tests(reference.BANK,guard),"run_counts":counts_module.exact_tests(join,guard),
        "actual":actual.exact_tests(reference.BANK,fitter,normalized,geometry,guard)}
    built=local.build_local_pairs(z,alpha,source,guard,deadline,PRECISION,progress=progress)
    v,core=built["verifier"],built["core"]; degrees=built["rational_coefficients"]
    original_arb=dict(built["coefficients"]); table_hash=support.coefficient_fingerprint(degrees,guard)
    if table_hash!=built["report"]["rational_coefficient_sha256"] or table_hash!=witness["local_coefficient_sha256"]:
        raise ArithmeticError("rebuilt exact local table differs from witness")
    proposed,fit=actual.fit_actual_markers(MappingProxyType(degrees),args.anchor,fitter,guard,progress)
    fitted,_,bank_hash=support.validate_bank(proposed,reference.BANK)
    if bank_hash!=witness["positive_rational_bank_sha256"]: raise ArithmeticError("anchor-fitted bank differs from witness")
    if support.coefficient_fingerprint(degrees,guard)!=table_hash or len(original_arb)!=len(built["coefficients"]):
        raise ArithmeticError("fitter changed local caps")
    if any(row is not built["coefficients"][key] for key,row in original_arb.items()):
        raise ArithmeticError("fitter changed original directed rows")
    regression=generator.old_reference(args.anchor,z,alpha)
    if regression and (table_hash,bank_hash)!=(base.TABLE_SHA,base.BANK_SHA):
        raise ArithmeticError("old3072 regression mismatch")
    exact,validation=validate_witness(witness,args.expect_witness_sha,fitted,table_hash,bank_hash,parser,guard)
    numerators=tuple(tuple(tuple(v.exact(x) for x in row) for row in region) for region in exact)
    cells=normalized.normalized_cells(fitted)
    marker_logs=tuple((v.exact(row[3]).log(),v.exact(row[4]).log()) for row in cells)
    den=actual.actual_denominators(v,core.BETA,alpha,band,guard)
    tests["directed_actual"]=actual.directed_tests(v,den,extraction,fitted,alpha,core.BETA,join,band,guard)
    del den
    result=generator.evaluate_band(v,core,run,actual,normalized,band,geometry,counts_module,source,join,
        cells,numerators,marker_logs,args.first,args.last,z,alpha,guard,progress)
    guard()
    print(json.dumps({"scope":"complete actual-size band conditional on trusted verified dynamic maximum receipt",
        "directed_arithmetic":True,"whole_code_certificate":False,
        "trusted_verified_generator_receipt_required":True,"regional_upper_bounds_independently_recomputed":False,
        "regional_tables_built":0,"regional_polynomial_products":0,"local_table_rebuilt_for_provenance":True,
        "parameters":{"anchor_actual":args.anchor,"actual_first":args.first,"actual_last":args.last,
            "z":str(z),"alpha":str(alpha),"recorded_integer_bits":bits,
            "working_precision":built["report"]["working_precision"],"seconds":args.seconds,"memory_budget_bytes":MEMORY},
        "requested_band":result,"input":input_report,"validation":validation,
        "local_coefficient_sha256":table_hash,"positive_rational_bank_sha256":bank_hash,
        "old3072_regression_hashes_enforced":regression,"unnormalized_exact_dyadics_restored":True,
        "actual_source_counts_output_replayed_once":True,"alpha_not_reapplied_to_numerators":True,
        "degree_reversals":0,"frozen_module_globals_changed":False,
        "fit":{"scope":"FLOAT proposal used only to verify the exact witness markers",
            "successful_proposal_evaluations_recorded":fit["successful_proposal_evaluations_recorded"],"fit_failures":fit["fit_failures"]},
        "tests":tests,"local_checks":{"tests":built["report"]["tests"],"field_tests":built["report"]["field_tests"]},
        "elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "nonnegative_intersections":v.intersections,
        "pins":{"generator":GENERATOR_PIN,"parser":PARSER_PIN,"base":generator.BASE_PIN,
            "maximum":generator.MAXIMUM_PIN,"run":generator.RUN_PIN,"counts":generator.COUNTS_PIN,
            "actual":run.ACTUAL_PIN,"support":actual.SUPPORT_PIN,"fitter":support.FITTER_PIN,
            "normalized":support.NORMALIZED_PIN,"band":normalized.BAND_PIN,"geometry":normalized.GEOMETRY_PIN,
            "reference":band.HOLE_PIN,"ordinary":reference.ORDINARY_PIN,"local":large.LOCAL_PIN,
            "source":large.ANCHOR_PIN,"join":large.JOIN_PIN,"normalizers":large.REGION_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},separators=(",",":")),flush=True)


if __name__=="__main__":
    main()
