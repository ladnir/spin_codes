"""Replay actual-size geometry/count bounds from a verified compact bank.

The input is ONLY the exact_maximum_witness JSON object exported by the
hash-pinned residue_max_band generator. Require its canonical SHA256 and
at most64KiB. Check exact headers,31 ordered marker rows, bounded integer
dyadics and separate regular/terminal scales. Regenerate the fixed local
table and rational marker bank, and require their recorded hashes and
the same exact markers in every witness row. Restore exact Fractions,
then convert them directly to directed Arb numerators.

This adapter trusts the previously verified generator receipt for the
UPPER-BOUND VALIDITY of the six saved values. A matching hash establishes
identity, not that arbitrary supplied values are valid coefficient caps.
No regional polynomial is rebuilt here. The adapter independently replays
the requested complete band geometry, actual source factors, exact run
counts, both terminal orders and final output factor using those caps.

The numerator bank remains unnormalized: no new alpha, source factor or
target denominator is absorbed into its six dyadics. Fixed z1/10 and
alpha119/1280, actual coordinates only, inclusive sizes641..4096.
Parent execution only:300/600s,512MiB,20GiB free,-B,stdout and no generated
files. Module import performs no work or proof-provider imports.
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
GENERATOR_PIN="b2761383e6a929978764baa1af9f777e22907c662d505f352015e16f3d8a8938"
MAXIMUM_PIN="d03a5dff0c607cac71f7b6da45c355c07654ca91bea879d726e58ce05c098036"
TABLE_SHA="636b295acf8ba93a20f76a983ab4723c290f9967a252a04b9562148a3a7c539d"
BANK_SHA="3e4843f4abd6be8730f72a5f77c02f4b4e554bf32a834464bda3df8643f551e7"
Z,ALPHA=F(1,10),F(119,1280)
ANCHOR,PRECISION,MEMORY=3072,384,512*(1<<20)
MAX_BYTES,MAX_EXPONENT=64*(1<<10),262144
HEADER_KEYS=frozenset(("coordinate","moduli","integer_bits","z","alpha",
    "local_coefficient_sha256","positive_rational_bank_sha256","rows"))


def load_generator():
    path=Path(__file__).with_name("cheap_lane_residue_max_band.py")
    if hashlib.sha256(path.read_bytes()).hexdigest()!=GENERATOR_PIN:
        raise RuntimeError("frozen compact-bank generator source mismatch")
    spec=importlib.util.spec_from_file_location("compact_maximum_generator",path)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    if (module.MAXIMUM_PIN,module.TABLE_SHA,module.BANK_SHA,module.Z,module.ALPHA)!=\
            (MAXIMUM_PIN,TABLE_SHA,BANK_SHA,Z,ALPHA):
        raise ArithmeticError("frozen generator witness contract changed")
    return module


def sha256_argument(text):
    if len(text)!=64 or any(c not in "0123456789abcdefABCDEF" for c in text):
        raise argparse.ArgumentTypeError("expected64 hexadecimal SHA256 characters")
    return text.lower()


def canonical_bytes(witness):
    return json.dumps(witness,sort_keys=True,separators=(",",":"),ensure_ascii=True,
        allow_nan=False).encode("ascii")


def witness_hash(witness):
    return hashlib.sha256(canonical_bytes(witness)).hexdigest()


def decode_witness(data):
    """Strict bounded JSON syntax; no floats, duplicate keys or huge integers."""
    if not isinstance(data,bytes) or not 1<=len(data)<=MAX_BYTES:
        raise ValueError("witness must contain1..65536 bytes")
    def integer(text):
        if len(text)>8: raise ValueError("JSON integer exceeds bounded header/exponent syntax")
        return int(text)
    def reject_number(text):
        raise ValueError("floating/nonfinite JSON numbers are forbidden")
    def unique(pairs):
        out={}
        for key,value in pairs:
            if key in out: raise ValueError("duplicate JSON key")
            out[key]=value
        return out
    try:
        witness=json.loads(data.decode("utf-8"),parse_int=integer,parse_float=reject_number,
            parse_constant=reject_number,object_pairs_hook=unique)
    except (UnicodeError,json.JSONDecodeError,RecursionError) as exc:
        raise ValueError("invalid bounded witness JSON") from exc
    if type(witness) is not dict or set(witness)!=HEADER_KEYS:
        raise ValueError("expected the raw exact_maximum_witness object, not a full receipt")
    return witness


def read_witness(path,expected_sha):
    path=Path(path)
    if not path.is_file(): raise ValueError("witness path must name an existing regular file")
    if not 1<=path.stat().st_size<=MAX_BYTES:
        raise ValueError("witness file exceeds64KiB or is empty")
    with path.open("rb") as stream: data=stream.read(MAX_BYTES+1)
    witness=decode_witness(data)
    actual=witness_hash(witness)
    if actual!=expected_sha:
        raise ValueError("canonical witness SHA256 differs from the required receipt")
    return witness,{"input_bytes":len(data),"input_file_sha256":hashlib.sha256(data).hexdigest(),
        "canonical_witness_sha256":actual,"canonical_hash_checked":True,
        "canonical_hash_format":"sorted-key compact JSON of exact_maximum_witness"}


def decode_dyadic(pair,bits):
    if type(pair) is not list or len(pair)!=2:
        raise ValueError("dyadic must be[mantissa_decimal_string,integer_exponent]")
    text,exponent=pair
    if type(text) is not str or not 1<=len(text)<=78 or not text.isascii() or not text.isdecimal():
        raise ValueError("bounded nonnegative decimal mantissa required")
    if len(text)>1 and text[0]=="0": raise ValueError("mantissa must be canonical decimal")
    mantissa=int(text)
    if mantissa.bit_length()>bits: raise ValueError("mantissa exceeds recorded shared bit budget")
    if type(exponent) is not int or abs(exponent)>MAX_EXPONENT:
        raise ValueError("dyadic exponent exceeds safe exact-conversion range")
    value=F(mantissa<<exponent) if exponent>=0 else F(mantissa,1<<(-exponent))
    return value,exponent


def validate_witness(witness,expected_sha,bank,guard):
    """Validate identity/layout/markers, NOT upper-bound validity of receipt."""
    if type(witness) is not dict or set(witness)!=HEADER_KEYS:
        raise ValueError("unexpected witness keys")
    if len(canonical_bytes(witness))>MAX_BYTES:
        raise ValueError("canonical witness exceeds64KiB")
    if witness_hash(witness)!=expected_sha:
        raise ValueError("canonical witness hash mismatch")
    if (witness["coordinate"],witness["moduli"],witness["z"],witness["alpha"],
            witness["local_coefficient_sha256"],witness["positive_rational_bank_sha256"])!=\
            ("actual",[256,256],str(Z),str(ALPHA),TABLE_SHA,BANK_SHA):
        raise ValueError("witness headers differ from fixed actual3072 construction")
    if type(witness["moduli"]) is not list or any(type(x) is not int for x in witness["moduli"]):
        raise ValueError("exact integer modulus layout required")
    bits=witness["integer_bits"]
    if type(bits) is not int or bits not in (192,256):
        raise ValueError("recorded mantissa bits must be192 or256")
    rows=witness["rows"]
    if type(rows) is not list or len(rows)!=31 or len(bank)!=31:
        raise ValueError("exactly31 ordered witness and reference rows required")
    numerators=[]; exponent_pairs=[]; zero_entries=0; maximum_mantissa_bits=0
    for index,(row,reference) in enumerate(zip(rows,bank)):
        guard()
        if type(row) is not list or len(row)!=5:
            raise ValueError("witness row must be[center,u,v,R4,V2]")
        center,u,w,regular,terminal=row
        if type(center) is not int or (center,u,w)!=(reference[0],str(reference[3]),str(reference[4])):
            raise ValueError("witness marker row is permuted, rescaled or inconsistent with regenerated bank")
        if type(u) is not str or type(w) is not str or max(len(u),len(w))>80:
            raise ValueError("bounded canonical rational marker strings required")
        if type(regular) is not list or type(terminal) is not list or (len(regular),len(terminal))!=(4,2):
            raise ValueError("six-channel R4/V2 shape required")
        values=[]; scales=[]
        for group in (regular,terminal):
            decoded=[decode_dyadic(pair,bits) for pair in group]
            if len({exponent for _,exponent in decoded})!=1:
                raise ValueError("each regular/terminal object must retain its own shared exponent")
            values.append(tuple(value for value,_ in decoded)); scales.append(decoded[0][1])
            zero_entries+=sum(value==0 for value,_ in decoded)
            maximum_mantissa_bits=max(maximum_mantissa_bits,max(int(pair[0]).bit_length() for pair in group))
        numerators.append(tuple(values)); exponent_pairs.append([center,*scales])
    if len(numerators)!=31: raise ArithmeticError("incomplete exact numerator bank")
    return tuple(numerators),{"rows":31,"exact_dyadics":186,"integer_bits":bits,
        "channel_order":["R00","R01","R10","R11","V0","V1"],
        "maximum_mantissa_bits":maximum_mantissa_bits,"zero_entries":zero_entries,
        "maximum_absolute_exponent":max(abs(e) for row in exponent_pairs for e in row[1:]),
        "separate_object_exponents":exponent_pairs,"regenerated_marker_rows_checked":31,
        "exact_dyadic_conversion":True,"trusted_generator_receipt_required":True,
        "upper_bound_values_independently_recomputed":False}


def exact_tests(bank,guard):
    """Stand-in validation tests; no regional arithmetic or external input."""
    checks=0
    def clone(value): return json.loads(json.dumps(value))
    rows=[]
    for i,row in enumerate(bank):
        rows.append([row[0],str(row[3]),str(row[4]),
            [[str(6*i+j),-11] for j in range(4)],[[str(6*i+j),7] for j in (4,5)]])
    witness={"coordinate":"actual","moduli":[256,256],"integer_bits":192,
        "z":str(Z),"alpha":str(ALPHA),"local_coefficient_sha256":TABLE_SHA,
        "positive_rational_bank_sha256":BANK_SHA,"rows":rows}
    digest=witness_hash(witness)
    decoded,report=validate_witness(witness,digest,bank,guard)
    assert len(decoded)==31 and report["exact_dyadics"]==186; checks+=1
    for i,(regular,terminal) in enumerate(decoded):
        assert regular==tuple(F(6*i+j,2048) for j in range(4)); checks+=1
        assert terminal==tuple(F((6*i+j)*128) for j in (4,5)); checks+=1
    assert witness_hash(decode_witness(json.dumps(witness,indent=2).encode("ascii")))==digest; checks+=1
    # These preserve the expected trusted hash: ANY changed value/order fails.
    mutated=[]
    bad=clone(witness); bad["rows"][0],bad["rows"][1]=bad["rows"][1],bad["rows"][0]; mutated.append(bad)
    bad=clone(witness); bad["rows"][0][3][1][0]="2"; mutated.append(bad)
    bad=clone(witness); bad["rows"][0][4][0][1]=-11; mutated.append(bad)
    bad=clone(witness); bad["rows"][0][3][0],bad["rows"][0][3][1]=bad["rows"][0][3][1],bad["rows"][0][3][0]; mutated.append(bad)
    for bad in mutated:
        try: validate_witness(bad,digest,bank,guard)
        except ValueError: checks+=1
        else: raise AssertionError("changed witness accepted under original trusted hash")
    # Even with a fresh identity hash, malformed headers/layout/markers fail.
    invalid=[]
    for key,value in (("coordinate","hole"),("moduli",[128,256]),("integer_bits",True),
            ("z","1/8"),("alpha","1/2"),("local_coefficient_sha256","0"*64),
            ("positive_rational_bank_sha256","0"*64)):
        bad=clone(witness); bad[key]=value; invalid.append(bad)
    bad=clone(witness); bad["rows"].pop(); invalid.append(bad)
    bad=clone(witness); bad["rows"].append(clone(bad["rows"][-1])); invalid.append(bad)
    bad=clone(witness); bad["rows"][0],bad["rows"][1]=bad["rows"][1],bad["rows"][0]; invalid.append(bad)
    bad=clone(witness); bad["rows"][0][1]="2"; invalid.append(bad)
    bad=clone(witness); bad["rows"][0][3].pop(); invalid.append(bad)
    bad=clone(witness); bad["rows"][0][4][1][1]=8; invalid.append(bad)
    bad=clone(witness); bad["rows"][0][3][0]=["-1",0]; invalid.append(bad)
    bad=clone(witness); bad["rows"][0][3][0]=["01",0]; invalid.append(bad)
    bad=clone(witness); bad["rows"][0][3][0]=[str(1<<192),0]; invalid.append(bad)
    bad=clone(witness); bad["rows"][0][3][0]=["0",MAX_EXPONENT+1]; invalid.append(bad)
    bad=clone(witness); bad["rows"][0][3][0]=["0",False]; invalid.append(bad)
    bad=clone(witness); bad["extra"]=0; invalid.append(bad)
    for bad in invalid:
        try: validate_witness(bad,witness_hash(bad),bank,guard)
        except ValueError: checks+=1
        else: raise AssertionError("malformed witness accepted")
    for data in (b'{"a":1,"a":2}',b'{"a":1.0}',b'{"a":NaN}',b'{"a":123456789}',
                 b'[]',b' '* (MAX_BYTES+1)):
        try: decode_witness(data)
        except ValueError: checks+=1
        else: raise AssertionError("invalid bounded JSON accepted")
    assert decode_dyadic(["0",-7],192)==(F(0),-7); checks+=1
    assert decode_dyadic(["3",-5],192)==(F(3,32),-5); checks+=1
    assert decode_dyadic(["3",7],192)==(F(384),7); checks+=1
    return {"exact_checks":checks,"mutation_rejection_uses_original_expected_hash":True,
        "header_layout_marker_checks_independent_of_identity_hash":True,
        "standin_values_are_not_claimed_upper_bounds":True,"no_regional_products":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--witness",type=Path,required=True)
    parser.add_argument("--expect-witness-sha",type=sha256_argument,required=True)
    parser.add_argument("--first",type=int,default=3008)
    parser.add_argument("--last",type=int,default=3136)
    parser.add_argument("--seconds",type=int,choices=(300,600),default=300)
    args=parser.parse_args()
    if not 641<=args.first<=args.last<=4096: raise ValueError("actual band must lie in641..4096")
    if shutil.disk_usage(Path(__file__).resolve().anchor).free<20*(1<<30):
        raise RuntimeError("requires20GiB free on source volume")
    started=time.monotonic(); deadline=started+args.seconds
    witness,input_report=read_witness(args.witness,args.expect_witness_sha)
    generator=load_generator(); maximum=generator.load_maximum()
    complete=maximum.load_complete(); selected=complete.load_selected(); run=selected.load_run()
    if (selected.RUN_PIN,run.COUNTS_PIN)!=(generator.RUN_PIN,generator.COUNTS_PIN):
        raise ArithmeticError("frozen run/count provenance changed")
    actual=run.load_actual(); support=actual.load_support()
    normalized=support.load_normalized(); band=normalized.load_band()
    counts_module=band.load("cheap_lane_cell_run_counts.py",generator.COUNTS_PIN,"max_replay_counts")
    fitter=band.load(support.FITTER_FILE,support.FITTER_PIN,"max_replay_fitter")
    geometry=band.load(normalized.GEOMETRY_FILE,normalized.GEOMETRY_PIN,"max_replay_geometry")
    reference=band.load("cheap_lane_ordinary_hole_directed.py",band.HOLE_PIN,"max_replay_reference")
    ordinary=reference.load("cheap_lane_ordered_ordinary_directed.py",reference.ORDINARY_PIN,"max_replay_ordinary")
    large=ordinary.load_large()
    anchor=large.load("cheap_lane_cone_cauchy_anchors.py",large.ANCHOR_PIN,"max_replay_source")
    local=large.load("cheap_lane_ordered_local_directed.py",large.LOCAL_PIN,"max_replay_local")
    extraction=large.load("cheap_lane_band_regions.py",large.REGION_PIN,"max_replay_normalizers")
    join=large.load("cheap_lane_directed_run_mass.py",large.JOIN_PIN,"max_replay_terminals")
    peak=anchor.peak_memory_reader()
    def guard():
        if time.monotonic()>deadline: raise TimeoutError("compact-bank replay deadline; no incomplete claim")
        if peak()>MEMORY: raise MemoryError("compact-bank replay512MiB limit exceeded")
    def progress(phase,**values):
        if phase=="actual_marker_fit" and values.get("center") not in (0,512,1024): return
        if phase=="run_band_geometry" and values.get("row",0)%8: return
        kept={key:values[key] for key in ("center","status","evaluations","row","vertices") if key in values}
        print(json.dumps({"phase":phase,"seconds":time.monotonic()-started,
            "peak_process_bytes":peak(),**kept}),flush=True)
    tests={"receipt_validation":exact_tests(reference.BANK,guard),
        "run_wrapper":run.wrapper_tests(guard),"bank":support.wrapper_tests(reference.BANK,guard),
        "actual":actual.exact_tests(reference.BANK,fitter,normalized,geometry,guard),
        "run_counts":counts_module.exact_tests(join,guard)}
    built=local.build_local_pairs(Z,ALPHA,anchor,guard,deadline,PRECISION,progress=progress)
    v,core=built["verifier"],built["core"]; degrees=built["rational_coefficients"]
    original_arb=dict(built["coefficients"])
    if (support.coefficient_fingerprint(degrees,guard),built["report"]["rational_coefficient_sha256"])!=\
            (TABLE_SHA,TABLE_SHA): raise ArithmeticError("local table differs from witness provenance")
    proposed,fit=actual.fit_actual_markers(MappingProxyType(degrees),ANCHOR,fitter,guard,progress)
    fitted,_,bank_hash=support.validate_bank(proposed,reference.BANK)
    if bank_hash!=BANK_SHA: raise ArithmeticError("regenerated bank differs from witness provenance")
    if support.coefficient_fingerprint(degrees,guard)!=TABLE_SHA or len(original_arb)!=len(built["coefficients"]):
        raise ArithmeticError("fitter changed local coefficients")
    if any(row is not built["coefficients"][key] for key,row in original_arb.items()):
        raise ArithmeticError("fitter changed original Arb rows")
    exact,validation=validate_witness(witness,args.expect_witness_sha,fitted,guard)
    numerators=tuple(tuple(tuple(v.exact(x) for x in row) for row in region) for region in exact)
    cells=normalized.normalized_cells(fitted)
    marker_logs=tuple((v.exact(cell[3]).log(),v.exact(cell[4]).log()) for cell in cells)
    den=actual.actual_denominators(v,core.BETA,ALPHA,band,guard)
    tests["directed_actual"]=actual.directed_tests(v,den,extraction,fitted,ALPHA,core.BETA,join,band,guard)
    del den
    result=generator.evaluate_band(v,core,run,actual,normalized,band,geometry,counts_module,anchor,join,
        cells,numerators,marker_logs,args.first,args.last,guard,progress)
    guard()
    print(json.dumps({"scope":"complete actual-size band replay conditional on a trusted verified maximum receipt",
        "directed_arithmetic":True,"whole_code_certificate":False,
        "trusted_verified_generator_receipt_required":True,
        "regional_upper_bounds_independently_recomputed":False,"regional_tables_built":0,
        "regional_polynomial_products":0,"local_table_rebuilt_for_provenance":True,
        "parameters":{"actual_first":args.first,"actual_last":args.last,"anchor_actual":ANCHOR,
            "z":str(Z),"alpha":str(ALPHA),"recorded_integer_bits":witness["integer_bits"],
            "working_precision":built["report"]["working_precision"],"seconds":args.seconds,
            "memory_budget_bytes":MEMORY},
        "requested_band":result,"input":input_report,"validation":validation,
        "local_coefficient_sha256":TABLE_SHA,"positive_rational_bank_sha256":bank_hash,
        "unnormalized_exact_dyadics_restored":True,"source_output_and_counts_replayed":True,
        "alpha_not_reapplied_to_numerators":True,"degree_reversals":0,
        "fit":{"scope":"frozen proposal procedure used to verify exact witness markers",
            "successful_proposal_evaluations_recorded":fit["successful_proposal_evaluations_recorded"],
            "fit_failures":fit["fit_failures"]},
        "tests":tests,"local_checks":{"tests":built["report"]["tests"],"field_tests":built["report"]["field_tests"]},
        "elapsed_seconds":time.monotonic()-started,"peak_process_bytes":peak(),"generated_bytes":0,
        "nonnegative_intersections":v.intersections,"pins":{"generator":GENERATOR_PIN,
            "maximum":MAXIMUM_PIN,"complete3072":maximum.COMPLETE_PIN,"selected":maximum.SELECTED_PIN,
            "run":generator.RUN_PIN,"counts":generator.COUNTS_PIN,"actual":run.ACTUAL_PIN,
            "support":actual.SUPPORT_PIN,"fitter":support.FITTER_PIN,"normalized":support.NORMALIZED_PIN,
            "band":normalized.BAND_PIN,"geometry":normalized.GEOMETRY_PIN,"reference":band.HOLE_PIN,
            "ordinary":reference.ORDINARY_PIN,"local":large.LOCAL_PIN,"source":large.ANCHOR_PIN,
            "join":large.JOIN_PIN,"normalizers":large.REGION_PIN},
        "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},separators=(",",":")),flush=True)


if __name__=="__main__":
    main()
