"""Certify every nonzero outer occupancy of the fixed K16 wider byte code.

Numerical display margins are diagnostics. Acceptance is an exact rational
comparison of the complete sum of dyadic upper endpoints with 2^-40.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import log2
from pathlib import Path
import re
from time import monotonic
import numpy as np
import outward_positive as op
import scaled_positive as sp
import global_outward as go
import macro_outward as mo
import q1_outward as q1

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]


def checked_pins(pins):
    return bool(pins) and all(Path(name).is_file() and hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest
                             for name,digest in pins.items())


def merge_pins(target, incoming):
    for name,digest in incoming.items():
        if name in target and target[name]!=digest: raise ArithmeticError('conflicting source pins')
        target[name]=digest


def source_pins():
    pins=q1.source_pins()
    merge_pins(pins,mo.source_pins())
    paths=[Path(__file__).resolve(),HERE/'test_certify_wider24.py',
           Path(go.__file__).resolve(),HERE/'test_global_outward.py']
    # Associate the proof receipt with the exact experimental implementation.
    experiment=REPO/'spin'/'experiments'/'packet8_wider24'
    paths += [experiment/name for name in ('Packet8Wide24.h','Setup.cpp','Scalar.cpp','Fast.cpp',
                                          'Outer.cpp','bench.cpp','CMakeLists.txt','run_ab.sh')]
    # Follow literal local includes; the certificate does not pin the compiler
    # or its system headers. The timed executable has a separate recorded hash.
    visited=set()
    pending=list(paths)
    while pending:
        path=pending.pop().resolve()
        if path in visited:continue
        visited.add(path)
        if path.suffix not in ('.h','.cpp'):continue
        for name in re.findall(r'^\s*#\s*include\s*"([^"]+)"',path.read_text(encoding='utf-8'),re.M):
            dependency=(path.parent/name).resolve()
            if not dependency.is_relative_to(REPO):raise ValueError('include outside repository')
            if not dependency.is_file():raise ValueError(f'missing local include: {dependency}')
            pending.append(dependency)
    paths=list(visited)
    for path in paths: pins[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    return pins


def load_macro(path):
    saved=json.loads(path.read_text(encoding='utf-8'))
    if (saved.get('schema')!='packet8-actual24-G4-outward-1'
            or not saved.get('source_pins_verified_at_finish')
            or not saved.get('outward_under_stated_arithmetic_contract')
            or not checked_pins(saved.get('source_sha256',{}))
            or saved['diagnostics']['group_steps']!=4
            or saved['diagnostics']['macro_windows']!=32
            or saved['diagnostics']['physical_windows']!=8
            or not saved['diagnostics']['tuple_weights_outside_fractional_power']):
        raise ValueError('authenticated outward G4 receipt required')
    active,z,local=mo.load_local(Path(saved['local_receipt']))
    alpha=Fraction(saved['alpha'])
    if (z!=Fraction(saved['z']) or saved['map_record']!=local['map_record']
            or Fraction(saved['a'])!=((1+z)/2)**8 or not 0<alpha<=1
            or Fraction(saved['diagnostics']['alpha'])!=alpha):
        raise ValueError('inconsistent local/macro witness')
    macro=np.asarray(saved['macro_operator_upper'],dtype=float)
    hexadecimal=np.asarray([[[float.fromhex(x) for x in row] for row in m]
                            for m in saved['macro_operator_upper_hex']])
    if (macro.shape!=(33,10,10) or not np.isfinite(macro).all() or np.any(macro<0)
            or not np.array_equal(macro,hexadecimal)):
        raise ValueError('finite hexadecimal-confirmed macro endpoints required')
    return macro,active,z,alpha,saved


def dyadic_record(value):
    value=Fraction(value)
    if value<0 or value.denominator&(value.denominator-1):
        raise ValueError('nonnegative exact dyadic required')
    return dict(numerator_hex=hex(value.numerator),exponent=-(value.denominator.bit_length()-1))


def read_dyadic(record):
    n=int(record['numerator_hex'],16)
    e=record['exponent']
    if n<0 or type(e) is not int: raise ValueError('nonnegative dyadic required')
    return Fraction(n)*Fraction(2)**e


def round_probability(value):
    """Weaken excessively small endpoints to keep receipts compact."""
    return max(Fraction(1,2**200),min(Fraction(1),Fraction(value)))


def bits(value):
    value=Fraction(value)
    if not value:return float('inf')
    return log2(value.denominator)-log2(value.numerator)


def exact_union(endpoints, *, groups=256, target=40):
    if set(endpoints)!=set(range(1,groups+1)):
        raise ValueError('every nonzero occupancy must be supplied')
    values=[Fraction(endpoints[q]) for q in range(1,groups+1)]
    if any(x<0 or x>1 for x in values):raise ValueError('occupancy probability uppers required')
    total=sum(values,Fraction(0))
    return total,total<=Fraction(1,2**target)


def run(args):
    if args.output.exists() or not args.macros:
        raise ValueError('at least one macro and fresh output required')
    op.check_runtime()
    pins,started=source_pins(),monotonic()
    best={q:Fraction(1) for q in range(2,257)}
    choices={q:'trivial probability cap' for q in range(2,257)}
    best_support=tuple(Fraction(1) for _ in range(65))
    seen_locals=set()
    receipt=dict(schema='packet8-wider24-whole-code-outward-1',
        geometry=dict(K=65536,N=131072,rate='1/2',outer_groups=256,group_dimension=256,
            group_output_bits=512,symbol_bits=32,regions=64,slots_per_region=256,
            packet_bits=8,physical_t=64,state_bits=24,steps_per_region=32,
            cutoff=13107,zero_initial_state=True,final_flush=False,
            continuous_state=True),target_margin_bits=40,source_sha256=pins,
        source_pins_verified_at_finish=False,whole_code_certificate=False,
        outward_under_stated_arithmetic_contract=True,
        positive_runtime_check_passed=True,selected_q_endpoint_floor='2^-200',
        subset_union_outside_fractional_power=True,trials=[])
    for path in args.macros:
        macro,active,z,alpha,saved=load_macro(path)
        merge_pins(pins,saved['source_sha256'])
        pins[str(path.resolve())]=hashlib.sha256(path.read_bytes()).hexdigest()
        regional=go.placement(macro,256)
        trial=dict(macro_receipt=str(path.resolve()),z=str(z),alpha=str(alpha),endpoints={})
        for q in range(2,257):
            upper=go.occupancy_upper(regional[q],q,z,alpha)
            bound=round_probability(sp.scalar_fraction(upper))
            trial['endpoints'][str(q)]=dyadic_record(bound)
            if bound<best[q]:best[q],choices[q]=bound,str(path.resolve())
        if z not in seen_locals:
            support=q1.q1_support_upper(active,z)
            best_support=q1.min_supports(best_support,support)
            trial['q1_support_uppers']=[dyadic_record(x) for x in support]
            seen_locals.add(z)
        receipt['trials'].append(trial)
        if not checked_pins(pins):raise ArithmeticError('source changed during global evaluation')
        print(f'global alpha={alpha}, z={float(z):.7g}: minimum q>=2 margin '
              f'{bits(max(best.values())):.8f} bits; elapsed {monotonic()-started:.2f}s',flush=True)
    best[1]=round_probability(q1.combine_supports(best_support))
    choices[1]='exact wider outer shells and per-support minimum of active-label witnesses'
    total,accepted=exact_union(best)
    worst=max(best,key=best.get)
    receipt.update(selected_endpoints={str(q):dyadic_record(best[q]) for q in range(1,257)},
        selected_witnesses={str(q):choices[q] for q in range(1,257)},
        q1_selected_support_uppers=[dyadic_record(x) for x in best_support],
        q1_exact_expected_shells=True,q1_margin_bits_display=bits(best[1]),
        exact_union=dyadic_record(total),exact_union_below_2_minus40=accepted,
        combined_margin_bits_display=bits(total),weakest_q=worst,
        weakest_margin_bits_display=bits(best[worst]),all_occupancies_checked=True,
        elapsed_seconds=monotonic()-started)
    if not checked_pins(pins):raise ArithmeticError('final source authentication failed')
    receipt['source_pins_verified_at_finish']=True
    receipt['whole_code_certificate']=accepted
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(f'EXACT UNION PASS={accepted}; margin={bits(total):.9f} bits; '
          f'weakest q={worst}; q1={bits(best[1]):.9f}',flush=True)
    return receipt


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--macros',type=Path,nargs='+',required=True)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
