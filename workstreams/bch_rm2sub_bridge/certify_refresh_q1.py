"""Outward Q1 refinement using exact uniform refresh after zero input.

Uses fresh nonzero multipliers, the existing fixed map, and fixed BCH dual.
Producer: 256-bit Arb, binary epoch powering. Replay: 512-bit Arb, linear
epoch iteration. Neither changes an existing receipt or assumes a spectrum.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import sys

from flint import arb, ctx
import bridge as base
import larger_state_maps as maps
sys.path.insert(0, str(base.BCH / 'code'))
from audit_bch_q1_full_arb import rational


def product(a, b):
    return tuple(a[4*i] * b[j] + a[4*i+1] * b[4+j]
                 + a[4*i+2] * b[8+j] + a[4*i+3] * b[12+j]
                 for i in range(4) for j in range(4))


def row_product(a, b):
    return (a[0]*b[0]+a[1]*b[4]+a[2]*b[8]+a[3]*b[12],
            a[0]*b[1]+a[1]*b[5]+a[2]*b[9]+a[3]*b[13],
            a[0]*b[2]+a[1]*b[6]+a[2]*b[10]+a[3]*b[14],
            a[0]*b[3]+a[1]*b[7]+a[2]*b[11]+a[3]*b[15])


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def epoch(t, s, spectrum, z, number):
    m=(1 << s)-1
    assert s >= 2 and sum(spectrum.values()) == m and min(spectrum) >= 1
    kappa=number(m)/(m-1)
    m0=sum((number(n)*z**w for w,n in spectrum.items()),number(0))/m
    m1=sum((number(n)*(w*z**(w-1)+(t-w)*z**(w+1))/t
            for w,n in spectrum.items()),number(0))/m
    zero=[number(0) for _ in range(16)];one=list(zero)
    zero[0]=number(1);zero[6]=z**min(spectrum)
    zero[10]=m0;zero[14]=kappa*m0
    one[1]=z
    for state,moment in ((1,z**(min(spectrum)-1)),(2,m1),(3,kappa*m1)):
        one[4*state]=moment/m
        one[4*state+3]=moment*(m-1)/m
    return tuple(zero),tuple(one)


def region(zero, one, count, number, linear=False):
    assert count >= 1
    rz=tuple(number(i==j) for i in range(4) for j in range(4))
    ra=tuple(number(0) for _ in range(16))
    if linear:
        for _ in range(count):
            ra=add(product(ra,zero),product(rz,one));rz=product(rz,zero)
    else:
        remaining=count;bz,ba=zero,one
        while remaining:
            if remaining&1:
                ra=add(product(ra,bz),product(rz,ba));rz=product(rz,bz)
            remaining>>=1
            if remaining:
                ba=add(product(ba,bz),product(bz,ba));bz=product(bz,bz)
    return rz,tuple(v/count for v in ra)


def coefficients(zero, one, block, number):
    empty=(number(0),)*4
    current=[(number(1),number(0),number(0),number(0))]
    for n in range(block):
        updated=[empty for _ in range(n+2)]
        for j,row in enumerate(current):
            a=row_product(row,zero);b=row_product(row,one)
            updated[j]=add(updated[j],a);updated[j+1]=add(updated[j+1],b)
        current=updated
    return [sum(row,number(0))/math.comb(block,w) for w,row in enumerate(current)]


def run(discovery, output, verify=False):
    discovery=discovery.resolve();output=output.resolve()
    from migrate_legacy_workspace import restore_or_verify
    restore_or_verify(base.ROOT)
    saved=base.read(output) if verify else None
    if not verify:
        assert not output.exists()
    screen=base.read(discovery)
    row=next(r for r in screen['rows'] if r['message_exponent']==20)
    t,s,spectrum,_=maps.load('t64_s20')
    witnesses=sorted({F.from_float(math.exp(row['witness_log_a'][str(w)]))/8192
                      for w in (38,40,42,44)}) if not verify else [base.decode(v) for v in saved['lambda_witnesses']]
    ctx.prec=512 if verify else 256
    best={w:F(8192) for w in base.WEIGHTS}
    for lam in witnesses:
        a=arb(lam.numerator)/lam.denominator;z=(-a).exp()
        moments=coefficients(*region(*epoch(t,s,spectrum,z,arb),8192//t,arb,linear=verify),256,arb)
        correction=8192*(base.CUTOFF*a).exp()
        for w in best:
            best[w]=min(best[w],rational((moments[w]*correction).upper()))
        print('replay' if verify else 'producer','scaled tilt',float(lam*8192),flush=True)
    upper, factor, rest=base.bch_bound(best)
    original_q1=base.decode(base.read(base.HERE/'generated/larger_t64_s20_q1_outward.json')['Q1_upper'])
    original_full=base.decode(base.read(base.HERE/'generated/larger_coverage_full.json')['covered_union_upper'])
    assert upper<original_q1
    combined=original_full-original_q1+upper
    assert combined<F(1,1<<50)
    if verify:
        for name,expected in saved['source_sha256'].items():
            assert base.sha(base.ROOT/name)==expected
        assert all(best[w]<=base.decode(saved['coefficient_upper'][str(w)]) for w in best)
        stored={int(w):base.decode(v) for w,v in saved['coefficient_upper'].items()}
        assert base.bch_bound(stored)==tuple(base.decode(saved[k]) for k in ('Q1_upper','factor','rest'))
        assert base.decode(saved['combined_full_upper'])==original_full-original_q1+base.decode(saved['Q1_upper'])
        base.write_new(output.with_name(output.stem+'_replay.json'),dict(
            status='REFRESH_Q1_512_BIT_LINEAR_EPOCH_REPLAY_PASSED',producer_sha256=base.sha(output),
            coefficients_checked=len(best),same_frozen_higher_occupancy_bounds=True))
        print('512-bit replay passed',flush=True)
        return
    dependencies=[Path(__file__),Path(base.__file__),Path(maps.__file__),discovery,
        base.HERE/'generated/larger_t64_s20_q1_outward.json',base.HERE/'generated/larger_coverage_full.json',
        base.HERE/'MIGRATION_MANIFEST.json',base.HERE/'migrate_legacy_workspace.py']
    base.write_new(output,dict(status='OUTWARD_REFRESH_Q1_IMPROVEMENT',configuration='t64_s20',
        message_bits=1<<20,cutoff=base.CUTOFF,lambda_witnesses=[base.encode(v) for v in witnesses],
        coefficient_upper={str(w):base.encode(v) for w,v in best.items()},
        Q1_upper=base.encode(upper),factor=base.encode(factor),rest=base.encode(rest),
        combined_full_upper=base.encode(combined),
        q1_margin_bits=math.log2(upper.denominator)-math.log2(upper.numerator),
        full_margin_bits=math.log2(combined.denominator)-math.log2(combined.numerator),
        original_full_margin_bits=math.log2(original_full.denominator)-math.log2(original_full.numerator),
        scope='Same fixed encoder and cutoff as T64_S20_FULL_CLOSURE.md; only Q1 replaced; no spectrum assumption.',
        source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in dependencies}))
    print('full margin',math.log2(combined.denominator)-math.log2(combined.numerator),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--discovery',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();run(args.discovery,args.output,args.verify)
