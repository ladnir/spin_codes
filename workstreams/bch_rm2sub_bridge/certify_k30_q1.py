"""Q1-only certificate at K=2^30; does not import smaller-K occupancy bounds."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

from flint import arb, arb_poly, ctx
import bridge as base
import larger_state_maps as maps
import certify_refresh_q1 as refresh
from audit_bch_q1_full_arb import rational
from migrate_legacy_workspace import restore_or_verify

ROWS = 1 << 23
CUTOFF = (256 * ROWS) // 10


def polynomial_region(zero, one, count):
    power = tuple(arb_poly([a,b]) for a,b in zip(zero,one))
    current = tuple(arb_poly([int(i == j)]) for i in range(4) for j in range(4))
    def multiply(a,b):
        return tuple(sum((a[4*i+k]*b[4*k+j] for k in range(4)),arb_poly()).truncate(2)
                     for i in range(4) for j in range(4))
    # Left-to-right powering, separate from the producer's paired matrices.
    for bit in bin(count)[2:]:
        current = multiply(current,current)
        if bit == '1':
            current = multiply(current,power)
    return tuple(v[0] for v in current),tuple(v[1]/count for v in current)


def run(output, verify=False):
    output = output.resolve()
    old = base.read(output) if verify else None
    if not verify:
        assert not output.exists()
    restore_or_verify(base.ROOT)
    t,s,spectrum,_ = maps.load('t64_s20')
    witnesses = [F(a,100*ROWS) for a in (264,295,332,376)]
    ctx.prec = 512 if verify else 256
    best = {w:F(ROWS) for w in base.WEIGHTS}
    for witness in witnesses:
        lam = arb(witness.numerator)/witness.denominator
        epoch = refresh.epoch(t,s,spectrum,(-lam).exp(),arb)
        region = polynomial_region(*epoch,ROWS//t) if verify else refresh.region(*epoch,ROWS//t,arb)
        co = refresh.coefficients(*region,256,arb)
        factor = ROWS*(CUTOFF*lam).exp()
        for w in best:
            best[w] = min(best[w],rational((co[w]*factor).upper()))
        print('replay' if verify else 'producer','scaled tilt',float(witness*ROWS),flush=True)
    bound,factor,rest = base.bch_bound(best)
    assert bound < F(1,1 << 40)
    if verify:
        for name,digest in old['source_sha256'].items():
            assert base.sha(base.ROOT/name) == digest
        assert old['parameters'] == dict(message_bits=1<<30,output_bits=1<<31,
            outer_rows=ROWS,step_bits=64,state_bits=20,cutoff=CUTOFF)
        assert all(best[w] <= base.decode(old['coefficient_upper'][str(w)]) for w in best)
        saved = {int(w):base.decode(v) for w,v in old['coefficient_upper'].items()}
        assert base.bch_bound(saved) == tuple(base.decode(old[k]) for k in ('Q1_upper','factor','rest'))
        assert F(1,1<<40)-base.decode(old['Q1_upper']) == base.decode(old['remaining_failure_budget'])
        base.write_new(output.with_name(output.stem+'_replay.json'),dict(
            status='K30_Q1_512_BIT_POLYNOMIAL_REPLAY_PASSED',producer_sha256=base.sha(output),
            coefficients_checked=92,full_distance_proved=False))
    else:
        paths = [Path(__file__),Path(refresh.__file__),Path(base.__file__),Path(maps.__file__),
                 base.HERE/'MIGRATION_MANIFEST.json',base.HERE/'migrate_legacy_workspace.py']
        base.write_new(output,dict(status='K30_Q1_OUTWARD_CERTIFICATE',configuration='t64_s20',
            parameters=dict(message_bits=1<<30,output_bits=1<<31,outer_rows=ROWS,step_bits=64,state_bits=20,cutoff=CUTOFF),
            Q1_upper=base.encode(bound),factor=base.encode(factor),rest=base.encode(rest),
            coefficient_upper={str(w):base.encode(v) for w,v in best.items()},
            lambda_witnesses=[base.encode(v) for v in witnesses],
            remaining_failure_budget=base.encode(F(1,1<<40)-bound),
            margin_bits=math.log2(bound.denominator)-math.log2(bound.numerator),
            full_distance_proved=False,covered_occupancies=[1],
            source_sha256={p.relative_to(base.ROOT).as_posix():base.sha(p) for p in paths}))
    print('Q1 margin',math.log2(bound.denominator)-math.log2(bound.numerator),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();run(args.output,args.verify)
