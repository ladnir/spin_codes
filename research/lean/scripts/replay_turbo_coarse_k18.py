"""Bounded full K18 proof-witness DIAGNOSTIC; not a Lean certificate.

Retains the current RS8 q1/q2 endpoints and recomputes every q>=3 using
compact count candidates, uncapped births, and the retained witness schedule.
"""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path
import subprocess
import sys
from time import monotonic

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'research/workstreams/k18_margin60'))
import replay_fractional as fractional
from flint import arb, arb_mat, ctx
from turbo_coarse_physical import CoarsePhysical, aq, up, binary_target_floor


def run(args):
    if args.output.exists():
        raise FileExistsError('fresh output required')
    started = monotonic()
    ctx.prec = 256
    base = args.receipt_root
    candidate = json.loads(args.candidate.read_text())
    current = json.loads((base/'k18-rs8-gf32-s22-v1.json').read_text())
    ordinary = json.loads((base/'rs16-s22-k262144-whole-fresh-p256-exact.json').read_text())
    dense = json.loads((base/'rs16-s22-k262144-whole-fresh-p256-dense.json').read_text())
    model = CoarsePhysical(candidate)
    chosen = {q:current['occupancy_uppers'][q] for q in ('1','2')}
    result = dict(schema='turbo-coarse-k18-full-diagnostic-1', status='partial',
        scope='Outward arithmetic from untrusted exact count candidates plus retained q1/q2; '
              'not Lean-authenticated local rows or a final certificate.',
        candidate_sha256=hashlib.sha256(args.candidate.read_bytes()).hexdigest(),
        groups=2048, cutoff=52428, precision=256, occupancy_uppers=chosen, witnesses=[])
    def save():
        result['elapsed_seconds'] = monotonic()-started
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    beta = Q(1<<256, 65535**8)
    schedules = {}
    for q in range(3,385):
        tilt, alpha = ((Q(17,200), Q(4,5)) if 62 <= q <= 143 else
                       (Q(ordinary['occupancy_choices'][str(q)]['tilt']), Q(1)))
        schedules.setdefault((tilt,alpha),[]).append(q)
    for (tilt,alpha), qs in schedules.items():
        local = model.family(tilt)
        blocks = fractional.fractional_blocks(local,alpha)
        regional = fractional.placement.placement_power(blocks,epochs=64,windows=32,maximum_groups=max(qs))
        for q in qs:
            chosen[str(q)] = fractional.endpoint(fractional.scalar_bound(regional[q],q,tilt,alpha))
        summary = dict(kind='exact_regional',tilt=str(tilt),alpha=str(alpha),
            count=len(qs),minimum=min(qs),maximum=max(qs),
            worst_margin_bits=min(fractional.bits(fractional.fraction(chosen[str(q)])) for q in qs))
        result['witnesses'].append(summary)
        save()
        print(json.dumps(summary),flush=True)
    schedules = {}
    for q in range(385,2049):
        witness = dense['occupancy_choices'][str(q)]
        schedules.setdefault((Q(witness['tilt']),Q(witness['marker_probability'])),[]).append(q)
    for (tilt,p), qs in schedules.items():
        local = model.family(tilt)
        theta = aq(15*p/16)
        matrix = arb_mat(11,11)
        for j in range(17):
            matrix += (comb(16,j)*theta**j*(1-theta)**(16-j))*local[j]
        power = matrix**8192
        moment = sum((power[0,k] for k in range(11)),arb(0))
        for q in qs:
            logchoose = arb(comb(2048,q)).log()
            loss = arb(0) if p == 1 else -logchoose-q*aq(p).log()-(2048-q)*aq(1-p).log()
            chosen[str(q)] = fractional.endpoint(up((logchoose+q*aq(beta).log()+aq(tilt)*52428+64*loss).exp()*moment))
        summary = dict(kind='iid_conditioning',tilt=str(tilt),marker_probability=str(p),
            count=len(qs),minimum=min(qs),maximum=max(qs),
            worst_margin_bits=min(fractional.bits(fractional.fraction(chosen[str(q)])) for q in qs))
        result['witnesses'].append(summary)
        save()
        print(json.dumps(summary),flush=True)
    if set(chosen) != {str(q) for q in range(1,2049)}:
        raise ArithmeticError('incomplete occupation coverage')
    total = sum((fractional.fraction(v) for v in chosen.values()),Q())
    tail = sum((fractional.fraction(v) for q,v in chosen.items() if int(q)>=3),Q())
    result.update(rounded_58_24_safe_threshold_met=total<binary_target_floor(Q(1456,25)),
        status='complete',exact_union=fractional.encode_exact(total),
        exact_tail_union=fractional.encode_exact(tail),
        margin_bits_display=fractional.bits(total),tail_margin_bits_display=fractional.bits(tail))
    save()
    print(json.dumps({k:result[k] for k in ('status','margin_bits_display','tail_margin_bits_display','elapsed_seconds')}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate',type=Path,required=True)
    parser.add_argument('--receipt-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--max-seconds',type=float,default=360)
    parser.add_argument('--_worker',action='store_true',help=argparse.SUPPRESS)
    args=parser.parse_args()
    if args._worker:
        run(args)
    else:
        try:
            cp=subprocess.run([sys.executable,str(Path(__file__).resolve()),*sys.argv[1:],'--_worker'],timeout=args.max_seconds)
            raise SystemExit(cp.returncode)
        except subprocess.TimeoutExpired:
            print('Bounded worker timed out; completed witnesses retained in output',file=sys.stderr)
            raise SystemExit(2)
