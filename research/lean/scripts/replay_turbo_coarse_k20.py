"""Full K20 cap-free proof-witness DIAGNOSTIC; no Lean certificate claim.

Use exact retained per-occupancy witnesses, coarse j>=3 birth, candidate
profiles and exact j1/j2 birth counts. Preserve rare/q1/q2-regular components.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
import hashlib
import json
from math import comb
from pathlib import Path
import subprocess
import sys
from time import monotonic

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'research/workstreams/k18_margin60'))
sys.path.insert(0,str(ROOT/'research/workstreams/k20_small_outer_screen'))
import replay_fractional as replay
import dense_complete as dense_arithmetic
from flint import arb,ctx
from turbo_coarse_physical import CoarsePhysical,aq,up,binary_target_floor


def run(args):
    if args.output.exists():raise FileExistsError('fresh output required')
    ctx.prec=256
    start=monotonic()
    current=json.loads((args.receipt_root/'k20-two-ms-rs8-transfer-v1.json').read_text())
    reference=json.loads((args.receipt_root/'k20-small-outer-whole-v1.json').read_text())
    components=[json.loads(Path(c['path']).read_text())for c in reference['components']]
    candidate=json.loads(args.candidate.read_text())
    model=CoarsePhysical(candidate)
    chosen={}
    result=dict(schema='turbo-coarse-k20-full-diagnostic-1',status='partial',
        scope='Outward replay from untrusted compact counts plus retained rare/q1/q2; not a Lean certificate.',
        groups=8192,cutoff=209715,precision=256,
        candidate_sha256=hashlib.sha256(args.candidate.read_bytes()).hexdigest(),
        retained_components={k:current['component_uppers'][k]for k in ('rare','q1','q2_regular')},
        occupancy_uppers=chosen,witnesses=[])
    def save():
        result['elapsed_seconds']=monotonic()-start
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    beta=Q(1<<256,65535**8)
    schedules=defaultdict(list)
    for q in range(3,616):
        c=components[reference['selected_components'][str(q)]]
        tilt=Q(c['tilt'])if 'tilt'in c else Q(c['occupancy_choices'][str(q)]['tilt'])
        alpha=Q(c.get('alpha','1'))
        schedules[tilt,alpha].append(q)
    for (tilt,alpha),qs in schedules.items():
        blocks=replay.fractional_blocks(model.family(tilt),alpha)
        regional=replay.placement.placement_power(blocks,epochs=256,windows=32,maximum_groups=max(qs))
        for q in qs:
            power=regional[q]**64
            moment=sum((power[0,j]for j in range(11)),arb(0))
            value=arb(comb(8192,q))*(aq(alpha)*(q*aq(beta).log()+aq(tilt)*209715)).exp()*moment
            chosen[str(q)]=replay.endpoint(value)
        summary=dict(kind='exact_regional',tilt=str(tilt),alpha=str(alpha),count=len(qs),
            minimum=min(qs),maximum=max(qs),
            worst_margin_bits=min(replay.bits(replay.fraction(chosen[str(q)]))for q in qs))
        result['witnesses'].append(summary);save();print(json.dumps(summary),flush=True)
    dense=components[10]
    schedules=defaultdict(list)
    for q in range(616,8193):
        w=dense['occupancy_choices'][str(q)]
        schedules[Q(w['tilt']),Q(w['alpha']),Q(w['marker_probability'])].append(q)
    block_cache={}
    logs=dense_arithmetic._logs(616,8192)
    for (tilt,alpha,p),qs in schedules.items():
        if (tilt,alpha)not in block_cache:
            block_cache[tilt,alpha]=replay.fractional_blocks(model.family(tilt),alpha)
        mixed=dense_arithmetic.marker_mixture(block_cache[tilt,alpha],p)
        power=mixed**16384
        moment=sum((power[0,j]for j in range(11)),arb(0))
        log_moment=moment.log()
        for q in qs:
            value=dense_arithmetic.coefficient_log_bound(q,p,log_moment,logs[q],tilt,alpha)
            chosen[str(q)]=replay.endpoint(value.exp())
        summary=dict(kind='iid_fractional_conditioning',tilt=str(tilt),alpha=str(alpha),
            marker_probability=str(p),count=len(qs),minimum=min(qs),maximum=max(qs),
            worst_margin_bits=min(replay.bits(replay.fraction(chosen[str(q)]))for q in qs))
        result['witnesses'].append(summary);save();print(json.dumps(summary),flush=True)
    if set(chosen)!={str(q)for q in range(3,8193)}:raise ArithmeticError('incomplete tail coverage')
    tail=sum((replay.fraction(v)for v in chosen.values()),Q())
    total=tail+sum((replay.fraction(v)for v in result['retained_components'].values()),Q())
    result.update(status='complete',exact_tail_union=replay.encode_exact(tail),
        exact_union=replay.encode_exact(total),margin_bits_display=replay.bits(total),
        tail_margin_bits_display=replay.bits(tail),
        rounded_51_82_safe_threshold_met=total<binary_target_floor(Q(2591,50)))
    save();print(json.dumps({k:result[k]for k in ('status','margin_bits_display','tail_margin_bits_display',
        'rounded_51_82_safe_threshold_met','elapsed_seconds')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--receipt-root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--max-seconds',type=float,default=600)
    p.add_argument('--_worker',action='store_true',help=argparse.SUPPRESS)
    a=p.parse_args()
    if a._worker:run(a)
    else:
        try:
            cp=subprocess.run([sys.executable,str(Path(__file__).resolve()),*sys.argv[1:],'--_worker'],timeout=a.max_seconds)
            raise SystemExit(cp.returncode)
        except subprocess.TimeoutExpired:
            print('Bounded worker timed out; completed witness checkpoint retained',file=sys.stderr)
            raise SystemExit(2)
