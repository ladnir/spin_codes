"""Refine the central BA spectrum supports without changing frozen inputs.

New affine supports are verified over their active intervals by the frozen
outward box formulas. Exact binary partition paths permit independent replay
of coverage. Reflections use p(b,w)=p(b,1-w) in the BA objective.
"""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
import sys
import time

from flint import ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import screen
import certify_imt_dense as dense
sys.path.insert(0,str(screen.FROZEN))
import certify_golay_ba_concave_majorant as old


def supports():
    ctx.prec = 256
    prior = screen.load_segments(screen.FROZEN/'golay_ba3_concave_majorant.json')
    lines = {(s.slope,s.intercept) for s in prior if s.slope>=0}
    added = set()
    for slope in (F(9,100),F(7,100),F(5,100),F(3,100),F(1,100)):
        # Entropy tangent plus 1e-5 is a proposal, not its justification.
        intercept = dense.upper((1+(-dense.number(slope)).exp()).log()
                                -dense.number(2).log()/2+dense.number(F(1,100000)))
        added.add((slope,intercept))
    lines |= added
    boundaries = {F(13,125),F(1,2)}
    for a,b in lines:
        for c,d in lines:
            if a != c:
                x = (d-b)/(a-c)
                if F(13,125)<x<F(1,2):boundaries.add(x)
    boundaries = sorted(boundaries)
    result = []
    for lo,hi in zip(boundaries,boundaries[1:]):
        slope,intercept = min(lines,key=lambda line:line[0]*(lo+hi)/2+line[1])
        if result and (result[-1].slope,result[-1].intercept)==(slope,intercept):
            prev=result.pop();lo=prev.omega_lo
        result.append(old.Segment(len(result),lo,hi,slope,intercept))
    assert result[-1].slope == 0
    assert all(a.slope>b.slope for a,b in zip(result,result[1:]))
    new = []
    for segment in result:
        if (segment.slope,segment.intercept) in added:
            new.append(segment)
        else:
            assert any(s.slope==segment.slope and s.intercept==segment.intercept
                       and s.lower<=segment.omega_lo<segment.omega_hi<=s.upper for s in prior)
    return result,new


def root(segment):
    return old.Box(segment,F(0),F(1),F(0),F(1),segment.omega_lo,segment.omega_hi)


def replay_partition(segment,leaves):
    """Traverse a complete binary partition; reject duplicate/unused paths."""
    paths = {row['path']:row for row in leaves}
    assert len(paths)==len(leaves) and all(set(p)<=set('01') for p in paths)
    prefixes={p[:i] for p in paths for i in range(len(p))}
    todo=[('',root(segment))];seen=set();maximum=-math.inf
    while todo:
        path,box=todo.pop()
        if path in paths:
            assert path not in prefixes,'overlapping ancestor leaf'
            bound=box.objective_upper()
            assert bound<0,(segment.index,path,bound)
            saved=paths[path]['upper']
            assert (saved is None and bound==-math.inf) or (saved is not None and bound<=saved)
            maximum=max(maximum,bound);seen.add(path)
        else:
            assert path in prefixes,'coverage gap'
            a,b=box.split();todo.extend([(path+'0',a),(path+'1',b)])
    assert seen==set(paths)
    return maximum


def main():
    if not __debug__:raise RuntimeError('Do not disable certificate assertions.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--verify',type=Path)
    parser.add_argument('--seconds',type=float,default=180)
    parser.add_argument('--max-boxes',type=int,default=200000)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    segments,new=supports()
    prior=screen.FROZEN/'golay_ba3_concave_majorant.json'
    saved=json.loads(args.verify.read_text()) if args.verify else None
    if saved:
        assert saved['status']=='proved'
        for name,digest in saved['source_sha256'].items():
            assert screen.sha(screen.ROOT/name)==digest,name
        old.iv.dps=100
        for segment in new:
            leaves=saved['new_support_checks'][str(segment.index)]
            maximum=replay_partition(segment,leaves)
            print('replayed support',segment.index,'leaves',len(leaves),'maximum',maximum,flush=True)
        payload=dict(saved,replay_interval_dps=old.iv.dps,replay_passed=True)
        payload['source_sha256']=dict(saved['source_sha256'],**{args.verify.resolve().relative_to(screen.ROOT).as_posix():screen.sha(args.verify)})
    else:
        old.iv.dps=70
        checks={str(s.index):[] for s in new}
        pending=[(s,'',root(s)) for s in reversed(new)]
        count=0;maximum=-math.inf;start=time.monotonic();last=start
        while pending and count<args.max_boxes and time.monotonic()-start<args.seconds:
            segment,path,box=pending.pop()
            upper=box.objective_upper();count+=1
            if upper<0:
                checks[str(segment.index)].append(dict(path=path,upper=upper if math.isfinite(upper) else None))
                maximum=max(maximum,upper)
            else:
                if len(path)>=90:raise RuntimeError(('depth limit',segment.index,path,upper))
                a,b=box.split();pending.extend([(segment,path+'1',b),(segment,path+'0',a)])
            now=time.monotonic()
            if now-last>=10:
                print('boxes',count,'accepted',sum(map(len,checks.values())),'pending',len(pending),'support',segment.index,'seconds',round(now-start,1),flush=True);last=now
        payload=dict(schema='golay-ba3-refined-concave-majorant-v1',status='proved' if not pending else 'unresolved',
                     interval_dps=70,weight_interval=['13/125','112/125'],
                     central_constant_upper=str(segments[-1].intercept),
                     left_segments=[dict(weight=[str(s.omega_lo),str(s.omega_hi)],slope=str(s.slope),intercept=str(s.intercept)) for s in segments[:-1]],
                     right_segments_are_reflections=True,new_support_checks=checks,
                     processed_boxes=count,remaining_boxes=len(pending),maximum_accepted_upper=maximum,
                     seconds=time.monotonic()-start,
                     source_sha256={p.relative_to(screen.ROOT).as_posix():screen.sha(p) for p in [Path(__file__),Path(old.__file__),Path(dense.__file__),Path(screen.__file__),prior]},
                     scope='New supports certified on active intervals; inherited supports unchanged; reflection is an exact symmetry of the asymptotic BA objective.')
    with args.output.open('x',encoding='utf-8') as stream:
        json.dump(payload,stream,indent=2);stream.write('\n')
    print(json.dumps({k:v for k,v in payload.items() if k not in ('new_support_checks','source_sha256','left_segments')}),flush=True)
    if payload['status']!='proved':raise SystemExit(1)


if __name__=='__main__':main()
