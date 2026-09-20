"""Outward before/after comparisons at saved difficult points; not a cover."""
import argparse
from fractions import Fraction as F
from pathlib import Path
from flint import arb,ctx
import overlap_transfer as sharpened
model=sharpened.model


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=sharpened.HERE/'POINTS.json')
    p.add_argument('--verify',action='store_true')
    args=p.parse_args()
    saved=model.base.read(args.output) if args.verify else None
    if saved:model.authenticate(saved)
    else:assert not args.output.exists()
    ctx.prec=512 if saved else 256
    gap=model.base.read(model.HERE/'DENSE_GAP_POINT.json')
    points=[gap]+model.base.read(model.HERE/'DENSE_RETUNED_POINTS.json')['points']
    results=[]
    cache={}
    for i,row in enumerate(points):
        exponent=row['message_exponent']
        if exponent not in cache:cache[exponent]=(sharpened.original.Checker(exponent),sharpened.Checker(exponent))
        old,new=cache[exponent]
        q,v,witness=row['occupation'],F(row['coordinate']),row['witness']
        before=old.bound(q,q,v,v,witness)
        after=new.bound(q,q,v,v,witness)
        bits=model.exact((after/arb(2).log()).upper())
        power=-(-bits.numerator//bits.denominator)
        result=dict(message_exponent=exponent,occupation=q,coordinate=str(v),witness=witness,
                    old_margin_bits=float(-before/arb(2).log()),margin_bits=float(-after/arb(2).log()),upper_power=power)
        if saved:
            previous=saved['results'][i]
            assert [previous[k] for k in ('message_exponent','occupation','coordinate','witness')]==[result[k] for k in ('message_exponent','occupation','coordinate','witness')]
            assert power<=previous['upper_power']
        results.append(result)
        print(ctx.prec,exponent,q,str(v),result['old_margin_bits'],'->',result['margin_bits'],flush=True)
    if saved:
        model.base.write_new(args.output.with_name(args.output.stem+'_replay.json'),
                             dict(status='OVERLAP_POINT_512_BIT_REPLAY_PASSED',producer_sha256=model.base.sha(args.output),full_distance_proved=False))
    else:
        sources=model.sources()
        for path in (sharpened.HERE/'PACKING_AUDIT.json',model.HERE/'DENSE_GAP_POINT.json',model.HERE/'DENSE_RETUNED_POINTS.json'):
            sources[path.relative_to(model.ROOT).as_posix()]=model.base.sha(path)
        model.base.write_new(args.output,dict(status='OUTWARD_OVERLAP_POINT_BOUNDS_NOT_COVER',precision_bits=256,
                                             results=results,full_distance_proved=False,source_sha256=sources))


if __name__=='__main__':main()
