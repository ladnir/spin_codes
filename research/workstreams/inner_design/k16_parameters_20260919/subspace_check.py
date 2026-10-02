"""Outward point checks for an explicitly saved expansion subspace."""
import argparse
from fractions import Fraction as F
from pathlib import Path
from flint import arb,ctx
import subspace
import model


def run(a):
    out=a.output.resolve()
    saved=model.base.base.read(out) if a.verify else None
    if saved:
        model.base.authenticate(saved)
        assert (saved['s'],saved['mode'])==(a.s,a.mode)
    elif out.exists(): raise FileExistsError(out)
    source=model.base.base.read(a.maps)
    record=next(r['inner'] for r in source['rows'] if r['s']==a.s and r['found'])
    ctx.prec=512 if saved else 256
    e=subspace.Engine(record)
    if saved: assert saved['instance']==e.identity()
    result=dict(status='K16_SUBSPACE_OUTWARD_POINTS',s=a.s,mode=a.mode,
        instance=e.identity(),full_distance_proved=False,precision_bits=ctx.prec)
    screen=model.study.core.grid.ladder.screen
    if a.mode=='q1':
        result['q1']=screen.q1(e,saved['q1'] if saved else None)
        print('outward subspace Q1',a.s,result['q1']['margin_bits'],flush=True)
    elif a.mode=='sparse':
        result['sparse']=[]
        for i,q in enumerate(a.q):
            row=screen.sparse(e,q,saved['sparse'][i] if saved else None)
            result['sparse'].append(row)
            e.region.cache_clear();e.epoch.cache_clear()
            print('outward subspace sparse',a.s,q,row['margin_bits'],flush=True)
    else:
        c=model.Checker(64,a.s)
        c.engine=e
        result['points']=[]
        for i,q in enumerate(a.q):
            v=F(a.coordinate)
            witness=saved['points'][i]['witness'] if saved else c.witness(q,q,v,v)
            upper=c.bound(q,q,v,v,witness).exp()
            bound=model.base.pack(upper)
            if saved:
                assert upper <= model.base.unpack(saved['points'][i]['upper'])
                bound=saved['points'][i]['upper']
            row=dict(q=q,coordinate=str(v),witness=witness,upper=bound,
                margin_bits=float(-model.base.unpack(bound).log()/arb(2).log()))
            result['points'].append(row)
            print('outward subspace dense',a.s,q,str(v),row['margin_bits'],flush=True)
    result['source_sha256']=model.base.sources()
    result['source_sha256'][a.maps.resolve().relative_to(model.base.ROOT).as_posix()]=model.base.base.sha(a.maps)
    if saved:
        result=dict(status='K16_SUBSPACE_POINTS_512_REPLAY_PASSED',
            producer_sha256=model.base.base.sha(out),full_distance_proved=False)
        out=out.with_name(out.stem+'_replay.json')
    model.base.base.write_new(out,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--maps',type=Path,required=True)
    p.add_argument('--s',type=int,required=True)
    p.add_argument('--mode',choices=['q1','sparse','dense'],required=True)
    p.add_argument('--q',type=int,nargs='+',default=[2,4,8,16,32,63])
    p.add_argument('--coordinate',default='21/128')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    run(p.parse_args())
