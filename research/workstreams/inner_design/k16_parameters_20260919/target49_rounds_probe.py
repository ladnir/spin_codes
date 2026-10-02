"""Reevaluate the old cover's largest leaves with the two-round engine."""
import argparse
from pathlib import Path
from flint import arb,ctx
import model
import target49_rounds_model as rounds
import subspace_cover


def run(a):
    saved=model.base.base.read(a.output) if a.verify else None
    if saved: model.base.authenticate(saved)
    elif a.output.exists(): raise FileExistsError(a.output)
    ctx.prec=512 if saved else 256
    maps=model.base.base.read(a.maps)
    record=next(r['inner'] for r in maps['rows'] if r['s']==12 and r['found'])
    c=rounds.Checker(record)
    seed=model.base.base.read(a.seed)
    model.base.authenticate(seed)
    control=subspace_cover.Checker(record)
    assert seed['instance']==control.engine.identity()
    keys=[r['key'] for r in saved['rows']] if saved else sorted(seed['leaves'],key=lambda k:seed['leaves'][k]['power'],reverse=True)[:a.count]
    rows=[]
    for i,key in enumerate(keys):
        node=seed['leaves'][key]
        geometry=subspace_cover.budget_dense.geometry.geometry(node)
        upper=c.bound(*geometry,node['witness']).exp()
        packed=model.base.pack(upper)
        if saved:
            assert upper<=model.base.unpack(saved['rows'][i]['upper'])
            packed=saved['rows'][i]['upper']
        row=dict(key=key,old_power=node['power'],upper=packed,
                 margin_bits=float(-model.base.unpack(packed).log()/arb(2).log()))
        rows.append(row)
        print('two-round dense',key,'old power',node['power'],'new margin',row['margin_bits'],flush=True)
    result=dict(status='K16_TARGET49_TWO_ROUND_DENSE_PROBE',instance=c.engine.identity(),
        rows=rows,full_distance_proved=False,source_sha256=model.base.sources())
    for path in (Path(__file__),Path(rounds.__file__),a.maps,a.seed,Path(subspace_cover.__file__)):
        result['source_sha256'][path.resolve().relative_to(model.base.ROOT).as_posix()]=model.base.base.sha(path)
    if saved:
        assert saved['instance']==result['instance']
        result=dict(status='K16_TARGET49_TWO_ROUND_DENSE_512_REPLAY_PASSED',
            producer_sha256=model.base.base.sha(a.output),full_distance_proved=False)
        output=a.output.with_name(a.output.stem+'_replay.json')
    else: output=a.output
    model.base.base.write_new(output,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--maps',type=Path,required=True)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--count',type=int,default=10)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    run(p.parse_args())
