"""Small-step dense feasibility checks using the already-derived routing bounds."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from flint import arb,ctx
import model
import refine
import subspace_cover


class Checker(model.Checker):
    def __init__(self,t,s):
        super().__init__(t,s)
        self.engine=refine.Engine(t,s)
        self.constant_comparison=lru_cache(maxsize=2048)(self._constant_comparison)

    _comparison=subspace_cover.Checker._comparison
    _constant_comparison=subspace_cover.Checker._constant_comparison
    split_bound=subspace_cover.Checker.split_bound
    bound=subspace_cover.Checker.bound


def run(a):
    out=a.output.resolve()
    saved=model.base.base.read(out) if a.verify else None
    if saved: model.base.authenticate(saved)
    elif out.exists(): raise FileExistsError(out)
    ctx.prec=512 if saved else 256
    c=Checker(a.t,a.s)
    if saved: assert saved['instance']==c.engine.identity()
    points=[]
    for i,q in enumerate(a.q):
        v=F(a.coordinate)
        witness=saved['points'][i]['witness'] if saved else c.witness(q,q,v,v)
        upper=c.bound(q,q,v,v,witness).exp()
        bound=model.base.pack(upper)
        if saved:
            assert (q,str(v))==(saved['points'][i]['q'],saved['points'][i]['coordinate'])
            assert upper<=model.base.unpack(saved['points'][i]['upper'])
            bound=saved['points'][i]['upper']
        margin=float(-model.base.unpack(bound).log()/arb(2).log())
        points.append(dict(q=q,coordinate=str(v),witness=witness,upper=bound,margin_bits=margin))
        print('dense refined',a.t,a.s,q,str(v),margin,flush=True)
    result=dict(status='K16_TARGET49_REFINED_DENSE_POINTS',instance=c.engine.identity(),
        points=points,full_distance_proved=False,source_sha256=model.base.sources())
    for module in (Path(__file__),Path(subspace_cover.__file__)):
        result['source_sha256'][module.resolve().relative_to(model.base.ROOT).as_posix()]=model.base.base.sha(module)
    if saved:
        result=dict(status='K16_TARGET49_DENSE_512_REPLAY_PASSED',
            producer_sha256=model.base.base.sha(out),full_distance_proved=False)
        out=out.with_name(out.stem+'_replay.json')
    model.base.base.write_new(out,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--t',type=int,default=32)
    p.add_argument('--s',type=int,default=15)
    p.add_argument('--q',type=int,nargs='+',default=[64,128,218,384,512])
    p.add_argument('--coordinate',default='21/128')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    run(p.parse_args())
