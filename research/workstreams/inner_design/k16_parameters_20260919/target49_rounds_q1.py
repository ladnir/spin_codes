"""Outward Q1 evaluation for two rounds on a selected t64 expansion subspace."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from flint import ctx
import model
import subspace
import target49_q1


class Engine(subspace.Engine):
    def __init__(self,record):
        super().__init__(record)
        self.n=2*len(self.levels)+2
        images=[sum(((a&q).bit_count()&1)<<p for p,a in enumerate(self.a_columns)) for q in self.columns]
        self.iw=[x.bit_count() for x in images]
        self.cw=[(x^(1<<p)).bit_count() for p,x in enumerate(images)]

    def identity(self):
        result=super().identity()
        result['inner']['transvection_rounds']=2
        result['setup']='independent row/region permutations; two independent transvections before feedback per update; zero start; output before update; persistent state; no flush'
        return result

    @lru_cache(maxsize=24)
    def epoch(self,tilt):
        z=(-model.number(tilt).exp()).exp()
        zero,one,n=model.study.core.mixing.transfers(self.spectrum,self.iw,self.cw,z,model.number(F(1,4)))
        assert n==self.n
        return tuple(map(model.up,zero)),tuple(map(model.up,one))


def run(a):
    out=a.output.resolve()
    saved=model.base.base.read(out) if a.verify else None
    if saved: model.base.authenticate(saved)
    elif out.exists(): raise FileExistsError(out)
    ctx.prec=512 if saved else 256
    source=model.base.base.read(a.maps)
    record=next(r['inner'] for r in source['rows'] if r['found'] and r['s']==a.s)
    e=Engine(record)
    if saved: assert saved['instance']==e.identity()
    audit=model.study.core.mixing.audit
    tilts=saved['tilts'] if saved else [str(F(i,40)) for i in (-280,-240,-224,-220,-216,-212,-208,-200,-160,-120)]
    best=None
    for tilt in tilts:
        row=audit.coefficients(e,F(tilt),linear=bool(saved))
        best=row if best is None else audit.minimum(best,row)
        print('two-round subspace',a.s,tilt,audit.bits(model.base.base.bch_bound(best)[0]),flush=True)
    if saved:
        accepted={int(w):model.base.base.decode(v) for w,v in saved['coefficients'].items()}
        assert set(best)==set(accepted) and all(best[w]<=accepted[w] for w in best)
        best=accepted
    upper,_,_=model.base.base.bch_bound(best)
    result=dict(status='K16_TWO_ROUND_SUBSPACE_Q1_OUTWARD',instance=e.identity(),tilts=tilts,
        coefficients={str(w):model.base.base.encode(v) for w,v in best.items()},
        upper=model.base.base.encode(upper),margin_bits=audit.bits(upper),
        full_distance_proved=False,source_sha256=model.base.sources())
    for path in (Path(__file__),Path(model.study.core.mixing.__file__),a.maps):
        result['source_sha256'][path.resolve().relative_to(model.base.ROOT).as_posix()]=model.base.base.sha(path)
    if saved:
        assert result['upper']==saved['upper']
        result=dict(status='K16_TWO_ROUND_SUBSPACE_Q1_512_LINEAR_REPLAY_PASSED',
            producer_sha256=model.base.base.sha(out),full_distance_proved=False)
        out=out.with_name(out.stem+'_replay.json')
    model.base.base.write_new(out,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--s',type=int,default=12)
    p.add_argument('--maps',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    run(p.parse_args())
