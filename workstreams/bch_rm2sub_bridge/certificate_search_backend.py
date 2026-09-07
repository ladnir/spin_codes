"""Fixed-instance numerical engine. No scheduling or acceptance policy here."""
import math
from fractions import Fraction as F
from pathlib import Path
from flint import arb,ctx

import certificate_search_core as core
import frontier_sparse as sparse
import certify_refresh_q1 as refresh

base=core.base


def sources():
    result=core.inputs.source_hashes();result.update(sparse.sources())
    result.update(core.outer_dependencies())
    for path in (Path(__file__),Path(core.__file__),Path(refresh.__file__)):
        result[path.relative_to(base.ROOT).as_posix()]=base.sha(path)
    return result


def compute(task,precision=256,saved=None):
    spec=task['instance']
    core.require(spec==core.instance(spec['configuration'],spec['message_exponent']),'Instance fingerprint mismatch')
    ctx.prec=precision;t,s,spectrum,kernel=core.inputs.load(spec['configuration'])
    rows,cutoff=spec['rows'],spec['cutoff']
    if task['kind']=='q1':
        witnesses=task.get('scaled_tilts',[264,295,332,376])
        core.require(witnesses and all(type(a) is int and 1<=a<=10000 for a in witnesses),'Invalid Q1 witnesses')
        best={w:F(rows) for w in base.WEIGHTS}
        for a in witnesses:
            lam=arb(a)/(100*rows)
            reg=refresh.region(*refresh.epoch(t,s,spectrum,(-lam).exp(),arb),rows//t,arb,linear=precision>=512)
            co=refresh.coefficients(*reg,256,arb);factor=rows*(cutoff*lam).exp()
            for w in best: best[w]=min(best[w],sparse.rational((co[w]*factor).upper()))
        if saved:
            recorded={int(w):base.decode(v) for w,v in saved['coefficient_upper'].items()}
            core.require(set(recorded)==set(best) and all(best[w]<=recorded[w] for w in best),'Q1 coefficient replay failed')
            best=recorded
        upper,_,_=base.bch_bound(best)
        return dict(bounds=[dict(occupation=1,upper=base.encode(upper))],
                    coefficient_upper={str(w):base.encode(v) for w,v in best.items()})
    core.require(task['kind']=='range','Unsupported certificate backend')
    lo,hi=task['interval'];tilt=task['tilt']
    core.require(type(lo) is int and type(hi) is int and 2<=lo<=hi<=rows,'Invalid occupancy interval')
    core.require(type(tilt) is int and -120<=tilt<=30,'Invalid range tilt')
    region=sparse.poly.regions(t,s,spectrum,kernel,(-(arb(tilt)/10).exp()).exp(),hi,length=rows)
    print('region ready',spec['configuration'],lo,hi,'precision',precision,flush=True)
    witness=saved['witness'] if saved else task.get('witness')
    if witness is None: witness=sparse.choose(region,lo,tilt,rows,cutoff,core.inputs.caps_module.caps())
    core.require(witness['tilt']==tilt,'Witness tilt mismatch')
    values=sparse.evaluate(region,witness,lo,hi,rows,cutoff,core.inputs.caps_module.caps())
    return dict(witness=witness,bounds=[dict(occupation=r['occupation'],upper=base.encode(F(2)**r['upper_power']))
                for r in values],diagnostic_margins=[r['margin_bits_diagnostic'] for r in values])


def produce(task,output):
    result=compute(task)
    base.write_new(output,dict(schema='bch-certificate-search-v1',task=task,result=result,
        status='OUTWARD_PRODUCER_NOT_YET_REPLAYED',source_sha256=sources()))


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--task',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();produce(base.read(a.task),a.output)
