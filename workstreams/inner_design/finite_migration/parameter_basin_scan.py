"""Scan separated IMT moment basins before another local dense refinement."""
import argparse
import copy
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
import parameter_scoped_cover as scoped


def run(seed_path,replay_path,output):
    bank=scoped.bank
    model,base=bank.model,bank.base
    saved,replay=[model.base.read(p) for p in (seed_path,replay_path)]
    model.authenticate(saved); model.authenticate(replay)
    assert saved['status']==scoped.STATUS and replay['status']=='VERIFIED_'+scoped.STATUS
    assert replay['producer_sha256']==model.base.sha(seed_path)
    with bank.previous.maps.use():
        moments=scoped.bridge.Moments()
        scoped.check(saved,moments)
        checker=moments.activation
        dense=copy.deepcopy(saved['dense'])
        target=max(dense['selected_boxes'],key=lambda b:b['own_log_bound'])
        target_before=copy.deepcopy(target)
        cache={}
        def evaluate(theta,tilt):
            theta,tilt=float(theta),float(tilt)
            if not .000001<=theta<=.999999 or not -8<=tilt<=1:
                return math.inf,None
            key=(theta,tilt)
            if key not in cache:
                if len(checker.epoch_cache)>4:
                    checker.epoch_cache.clear()
                term=checker.moment(theta,tilt)
                value,w=bank.proposal_bound(checker,target,theta,tilt,term)
                cache[key]=dict(theta=theta,tilt=tilt,moment_log=term,log_bound=value,witness=w)
                if len(cache)%100==0:
                    print('basin moments',len(cache),'best point bits',
                          -min(row['log_bound'] for row in cache.values())/math.log(2),flush=True)
            return cache[key]['log_bound'],cache[key]
        for tilt in np.arange(-4.,1.001,.25):
            for theta in sorted(set(np.linspace(.025,.975,39))|{.5}):
                evaluate(theta,tilt)
        coarse_best=min(cache.values(),key=lambda row:row['log_bound'])
        seeds=[]
        for row in sorted(cache.values(),key=lambda r:r['log_bound']):
            if all(abs(row['tilt']-other['tilt'])>=.35 or abs(row['theta']-other['theta'])>=.075 for other in seeds):
                seeds.append(row)
            if len(seeds)==4:
                break
        refined=[]
        for initial in seeds:
            def objective(x):
                return evaluate(float(expit(x[0])),x[1])[0]/(128*16384)
            answer=minimize(objective,[math.log(initial['theta']/(1-initial['theta'])),initial['tilt']],
                method='Nelder-Mead',options=dict(maxiter=130,xatol=1e-6,fatol=1e-10))
            _,row=evaluate(float(expit(answer.x[0])),answer.x[1])
            assert row is not None
            refined.append(dict(**row,success=bool(answer.success),evaluations=answer.nfev))
            print('refined basin point bits',-row['log_bound']/math.log(2),flush=True)
        entries=copy.deepcopy(saved.get('entries',[]))
        for row in refined:
            changed=0
            for box in dense['selected_boxes']:
                value,w=bank.proposal_bound(checker,box,row['theta'],row['tilt'],row['moment_log'])
                if value<box['own_log_bound']:
                    box.update(own_log_bound=value,witness=w)
                    changed+=1
            entries.append(dict(theta=row['theta'],tilt=row['tilt'],moment_log=row['moment_log'],changed=changed))
        dense['log_union_upper']=float(np.logaddexp.reduce([b['own_log_bound'] for b in dense['selected_boxes']]))
        dense['margin_bits']=-dense['log_union_upper']/math.log(2)
        sources=base.grid.ladder.candidate.sources()
        for path in (seed_path,replay_path):
            sources[path.relative_to(base.grid.ROOT).as_posix()]=model.base.sha(path)
        model.base.write_new(output,dict(status=scoped.STATUS,geometry=saved['geometry'],inner=saved['inner'],
            dense=dense,entries=entries,scan=dict(target=target_before,coarse_best=coarse_best,
                refined=refined,evaluations=len(cache)),full_distance_proved=False,source_sha256=sources))
        print('basin cover margin',dense['margin_bits'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('seed','replay','output'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    run(a.seed.resolve(),a.replay.resolve(),a.output.resolve())
