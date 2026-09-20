"""Add the exact fair-input causal anchor to a replayed IMT dense suffix."""
import argparse
import copy
import math
from pathlib import Path

import numpy as np
import parameter_scoped_cover as scoped
import parameter_causal_dense as causal


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
        total=checker.block*checker.length
        lam=math.log((total-checker.cutoff)/checker.cutoff)
        tilt=math.log(lam)
        term=checker.moment(.5,tilt)
        exact_causal=causal.causal_moment(.5,lam,total)
        assert abs(term-exact_causal)<1e-7
        dense=copy.deepcopy(saved['dense'])
        changed=0
        for box in dense['selected_boxes']:
            value,witness=bank.proposal_bound(checker,box,.5,tilt,term)
            if value<box['own_log_bound']:
                box.update(own_log_bound=value,witness=witness)
                changed+=1
        dense['log_union_upper']=float(np.logaddexp.reduce([b['own_log_bound'] for b in dense['selected_boxes']]))
        dense['margin_bits']=-dense['log_union_upper']/math.log(2)
        sources=base.grid.ladder.candidate.sources()
        for path in (seed_path,replay_path):
            sources[path.relative_to(base.grid.ROOT).as_posix()]=model.base.sha(path)
        model.base.write_new(output,dict(status=scoped.STATUS,geometry=saved['geometry'],inner=saved['inner'],
            dense=dense,entries=[*saved.get('entries',[]),dict(theta=.5,tilt=tilt,moment_log=term,changed=changed)],
            fair_input_causal_log=exact_causal,full_distance_proved=False,source_sha256=sources))
        print('uniform anchor improved',changed,'boxes, dense bits',dense['margin_bits'],
              'weak boxes',sum(b['own_log_bound']>-60*math.log(2) for b in dense['selected_boxes']),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('seed','replay','output'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    run(a.seed.resolve(),a.replay.resolve(),a.output.resolve())
