"""Bisect weak dense-suffix boxes using retained complete-moment witnesses."""
import argparse
import copy
import heapq
import math
from pathlib import Path

import numpy as np
import parameter_scoped_cover as scoped
import type_box_coverage


def run(seed_path,replay_path,output,maximum):
    bank=scoped.bank
    base,model=bank.base,bank.model
    saved,replay=[model.base.read(p) for p in (seed_path,replay_path)]
    model.authenticate(saved); model.authenticate(replay)
    assert saved['status']==scoped.STATUS and replay['status']=='VERIFIED_'+scoped.STATUS
    assert replay['producer_sha256']==model.base.sha(seed_path)
    with bank.previous.maps.use():
        moments=scoped.bridge.Moments()
        scoped.check(saved,moments)
        dense=copy.deepcopy(saved['dense'])
        boxes={i:b for i,b in enumerate(dense['selected_boxes'])}
        queue=[(-b['own_log_bound'],i) for i,b in boxes.items()]
        heapq.heapify(queue)
        witnesses=[]
        for row in saved.get('entries',[]):
            theta,tilt=row['theta'],row['tilt']
            if len(moments.activation.epoch_cache)>4:
                moments.activation.epoch_cache.clear()
            term=moments.activation.moment(theta,tilt)
            assert abs(term-row['moment_log'])<1e-7
            witnesses.append((theta,tilt,term))
        next_id=len(boxes)
        singleton_count=0
        used=0
        while used<maximum and queue:
            _,index=heapq.heappop(queue)
            parent=boxes[index]
            if parent['own_log_bound'] < -60*math.log(2):
                break
            children=scoped.partition.split_box(parent,16384)
            if not children:
                singleton_count+=1
                continue
            del boxes[index]
            for child in children:
                child['witness']=copy.deepcopy(parent['witness'])
                value=moments.bound(child)
                for theta,tilt,term in witnesses:
                    trial,w=bank.proposal_bound(moments.original,child,theta,tilt,term)
                    if trial<value:
                        value=trial
                        child['witness']=w
                child['own_log_bound']=value
                boxes[next_id]=child
                heapq.heappush(queue,(-value,next_id))
                next_id+=1
            used+=1
            if used%64==0:
                log_union=float(np.logaddexp.reduce([b['own_log_bound'] for b in boxes.values()]))
                print('scoped splits',used,'margin',-log_union/math.log(2),'weak singletons',singleton_count,flush=True)
        dense['selected_boxes']=list(boxes.values())
        count=type_box_coverage.check(dense['selected_boxes'],16384,dense['occupation_min'],3)
        dense['integer_types_checked']=str(count)
        dense['log_union_upper']=float(np.logaddexp.reduce([b['own_log_bound'] for b in boxes.values()]))
        dense['margin_bits']=-dense['log_union_upper']/math.log(2)
        sources=base.grid.ladder.candidate.sources()
        for path in (seed_path,replay_path):
            sources[path.relative_to(base.grid.ROOT).as_posix()]=model.base.sha(path)
        model.base.write_new(output,dict(status=scoped.STATUS,geometry=saved['geometry'],inner=saved['inner'],
            dense=dense,entries=saved.get('entries',[]),refinement_splits=used,weak_singletons=singleton_count,
            full_distance_proved=False,source_sha256=sources))
        print('scoped partition',len(boxes),'boxes, margin',dense['margin_bits'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('seed','replay','output'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--splits',type=int,default=256)
    a=p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    run(a.seed.resolve(),a.replay.resolve(),a.output.resolve(),a.splits)
