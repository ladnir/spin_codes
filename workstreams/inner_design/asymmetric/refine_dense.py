"""Bounded witness retuning for independent A/B; always retain a full cover."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import screen_dense as screen
g=screen.g
HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--name',default='greedy3_2')
    parser.add_argument('--distance',choices=['33/200','19/100'],default='19/100')
    parser.add_argument('--nodes',type=int,default=40);parser.add_argument('--resume',action='store_true')
    args=parser.parse_args();suffix=args.distance.replace('/','_')
    a_path=HERE.parent/'NO_CONSTANT_MAP.json';a_record=json.loads(a_path.read_text())
    b_path=HERE/'FEEDBACK_SCREEN.json';bank=json.loads(b_path.read_text())
    for name,digest in bank['source_sha256'].items():assert hashlib.sha256((g.tv.ROOT/name).read_bytes()).hexdigest()==digest,name
    b_record=next(r for r in bank['candidates'] if r['name']==args.name)
    original_path=HERE.parent/f'NO_CONSTANT_DENSE_COVER_{suffix}.json'
    output=HERE/f'DENSE_{args.name}_{suffix}.json'
    prior_path=output if args.resume else original_path
    cover=json.loads(prior_path.read_text());previous_hash=hashlib.sha256(prior_path.read_bytes()).hexdigest()
    if args.resume:
        assert cover['candidate']==args.name and cover['distance_target']==args.distance
        for name,digest in cover['source_sha256'].items():assert hashlib.sha256((g.tv.ROOT/name).read_bytes()).hexdigest()==digest,name
    g.check_coverage(cover['selected_boxes'],32768,cover['occupation_min'])
    model=screen.Model(a_record,b_record,cover)
    leaves=[model.seed(b) for b in cover['selected_boxes']]
    threshold=-80*math.log(2);processed=0
    while processed<args.nodes:
        index=max(range(len(leaves)),key=lambda i:leaves[i]['log_bound'])
        if leaves[index]['log_bound']<threshold:break
        node=model.optimize(leaves[index]);processed+=1
        print(args.name,args.distance,'node',processed,'margin',-node['log_bound']/math.log(2),'depth',node['depth'],flush=True)
        if node['log_bound']<threshold:leaves[index]=node;continue
        corners=g.typed.vertices(node['lower'],node['upper'],32768)
        children=g.typed.split_box(np.array(node['lower']),np.array(node['upper']),corners)
        if not children or node['depth']>=16:leaves[index]=node;break
        leaves.pop(index)
        for lower,upper in children:
            leaves.append(model.seed(dict(lower=lower.tolist(),upper=upper.tolist(),witness=node['witness'],depth=node['depth']+1)))
    g.check_coverage(leaves,32768,cover['occupation_min'])
    sources=[Path(__file__),a_path,b_path,original_path,Path(screen.__file__),Path(screen.dense_cover.__file__),Path(g.__file__)]
    result=dict(status='BINARY64_ASYMMETRIC_DENSE_COVER_NOT_CERTIFICATE',candidate=args.name,
        distance_target=args.distance,bad_weight=cover['bad_weight'],occupation_min=cover['occupation_min'],occupation_max=32768,
        processed_nodes=processed+(cover['processed_nodes'] if args.resume else 0),previous_cover_sha256=previous_hash,
        leaf_count=len(leaves),unresolved_leaves=sum(v['log_bound']>=threshold for v in leaves),
        dense_margin_bits=-float(np.logaddexp.reduce([v['log_bound'] for v in leaves]))/math.log(2),
        exact_partition_checked=True,bands=cover['bands'],probability_banks=cover['probability_banks'],selected_boxes=leaves,
        source_sha256={p.relative_to(g.tv.ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('candidate','distance_target','processed_nodes','leaf_count','unresolved_leaves','dense_margin_bits')}),flush=True)


if __name__=='__main__':main()
