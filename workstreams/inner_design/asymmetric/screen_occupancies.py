"""Reuse the generic transfer with independent A/B data; binary64 only."""
import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import search_feedback as search
g=search.g
HERE=Path(__file__).resolve().parent


def run(record,a_columns,witnesses,maximum):
    columns=record['b_columns'];s=record['s'];t=record['t']
    a_weights=search.wm.all_weights(a_columns,s)
    a_spectrum={int(w):int(n) for w,n in zip(*np.unique(a_weights,return_counts=True)) if w}
    dual={int(w):int(n) for w,n in record['dual_spectrum'].items() if int(w)}
    kernel=[record['kernel_spectrum'].get(str(j),0) for j in range(t+1)]
    caps=g.fiber_caps(t,s,dual,kernel)
    low=search.low_cancellation(a_columns,columns,s)
    @lru_cache(maxsize=128)
    def at(z,q):
        epoch=g.epochs(a_spectrum,kernel,columns,caps,low,math.exp(z),1,q)
        return g.regions(epoch,t,32768,q)
    counts={w:n for w,n in enumerate(g.tv.smaller_outer.spectrum()) if w and n}
    results=[]
    for row in witnesses['results']:
        pairs=[];higher=[]
        for witness in row['q2']:
            z=float(g.F(witness['log_tilt']));ids=witness['band_indices']
            roots,active,inactive=g.tv.fixed.general.density_roots(counts,128,witnesses['q2_bands'],float(g.F(witness['shift'])))
            matrix=g.mixture(at(z,2),active[ids],inactive[ids])
            value=g.terminal(matrix,128)+128*sum(roots[ids])+math.log(math.comb(32768,2))
            if ids[0]!=ids[1]:value+=math.log(2)
            pairs.append(value+row['bad_weight']*math.exp(z))
        result=dict(distance_target=row['distance_target'],q2_margin_bits=-float(np.logaddexp.reduce(pairs))/math.log(2),higher=[])
        for witness in row['adaptive']:
            q=witness['occupation']
            if q>maximum:break
            z=float(g.F(witness['log_tilt']))
            roots,active,inactive=g.tv.fixed.general.density_roots(counts,128,witnesses['adaptive_bands'],float(g.F(witness['shift'])))
            value=g.adaptive(at(z,q),roots,active,inactive,q)+row['bad_weight']*math.exp(z)
            higher.append(value);result['higher'].append(dict(occupation=q,margin_bits=-value/math.log(2)))
            if q%16==0:print(record['name'],row['distance_target'],'sparse',q,flush=True)
        result['q2_through_maximum_margin_bits']=-float(np.logaddexp.reduce(pairs+higher))/math.log(2)
        q1=next(v['margin_bits'] for v in record['q1'] if v['distance_target']==row['distance_target'])
        result['partial_union_margin_bits']=-float(np.logaddexp.reduce([-q1*math.log(2)]+pairs+higher))/math.log(2)
        result['covered_occupations']=[1,min(maximum,row['maximum_sparse'])]
        results.append(result)
        print(record['name'],json.dumps(result),flush=True)
    return dict(name=record['name'],results=results)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--maximum',type=int,default=8)
    parser.add_argument('--names',nargs='+',default=['symmetric_control','mixed2_3','greedy3_2'])
    args=parser.parse_args();assert 3<=args.maximum<=128
    path=HERE/'FEEDBACK_SCREEN.json';screen=json.loads(path.read_text())
    for name,digest in screen['source_sha256'].items():assert hashlib.sha256((g.tv.ROOT/name).read_bytes()).hexdigest()==digest,name
    map_path=HERE.parent/'NO_CONSTANT_MAP.json';a_columns=json.loads(map_path.read_text())['columns']
    witness_path=g.tv.fixed.HERE/'SMALLER_OUTWARD_WITNESSES.json';witnesses=json.loads(witness_path.read_text())
    candidates={r['name']:r for r in screen['candidates']}
    results=[run(candidates[name],a_columns,witnesses,args.maximum) for name in args.names]
    sources=[Path(__file__),path,map_path,witness_path,Path(search.__file__),Path(g.__file__),Path(g.tv.fixed.general.__file__)]
    payload=dict(status='BINARY64_PARTIAL_OCCUPANCIES_NOT_CERTIFICATE',maximum_occupation=args.maximum,
                 original_witnesses_not_retuned=True,results=results,
                 source_sha256={p.relative_to(g.tv.ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (HERE/f'OCCUPANCIES_Q{args.maximum}.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
