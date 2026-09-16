"""Adaptive dense witness search; all output is binary64, not a certificate.

Retunes both valid moment envelopes and subdivides count boxes when needed.
Always retains an exact partition, including any unresolved leaves at the cap.
"""
import argparse
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import gammaln
import general_occupancies as g

HERE=Path(__file__).resolve().parent


class Model:
    def __init__(self,record,dense,cutoff):
        self.columns=record['columns'];self.spectrum={int(w):int(n) for w,n in record['spectrum'].items() if int(w)}
        self.kernel=[int(record['kernel'].get(str(j),0)) for j in range(129)]
        self.caps=g.fiber_caps(128,19,self.spectrum,self.kernel);self.low=g.low_cancellation(self.columns,19)
        self.dense=dense;self.cutoff=cutoff
        self.counts={w:n for w,n in enumerate(g.tv.smaller_outer.spectrum()) if w and n}
        self.combinations=np.array([math.log(math.comb(128,j)) for j in range(129)])
        self.indices=np.arange(129)

    @lru_cache(maxsize=256)
    def epoch(self,z):return g.epochs(self.spectrum,self.kernel,self.columns,self.caps,self.low,math.exp(z),1,128)

    def geometry(self,box):
        corners=g.typed.vertices(box['lower'],box['upper'],32768)
        assert len(corners)>0
        bank=box['witness']['probability_bank'];ps=np.array(self.dense['probability_banks'][bank]['probabilities'])
        costs=[0.]
        for band,p in zip(self.dense['bands'],ps[1:]):
            costs.append(math.log(self.counts[128]) if band==[128] else max(
                math.log(self.counts[w])-math.log(math.comb(128,w))-w*math.log(p)-(128-w)*math.log1p(-p) for w in band))
        offsets=corners@np.array(costs)-127*(gammaln(32769)-gammaln(corners+1).sum(axis=1))
        lattice=g.typed.lattice_log_count(box['lower'],box['upper'])
        return corners,ps,offsets,lattice

    def evaluate(self,geometry,z,proposal):
        corners,ps,offsets,lattice=geometry
        p=np.asarray(proposal,dtype=float);p=p/p.sum();theta=float(ps@p)
        assert 0<theta<1 and np.all(p>0)
        probabilities=self.combinations+self.indices*math.log(theta)+(128-self.indices)*math.log1p(-theta)
        matrix=np.logaddexp.reduce(self.epoch(z)+probabilities[:,None,None],axis=0)
        ordinary=g.terminal(matrix,32768)
        direct=g.terminal(g.bernoulli_epoch(self.spectrum,self.kernel,theta,math.exp(z)),32768)
        family='occupation' if ordinary<=direct else 'bernoulli'
        value=min(ordinary,direct)+self.cutoff*math.exp(z)+max(offsets-128*(corners@np.log(p)))+lattice
        return float(value),family

    def seed(self,box):
        geometry=self.geometry(box);w=box['witness']
        value,family=self.evaluate(geometry,float(w['log_surprisal']),w['proposal'])
        return dict(lower=box['lower'],upper=box['upper'],witness=dict(w,transfer_family=family),log_bound=value,depth=box.get('depth',0))

    def optimize(self,box):
        geometry=self.geometry(box);corners=geometry[0];w=box['witness'];bank=w['probability_bank']
        best=(box['log_bound'],float(w['log_surprisal']),np.array(w['proposal']),w['transfer_family'])
        centroid=np.maximum(corners.mean(axis=0)/32768,1e-9);centroid/=centroid.sum()
        starts=[]
        grid=sorted(set(list(np.arange(-9.,.01,.5))+[float(w['log_surprisal'])]))
        for p in (best[2],centroid):
            for z in grid:
                value,family=self.evaluate(geometry,z,p);starts.append((value,z,p,family))
                if value<best[0]:best=(value,z,p,family)
        # Optimize the proposal with tilt fixed: epoch moments are cached and
        # no producer parameter changes invisibly within a final box witness.
        selected=[]
        for entry in sorted(starts,key=lambda item:item[0]):
            if entry[1] not in [v[1] for v in selected]:selected.append(entry)
            if len(selected)==2:break
        for _,z,p,_ in selected:
            def objective(x):
                logs=np.r_[x,0.];p=np.exp(logs-float(np.logaddexp.reduce(logs)))
                return self.evaluate(geometry,z,p)[0]/4194304
            optimum=minimize(objective,np.log(p[:2]/p[2]),method='L-BFGS-B',
                             bounds=[(-25.,25.),(-25.,25.)],options=dict(maxiter=60,ftol=1e-13,gtol=1e-9))
            logs=np.r_[optimum.x,0.];p=np.exp(logs-float(np.logaddexp.reduce(logs)))
            for tilt in np.arange(z-.375,z+.376,.125):
                value,family=self.evaluate(geometry,float(tilt),p)
                if value<best[0]:best=(value,float(tilt),p,family)
        value,z,p,family=best
        return dict(lower=box['lower'],upper=box['upper'],depth=box['depth'],log_bound=value,
                    witness=dict(probability_bank=bank,log_surprisal=z,proposal=p.tolist(),transfer_family=family))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nodes',type=int,default=80)
    parser.add_argument('--map',type=Path,default=HERE/'NO_CONSTANT_MAP.json')
    parser.add_argument('--distance',choices=['33/200','19/100'],default='33/200')
    parser.add_argument('--resume',action='store_true',help='continue the saved verified partition')
    args=parser.parse_args()
    record=json.loads(args.map.read_text());rows=g.tv.fixed.maps.generators(record['columns'],19)
    exact=g.tv.fixed.maps.spectrum(rows);dual=g.tv.fixed.maps.dual_spectrum(exact,128,19)
    assert exact=={int(w):int(n) for w,n in record['spectrum'].items()}
    assert dual=={int(w):int(n) for w,n in record['kernel'].items()}
    path=g.tv.fixed.HERE/'SMALLER_OUTWARD_WITNESSES.json';old=json.loads(path.read_text())
    row=next(r for r in old['results'] if r['distance_target']==args.distance)
    model=Model(record,row['dense'],row['bad_weight'])
    output=HERE/f'NO_CONSTANT_DENSE_COVER_{args.distance.replace("/","_")}.json'
    previous_nodes=0;resume_hash=None
    if args.resume:
        previous=json.loads(output.read_text())
        assert previous['distance_target']==args.distance and previous['bad_weight']==row['bad_weight']
        assert previous['map_file']==args.map.resolve().relative_to(g.tv.ROOT).as_posix()
        g.check_coverage(previous['selected_boxes'],32768,row['maximum_sparse']+1)
        leaves=[model.seed(b) for b in previous['selected_boxes']]
        previous_nodes=previous['processed_nodes'];resume_hash=hashlib.sha256(output.read_bytes()).hexdigest()
    else:leaves=[model.seed(b) for b in row['dense']['selected_boxes']]
    threshold=-80*math.log(2);processed=0
    while processed<args.nodes:
        index=max(range(len(leaves)),key=lambda i:leaves[i]['log_bound'])
        if leaves[index]['log_bound']<threshold:break
        node=model.optimize(leaves[index]);processed+=1
        print(args.distance,'node',processed,'margin',-node['log_bound']/math.log(2),'depth',node['depth'],flush=True)
        if node['log_bound']<threshold:
            leaves[index]=node;continue
        corners=g.typed.vertices(node['lower'],node['upper'],32768)
        children=g.typed.split_box(np.array(node['lower']),np.array(node['upper']),corners)
        if not children or node['depth']>=12:
            leaves[index]=node;break
        leaves.pop(index)
        for lower,upper in children:
            child=dict(lower=lower.tolist(),upper=upper.tolist(),witness=node['witness'],depth=node['depth']+1)
            leaves.append(model.seed(child))
    g.check_coverage(leaves,32768,row['maximum_sparse']+1)
    result=dict(status='BINARY64_DENSE_COVER_NOT_CERTIFICATE',distance_target=args.distance,
                bad_weight=row['bad_weight'],map_file=args.map.resolve().relative_to(g.tv.ROOT).as_posix(),
                processed_nodes=previous_nodes+processed,resumed_from_sha256=resume_hash,leaf_count=len(leaves),
                unresolved_leaves=sum(b['log_bound']>=threshold for b in leaves),
                dense_margin_bits=-float(np.logaddexp.reduce([b['log_bound'] for b in leaves]))/math.log(2),
                exact_partition_checked=True,occupation_min=row['maximum_sparse']+1,occupation_max=32768,
                bands=row['dense']['bands'],probability_banks=row['dense']['probability_banks'],selected_boxes=leaves,
                source_sha256={p.resolve().relative_to(g.tv.ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in [Path(__file__),Path(g.__file__),path,args.map]})
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('selected_boxes','source_sha256','probability_banks','bands')}),flush=True)


if __name__=='__main__':main()
