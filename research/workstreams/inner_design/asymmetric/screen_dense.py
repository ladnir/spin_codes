"""Full dense coverage for independent A/B, using fixed saved witnesses.

Both candidate transfers are valid for independent maps. In particular the
Fourier cap below never uses A(q+a) to constrain B^T a intersections.
Binary64 diagnostics only; no outward certificate is produced here.
"""
import argparse
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import search_feedback as search
import dense_cover
g=search.g
HERE=Path(__file__).resolve().parent


def bernoulli(a_spectrum,b_spectrum,kernel,theta,lam):
    assert 0<theta<1 and lam>0
    t=len(kernel)-1;m=sum(a_spectrum.values());assert sum(b_spectrum.values())==m
    levels=sorted(a_spectrum);z=math.exp(-lam)
    g0=1-theta+theta*z;g1=theta+(1-theta)*z
    rho0=abs(1-2*theta*z/g0);rho1=abs(1-2*(1-theta)*z/g1)
    entries=[]
    for v in levels:
        cap=1.
        for w,count in b_spectrum.items():
            overlap=min(v,w) if rho1>=rho0 else max(0,v+w-t)
            cap+=count*rho0**(w-overlap)*rho1**overlap
        cap/=m+1
        moment=(t-v)*math.log(g0)+v*math.log(g1)
        entries.append(moment+math.log(cap/2+1/(2*m)))
    n=len(levels)+2;matrix=np.full((n,n),-math.inf)
    zero=[];live=[]
    for j,k in enumerate(kernel):
        log_weight=j*(math.log(theta)-lam)+(t-j)*math.log1p(-theta)
        if k:zero.append(math.log(k)+log_weight)
        if math.comb(t,j)>k:live.append(math.log(math.comb(t,j)-k)+log_weight)
    matrix[0,0]=np.logaddexp.reduce(zero)
    if live:matrix[0,1]=np.logaddexp.reduce(live)
    for i,value in [(1,max(entries))]+[(i+2,v) for i,v in enumerate(entries)]:
        matrix[i,0]=value
        for k,w in enumerate(levels):matrix[i,k+2]=value+math.log(a_spectrum[w])
    return matrix


class Model(dense_cover.Model):
    def __init__(self,a_record,b_record,cover):
        super().__init__(a_record,cover,cover['bad_weight'])
        self.columns=b_record['b_columns']
        self.kernel=[b_record['kernel_spectrum'].get(str(j),0) for j in range(129)]
        self.b_spectrum={int(w):int(n) for w,n in b_record['dual_spectrum'].items() if int(w)}
        self.caps=g.fiber_caps(128,19,self.b_spectrum,self.kernel)
        self.low=search.low_cancellation(a_record['columns'],self.columns,19)

    def evaluate(self,geometry,z,proposal):
        corners,ps,offsets,lattice=geometry
        p=np.array(proposal,dtype=float);p/=p.sum();theta=float(ps@p)
        assert 0<theta<1 and np.all(p>0)
        probabilities=self.combinations+self.indices*math.log(theta)+(128-self.indices)*math.log1p(-theta)
        ordinary=g.terminal(np.logaddexp.reduce(self.epoch(z)+probabilities[:,None,None],axis=0),32768)
        direct=g.terminal(bernoulli(self.spectrum,self.b_spectrum,self.kernel,theta,math.exp(z)),32768)
        family='occupation' if ordinary<=direct else 'asymmetric_bernoulli'
        value=min(ordinary,direct)+self.cutoff*math.exp(z)+max(offsets-128*(corners@np.log(p)))+lattice
        return float(value),family


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--names',nargs='+',default=['mixed2_3','greedy3_2'])
    args=parser.parse_args()
    a_path=HERE.parent/'NO_CONSTANT_MAP.json';a_record=json.loads(a_path.read_text())
    b_path=HERE/'FEEDBACK_SCREEN.json';screen=json.loads(b_path.read_text())
    for name,digest in screen['source_sha256'].items():assert hashlib.sha256((g.tv.ROOT/name).read_bytes()).hexdigest()==digest,name
    candidates={r['name']:r for r in screen['candidates']};results=[]
    sources=[Path(__file__),a_path,b_path,Path(search.__file__),Path(g.__file__),Path(dense_cover.__file__)]
    for name in args.names:
        row=dict(name=name,results=[])
        for suffix in ('33_200','19_100'):
            cover_path=HERE.parent/f'NO_CONSTANT_DENSE_COVER_{suffix}.json';sources.append(cover_path)
            cover=json.loads(cover_path.read_text());g.check_coverage(cover['selected_boxes'],32768,cover['occupation_min'])
            model=Model(a_record,candidates[name],cover);logs=[];families=[]
            for i,box in enumerate(cover['selected_boxes']):
                w=box['witness'];value,family=model.evaluate(model.geometry(box),w['log_surprisal'],w['proposal'])
                logs.append(value);families.append(family)
                if i%64==0:print(name,suffix,'dense',i+1,'/',len(cover['selected_boxes']),flush=True)
            result=dict(distance_target=cover['distance_target'],covered_occupations=[cover['occupation_min'],32768],
                        exact_partition_checked=True,dense_margin_bits=-float(np.logaddexp.reduce(logs))/math.log(2),
                        box_log_bounds=logs,transfer_families=families,boxes_above_unit_bound=sum(v>=0 for v in logs))
            row['results'].append(result)
            print(name,suffix,'dense margin',result['dense_margin_bits'],'unit failures',result['boxes_above_unit_bound'],flush=True)
        results.append(row)
    payload=dict(status='BINARY64_DENSE_ONLY_NOT_CERTIFICATE',original_witnesses_not_retuned=True,results=results,
        source_sha256={p.relative_to(g.tv.ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    (HERE/'DENSE_SCREEN.json').write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
