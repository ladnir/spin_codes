"""Outward all-occupancy certificate for independent expansion and feedback.

Fixed instance: BCH [128,32,32], K=2^20, balanced A, greedy3_2 B, one
transvection per epoch. --verify replays at 512 bits against 256-bit bounds.
"""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path

import search_feedback as search
import certify_no_constant as base
from flint import arb,ctx

g=search.g
old=base.old
single=base.single
up=base.up
terminal=base.terminal
T,S,B,L,N=128,19,128,32768,4194304
HERE=Path(__file__).resolve().parent
INNER=dict(t=T,s=S,transvection_rounds=1,expansion_map='../NO_CONSTANT_MAP.json',
           feedback_map='FEEDBACK_SCREEN.json',feedback_name='greedy3_2')


class Engine(base.Engine):
    def __init__(self,a_record):
        assert (a_record['t'],a_record['s'])==(T,S)
        self.a_columns=list(a_record['columns'])
        b_record=next(v for v in json.loads((HERE/'FEEDBACK_SCREEN.json').read_text())['candidates'] if v['name']=='greedy3_2')
        self.columns=list(b_record['b_columns'])
        assert len(self.a_columns)==len(self.columns)==T
        assert search.rank(self.a_columns)==search.rank(self.columns)==S
        assert all(self.columns) and len(set(self.columns))==T
        assert all(v.bit_count()==3 for v in self.columns)
        a_exact=g.tv.fixed.maps.spectrum(g.tv.fixed.maps.generators(self.a_columns,S))
        b_exact=g.tv.fixed.maps.spectrum(g.tv.fixed.maps.generators(self.columns,S))
        dual=g.tv.fixed.maps.dual_spectrum(b_exact,T,S)
        assert a_exact=={int(w):int(v) for w,v in a_record['spectrum'].items()}
        assert b_exact=={int(w):int(v) for w,v in b_record['dual_spectrum'].items()}
        assert dual=={int(w):int(v) for w,v in b_record['kernel_spectrum'].items()}
        assert a_exact[0]==b_exact[0]==1
        self.spectrum={w:v for w,v in a_exact.items() if w}
        self.b_spectrum={w:v for w,v in b_exact.items() if w}
        self.kernel=[dual.get(j,0) for j in range(T+1)]
        assert self.kernel[1:5]==[0,0,0,291]
        self.levels=sorted(self.spectrum)
        self.n=len(self.levels)+2
        self.m=(1<<S)-1
        self.caps=g.fiber_caps(T,S,self.b_spectrum,self.kernel)
        self.low=search.low_cancellation(self.a_columns,self.columns,S)
        for j in (1,2):assert self.low[j]['nonzero_input_count']==math.comb(T,j)

    def bernoulli(self,theta,lam):
        assert theta>0 and theta<1 and lam>0
        z=(-lam).exp();g0=1-theta+theta*z;g1=theta+(1-theta)*z
        rho0=min(arb(1),up(abs(1-2*theta*z/g0)))
        rho1=min(arb(1),up(abs(1-2*(1-theta)*z/g1)))
        n=self.n;m=self.m;entries=[]
        for v in self.levels:
            cap=arb(1)
            for w,count in self.b_spectrum.items():
                # No A(q+a) identity: evaluate both feasible overlap endpoints.
                ends=(max(0,v+w-T),min(v,w))
                cap+=count*max(up(rho0**(w-h)*rho1**h) for h in ends)
            moment=g0**(T-v)*g1**v
            entries.append(up(moment*(cap/(2*(m+1))+arb(1)/(2*m))))
        result=[arb(0) for _ in range(n*n)]
        result[0]=sum((self.kernel[j]*(theta*z)**j*(1-theta)**(T-j) for j in range(T+1)),arb(0))
        result[1]=sum(((math.comb(T,j)-self.kernel[j])*(theta*z)**j*(1-theta)**(T-j) for j in range(T+1)),arb(0))
        for i,entry in [(1,max(entries))]+[(i+2,entry) for i,entry in enumerate(entries)]:
            result[i*n]=entry
            for k,w in enumerate(self.levels):result[i*n+k+2]=entry*self.spectrum[w]
        return tuple(map(up,result))

    def dense_box(self,box,cover,counts,cutoff):
        # The inherited routine selects one complete transfer per witness.
        family=box['witness']['transfer_family']
        assert family in ('occupation','asymmetric_bernoulli')
        if family=='asymmetric_bernoulli':
            box=dict(box,witness=dict(box['witness'],transfer_family='bernoulli'))
        return super().dense_box(box,cover,counts,cutoff)

    def q1(self,log_lam):
        # Generic fixed-j bounds keep A and B separate, including exact j=1
        # cancellation. Never invoke the symmetric single.transfers shortcut.
        zero,one=self.epoch(log_lam)[:2]
        return single.moments(*single.regions(zero,one,self.n,L//T),self.n,B)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();ctx.prec=512 if args.verify else 256
    receipt=HERE/'ASYMMETRIC_MARGIN_CERTIFICATE.json'
    saved=json.loads(receipt.read_text()) if args.verify else None
    if saved:
        assert saved['status']=='OUTWARD_ASYMMETRIC_ALL_OCCUPANCY_CERTIFICATE' and saved['precision_bits']==256
        assert saved['message_bits']==1<<20 and saved['output_bits']==N and len(saved['results'])==2
        assert saved['inner']==INNER
        for name,digest in saved['source_sha256'].items():assert g.tv.fixed.sha(g.tv.ROOT/name)==digest,name
    map_path=HERE.parent/'NO_CONSTANT_MAP.json';map_record=json.loads(map_path.read_text());engine=Engine(map_record)
    original_path=g.tv.fixed.HERE/'SMALLER_OUTWARD_WITNESSES.json';original=json.loads(original_path.read_text())
    counts={w:n for w,n in enumerate(g.tv.smaller_outer.spectrum()) if w and n}
    construction=g.tv.smaller_outer.construction();assert construction['dimension']==32 and construction['length']==B
    if saved:assert saved['outer']==construction
    sources=[map_path,original_path,Path(__file__),Path(single.__file__),Path(old.__file__),Path(g.__file__),Path(base.__file__),
             Path(search.__file__),Path(search.wm.__file__),HERE/'FEEDBACK_SCREEN.json',
             Path(g.tv.__file__),Path(g.tv.fixed.__file__),Path(g.tv.fixed.maps.__file__),Path(g.tv.fixed.outer.__file__),
             Path(g.tv.fixed.outer.bch.__file__),Path(g.tv.smaller_outer.__file__),Path(g.typed.__file__),
             g.tv.fixed.HERE/'verify_fixed_inner_results.py',g.tv.fixed.HERE/'sources/EBCH128_29.wd',g.tv.fixed.HERE/'sources/EBCH128_36.wd']
    results=[]
    for index,row in enumerate(original['results']):
        prior=saved['results'][index] if saved else None
        delta=F(row['distance_target']);target=40 if delta==F(33,200) else 30;cutoff=N*delta.numerator//delta.denominator
        assert cutoff==row['bad_weight']
        cover_path=HERE/f'DENSE_greedy3_2_{str(delta).replace("/","_")}.json'
        cover=json.loads(cover_path.read_text());sources.append(cover_path)
        assert cover['distance_target']==str(delta) and cover['bad_weight']==cutoff
        assert cover['occupation_min']==row['maximum_sparse']+1 and cover['occupation_max']==L
        g.check_coverage(cover['selected_boxes'],L,row['maximum_sparse']+1)
        if prior:
            assert prior['distance_target']==str(delta) and prior['bad_weight']==cutoff and prior['target_bits']==target
            assert prior['maximum_sparse']==row['maximum_sparse'] and prior['exact_occupation_coverage_checked'] is True
            assert len(prior['q2_terms'])==len(row['q2']) and len(prior['adaptive_terms'])==len(row['adaptive'])
            assert len(prior['dense_terms'])==len(cover['selected_boxes'])
        def record(value,previous=None):
            if previous is not None:
                assert up(value)<=old.unpack(previous),'512-bit replay exceeds saved bound'
                return previous
            return old.pack(value)
        z=F(-939,100) if target==40 else F(-966,100);lam=old.number(z).exp()
        n=engine.n
        q1=engine.q1(z)
        value=L*(cutoff*lam).exp()*sum((num*q1[w] for w,num in counts.items()),arb(0))
        q1_bound=record(value,prior['q1_upper'] if prior else None)
        pairs=[]
        expected=[list(map(int,c)) for c in g.tv.fixed.composition.compositions(len(original['q2_bands']),2)]
        assert [w['band_indices'] for w in row['q2']]==expected
        for i,witness in enumerate(row['q2']):
            z=F(witness['log_tilt']);ps,gs=old.parameters(counts,original['q2_bands'],old.number(witness['shift']))
            ids=witness['band_indices'];matrix=engine.mixture(engine.region(z,2),[ps[i] for i in ids])
            value=terminal(matrix,n,B)*gs[ids[0]]*gs[ids[1]]*(2 if ids[0]!=ids[1] else 1)*math.comb(L,2)*(cutoff*old.number(z).exp()).exp()
            pairs.append(record(value,prior['q2_terms'][i] if prior else None))
        sparse=[]
        assert [w['occupation'] for w in row['adaptive']]==list(range(3,row['maximum_sparse']+1))
        for i,witness in enumerate(row['adaptive']):
            z=F(witness['log_tilt']);q=witness['occupation']
            ps,gs=old.parameters(counts,original['adaptive_bands'],old.number(witness['shift']))
            value=engine.adaptive(engine.region(z,q),ps,gs,q)*(cutoff*old.number(z).exp()).exp()
            sparse.append(record(value,prior['adaptive_terms'][i] if prior else None))
            if q%16==0:print(ctx.prec,str(delta),'sparse',q,flush=True)
        dense=[]
        for i,box in enumerate(cover['selected_boxes']):
            value=engine.dense_box(box,cover,counts,cutoff)
            dense.append(record(value,prior['dense_terms'][i] if prior else None))
            if i%32==0:print(ctx.prec,str(delta),'dense',i+1,'/',len(cover['selected_boxes']),flush=True)
        union=old.unpack(q1_bound)+sum((old.unpack(v) for v in pairs+sparse+dense),arb(0))
        union_bound=record(union,prior['union_upper'] if prior else None)
        assert old.unpack(union_bound)<arb(2)**(-target),'target not certified'
        result=dict(distance_target=str(delta),bad_weight=cutoff,target_bits=target,q1_upper=q1_bound,
                    q2_terms=pairs,adaptive_terms=sparse,dense_terms=dense,maximum_sparse=row['maximum_sparse'],
                    exact_occupation_coverage_checked=True,union_upper=union_bound,
                    margin_bits_diagnostic=old.margin(old.unpack(union_bound)))
        results.append(result);print('PASS',ctx.prec,str(delta),result['margin_bits_diagnostic'],flush=True)
    payload=dict(status='OUTWARD_ASYMMETRIC_ALL_OCCUPANCY_CERTIFICATE',precision_bits=ctx.prec,outer=construction,
                 message_bits=1<<20,output_bits=N,inner=INNER,
                 results=results,source_sha256={p.relative_to(g.tv.ROOT).as_posix():g.tv.fixed.sha(p) for p in sources})
    if saved:
        (HERE/'ASYMMETRIC_MARGIN_CERTIFICATE_REPLAY.json').write_text(json.dumps(dict(
            status='HIGHER_PRECISION_REPLAY_PASSED',precision_bits=ctx.prec,certificate_sha256=g.tv.fixed.sha(receipt),
            margins_bits=[r['margin_bits_diagnostic'] for r in results]),indent=2)+'\n',encoding='utf-8')
    else:receipt.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
