"""Outward all-occupancy certificate for the balanced-image one-round mixer.

Only fixed witnesses are imported from binary64. Reconstructs map spectra,
fiber caps, all transfers, exact cover geometry, and unions in Arb. No
supported implementation is changed. --verify replays saved bounds at 512 bits.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

from flint import arb,arb_mat,arb_poly,ctx
import general_occupancies as g
import certify_q1 as single
import certify_smaller_margins as old

HERE=Path(__file__).resolve().parent
T,S,B,L,N=128,19,128,32768,4194304
up=old.upper


def terminal(flat,n,power):
    matrix=arb_mat([list(flat[n*i:n*(i+1)]) for i in range(n)])**power
    return sum((matrix[0,j] for j in range(n)),arb(0))


class Engine:
    def __init__(self,record):
        assert (record['t'],record['s'])==(T,S) and len(record['columns'])==T
        self.columns=record['columns'];rows=g.tv.fixed.maps.generators(self.columns,S)
        exact=g.tv.fixed.maps.spectrum(rows);dual=g.tv.fixed.maps.dual_spectrum(exact,T,S)
        assert exact=={int(w):int(v) for w,v in record['spectrum'].items()}
        assert dual=={int(w):int(v) for w,v in record['kernel'].items()}
        assert len(set(self.columns))==T and all(self.columns)
        assert not any((a&b).bit_count()&1 for a in rows for b in rows)
        assert T not in exact and min(w for w in exact if w)==48 and max(exact)==80
        self.spectrum={w:v for w,v in exact.items() if w}
        self.kernel=[dual.get(j,0) for j in range(T+1)]
        assert min(j for j,v in dual.items() if j and v)==5
        self.levels=sorted(self.spectrum);self.n=len(self.levels)+2;self.m=(1<<S)-1
        self.caps=g.fiber_caps(T,S,self.spectrum,self.kernel)
        self.low=g.low_cancellation(self.columns,S)

    @lru_cache(maxsize=192)
    def epoch(self,log_lam):
        lam=old.number(log_lam).exp();z=(-lam).exp();n=self.n;m=self.m
        output=[];powers=[z**j for j in range(2*T+1)]
        def hist(items):return sum((count*powers[w] for w,count in items),arb(0))
        for j in range(T+1):
            total=math.comb(T,j);nonzero=total-self.kernel[j]
            nk=arb(nonzero)/total;r=arb(int(self.caps[j]['cap']))/total
            moments=[up(sum((math.comb(w,h)*math.comb(T-w,j-h)*powers[w+j-2*h]
                         for h in range(max(0,j-T+w),min(w,j)+1)),arb(0))/total) for w in self.levels]
            d=max(moments)
            cd=min(up(d),up(nk),up(r*powers[min(abs(w-j) for w in self.levels)]))
            if j in self.low:
                cd=min(cd,max((up(hist(pattern)/total) for pattern in self.low[j]['patterns']),default=arb(0)))
            matrix=[arb(0) for _ in range(n*n)]
            matrix[0]=arb(self.kernel[j])/total*powers[j];matrix[1]=nk*powers[j]
            matrix[n]=cd/2+min(up(d),up(nk))/(2*m);matrix[n+1]=d/2
            for i,v in enumerate(self.levels):
                matrix[n+i+2]=d*self.spectrum[v]/(2*m)
                mass=min(nonzero,self.spectrum[v]*int(self.caps[j]['cap']))
                cancel=min(up(moments[i]),up(arb(mass)*powers[abs(v-j)]/(self.spectrum[v]*total)))
                if j in self.low:
                    cancel=min(cancel,up(hist(self.low[j]['by_weight'].get(v,{}).items())/(self.spectrum[v]*total)))
                matrix[(i+2)*n]=cancel/2+min(up(moments[i]),up(nk))/(2*m)
                for k,w in enumerate(self.levels):matrix[(i+2)*n+k+2]=moments[i]*self.spectrum[w]/(2*m)
                if nonzero:matrix[(i+2)*n+1]=moments[i]/2
                else:matrix[(i+2)*n+i+2]+=moments[i]/2
            output.append(tuple(map(up,matrix)))
        return output

    def poly_mul(self,a,b,maximum):
        n=self.n
        return tuple(sum((a[n*i+k]*b[n*k+j] for k in range(n)),arb_poly()).truncate(maximum+1)
                     for i in range(n) for j in range(n))

    @lru_cache(maxsize=192)
    def region(self,log_lam,maximum):
        n=self.n;rows=self.epoch(log_lam)
        a=tuple(arb_poly([row[k]*math.comb(T,j) for j,row in enumerate(rows[:min(T,maximum)+1])]) for k in range(n*n))
        result=tuple(arb_poly([int(i==j)]) for i in range(n) for j in range(n))
        remaining=L//T
        while remaining:
            if remaining&1:result=self.poly_mul(result,a,maximum)
            remaining>>=1
            if remaining:a=self.poly_mul(a,a,maximum)
        return [tuple(max(arb(0),up(p[j]/math.comb(L,j))) for p in result) for j in range(maximum+1)]

    def bernoulli(self,theta,lam):
        z=(-lam).exp();g0=1-theta+theta*z;g1=theta+(1-theta)*z
        p0=theta*z/g0;p1=(1-theta)*z/g1
        rho0=min(arb(1),up(abs(1-2*p0)));rho1=min(arb(1),up(abs(1-2*p1)))
        n=self.n;m=self.m;entries=[]
        for v in self.levels:
            cap=1+rho1**v
            for w,count in self.spectrum.items():
                remaining=count-int(w==v)
                if not remaining:continue
                overlaps=[(v+w-d)//2 for d in self.levels if (v+w-d)%2==0
                          and max(0,v+w-T)<=(v+w-d)//2<=min(v,w)]
                assert overlaps
                cap+=remaining*max(up(rho0**(w-b)*rho1**b) for b in overlaps)
            cap/=m+1
            moment=g0**(T-v)*g1**v
            entries.append(up(moment*(cap/2+arb(1)/(2*m))))
        result=[arb(0) for _ in range(n*n)]
        result[0]=sum((self.kernel[j]*(theta*z)**j*(1-theta)**(T-j) for j in range(T+1)),arb(0))
        result[1]=sum(((math.comb(T,j)-self.kernel[j])*(theta*z)**j*(1-theta)**(T-j) for j in range(T+1)),arb(0))
        for i,entry in [(1,max(entries))]+[(i+2,entry) for i,entry in enumerate(entries)]:
            result[i*n]=entry
            for k,w in enumerate(self.levels):result[i*n+k+2]=entry*self.spectrum[w]
        return tuple(map(up,result))

    def mixture(self,rows,ps):
        probabilities=[arb(1)]
        for p in ps:
            following=[arb(0)]*(len(probabilities)+1)
            for j,v in enumerate(probabilities):following[j]+=v*(1-p);following[j+1]+=v*p
            probabilities=following
        return tuple(sum((p*rows[j][k] for j,p in enumerate(probabilities)),arb(0)) for k in range(self.n*self.n))

    def adaptive(self,rows,ps,gammas,q):
        roots=[up((gamma.log()/B).exp()) for gamma in gammas];current=rows[:q+1]
        for _ in range(q):
            current=[tuple(max(up(root*((1-p)*a[k]+p*b[k])) for root,p in zip(roots,ps))
                           for k in range(self.n*self.n)) for a,b in zip(current[:-1],current[1:])]
        return terminal(current[0],self.n,B)*math.comb(L,q)*len(ps)**q

    def dense_box(self,box,cover,counts,cutoff):
        witness=box['witness'];bank=cover['probability_banks'][witness['probability_bank']]
        ps=[old.number(v) for v in bank['probabilities']]
        assert ps[0]==0 and ps[-1]==1
        gammas=[arb(1)]+old.costs(counts,cover['bands'],ps[1:])
        raw=[F.from_float(v) for v in witness['proposal']];assert all(v>0 for v in raw)
        total=sum(raw);proposal=[old.number(v/total) for v in raw]
        theta=sum((p*v for p,v in zip(ps,proposal)),arb(0))
        log_lam=F.from_float(witness['log_surprisal']);lam=old.number(log_lam).exp()
        family=witness['transfer_family']
        if family=='occupation':
            rows=self.epoch(log_lam)
            probabilities=[math.comb(T,j)*theta**j*(1-theta)**(T-j) for j in range(T+1)]
            matrix=tuple(sum((p*rows[j][k] for j,p in enumerate(probabilities)),arb(0)) for k in range(self.n*self.n))
        else:
            assert family=='bernoulli';matrix=self.bernoulli(theta,lam)
        moment=terminal(matrix,self.n,N//T).log()
        corners=g.typed.vertices(box['lower'],box['upper'],L);values=[]
        for corner in corners:
            c=list(map(int,corner));mult=math.comb(L,c[0])*math.comb(L-c[0],c[1])
            value=cutoff*lam+moment-(B-1)*arb(mult).log()
            value+=sum((num*(gamma.log()-B*p.log()) for num,gamma,p in zip(c,gammas,proposal)),arb(0))
            values.append(up(value))
        widths=[hi-lo+1 for lo,hi in zip(box['lower'],box['upper'])]
        return max(values).exp()*(math.prod(widths)//max(widths))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();ctx.prec=512 if args.verify else 256
    receipt=HERE/'NO_CONSTANT_MARGIN_CERTIFICATE.json'
    saved=json.loads(receipt.read_text()) if args.verify else None
    if saved:
        assert saved['status']=='OUTWARD_ALL_OCCUPANCY_CERTIFICATE' and saved['precision_bits']==256
        assert saved['message_bits']==1<<20 and saved['output_bits']==N and len(saved['results'])==2
        assert saved['inner']==dict(t=T,s=S,transvection_rounds=1,map='NO_CONSTANT_MAP.json')
        for name,digest in saved['source_sha256'].items():assert g.tv.fixed.sha(g.tv.ROOT/name)==digest,name
    map_path=HERE/'NO_CONSTANT_MAP.json';map_record=json.loads(map_path.read_text());engine=Engine(map_record)
    original_path=g.tv.fixed.HERE/'SMALLER_OUTWARD_WITNESSES.json';original=json.loads(original_path.read_text())
    counts={w:n for w,n in enumerate(g.tv.smaller_outer.spectrum()) if w and n}
    construction=g.tv.smaller_outer.construction();assert construction['dimension']==32 and construction['length']==B
    if saved:assert saved['outer']==construction
    sources=[map_path,original_path,Path(__file__),Path(single.__file__),Path(old.__file__),Path(g.__file__),
             Path(g.tv.__file__),Path(g.tv.fixed.__file__),Path(g.tv.fixed.maps.__file__),Path(g.tv.fixed.outer.__file__),
             Path(g.tv.fixed.outer.bch.__file__),Path(g.tv.smaller_outer.__file__),Path(g.typed.__file__),
             g.tv.fixed.HERE/'verify_fixed_inner_results.py',g.tv.fixed.HERE/'sources/EBCH128_29.wd',g.tv.fixed.HERE/'sources/EBCH128_36.wd']
    results=[]
    for index,row in enumerate(original['results']):
        prior=saved['results'][index] if saved else None
        delta=F(row['distance_target']);target=40 if delta==F(33,200) else 30;cutoff=N*delta.numerator//delta.denominator
        assert cutoff==row['bad_weight']
        cover_path=HERE/f'NO_CONSTANT_DENSE_COVER_{str(delta).replace("/","_")}.json'
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
        zero,one,n=single.transfers(engine.spectrum,engine.columns,lam,1)
        q1=single.moments(*single.regions(zero,one,n,L//T),n,B)
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
    payload=dict(status='OUTWARD_ALL_OCCUPANCY_CERTIFICATE',precision_bits=ctx.prec,outer=construction,
                 message_bits=1<<20,output_bits=N,inner=dict(t=T,s=S,transvection_rounds=1,map='NO_CONSTANT_MAP.json'),
                 results=results,source_sha256={p.relative_to(g.tv.ROOT).as_posix():g.tv.fixed.sha(p) for p in sources})
    if saved:
        (HERE/'NO_CONSTANT_MARGIN_CERTIFICATE_REPLAY.json').write_text(json.dumps(dict(
            status='HIGHER_PRECISION_REPLAY_PASSED',precision_bits=ctx.prec,certificate_sha256=g.tv.fixed.sha(receipt),
            margins_bits=[r['margin_bits_diagnostic'] for r in results]),indent=2)+'\n',encoding='utf-8')
    else:receipt.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
