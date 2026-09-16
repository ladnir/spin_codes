"""All-epoch-occupancy envelope and sparse binary64 witness replay.

Exact integer syndrome-fiber caps come from packing, Fourier triangle and
Parseval. This file is diagnostic: no binary64 result is a certificate.
"""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import transvection as tv
import weight_memory as wm
import typed_dense_boxes as typed
from verify_fixed_inner_results import check_coverage

HERE=Path(__file__).resolve().parent


def krawtchouk(t,w,j):
    return sum((-1)**i*math.comb(w,i)*math.comb(t-w,j-i)
               for i in range(max(0,j-t+w),min(w,j)+1))


def fiber_caps(t,s,spectrum,kernel):
    """Integer upper bounds on each nonzero syndrome fiber at fixed weight.

    Parseval is centered on the NONZERO syndromes; the known zero fiber is
    subtracted exactly before applying the finite-population maximum bound.
    """
    m=(1<<s)-1
    assert sum(spectrum.values())==m and len(kernel)==t+1
    assert kernel[0]==1 and sum(kernel)==1<<(t-s)
    distance=next((j for j in range(1,t+1) if kernel[j]),t+1)
    radius=(distance-1)//2
    result=[]
    for j in range(t+1):
        total=math.comb(t,j);zero=kernel[j];nonzero=total-zero
        transforms=[(n,krawtchouk(t,w,j)) for w,n in spectrum.items()]
        numerator=total*total+sum(n*k*k for n,k in transforms)
        assert numerator%(m+1)==0
        squares=numerator//(m+1)
        radicand=(m-1)*(m*(squares-zero*zero)-nonzero*nonzero)
        assert radicand>=0
        root=math.isqrt(radicand)
        if root*root<radicand:root+=1
        parseval=(nonzero+root)//m
        triangle=(total+sum(n*abs(k) for n,k in transforms))//(m+1)
        h=min(j,t-j)
        packing=1 if h<radius else math.comb(t,h-radius)//math.comb(h,radius)
        cap=min(nonzero,packing,triangle,parseval)
        assert cap>=0
        result.append(dict(weight=j,total=str(total),zero=str(zero),cap=str(cap),
                           packing=str(packing),fourier_triangle=str(triangle),parseval=str(parseval)))
    return result


def low_cancellation(columns,s,maximum=2):
    """Exact input/syndrome/output histograms, with no injectivity assumption."""
    t=len(columns);weights=wm.all_weights(columns,s)
    images=[sum(tv.dot(col,q)<<p for p,col in enumerate(columns)) for q in columns]
    result={}
    for j in range(1,maximum+1):
        by_syndrome=defaultdict(Counter);by_weight=defaultdict(Counter)
        for support in itertools.combinations(range(t),j):
            syndrome=0;emitted=0
            for p in support:syndrome^=columns[p];emitted^=images[p]^(1<<p)
            if not syndrome:continue
            out_weight=emitted.bit_count()
            by_syndrome[syndrome][out_weight]+=1
            by_weight[int(weights[syndrome])][out_weight]+=1
        patterns=sorted({tuple(sorted(hist.items())) for hist in by_syndrome.values()})
        result[j]=dict(patterns=patterns,by_weight=dict(by_weight),
                       nonzero_input_count=sum(sum(v.values()) for v in by_syndrome.values()))
    return result


def exceptional_coset(columns,s,spectrum):
    """Exact B-syndrome enumerator for the unique all-one A-image state."""
    t=len(columns);weights=wm.all_weights(columns,s)
    candidates=np.flatnonzero(weights==t)
    assert len(candidates)==1
    h=int(candidates[0]);states=np.arange(1<<s,dtype=np.uint32)
    signs=1-2*(np.bitwise_count(states & h).astype(np.int32)&1)
    signed={w:int(signs[weights==w].sum()) for w in range(t+1) if np.any(weights==w)}
    coset=[]
    for j in range(t+1):
        numerator=sum(n*krawtchouk(t,w,j) for w,n in signed.items())
        assert numerator%(1<<s)==0 and numerator>=0
        coset.append(numerator//(1<<s))
    assert sum(coset)==1<<(t-s)
    return h,coset


def logsum_hist(hist,lam):
    terms=[math.log(n)-lam*w for w,n in hist]
    return float(np.logaddexp.reduce(terms)) if terms else -math.inf


def epochs(spectrum,kernel,columns,caps,low,lam,rounds=1,maximum=128):
    """Row/column order Z,D,C_48,...,C_128; normalized shell mass caps."""
    levels=sorted(spectrum);s=(sum(spectrum.values())+1).bit_length()-1
    m=(1<<s)-1;t=len(columns);n=len(levels)+2
    assert 0<lam and (rounds is None or 1<=rounds) and 0<=maximum<=t
    # None selects the original exact-refresh mixer, not a zero-round product.
    epsilon=0. if rounds is None else 2.**-rounds
    le=math.log(epsilon) if epsilon else -math.inf;lf=math.log1p(-epsilon)
    result=np.full((maximum+1,n,n),-np.inf)
    for j in range(maximum+1):
        total=math.comb(t,j);log_total=math.log(total)
        nonzero=total-kernel[j]
        nk=math.log(nonzero)-log_total if nonzero else -math.inf
        beta=math.log(kernel[j])-log_total if kernel[j] else -math.inf
        moments=[]
        for w in levels:
            terms=[math.log(math.comb(w,h)*math.comb(t-w,j-h))-log_total-lam*(w+j-2*h)
                   for h in range(max(0,j-t+w),min(w,j)+1)]
            moments.append(float(np.logaddexp.reduce(terms)))
        arbitrary=max(moments);matrix=result[j]
        matrix[0,0]=beta-lam*j;matrix[0,1]=nk-lam*j
        cap=int(caps[j]['cap'])
        r=math.log(cap)-log_total if cap else -math.inf
        cd=min(arbitrary,nk,r-lam*min(abs(w-j) for w in levels))
        if j in low:
            cd=min(cd,max((logsum_hist(hist,lam)-log_total for hist in low[j]['patterns']),default=-math.inf))
        matrix[1,0]=np.logaddexp(le+cd,lf+min(arbitrary,nk)-math.log(m))
        matrix[1,1]=le+arbitrary
        for i,v in enumerate(levels):
            matrix[1,i+2]=lf+arbitrary+math.log(spectrum[v]/m)
            moment=moments[i]
            mass=min(nonzero,spectrum[v]*cap)
            c=min(moment,math.log(mass)-math.log(spectrum[v])-log_total-lam*abs(v-j)) if mass else -math.inf
            if j in low:
                hist=low[j]['by_weight'].get(v,{})
                c=min(c,logsum_hist(hist.items(),lam)-math.log(spectrum[v])-log_total)
            matrix[i+2,0]=np.logaddexp(le+c,lf+min(moment,nk)-math.log(m))
            for k,w in enumerate(levels):matrix[i+2,k+2]=lf+moment+math.log(spectrum[w]/m)
            if nonzero:
                matrix[i+2,1]=le+moment
            else:
                # Empty and full kernel inputs preserve the lazy state's shell.
                matrix[i+2,i+2]=np.logaddexp(matrix[i+2,i+2],le+moment)
    return result


def poly_mul(a,b,maximum):
    n=a.shape[1];limit=min(maximum,len(a)+len(b)-2)
    result=np.full((limit+1,n,n),-np.inf)
    for degree in range(limit+1):
        lo=max(0,degree-len(b)+1);hi=min(degree,len(a)-1)
        result[degree]=np.logaddexp.reduce(wm.product(a[lo:hi+1],b[degree-hi:degree-lo+1][::-1]),axis=0)
    return result


def regions(epoch,t,length,maximum):
    n=epoch.shape[1]
    current=np.full((1,n,n),-np.inf)
    for i in range(n):current[0,i,i]=0
    power=epoch[:maximum+1]+np.array([math.log(math.comb(t,j)) for j in range(min(maximum,t)+1)])[:,None,None]
    exponent=length//t
    assert length%t==0
    while exponent:
        if exponent&1:current=poly_mul(current,power,maximum)
        exponent>>=1
        if exponent:power=poly_mul(power,power,maximum)
    return current-np.array([math.log(math.comb(length,j)) for j in range(maximum+1)])[:,None,None]


def terminal(matrix,length):
    n=matrix.shape[0];current=np.full((1,n,n),-np.inf)
    for i in range(n):current[0,i,i]=0
    power=matrix[None]
    while length:
        if length&1:current=wm.product(current,power)
        length>>=1
        if length:power=wm.product(power,power)
    return float(np.logaddexp.reduce(current[0,0]))


def bernoulli_epoch(spectrum,kernel,theta,lam,rounds=1):
    """Direct tilted-syndrome smoothing for the dense Bernoulli surrogate.

    Weighting iid input X by z^wt(X+Aq) preserves independent output bits.
    Their parity biases are bounded by rho; Fourier inversion of B bounds
    every weighted syndrome atom. This handles lazy feedback without first
    forgetting its input-induced smoothing.
    """
    assert 0<theta<1 and lam>0
    z=math.exp(-lam);t=len(kernel)-1;m=sum(spectrum.values());levels=sorted(spectrum)
    g0=1-theta+theta*z;g1=theta+(1-theta)*z
    p0=theta*z/g0;p1=(1-theta)*z/g1
    rho0,rho1=abs(1-2*p0),abs(1-2*p1)
    eps=2.**-rounds
    caps=[]
    for v in levels:
        # A(q+a) constrains the overlap of A(q) and A(a). Account separately
        # for a=q, and for a=q+h if the all-one image h exists. All remaining
        # differences have a nonzero, non-all-one weight from this spectrum.
        total=1.+rho1**v
        has_complement=t in spectrum and v!=t
        if has_complement:total+=rho0**(t-v)
        distances=[d for d in levels if d!=t]
        for w,count in spectrum.items():
            remaining=count-int(w==v)-int(has_complement and w==t-v)
            if not remaining:continue
            overlaps=[(v+w-d)//2 for d in distances if (v+w-d)%2==0 and 0<=(v+w-d)//2<=min(v,w)
                      and (v+w-d)//2>=v+w-t]
            assert overlaps
            overlap=max(overlaps) if rho1>=rho0 else min(overlaps)
            total+=remaining*rho0**(w-overlap)*rho1**overlap
        caps.append(total/(m+1))
    moments=np.array([(t-w)*math.log(g0)+w*math.log(g1) for w in levels])
    n=len(levels)+2;matrix=np.full((n,n),-np.inf)
    zero_terms=[];active_terms=[]
    for j,count in enumerate(kernel):
        weight=j*(math.log(theta)-lam)+(t-j)*math.log1p(-theta)
        if count:zero_terms.append(math.log(count)+weight)
        active=math.comb(t,j)-count
        if active:active_terms.append(math.log(active)+weight)
    matrix[0,0]=float(np.logaddexp.reduce(zero_terms))
    matrix[0,1]=float(np.logaddexp.reduce(active_terms))
    entries=moments+np.log(eps*np.array(caps)+(1-eps)/m)
    for i,entry in [(1,float(max(entries)))]+[(i+2,float(v)) for i,v in enumerate(entries)]:
        matrix[i,0]=entry
        for k,w in enumerate(levels):matrix[i,k+2]=entry+math.log(spectrum[w])
    return matrix


def exceptional_epochs(spectrum,kernel,coset,caps,lam,rounds=1):
    """Alternative Z/G/H envelope: H is the all-one image, G excludes H.

    Keeps the rare exceptional state explicit instead of allowing every lazy
    transition to adopt it for free. No uniformity inside G is assumed.
    """
    t=len(kernel)-1;m=sum(spectrum.values());levels=[w for w in spectrum if w!=t]
    eps=2.**-rounds;eta=1-eps
    out=np.full((t+1,3,3),-np.inf)
    for j in range(t+1):
        total=math.comb(t,j);beta=kernel[j]/total;gamma=coset[j]/total
        assert kernel[j]+coset[j]<=total
        moments=[]
        for w in levels:
            moments.append(float(np.logaddexp.reduce([
                math.log(math.comb(w,v)*math.comb(t-w,j-v))-math.log(total)-lam*(w+j-2*v)
                for v in range(max(0,j-t+w),min(w,j)+1)])))
        d=max(moments);r=int(caps[j]['cap'])/total
        cancel=r*math.exp(-lam*min(abs(w-j) for w in levels))
        # A transition G->0 or G->H requires syndrome outside {0,H}.
        cancel=min(cancel,(total-kernel[j]-coset[j])/total,math.exp(d))
        matrix=np.zeros((3,3))
        other=(total-kernel[j]-coset[j])/total
        matrix[0]=[beta,other,gamma]
        matrix[2]=[eps*gamma+eta*(1-beta)/m,
                   eps*other+eta*(m-2+beta+gamma)/m,
                   eps*beta+eta*(1-gamma)/m]
        with np.errstate(divide='ignore'):
            out[j,0]=np.log(matrix[0])-lam*j
            out[j,2]=np.log(matrix[2])-lam*(t-j)
        lc=math.log(cancel) if cancel else -math.inf
        # Fresh zero/H mass can be intersected with their necessary syndromes.
        out[j,1,0]=np.logaddexp(math.log(eps)+lc,math.log(eta/m)+min(d,math.log1p(-beta) if beta<1 else -math.inf))
        out[j,1,2]=np.logaddexp(math.log(eps)+lc,math.log(eta/m)+min(d,math.log1p(-gamma) if gamma<1 else -math.inf))
        out[j,1,1]=d
    return out


def mixture(rows,active,inactive):
    probabilities=np.array([0.])
    for p,q in zip(active,inactive):
        out=np.full(len(probabilities)+1,-np.inf)
        out[:-1]=probabilities+q
        out[1:]=np.logaddexp(out[1:],probabilities+p)
        probabilities=out
    return np.logaddexp.reduce(rows[:len(probabilities)]+probabilities[:,None,None],axis=0)


def adaptive(rows,roots,active,inactive,q):
    current=rows[:q+1]
    for _ in range(q):
        updated=np.full_like(current[:-1],-np.inf)
        for rho,p,n in zip(roots,active,inactive):
            np.maximum(updated,rho+np.logaddexp(n+current[:-1],p+current[1:]),out=updated)
        current=updated
    return terminal(current[0],128)+math.log(math.comb(32768,q))+q*math.log(len(roots))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--maximum',type=int,default=16);parser.add_argument('--rounds',type=int,default=1)
    parser.add_argument('--dense',action='store_true')
    parser.add_argument('--map',type=Path,help='exact separately audited RM2 map, preserving t=128,s=19')
    args=parser.parse_args();assert 3<=args.maximum<=128
    spectrum,kernel,sources=tv.fixed.load_inner()
    map_path=tv.ROOT/'workstreams/bare_bch_rm2sub/generated/manifest.json'
    columns=json.loads(map_path.read_text())['t128_s19']['columns']
    if args.map:
        map_path=args.map.resolve();record=json.loads(map_path.read_text())
        assert (record['t'],record['s'])==(128,19)
        columns=record['columns'];rows=tv.fixed.maps.generators(columns,19)
        exact=tv.fixed.maps.spectrum(rows)
        assert exact=={int(w):int(n) for w,n in record['spectrum'].items()}
        dual=tv.fixed.maps.dual_spectrum(exact,128,19)
        assert dual=={int(w):int(n) for w,n in record['kernel'].items()}
        assert len(set(columns))==128 and all(columns) and not any((a&b).bit_count()&1 for a in rows for b in rows)
        spectrum={w:n for w,n in exact.items() if w};kernel=[dual.get(j,0) for j in range(129)]
    caps=fiber_caps(128,19,spectrum,kernel);low=low_cancellation(columns,19)
    (HERE/'SYNDROME_FIBER_CAPS.json').write_text(json.dumps(caps,indent=2)+'\n',encoding='utf-8')
    witness_path=tv.fixed.HERE/'SMALLER_OUTWARD_WITNESSES.json'
    witnesses=json.loads(witness_path.read_text())
    counts={w:n for w,n in enumerate(tv.smaller_outer.spectrum()) if w and n}
    @lru_cache(maxsize=128)
    def at(z,q):return regions(epochs(spectrum,kernel,columns,caps,low,math.exp(z),args.rounds,min(q,128)),128,32768,q)
    @lru_cache(maxsize=64)
    def epoch_at(z):return epochs(spectrum,kernel,columns,caps,low,math.exp(z),args.rounds,128)
    results=[]
    # Q1 screen uses the same pointwise shells and exact singleton cancellation.
    weights=wm.all_weights(columns,19);cw=tv.cancellation_weights(columns)
    q1_record=dict(s=19,t=128,spectrum=spectrum,levels=sorted(spectrum),
                   cancellation=[[int(weights[col]),y] for col,y in zip(columns,cw)])
    tilts=np.arange(-1050,-749,dtype=float)/100
    rz,ra=wm.regions(*wm.transfers(q1_record,np.exp(tilts),args.rounds),256)
    q1_moments=wm.coefficients(rz,ra-math.log(256),128)
    for row in witnesses['results']:
        pairs=[];others=[]
        for witness in row['q2']:
            z=float(F(witness['log_tilt']));ids=witness['band_indices']
            roots,active,inactive=tv.fixed.general.density_roots(counts,128,witnesses['q2_bands'],float(F(witness['shift'])))
            matrix=mixture(at(z,2),active[ids],inactive[ids])
            value=terminal(matrix,128)+128*sum(roots[ids])+math.log(math.comb(32768,2))
            if ids[0]!=ids[1]:value+=math.log(2)
            value+=row['bad_weight']*math.exp(z);pairs.append(value)
        result=dict(distance_target=row['distance_target'],rounds=args.rounds,
                    q2_margin_bits=-float(np.logaddexp.reduce(pairs))/math.log(2),higher=[])
        q1_values=np.minimum(0.,q1_moments+row['bad_weight']*np.exp(tilts)[:,None])
        q1_best=np.min(q1_values,axis=0)
        q1_log=float(np.logaddexp.reduce([math.log(32768*n)+q1_best[w] for w,n in counts.items()]))
        result['q1_margin_bits']=-q1_log/math.log(2)
        print(json.dumps(result),flush=True)
        for witness in row['adaptive']:
            q=witness['occupation']
            if q>args.maximum:break
            z=float(F(witness['log_tilt']))
            roots,active,inactive=tv.fixed.general.density_roots(counts,128,witnesses['adaptive_bands'],float(F(witness['shift'])))
            value=adaptive(at(z,q),roots,active,inactive,q)+row['bad_weight']*math.exp(z)
            others.append(value);result['higher'].append(dict(occupation=q,margin_bits=-value/math.log(2)))
            print(row['distance_target'],q,-value/math.log(2),flush=True)
        result['q2_through_maximum_margin_bits']=-float(np.logaddexp.reduce(pairs+others))/math.log(2)
        if args.dense:
            dense=row['dense'];check_coverage(dense['selected_boxes'],32768,row['maximum_sparse']+1)
            logs=[]
            for i,box in enumerate(dense['selected_boxes']):
                witness=box['witness'];z=witness['log_surprisal']
                bank=dense['probability_banks'][witness['probability_bank']]
                ps=np.array(bank['probabilities']);proposal=np.array(witness['proposal']);proposal/=proposal.sum()
                costs=[0.]
                for band,p in zip(dense['bands'],ps[1:]):
                    costs.append(math.log(counts[128]) if band==[128] else max(
                        math.log(counts[w])-math.log(math.comb(128,w))-w*math.log(p)-(128-w)*math.log1p(-p) for w in band))
                theta=float(ps@proposal)
                probabilities=np.array([math.log(math.comb(128,j))+j*math.log(theta)+(128-j)*math.log1p(-theta) for j in range(129)])
                matrix=np.logaddexp.reduce(epoch_at(z)+probabilities[:,None,None],axis=0)
                moment=min(terminal(matrix,32768),terminal(bernoulli_epoch(spectrum,kernel,theta,math.exp(z),args.rounds),32768))
                corners=typed.vertices(box['lower'],box['upper'],32768)
                value=max(typed.point_logs(corners,32768,128,costs,proposal,moment,row['bad_weight'],math.exp(z)))
                value+=typed.lattice_log_count(box['lower'],box['upper']);logs.append(float(value))
                if i%64==0:print(row['distance_target'],'dense',i+1,'margin',-value/math.log(2),flush=True)
            result['dense_margin_bits']=-float(np.logaddexp.reduce(logs))/math.log(2)
            result['dense_coverage']=[row['maximum_sparse']+1,32768]
            result['dense_box_log_bounds']=logs
            print(row['distance_target'],'DENSE UNION',result['dense_margin_bits'],flush=True)
            if args.maximum>=row['maximum_sparse']:
                result['full_union_margin_bits']=-float(np.logaddexp.reduce([q1_log]+pairs+others+logs))/math.log(2)
                print(row['distance_target'],'FULL DIAGNOSTIC',result['full_union_margin_bits'],flush=True)
        result['covered_sparse_occupations']=[1,min(args.maximum,row['maximum_sparse'])]
        results.append(result)
    sources += [Path(__file__),Path(wm.__file__),Path(tv.__file__),map_path,witness_path,Path(tv.fixed.general.__file__)]
    sources += [Path(typed.__file__),tv.fixed.HERE/'verify_fixed_inner_results.py']
    payload=dict(status='BINARY64_OCCUPATION_DIAGNOSTIC_NOT_CERTIFICATE',
                 full_occupation_coverage=bool(args.dense and args.maximum>=max(r['maximum_sparse'] for r in witnesses['results'])),results=results,
                 source_sha256={p.relative_to(tv.ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    output=HERE/f'{"NO_CONSTANT_" if args.map else ""}SPARSE_MIXER_R{args.rounds}_Q{args.maximum}{"_DENSE" if args.dense else ""}.json'
    output.write_text(json.dumps(payload,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
