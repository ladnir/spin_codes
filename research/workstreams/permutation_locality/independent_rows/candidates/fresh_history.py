"""Remember a fresh packet's weight and optional expansion-weight class.

The added coordinates are masses with uniform-density envelopes on exact
subsets of single-packet feedbacks. Local calculations use Arb; the global
search is binary64 and is not a certificate.
"""
import argparse
from math import comb,log

import numpy as np
from flint import arb,arb_mat
import attack_cache
import shape_potential as potential
import mass_density_screen as screen
from mixing_attack import retarget_matrix
from occupancy_memory import Z,F,M,C,U
from group_rank_one_verify import up
from group_moment import maps


def classes(details,split=True):
    histograms,_,fresh,_=details['histograms']
    levels=[sum(r*n for r,n in enumerate(h)) for h in histograms]
    result=[]
    for a,distribution in enumerate(fresh,1):
        assert sum(distribution.values())==32*comb(4,a)
        for v in sorted({levels[h] for h,n in distribution.items() if n}) if split else [None]:
            members={h:n for h,n in distribution.items() if n and (v is None or levels[h]==v)}
            result.append(dict(a=a,v=v,hist=members,n=sum(members.values()),
                               lower=min(levels[h] for h in members)))
    return result


def return_moment(details,shape,record,tilt):
    """Exact moment of output times lazy return, without the lazy probability."""
    if shape not in details['fresh_classes']:return None
    rows=details['fresh_classes'][shape]
    chosen=[hist for (a,v),hist in rows.items() if a==record['a']
            and (record['v'] is None or record['v']==v)]
    denominator=details['cancellations'][0][shape][2]*record['n']
    return sum((count*(-arb(tilt)*w).exp() for hist in chosen
                for w,count in enumerate(hist) if count),arb(0))/denominator


def single_density(records,tilt,bits=40):
    """Integer-weighted exact single-packet pairs; round exponentials up."""
    images,columns,_=maps();mask=(1<<64)-1
    low=np.array([x&mask for x in images],dtype=np.uint64)
    high=np.array([x>>64 for x in images],dtype=np.uint64)
    weights=np.bitwise_count(low)+np.bitwise_count(high)
    packets={a:[] for a in range(1,5)};seen=set()
    for window in range(32):
        for lane in range(1,16):
            state=0
            for bit in range(4):
                if lane>>bit&1:state^=columns[4*window+bit]
            if not state or state in seen:raise ValueError('feedback injectivity required')
            seen.add(state);word=lane<<(4*window)
            packets[lane.bit_count()].append((state,word&mask,word>>64))
    scale=1<<bits
    powers=np.array([int(((-arb(tilt)*w).exp()*scale).upper().ceil().unique_fmpz())
                     for w in range(129)],dtype=np.int64)
    result={}
    for r in records:
        sources=np.array([s for s,_,_ in packets[r['a']] if r['v'] is None or weights[s]==r['v']],dtype=np.int64)
        if len(sources)!=r['n']:raise ValueError('fresh class size mismatch')
        for b,choices in packets.items():
            sy=np.array([s for s,_,_ in choices],dtype=np.int64)
            lo=np.array([x for _,x,_ in choices],dtype=np.uint64)
            hi=np.array([x for _,_,x in choices],dtype=np.uint64)
            target=sources[:,None]^sy
            emitted=np.bitwise_count(low[sources,None]^lo)+np.bitwise_count(high[sources,None]^hi)
            den=len(sources)*len(choices)
            if den*int(powers.max())>=1<<63:raise ValueError('integer pair-sum overflow')
            values=np.zeros(len(images),dtype=np.int64)
            np.add.at(values,target.ravel(),powers[emitted].ravel())
            assert int(values.sum())==int(powers[emitted].sum())
            result[r['a'],r['v'],b]={'density':up(arb(int(values[1:].max()))/(den*scale)),
                **{cut:up(arb(int(values[(weights>0)&(weights<=cut)].sum()))/(den*scale)) for cut in (48,56)}}
    return result


def family_matrix(old,occupancy,shape,details,feedback,records,tilt,penalty,rounds,single=None):
    size=11+len(records);matrix=arb_mat(size,size)
    for i in range(11):
        for j in range(11):matrix[i,j]=old[i,j]
    spectrum=details['spectrum'];states=sum(spectrum.values());levels=sorted(spectrum)
    alpha=arb(0) if rounds is None else arb(2)**(-rounds);beta=1-alpha
    if occupancy==1:
        if shape is None:raise ValueError('single-packet shape must be retained')
        a=shape[0]
        for k,record in enumerate(records,11):
            if record['a']==a:matrix[Z,k]=up(old[Z,F]*record['n']/(32*comb(4,a)))
        matrix[Z,F]=0
    for k,record in enumerate(records,11):
        n=record['n'];v=record['lower'];a=record['a']
        if occupancy==0:
            moment=(-arb(tilt)*v).exp()
            matrix[k,k]=up(alpha*moment)
            for ell,w in enumerate(levels):matrix[k,U+ell]=up(beta*moment*spectrum[w]/states)
            continue
        # Two representations of the same incoming measure. Both bound the
        # same lazy-mature and refreshed-uniform output components.
        for target in range(11):
            via_fresh=up(old[F,target]*(32*comb(4,a))/n)
            via_mature=up(old[M,target]+old[C,target]/n
                          +int(v<=48)*old[9,target]+int(v<=56)*old[10,target])
            matrix[k,target]=min(via_fresh,via_mature)
        if shape is None:continue
        counts=tuple(shape.count(b) for b in range(1,5))
        moment=up(sum((number*details['moments'][counts][h]
                       for h,number in record['hist'].items()),arb(0))/n)
        scale=arb(penalty)**shape.count(4)
        closed=return_moment(details,shape,record,tilt)
        coarse=(-arb(tilt)*max(0,v-sum(shape))).exp()
        if shape in feedback:
            zero,peak,den,_=feedback[shape]
            atom=arb(max(zero,peak))/den;nonzero=arb(peak)/den
        else:
            atom=nonzero=arb(1)/((33-len(shape))*max(comb(4,b) for b in shape))
        density=up(alpha*scale*coarse*min(arb(1)/n,atom))
        lazy_zero=up(coarse*min(arb(1)/n,nonzero))
        if closed is not None:lazy_zero=min(lazy_zero,up(closed))
        matrix[k,Z]=min(matrix[k,Z],up(scale*(alpha*lazy_zero+beta*moment/states)))
        # Exact closed mass can be subtracted using its LOWER endpoint.
        remaining=moment if closed is None else up(moment-closed.lower())
        if remaining<0:raise ValueError('return mass exceeds output moment')
        matrix[k,M]=min(matrix[k,M],up(alpha*scale*remaining))
        matrix[k,C]=min(matrix[k,C],density)
        for ell,w in enumerate(levels):
            matrix[k,U+ell]=min(matrix[k,U+ell],up(beta*scale*moment*spectrum[w]/states))
        for target in (9,10):matrix[k,target]=min(matrix[k,target],up(alpha*scale*moment))
        if occupancy==1 and single is not None:
            exact=single[a,record['v'],shape[0]]
            matrix[k,C]=min(matrix[k,C],up(alpha*scale*exact['density']))
            for target,cut in ((9,48),(10,56)):
                matrix[k,target]=min(matrix[k,target],up(alpha*scale*exact[cut]))
    if any(matrix[i,j]<0 for i in range(size) for j in range(size)):
        raise ValueError('negative history coefficient')
    return matrix


def build(base,feedback,details,tilt,penalty='.9',rounds=2,split=True):
    old_tilt,old_penalty,old_rounds=details['parameters']
    if screen.Q(tilt)!=old_tilt:raise ValueError('mismatched local tilt')
    records=classes(details,split)
    single=single_density(records,tilt)
    shaped=potential.shape_matrices(base,feedback,details,tilt,str(old_penalty),old_rounds)
    result=[]
    for j in range(33):
        family=[]
        for shape,old in shaped.get(j,{None:base[j]}).items():
            old=retarget_matrix(old,j,details['spectrum'],tilt,old_rounds,rounds)
            ratio=arb(penalty)/arb(str(old_penalty))
            scale=ratio**shape.count(4) if shape is not None else max(arb(1),up(ratio**j))
            old=old*scale
            matrix=family_matrix(old,j,shape,details,feedback,records,tilt,penalty,rounds,single)
            family.append(np.array([[float(matrix[i,k]) for k in range(matrix.ncols())]
                                    for i in range(matrix.nrows())]))
        result.append(np.array(family))
    terminal=np.concatenate((screen.baseline.TAIL_TERMINAL,np.ones(len(records))))
    return result,terminal,records


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilt',default='.056')
    parser.add_argument('--penalty',default='.9')
    parser.add_argument('--rounds',type=int,default=2)
    parser.add_argument('--groups',type=int,nargs='+',default=[96,128])
    parser.add_argument('--iterations',type=int,default=10)
    parser.add_argument('--weight-only',action='store_true')
    parser.add_argument('--directory',default='tmp/four-bit-attack-cache')
    args=parser.parse_args()
    data=attack_cache.build([args.tilt],directory=args.directory)[args.tilt]
    base,feedback,details=data
    families,terminal,records=build(*data,args.tilt,args.penalty,args.rounds,not args.weight_only)
    print('FRESH HISTORY classes',[(r['a'],r['v'],r['n']) for r in records],flush=True)
    caps=screen.baseline.authenticated_caps()
    full=1/screen.Q(args.penalty)
    count=min(screen.Q(screen.baseline.integer_cdf(screen.baseline.weighted_cdf_upper(caps,1<<128,full_weight=full))[200]),
              screen.baseline.weighted_union_shells(caps,full_weight=full)[200])
    baseline_region=screen.float_placement(screen.as_array(base),max(args.groups))
    for q in args.groups:
        original,p=screen.score(baseline_region,q,200,count,args.tilt,cutoff=193986)
        masses=np.array([1.])
        for _ in range(q):masses=np.convolve(masses,[1-p,p])
        matrix=np.einsum('r,rij->ij',masses,baseline_region[:q+1])
        old=potential.matrix_potential(matrix)
        initial=np.concatenate((old,[min(old[F]*32*comb(4,r['a'])/r['n'],
            old[M]+old[C]/r['n']+int(r['lower']<=48)*old[9]+int(r['lower']<=56)*old[10]) for r in records]))
        witness=potential.common_potential(families,potential.placement_weights(q),masses,initial,
                                          iterations=args.iterations,terminal=terminal)
        outer=(float(args.tilt)*193986-q*screen.baseline.log_binomial_mass(256,200,p)
               +q*(log(count.numerator)-log(count.denominator))+log(comb(2048,q)))
        print('FRESH HISTORY q/tilt/rho/R',q,args.tilt,args.penalty,args.rounds,
              'log2',(witness['log_moment']+outer)/log(2),'control',original,'p',p,
              'witness',witness,flush=True)
    print('Binary64 selected-support screen only; no new full-code certificate.')


if __name__=='__main__':main()
