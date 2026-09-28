"""Independent-row SPIN with two-bit packets: isolated proof model.

Same BCH and IMT maps as the four-bit experiment. Only packet geometry
changes: 4096 row pairs, 64 slots/epoch, 64 epochs/region, 256 regions.
No production defaults, four-bit verifier, or old certificate are changed.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import comb
from pathlib import Path
import sys

import numpy as np
from flint import arb,arb_mat,ctx

ROOT=Path(__file__).resolve().parent.parent
for directory in (ROOT,ROOT/'independent_rows',ROOT/'independent_rows/candidates',
                  ROOT/'independent_rows/candidates/overlap_second'):
    sys.path.insert(0,str(directory))
from group_moment import maps
from group_rank_one_verify import up
from occupancy_model import placement
from occupancy_memory import rounded
from exact_feedback import inverse_counts
from second import moments as second_moments
from quadratic import witness
from density_second import shape_bounds as density_bounds
from feedback_density import class_transforms

Z,F,M,C,U=0,1,2,3,4
GROUPS=4096
REGIONS=256
WINDOWS=64
EPOCHS=64
THRESHOLD=209715


def checked_updates(value):
    if type(value) is not int or not 1<=value<=32:
        raise ValueError('integer update count in 1..32 required')
    return value


def resolve_updates(requested,record=None,*,retarget=False):
    """Legacy witnesses use two updates; changing the code requires retargeting."""
    saved=checked_updates(2 if record is None else record.get('updates',2))
    result=saved if requested is None else checked_updates(requested)
    if retarget and record is None:raise ValueError('retarget requires a saved partition')
    if record is not None and result!=saved and not retarget:
        raise ValueError('replay and resume must preserve the update count; use retarget')
    return result


def shapes(maximum):
    return sorted(((a,b) for a in range(maximum+1) for b in range(maximum+1-a)),
                  key=lambda s:(sum(s),s))


def characters(columns,state_bits,maximum):
    if (len(columns)%2 or not 1<=maximum<=len(columns)//2
            or not 1<=state_bits<=20 or any(not 0<=c<1<<state_bits for c in columns)):
        raise ValueError('complete two-bit geometry and bounded state dimension required')
    windows=len(columns)//2;states=np.arange(1<<state_bits,dtype=np.uint32)
    hist=np.zeros((len(states),3),dtype=np.uint8)
    for w in range(windows):
        h=((np.bitwise_count(states&columns[2*w])&1)
           +(np.bitwise_count(states&columns[2*w+1])&1))
        for r in range(3):hist[:,r]+=h==r
    keys=sum(hist[:,r].astype(np.uint64)*(windows+1)**r for r in range(3))
    unique,inverse,multiplicities=np.unique(keys,return_inverse=True,return_counts=True)
    records=np.array([[int(k//(windows+1)**r%(windows+1)) for r in range(3)] for k in unique])
    assert np.all(records.sum(axis=1)==windows)
    for a in sorted({0,1,min(17,len(states)-1),len(states)-1}):
        direct=Counter(sum((a&c).bit_count()%2 for c in columns[2*w:2*w+2]) for w in range(windows))
        assert tuple(records[inverse[a]])==tuple(direct[r] for r in range(3))
    selected=shapes(maximum);index={s:i for i,s in enumerate(selected)}
    denominator={s:comb(windows,sum(s))*comb(sum(s),s[0])*2**s[0] for s in selected}
    if max(denominator.values())>=1<<63:raise ValueError('shape counts exceed signed-int64 range')
    coefficients=np.zeros((len(selected),len(records)),dtype=np.int64);coefficients[0]=1
    cumulative=np.cumsum(records,axis=1)
    for w in range(windows):
        h=np.sum(cumulative<=w,axis=1)
        for a,b in reversed(selected[1:]):
            if a+b>w+1:continue
            if a:coefficients[index[a,b]]+=coefficients[index[a-1,b]]*(2-2*h)
            if b:coefficients[index[a,b]]+=coefficients[index[a,b-1]]*(1-2*(h%2))
    for s,row in zip(selected,coefficients):
        assert int(row[inverse[0]])==denominator[s]
        assert int(row.min())>=-denominator[s] and int(row.max())<=denominator[s]
    print('TWO-BIT character histograms',len(records),'nonempty shapes',len(selected)-1,flush=True)
    return selected,coefficients,denominator,multiplicities,inverse


def pattern_histogram(word,windows=64):
    weights=Counter((word>>(2*w)&3).bit_count() for w in range(windows))
    return tuple(weights[h] for h in range(3))


def output_coefficients(histogram,maximum,tilt):
    """Exact-law positive polynomial evaluated outward, without replacement."""
    if len(histogram)!=3 or min(histogram)<0 or Q(tilt)<0:
        raise ValueError('nonnegative paired-window histogram and tilt required')
    windows=sum(histogram)
    if not 0<=maximum<=windows:raise ValueError('invalid local degree')
    selected=shapes(maximum);coefficients={s:arb(int(s==(0,0))) for s in selected}
    lam=arb(Q(tilt).numerator)/Q(tilt).denominator;done=0
    for h,copies in enumerate(histogram):
        one=(2-h)*(-lam).exp()+h*lam.exp()
        two=(-lam*(2-2*h)).exp()
        for _ in range(copies):
            done+=1
            for a,b in reversed(selected[1:]):
                if a+b>done:continue
                if a:coefficients[a,b]+=coefficients[a-1,b]*one
                if b:coefficients[a,b]+=coefficients[a,b-1]*two
    v=histogram[1]+2*histogram[2]
    return {s:up(coefficients[s]*(-lam*v).exp()/
                 (comb(windows,sum(s))*comb(sum(s),s[0])*2**s[0])) for s in selected}


def census(maximum=8):
    if not 1<=maximum<=10:raise ValueError('complete census cutoff must be in [1,10]')
    images,columns,spectrum=maps();levels=sorted(spectrum)
    data=characters(columns,19,maximum)
    selected,table,denominators,_,inverse=data
    expansion=np.array([x.bit_count() for x in images],dtype=np.int64)
    prepared=class_transforms(expansion)
    feedback={}
    for shape,row in zip(selected[1:],table[1:]):
        D=denominators[shape];counts=inverse_counts(row[inverse],D)
        classes={v:int(counts[expansion==v].sum()) for v in levels}
        feedback[shape]=(counts,int(counts[0]),int(counts[1:].max()),D,classes)
        assert int(counts[0])+sum(classes.values())==D
    single={};inputs={};cancellation={};fresh_hist={}
    for b in (1,2):
        words=[mask<<(2*w) for w in range(64) for mask in (1,2,3) if mask.bit_count()==b]
        atoms={};cancel=Counter()
        for word in words:
            state=0
            for p in range(128):
                if word>>p&1:state^=columns[p]
            assert state and state not in atoms
            atoms[state]=word;cancel[(images[state].bit_count(),(images[state]^word).bit_count())]+=1
        single[b]=atoms;inputs[b]=words;cancellation[b]=cancel
        fresh_hist[b]=Counter(pattern_histogram(images[s]) for s in atoms)
        record=feedback[(int(b==1),int(b==2))]
        assert record[1:4]==(0,1,len(words))
    assert set(single[1]).isdisjoint(single[2])
    convolutions={(a,b):Counter(s^t for s in single[a] for t in single[b]) for a,b in product((1,2),repeat=2)}
    peaks={b:max(Q(max(n for s,n in convolutions[a,b].items() if s),len(single[a])*len(single[b]))
                 for a in (1,2)) for b in (1,2)}
    lo=np.array([x&((1<<64)-1) for x in images],dtype=np.uint64)
    hi=np.array([x>>64 for x in images],dtype=np.uint64)
    paired=(np.bitwise_count(lo&(lo>>np.uint64(1))&np.uint64(0x5555555555555555))
            +np.bitwise_count(hi&(hi>>np.uint64(1))&np.uint64(0x5555555555555555)))
    histograms={v:Counter() for v in levels}
    keys=expansion*65+paired
    keys,counts=np.unique(keys[1:],return_counts=True)
    for key,count in zip(keys,counts):
        v,h2=divmod(int(key),65);h1=v-2*h2;h0=64-h1-h2
        assert min(h0,h1,h2)>=0
        histograms[v][h0,h1,h2]=int(count)
    assert all(sum(histograms[v].values())==spectrum[v] for v in levels)
    for s in (1,17,524287):assert pattern_histogram(images[s]) in histograms[images[s].bit_count()]
    expansion_basis=[images[1<<b] for b in range(19)]
    second=second_moments(expansion_basis,columns,maximum,width=2,character_data=(selected,table,inverse))
    for shape,(_,zero,peak,D,_) in feedback.items():
        if sum(shape)<=4:print('TWO-BIT feedback',shape,'zero',Q(zero,D),'peak',Q(peak,D),flush=True)
    print('TWO-BIT fresh nonzero convolution peaks',peaks,flush=True)
    print('TWO-BIT zero atom exceeds nonzero peak in',
          sum(zero>peak for _,zero,peak,_,_ in feedback.values()),'censused shapes',flush=True)
    return dict(maximum=maximum,spectrum=spectrum,levels=levels,feedback=feedback,second=second,
                single=single,cancellation=cancellation,fresh_hist=fresh_hist,peaks=peaks,
                histograms=histograms,expansion=expansion,prepared=prepared)


def prepare_density(data,tilts):
    """Share exact translated convolutions across several output tilts."""
    tilts=tuple(dict.fromkeys(str(t) for t in tilts));cache=data.setdefault('density_cache',{})
    if any(Q(t)<=0 for t in tilts):raise ValueError('positive density tilts required')
    for shape,(counts,_,_,D,_) in data['feedback'].items():
        weights=(1,)*shape[0]+(2,)*shape[1]
        if sum(weights)>min(data['levels']):continue
        missing=[t for t in tilts if (ctx.prec,t,shape) not in cache]
        if not missing:continue
        values=density_bounds(counts,D,weights,data['expansion'],data['prepared'],
                              data['second'][weights][2],missing)
        for t,record in values.items():cache[ctx.prec,t,shape]=record
    print('TWO-BIT shared density preparation:',len(data['feedback']),'shapes,',len(tilts),'tilts',flush=True)


def conditional_atom(shape,feedback,windows=64):
    """Bound any syndrome by leaving a censused subset of packets unexposed.

    Conditional on the exposed slots, the remaining d packets are a uniform
    input of their shape conditioned to avoid those slots. Dividing its
    unconditional atom cap by that avoidance probability is valid for every
    exposed feedback value, including shifts to feedback zero.
    """
    a,b=shape;j=a+b
    if min(a,b)<0 or not 1<=j<=windows:raise ValueError('valid nonempty packet shape required')
    best=Q(1)
    for (x,y),(_,zero,peak,D,_) in feedback.items():
        d=x+y
        if not 1<=d<=j or x>a or y>b:continue
        cap=Q(max(zero,peak),D)*Q(comb(windows,d),comb(windows-j+d,d))
        best=min(best,cap)
    return min(Q(1),best)


def epoch_operators(data,tilt,penalty='1',*,maximum=64,translated_density=True,output_degree=None,mass_columns=False,updates=2):
    if type(updates) is not int or not 1<=updates<=32:raise ValueError('integer update count in 1..32 required')
    if not 1<=maximum<=64 or not 0<Q(penalty)<=1 or Q(tilt)<0:
        raise ValueError('valid local degree, penalty and nonnegative tilt required')
    levels=data['levels'];spectrum=data['spectrum'];m=(1<<19)-1;n=4+len(levels)
    lam=arb(Q(tilt).numerator)/Q(tilt).denominator
    powers=[up((-lam*w).exp()) for w in range(129)]
    alpha=arb(2)**-updates;beta=1-alpha
    histograms=set(h for by_level in data['histograms'].values() for h in by_level)
    moment_degree=data['maximum'] if output_degree is None else output_degree
    if not data['maximum']<=moment_degree<=64:raise ValueError('output degree must cover the census and be <=64')
    values={h:output_coefficients(h,moment_degree,tilt) for h in histograms}
    result=[]
    def average(hist,shape):return up(sum((values[h][shape]*count for h,count in hist.items()),arb(0))/sum(hist.values()))
    empty=arb_mat(n,n);empty[Z,Z]=1
    for i in (F,M):
        empty[i,i]=alpha*powers[min(levels)]
        moment=(max(average(h,(0,0)) for h in data['fresh_hist'].values()) if i==F else powers[min(levels)])
        for k,v in enumerate(levels):empty[i,U+k]=beta*moment*spectrum[v]/m
    empty[C,C]=alpha*powers[min(levels)]
    for i,v in enumerate(levels):
        empty[U+i,U+i]=alpha*powers[v]
        for k,w in enumerate(levels):empty[U+i,U+k]+=beta*powers[v]*spectrum[w]/m
    result.append(rounded(empty))
    for j in range(1,maximum+1):
        candidates=[]
        for b in range(j+1):
            shape=(j-b,b);weights=(1,)*(j-b)+(2,)*b;W=j+b
            raw=powers[min(abs(v-W) for v in levels)]
            local=data['feedback'].get(shape)
            atom=conditional_atom(shape,data['feedback'])
            if local:
                counts,zero,peak,D,classes=local;p0=Q(zero,D);atom=Q(peak,D);live=1-p0
            else:
                # Unknown zero probability: never subtract its upper bound.
                p0=atom if W%2==0 else Q(0);live=Q(1);point=raw
                uniform={v:powers[abs(v-W)] for v in levels}
                fresh=raw
            if j<=moment_degree:
                point=max(row[shape] for row in values.values())
                uniform={v:average(data['histograms'][v],shape) for v in levels}
                fresh=max(average(h,shape) for h in data['fresh_hist'].values())
            def aq(value):return arb(value.numerator)/value.denominator
            matrix=arb_mat(n,n);matrix[Z,Z]=powers[W]*aq(p0)
            if j==1:matrix[Z,F]=powers[W]
            else:
                matrix[Z,M]=powers[W]*aq(live);matrix[Z,C]=powers[W]*aq(atom)
            lazy_zero=raw*aq(min(atom,Q(1,64)))
            lazy_density=raw*aq(min(max(atom,p0),Q(1,64)))
            if j==1:
                weight=W;hist=data['cancellation'][weight];choices=len(data['single'][weight])
                lazy_zero=sum((count*powers[w] for (_,w),count in hist.items()),arb(0))/(choices*choices)
                lazy_density=raw*aq(data['peaks'][weight])
            matrix[F,Z]=alpha*lazy_zero+beta*fresh/m
            matrix[F,M]=alpha*fresh;matrix[F,C]=alpha*lazy_density
            matrix[M,Z]=beta*point/m;matrix[M,M]=alpha*point
            matrix[C,Z]=alpha*raw*aq(live);matrix[C,C]=alpha*raw
            if local:
                matrix[C,Z]=alpha*sum((powers[abs(v-W)]*count for v,count in classes.items()),arb(0))/D
                if Q(tilt)>0 and W<=min(levels):
                    second=data['second'][weights]
                    assert second[0]==D
                    restricted=second[1]-Q(W*W*zero,D)
                    quadratic,_,_=witness({v:Q(count,D) for v,count in classes.items()},W,restricted,tilt)
                    matrix[C,Z]=min(matrix[C,Z],alpha*quadratic)
            for i,moment in ((F,fresh),(M,point)):
                for k,v in enumerate(levels):matrix[i,U+k]=beta*moment*spectrum[v]/m
            for i,v in enumerate(levels):
                peak=powers[abs(v-W)];moment=uniform[v]
                cancellation=peak*(aq(Q(classes[v],D)) if local else aq(live))/spectrum[v]
                matrix[U+i,Z]=alpha*cancellation+beta*moment/m
                matrix[U+i,M]=alpha*moment;matrix[U+i,C]=alpha*peak/spectrum[v]
                for k,w in enumerate(levels):matrix[U+i,U+k]=beta*moment*spectrum[w]/m
            if local and translated_density and not mass_columns and Q(tilt)>0 and W<=min(levels):
                key=ctx.prec,str(tilt),shape
                cache=data.setdefault('density_cache',{})
                if key not in cache:
                    cache[key]=density_bounds(counts,D,weights,data['expansion'],data['prepared'],
                                              data['second'][weights][2],[tilt])[str(tilt)]
                bounds=cache[key]
                # The cached translated bounds include the reference lazy
                # factor 1/4 (density_second.shape_bounds rounds=2).
                density_scale=arb(2)**(2-updates)
                matrix[C,C]=min(matrix[C,C],density_scale*bounds['density'])
                for i,v in enumerate(levels):matrix[U+i,C]=min(matrix[U+i,C],density_scale*bounds['uniform'][v])
            if mass_columns:
                # A lazy atom is bounded both by total tilted mass and by the
                # feedback atom cap times the pointwise output tilt. Replace
                # each entire coupled column; never take unrelated entrywise
                # minima of the mass and density alternatives.
                matrix[M,Z]+=alpha*min(point,raw*aq(atom));matrix[C,Z]=0
                matrix[M,C]=alpha*min(point,raw*aq(max(atom,p0)));matrix[C,C]=0
            candidates.append(rounded(matrix*aq(Q(penalty)**b)))
        result.append(arb_mat([[max(t[i,k] for t in candidates) for k in range(n)] for i in range(n)]))
    print('TWO-BIT local operators',tilt,penalty,'updates',updates,'complete occupancies',maximum,flush=True)
    return result


def terminal():return np.array([1.,1.,1.,0.,1.,1.,1.,1.,1.])


def log_power_moment(matrix,length=256,terminal=None):
    """Normalize before the first squaring in binary64 witness selection."""
    from math import log
    from two_group_screen import log_power_moment as shared
    scale=float(np.max(matrix))
    if not np.isfinite(scale) or scale<=0:
        return float('inf')  # Lost numerical mass is not a successful witness.
    value=shared(matrix/scale,length,terminal)+length*log(scale)
    return value if np.isfinite(value) else float('inf')


def region(operators,q):
    return placement(operators,epochs=EPOCHS,windows=WINDOWS,maximum_groups=q,rounding=rounded)


def float_region(operators,q):
    degree=min(64,len(operators)-1);n=operators[0].nrows()
    if degree<min(q,64):raise ValueError('missing local occupancies')
    matrices=np.array([[[float(t[i,j]) for j in range(n)] for i in range(n)] for t in operators])
    current=np.zeros((q+1,n,n));current[0]=np.eye(n)
    for e in range(64):
        following=np.zeros_like(current)
        previous_slots=64*e;total_slots=previous_slots+64
        for k in range(min(q,degree)+1):
            stop=min(q,previous_slots+k)
            # Hypergeometric weights in r, via a positive recurrence. This is
            # a binary64 proposal only; model.region remains the exact-law
            # outward replay. Avoid millions of giant integer binomials here.
            weights=np.empty(stop-k+1)
            weights[0]=float(Q(comb(64,k),comb(total_slots,k)))
            if stop>k:
                r=np.arange(k,stop,dtype=float)
                ratios=(previous_slots-r+k)/(r+1-k)*(r+1)/(total_slots-r)
                weights[1:]=weights[0]*np.cumprod(ratios)
            following[k:stop+1]+=weights[:,None,None]*(current[:stop-k+1]@matrices[k])
        current=following
    return current
