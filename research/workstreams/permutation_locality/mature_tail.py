"""Retain mature mass below expansion weights 48 and 56.

The two added coordinates bound subsets of the existing mature component;
they are not extra probability mass and have terminal coefficient zero.
Selected-support diagnostic only until a complete outward cover is replayed.
"""
import argparse
from itertools import combinations
from math import comb,log,prod
from fractions import Fraction
import numpy as np
from flint import arb,arb_mat,ctx
from occupancy_memory import Z,F,M,C,U,TERMINAL,prepare,rounded
from occupancy_window_average import prepare_inputs,averages,refine as window_refine
from occupancy_fresh_moment import fresh_census,refine as fresh_refine
from occupancy_multi_average import refine as multi_refine,multiplicities
from fresh_collision import refine as collision_refine
from zero_moment import census as zero_census,refine as zero_refine
from occupancy_model import local_data,placement
from occupancy_allones import weighted_cdf
from occupancy_screen import optimize
from bch_joint_support import authenticated_caps,support_caps
from shortened_bound import dimension_caps
from basis_lattice import improve_caps
from group_rank_one_verify import up

L48,L56=9,10
TAIL_TERMINAL=np.append(TERMINAL,[0.,0.])


def census(inputs,prepared):
    low,high,words,spectrum=inputs
    weight=np.bitwise_count(low)+np.bitwise_count(high)
    states=np.arange(len(low),dtype=np.uint64)
    atoms=prepared[0][0][3]
    result={}
    for b,entries in words.items():
        hits=np.zeros((2,len(low)),dtype=np.uint16)
        for syndrome,_,_ in entries:
            new=weight[states^np.uint64(syndrome)]
            for i,cut in enumerate((48,56)):
                hits[i]+=(new>0)&(new<=cut)
        records={}
        for i,cut in enumerate((48,56)):
            arbitrary=Fraction(int(hits[i,1:].max()),len(entries))
            fresh=max(Fraction(sum(int(hits[i,s])*c for s,c in dist.items()),len(entries)*sum(dist.values()))
                      for dist in atoms.values())
            uniform={v:Fraction(int(hits[i,weight==v].sum()),len(entries)*n) for v,n in spectrum.items()}
            by_weight={v:Fraction(int(hits[i,weight==v].max()),len(entries)) for v in spectrum}
            records[cut]=arbitrary,fresh,uniform,by_weight
            winner=int(hits[i,1:].argmax())+1
            for state in (0,1,17,winner):
                actual=sum(0<int(weight[state^syndrome])<=cut for syndrome,_,_ in entries)
                assert actual==int(hits[i,state])
        result[b]=records
        print('Mature-tail census weight',b,'max probabilities',
              [(cut,str(row[0])) for cut,row in records.items()],flush=True)
    return result


def lift(base,spectrum,probabilities,tilt,penalty,cut=64,pairs=None,class_tail=False,input_penalty=1,odd_penalty=1):
    """Keep the old density bounds; refine mass with two overlapping tail caps."""
    assert cut in (56,64)
    levels=sorted(spectrum);m=(1<<19)-1
    result=[]
    for j,old in enumerate(base):
        t=arb_mat([[old[a,b] if a<9 and b<9 else arb(0) for b in range(11)] for a in range(11)])
        # Incoming fresh/zero/uniform mass sent into M also bounds either tail.
        for source in [Z,F]+list(range(U,U+5)):
            t[source,L48]=old[source,M];t[source,L56]=old[source,M]
        candidates=[]
        for shape in multiplicities(j):
            total=sum(a*n for a,n in enumerate(shape,1))
            scale=arb(penalty)**shape[3]*arb(input_penalty)**total*arb(odd_penalty)**(shape[0]+shape[2])
            moment={v:min(arb(1),up((-arb(tilt)*v).exp()*prod(
                ((arb(128-v)*(-arb(tilt)*a).exp()+arb(v)*(arb(tilt)*a).exp())/128)**n
                for a,n in enumerate(shape,1)))) for v in levels}
            a=max(moment[v] for v in levels if v>=cut)
            b=max(arb(0),up(moment[56]-a)) if cut==64 else arb(0)
            c=max(arb(0),up(moment[48]-a-b))
            candidates.append((up(a*scale),up(b*scale),up(c*scale)))
        a,b,c=(max(row[i] for row in candidates) for i in range(3))
        for source,factor in ((M,a),(L56,b),(L48,c)):
            t[source,M]=up(factor/2)
            if j:t[source,Z]=up(factor/(2*m))
            for k,v in enumerate(levels):t[source,U+k]=up(factor*spectrum[v]/(2*m))
        if not j:
            # Empty lazy steps preserve the state and its expansion weight.
            f48=(-arb(tilt)*48).exp();f56=(-arb(tilt)*56).exp()
            t[L48,L48]=up(f48/2)
            t[L56,L56]=up(f56/2)
            t[L48,L56]=up((f48-f56)/2)
        elif j==1:
            for target,cutoff in ((L48,48),(L56,56)):
                for source in [F,M]+list(range(U,U+5)):
                    terms=[]
                    for column_weight,record in probabilities.items():
                        arbitrary,fresh,uniform=record[cutoff][:3]
                        p=arbitrary if source==M else fresh if source==F else uniform[levels[source-U]]
                        v=48 if source in (F,M) else levels[source-U]
                        f=(-arb(tilt)*max(0,v-column_weight)).exp()*arb(penalty)**int(column_weight==4)*arb(input_penalty)**column_weight*arb(odd_penalty)**(column_weight%2)
                        terms.append(up(f*arb(p.numerator)/p.denominator/2))
                    t[source,target]=min(old[source,M],max(terms))
        else:
            # Higher-occupancy steps may put all outgoing mature mass in a tail.
            for source in (M,L48,L56):
                t[source,L48]=t[source,M];t[source,L56]=t[source,M]
            if pairs is not None:
                for target,cutoff in ((L48,48),(L56,56)):
                    terms={source:[] for source in [Z,F,M]+list(range(U,U+5))}
                    for shape in multiplicities(j):
                        expanded=[a for a,n in enumerate(shape,1) for _ in range(n)]
                        choices=set(combinations(expanded,2))
                        probability=min(Fraction(1),min(pairs[p][cutoff][0] for p in choices)*Fraction(32*31,(34-j)*(33-j)))
                        zero_probability=pairs[tuple(expanded)][cutoff][1] if j==2 else probability
                        total=sum(expanded);scale=arb(penalty)**shape[3]*arb(input_penalty)**total*arb(odd_penalty)**(shape[0]+shape[2])
                        for source in terms:
                            if source==Z:
                                f=(-arb(tilt)*total).exp();p=zero_probability;half=1
                            else:
                                v=48 if source in (F,M) else levels[source-U]
                                f=(-arb(tilt)*max(0,v-total)).exp();p=probability;half=2
                            terms[source].append(up(f*scale*arb(p.numerator)/(half*p.denominator)))
                    for source,values in terms.items():
                        t[source,target]=min(old[source,M],max(values))
                    t[L48,target]=t[L56,target]=arb(0)
        if class_tail and j and (j==1 or pairs is not None):
            for target,cutoff in ((L48,48),(L56,56)):
                coefficients=[]
                for shape in multiplicities(j):
                    expanded=[w for w,n in enumerate(shape,1) for _ in range(n)]
                    total=sum(expanded);scale=arb(penalty)**shape[3]*arb(input_penalty)**total*arb(odd_penalty)**(shape[0]+shape[2])
                    if j==1:
                        ps=probabilities[expanded[0]][cutoff][3]
                    elif j==2:
                        ps=pairs[tuple(expanded)][cutoff][2]
                    else:
                        choices=set(combinations(expanded,2))
                        probability=min(Fraction(1),min(pairs[p][cutoff][0] for p in choices)*Fraction(32*31,(34-j)*(33-j)))
                        ps={v:probability for v in levels}
                    moment={v:up((-arb(tilt)*max(0,v-total)).exp()*arb(p.numerator)/p.denominator)
                            for v,p in ps.items()}
                    a=max(moment[v] for v in levels if v>=cut)
                    b=max(arb(0),up(moment[56]-a)) if cut==64 else arb(0)
                    c=max(arb(0),up(moment[48]-a-b))
                    coefficients.append(tuple(up(x*scale/2) for x in (a,b,c)))
                for index,source in enumerate((M,L56,L48)):
                    t[source,target]=max(row[index] for row in coefficients)
        result.append(rounded(t))
    return result


def self_test(inputs,operators,tilt,penalty,rounds=1,input_penalty=1,odd_penalty=1):
    low,high,words,spectrum=inputs
    weights=np.bitwise_count(low)+np.bitwise_count(high)
    states=sorted(set([1,17]+[int(np.flatnonzero(weights==v)[0]) for v in sorted(spectrum)]))
    old_precision=ctx.prec
    n=operators[0].nrows()
    ctx.prec=old_precision+128
    powers=[(-arb(tilt)*w).exp() for w in range(129)]
    lazy_probability=arb(2)**(-rounds)
    checks=0
    try:
        for j in (0,1):
            cases={0:[(0,0,0)]} if not j else words
            for b,entries in cases.items():
                rho=arb(penalty)**int(b==4)*arb(input_penalty)**b*arb(odd_penalty)**(b%2)
                for state in states:
                    actual=[arb(0)]*n;lazy={}
                    for syndrome,lo,hi in entries:
                        weight=(int(low[state])^lo).bit_count()+(int(high[state])^hi).bit_count()
                        output=powers[weight]*rho/len(entries)
                        f=output*lazy_probability
                        refresh=output*(1-lazy_probability)/((1<<19)-1)
                        target=state^syndrome
                        if target:
                            actual[M]+=f
                            lazy[target]=lazy.get(target,arb(0))+f
                            if weights[target]<=48:actual[L48]+=f
                            if weights[target]<=56:actual[L56]+=f
                        else:actual[Z]+=f
                        if syndrome:actual[Z]+=refresh
                        for k,v in enumerate(sorted(spectrum)):
                            actual[U+k]+=refresh*spectrum[v]
                    actual[C]=max(lazy.values(),default=arb(0))
                    if n==12:actual[11]=max((value for state,value in lazy.items() if weights[state]==48),default=arb(0))
                    incoming=[arb(0)]*n
                    incoming[M]=incoming[C]=arb(1)
                    incoming[L48]=arb(int(weights[state]<=48))
                    incoming[L56]=arb(int(weights[state]<=56))
                    if n==12:incoming[11]=arb(int(weights[state]==48))
                    for k,value in enumerate(actual):
                        bound=sum((incoming[i]*operators[j][i,k] for i in range(n)),arb(0))
                        assert value<=bound or (value-bound).contains(0),(j,b,state,k,value,bound)
                        checks+=1
        assert all(TAIL_TERMINAL[i]==0 for i in (C,L48,L56))
    finally:
        ctx.prec=old_precision
    print('Mature-tail direct state/empty/single-window component checks:',checks,'rounds',rounds,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=13)
    parser.add_argument('--support',type=int,default=152)
    parser.add_argument('--tilt',default='.008')
    parser.add_argument('--penalty',default='.6875')
    parser.add_argument('--pair-tail',action='store_true')
    parser.add_argument('--class-tail',action='store_true')
    parser.add_argument('--density48',action='store_true')
    args=parser.parse_args()
    spectrum=authenticated_caps();dimensions=dimension_caps()
    caps=improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum)
    counts=weighted_cdf(spectrum,caps,dimensions,args.penalty,prefix_flags=True)
    prepared=prepare(local_data(4));fresh=fresh_census(prepared)
    inputs=prepare_inputs();tails=census(inputs,prepared);zeros=zero_census();ctx.prec=192
    pairs=None
    if args.pair_tail:
        from pair_tail import census as pair_census
        pairs=pair_census(inputs,prepared)
    values=averages(inputs,args.tilt)
    base=fresh_refine(prepared,fresh,args.tilt,args.groups,args.penalty)
    base=window_refine(prepared,values,args.tilt,args.groups,args.penalty,base)
    base=multi_refine(base,fresh,args.tilt,args.penalty)
    base=collision_refine(base,prepared,args.tilt,args.penalty)
    base=zero_refine(base,zeros,args.tilt,args.penalty)
    density_values=None
    if args.density48:
        import density48
        density_values=density48.averages(inputs,density48.census(inputs,prepared),args.tilt)
    for cut in (None,56,64):
        ops=base if cut is None else lift(base,inputs[3],tails,args.tilt,args.penalty,cut,pairs,args.class_tail)
        if cut is not None and args.density48:ops=density48.lift(ops,density_values,args.tilt,args.penalty,pairs)
        if cut is not None:self_test(inputs,ops,args.tilt,args.penalty)
        terminal=TERMINAL if cut is None else TAIL_TERMINAL
        if cut is not None and args.density48:terminal=np.append(terminal,0.)
        regions=placement(ops,rounding=rounded);n=len(terminal)
        arrays=[np.array([[float(t[i,j]) for j in range(n)] for i in range(n)]) for t in regions]
        value,_=optimize(arrays,(args.support,)*args.groups,float(args.tilt),terminal)
        score=(value+args.groups*log(counts[args.support])+log(comb(2048,args.groups)))/log(2)
        print('MATURE-TAIL selected point:',cut,'log2 upper',score,flush=True)
    print('Selected-point binary64 screen, not a complete outward certificate.')


if __name__=='__main__':main()
