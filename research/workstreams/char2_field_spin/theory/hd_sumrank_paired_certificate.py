"""Directed paired-profile regional bounds for the feedback-refresh candidate.

This proves only the row classes explicitly returned. It is not a whole-code
certificate until a coverage driver verifies their union has no gaps.
"""

from fractions import Fraction as F
from functools import lru_cache
from math import comb

from flint import arb,arb_poly,ctx

from hd_sumrank_transfer import profiles,minimum_live_weight
from hd_randomized_outer_quiet_exact import _arb_fraction


def round_up(value,bits=12):
    """Exactly round a nonnegative rational upward to a dyadic mantissa."""
    value=F(value)
    if not value:
        return value
    if value<0:
        raise ValueError('Negative positive-operator entry')
    exponent=value.numerator.bit_length()-value.denominator.bit_length()
    quantum=F(2)**(exponent-bits)
    scaled=value/quantum
    return ((scaled.numerator+scaled.denominator-1)//scaled.denominator)*quantum


def grouped_packets(z,state=16,cutoff=20,bits=12,maximum=None,feedthrough='cauchy',width=64):
    if width not in (8,16,32,64,128,256) or (width!=64 and feedthrough not in ('systematic_mds','bucket4','bucket8')):
        raise ValueError('Unsupported packet width for this construction')
    if maximum is None:maximum=width
    if not 0<=maximum<=width:raise ValueError('Invalid packet degree cutoff')
    if feedthrough in ('common1','common2','common4'):
        if state!=16 or width!=64:raise ValueError('Common refresh requires state16, width64')
        from hd_common_refresh_packets import grouped_packets as common_groups
        return common_groups(z,bits,maximum,int(feedthrough[-1]))
    if feedthrough in ('block4-rank2-fb4','block4-rank2-fb8'):
        if state!=16 or width!=64:raise ValueError('Rank-two block profile requires state16, width64')
        from hd_rank_two_block_packets import grouped_packets as rank_two_groups
        return rank_two_groups(z,int(feedthrough[-1]),bits,maximum,cutoff)
    if feedthrough in ('bucket4','bucket8'):
        from hd_bucket_feedback import grouped_packets as bucket_groups
        return bucket_groups(z,state,width,int(feedthrough[6:]),bits,maximum,cutoff)
    if feedthrough=='systematic_mds':
        from hd_systematic_mds import matrix
        return tuple((j,comb(width,j),tuple(round_up(x,bits) for x in matrix(j,z,state,width)))
                     for j in range(maximum+1))
    q,d=65536,65536**state-1
    groups={}; tails={}
    for profile,j,multiplicity,k,t,v in profiles(F(z),state=state,maximum=maximum,feedthrough=feedthrough):
        rank=min(state,sum(min(4,a) for a in profile))
        gamma=min(k,F(z)**minimum_live_weight(profile,state)/q**rank,d*v)
        matrix=(k,t-k,v/d if j else F(0),(v if j else F(q-1,q)*v)+gamma/q)
        if j>cutoff:
            if j not in tails: tails[j]=list(matrix)
            else: tails[j]=[max(x,y) for x,y in zip(tails[j],matrix)]
        else:
            key=(j,tuple(round_up(x,bits) for x in matrix))
            groups[key]=groups.get(key,0)+multiplicity
    result=[(j,count,matrix) for (j,matrix),count in groups.items()]
    result.extend((j,comb(64,j),tuple(round_up(x,bits) for x in matrix)) for j,matrix in tails.items())
    for j in range(maximum+1):
        assert sum(n for h,n,m in result if h==j)==comb(64,j)
    return tuple(result)


def pair_coefficients(z,alpha,state=16,cutoff=20,bits=12,progress=None,maximum=None,feedthrough='cauchy',width=64,right_feedthrough=None):
    groups=grouped_packets(z,state,cutoff,bits,maximum,feedthrough,width)
    exponent=_arb_fraction(F(alpha))
    groups=tuple((j,count,tuple(_arb_fraction(x) for x in matrix)) for j,count,matrix in groups)
    right=groups if right_feedthrough is None else tuple((j,count,tuple(_arb_fraction(x) for x in matrix))
             for j,count,matrix in grouped_packets(z,state,cutoff,bits,maximum,right_feedthrough,width))
    coefficients=[[arb(0) for _ in range(2*width+1)] for _ in range(4)]
    for index,(j,n,a) in enumerate(groups):
        for k,m,b in right:
            for r in range(2):
                for s in range(2):
                    value=a[2*r]*b[s]+a[2*r+1]*b[2+s]
                    if not value.is_zero():
                        coefficients[2*r+s][j+k]+=n*m*value**exponent
        if progress is not None and (index+1)%100==0:
            progress(index+1,len(groups))
    return coefficients,len(groups)


def poly_product(a,b,degree):
    return tuple(((a[2*r]*b[s]+a[2*r+1]*b[2+s]).truncate(degree+1)) for r in range(2) for s in range(2))


def scalar_power(matrix,exponent):
    result=(arb(1),arb(0),arb(0),arb(1))
    while exponent:
        if exponent&1:
            result=tuple(result[2*r]*matrix[s]+result[2*r+1]*matrix[2+s] for r in range(2) for s in range(2))
        exponent>>=1
        if exponent:
            matrix=tuple(matrix[2*r]*matrix[s]+matrix[2*r+1]*matrix[2+s] for r in range(2) for s in range(2))
    return result


def receipt(first,last,z,alpha,marker,*,state=16,cutoff=20,bits=12,precision=256,threshold=-64,progress=None,feedthrough='cauchy',width=64,right_feedthrough=None,outer_length=256,forbidden_weight=838860,paired=True,field_order=65536,packet_source=None):
    z,alpha,marker=F(z),F(alpha),F(marker)
    if outer_length not in (64,128,256) or not isinstance(forbidden_weight,int) or not 0<=forbidden_weight<2**21:
        raise ValueError('Invalid outer length or forbidden output weight')
    rows=2**21//outer_length
    packet_width=width*(2 if paired else 1)
    if rows%packet_width or not 1<=first<=last<=rows or not 0<z<1 or not 0<alpha<=1 or not 0<marker<1:
        raise ValueError('Invalid witness or interval')
    if right_feedthrough is not None and not paired:
        raise ValueError('Alternating feedback requires paired composition')
    if field_order<2 or field_order&(field_order-1) or (field_order!=65536 and packet_source is None):
        raise ValueError('A different field requires its own exact packet source')
    if packet_source is not None and paired:raise ValueError('Custom exact source currently requires single-packet composition')
    with ctx.workprec(precision):
        if paired:
            coefficients,count=pair_coefficients(z,alpha,state,cutoff,bits,progress,min(last,width),feedthrough,width,right_feedthrough)
        else:
            groups=packet_source(z,min(last,width)) if packet_source else grouped_packets(z,state,cutoff,bits,min(last,width),feedthrough,width)
            count=len(groups)
            coefficients=[[arb(0) for _ in range(width+1)] for _ in range(4)]
            aa=_arb_fraction(alpha)
            for j,multiplicity,entries in groups:
                for i,value in enumerate(entries):
                    if value:coefficients[i][j]+=multiplicity*_arb_fraction(value)**aa
        p=_arb_fraction(marker)
        base=tuple(arb_poly([v*p**j*(1-p)**(packet_width-j) for j,v in enumerate(row[:min(last,packet_width)+1])]) for row in coefficients)
        result=(arb_poly([1]),arb_poly([]),arb_poly([]),arb_poly([1]))
        exponent=rows//packet_width
        while exponent:
            if exponent&1: result=poly_product(result,base,last)
            exponent>>=1
            if exponent: base=poly_product(base,base,last)
        log2=lambda x:x.log()/arb(2).log()
        az,aa=_arb_fraction(z),_arb_fraction(alpha)
        log_beta=(field_order.bit_length()-1)*outer_length-(outer_length//2)*log2(arb(field_order-1))
        accepted=[]; union=arb(0); rejected=[]; worst=None
        for r in range(first,last+1):
            choose=arb(comb(rows,r))
            mass=choose*p**r*(1-p)**(rows-r)
            # All exact coefficients are nonnegative. Upper endpoints give
            # a positive majorant even if an irrelevant tiny ball crosses 0.
            region=tuple(arb(max(arb(0),poly[r].upper()))/mass for poly in result)
            matrix=scalar_power(region,outer_length)
            prefactor=arb(2)**(log2(choose)+aa*(r*log_beta-log2(arb(field_order-1))-forbidden_weight*log2(az)))
            value=prefactor*(matrix[0]+matrix[1])
            bound=log2(value)
            if bound<threshold:
                accepted.append(r); union+=value
                if worst is None or bound>worst[1]: worst=(r,bound)
            else:
                rejected.append((r,str(bound)))
        intervals=[]
        for r in accepted:
            if intervals and r==intervals[-1][1]+1: intervals[-1][1]=r
            else: intervals.append([r,r])
        return dict(construction='sumrank-feedback-refresh',feedthrough=feedthrough,right_feedthrough=right_feedthrough,state=state,packet_width=width,paired=paired,mark=str(z),alpha=str(alpha),
                    outer_length=outer_length,forbidden_weight=forbidden_weight,field_order=field_order,
                    marker=str(marker),cutoff=cutoff,mantissa_bits=bits,precision=precision,packet_groups=count,
                    requested_rows=[first,last],certified_intervals=intervals,threshold=threshold,
                    union_log2_upper=str(log2(union)) if accepted else None,
                    weakest_row=worst[0] if worst else None,weakest_log2_upper=str(worst[1]) if worst else None,
                    rejected_rows=rejected,directed=True,whole_code_certificate=False)


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first',type=int,required=True)
    parser.add_argument('--last',type=int,required=True)
    parser.add_argument('--mark',type=F,required=True)
    parser.add_argument('--alpha',type=F,required=True)
    parser.add_argument('--marker',type=F,required=True)
    parser.add_argument('--state',type=int,default=16)
    parser.add_argument('--cutoff',type=int,default=20)
    parser.add_argument('--precision',type=int,default=256)
    args=parser.parse_args()
    print(json.dumps(receipt(args.first,args.last,args.mark,args.alpha,args.marker,state=args.state,cutoff=args.cutoff,precision=args.precision,
                             progress=lambda done,total:print(json.dumps(dict(groups_done=done,groups_total=total)),flush=True))),flush=True)
