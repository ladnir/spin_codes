"""Directed packetwise witnesses for the larger-state feedback candidate."""

from fractions import Fraction as F
from math import comb

from flint import arb,ctx

from hd_sumrank_paired_certificate import grouped_packets,scalar_power,pair_coefficients
from hd_randomized_outer_quiet_exact import _arb_fraction


def witness(z,alpha,marker,state=16,precision=256,feedthrough='cauchy',width=64,paired=False,right_feedthrough=None,outer_length=256,forbidden_weight=838860,field_order=65536,packet_source=None):
    z,alpha,marker=F(z),F(alpha),F(marker)
    if outer_length not in (64,128,256) or not isinstance(forbidden_weight,int) or not 0<=forbidden_weight<2**21:
        raise ValueError('Invalid outer length or forbidden output weight')
    if not 0<z<1 or not 0<alpha<=1 or not 0<marker<1 or (2**21//outer_length)%(width*(2 if paired else 1)):
        raise ValueError('Invalid witness or packet geometry')
    if right_feedthrough is not None and not paired:raise ValueError('Alternating feedback requires paired composition')
    if field_order<2 or field_order&(field_order-1) or (field_order!=65536 and packet_source is None):
        raise ValueError('A different field requires its own exact packet source')
    if packet_source is not None and paired:raise ValueError('Custom exact source currently requires single-packet composition')
    with ctx.workprec(precision):
        p,a=_arb_fraction(marker),_arb_fraction(alpha)
        matrix=[arb(0) for _ in range(4)]
        # All profiles are kept; dyadic rounding is explicitly upward.
        packet_width=2*width if paired else width
        if paired:
            coefficients,_=pair_coefficients(z,alpha,state,cutoff=width,feedthrough=feedthrough,width=width,right_feedthrough=right_feedthrough)
            for i,row in enumerate(coefficients):
                matrix[i]=sum(value*p**j*(1-p)**(packet_width-j) for j,value in enumerate(row))
        else:
            groups=packet_source(z,width) if packet_source else grouped_packets(z,state,width,feedthrough=feedthrough,width=width)
            for j,count,entries in groups:
                weight=count*p**j*(1-p)**(width-j)
                for i,value in enumerate(entries):
                    if value: matrix[i]+=weight*_arb_fraction(value)**a
        run=scalar_power(tuple(matrix),2097152//packet_width)
        log2=lambda x:x.log()/arb(2).log()
        intercept=(log2(run[0]+run[1])-a*(log2(arb(field_order-1))+forbidden_weight*log2(_arb_fraction(z)))
                   -2097152*log2(1-p))
        slope=a*((field_order.bit_length()-1)*outer_length-(outer_length//2)*log2(arb(field_order-1)))-outer_length*log2(p)+outer_length*log2(1-p)
        return slope,intercept


def receipt(center,z,alpha,marker,state=16,precision=256,threshold=-64,feedthrough='cauchy',width=64,paired=False,right_feedthrough=None,outer_length=256,forbidden_weight=838860,field_order=65536,packet_source=None):
    rows=2097152//outer_length
    if not 1<=center<=rows:raise ValueError('Invalid center row count')
    with ctx.workprec(precision):
        slope,intercept=witness(z,alpha,marker,state,precision,feedthrough,width,paired,right_feedthrough,outer_length,forbidden_weight,field_order,packet_source)
        def bound(r):
            return r*slope+intercept-(outer_length-1)*arb(comb(rows,r)).log()/arb(2).log()
        if not bound(center)<threshold:
            return dict(center=center,accepted=False,log2_upper=str(bound(center)),whole_code_certificate=False)
        left,right=center,center
        while left>1 and bound(left-1)<threshold: left-=1
        while right<rows and bound(right+1)<threshold: right+=1
        worst=max(bound(left),bound(right))
        union=worst+arb(right-left+1).log()/arb(2).log()
        return dict(construction='sumrank-feedback-refresh',feedthrough=feedthrough,right_feedthrough=right_feedthrough,state=state,packet_width=width,paired=paired,mark=str(F(z)),alpha=str(F(alpha)),
                    outer_length=outer_length,forbidden_weight=forbidden_weight,field_order=field_order,
                    marker=str(F(marker)),certified_intervals=[[left,right]],center=center,
                    union_log2_upper=str(union),weakest_log2_upper=str(worst),accepted=True,
                    directed=True,whole_code_certificate=False)


if __name__=='__main__':
    import argparse,json
    from hd_sumrank_fast_screen import optimize
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--centers',type=int,nargs='+',default=(1024,1536,2048,3072,4096,6144))
    args=parser.parse_args()
    for center in args.centers:
        screen=optimize(center,state=16,iterations=60)
        z=F(screen['mark']).limit_denominator(1000000)
        alpha=F(screen['alpha']).limit_denominator(1000000)
        p=F(1/(1+2**(-screen['marker_log_odds']))).limit_denominator(1000000)
        print(json.dumps(receipt(center,z,alpha,p)),flush=True)
