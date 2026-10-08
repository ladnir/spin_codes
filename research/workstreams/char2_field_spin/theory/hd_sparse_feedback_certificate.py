"""Directed regional bounds for periodic 4/8-check feedback schedules.

These receipts certify only the explicitly returned active-row classes.
A complete certificate additionally needs gap-free coverage, a union bound,
and the shared-refresh density loss. No complete certificate is asserted here.
"""
from fractions import Fraction as F
from math import comb
from flint import arb,arb_poly,ctx
from hd_randomized_outer_quiet_exact import _arb_fraction
from hd_sumrank_paired_certificate import pair_coefficients,poly_product,scalar_power


def period_coefficients(z,alpha,cutoff=64,maximum=64,strong_pairs=1):
    if strong_pairs not in (1,3):raise ValueError('Use two or four pairs per period')
    common=dict(state=16,cutoff=cutoff,maximum=maximum,feedthrough='bucket4')
    left,n=pair_coefficients(z,alpha,**common)
    right,m=pair_coefficients(z,alpha,right_feedthrough='bucket8',**common)
    result=tuple(arb_poly(row) for row in left)
    factor=tuple(arb_poly(row) for row in right)
    for i in range(strong_pairs):result=poly_product(result,factor,128*(i+2))
    return result,(n,m)


def receipt(first,last,z,alpha,marker,*,precision=256,threshold=-72,cutoff=64,strong_pairs=1):
    z,alpha,marker=F(z),F(alpha),F(marker)
    if not 1<=first<=last<=8192 or not 0<z<1 or not 0<alpha<=1 or not 0<marker<1:
        raise ValueError('Invalid witness')
    with ctx.workprec(precision):
        coefficients,counts=period_coefficients(z,alpha,cutoff,min(last,64),strong_pairs)
        period=128*(strong_pairs+1);p=_arb_fraction(marker)
        base=tuple(arb_poly([row[j]*p**j*(1-p)**(period-j)
                            for j in range(min(last,period)+1)]) for row in coefficients)
        result=(arb_poly([1]),arb_poly([]),arb_poly([]),arb_poly([1]))
        exponent=8192//period
        while exponent:
            if exponent&1:result=poly_product(result,base,last)
            exponent>>=1
            if exponent:base=poly_product(base,base,last)
        log2=lambda x:x.log()/arb(2).log()
        aa,az=_arb_fraction(alpha),_arb_fraction(z)
        beta=4096-128*log2(arb(65535))
        accepted=[];rejected=[];union=arb(0);worst=None
        for r in range(first,last+1):
            choose=arb(comb(8192,r));mass=choose*p**r*(1-p)**(8192-r)
            region=tuple(arb(max(arb(0),row[r].upper()))/mass for row in result)
            run=scalar_power(region,256)
            prefactor=arb(2)**(log2(choose)+aa*(r*beta-log2(arb(65535))-838860*log2(az)))
            value=prefactor*(run[0]+run[1]);bound=log2(value)
            if bound<threshold:
                accepted.append(r);union+=value
                if worst is None or bound>worst[1]:worst=(r,bound)
            else:rejected.append([r,str(bound)])
        intervals=[]
        for r in accepted:
            if intervals and r==intervals[-1][1]+1:intervals[-1][1]=r
            else:intervals.append([r,r])
        return dict(construction='global-emission-sparser-eight-check-feedback',
                    feedback_schedule=[4,4]+[4,8]*strong_pairs,
                    packet_width=64,composition_width=period,state=16,
                    requested_rows=[first,last],certified_intervals=intervals,
                    mark=str(z),alpha=str(alpha),marker=str(marker),
                    precision=precision,cutoff=cutoff,packet_groups=counts,
                    threshold=threshold,union_log2_upper=str(log2(union)) if accepted else None,
                    weakest_row=worst[0] if worst else None,
                    weakest_log2_upper=str(worst[1]) if worst else None,
                    rejected_rows=rejected,directed=True,whole_code_certificate=False)


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first',type=int,required=True)
    parser.add_argument('--last',type=int,required=True)
    parser.add_argument('--mark',type=F,required=True)
    parser.add_argument('--alpha',type=F,required=True)
    parser.add_argument('--marker',type=F,required=True)
    parser.add_argument('--precision',type=int,default=256)
    parser.add_argument('--strong-pairs',type=int,choices=(1,3),default=1)
    args=parser.parse_args()
    print(json.dumps(receipt(args.first,args.last,args.mark,args.alpha,args.marker,
                             precision=args.precision,strong_pairs=args.strong_pairs)),flush=True)
