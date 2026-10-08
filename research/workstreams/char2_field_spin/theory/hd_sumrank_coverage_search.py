"""Select and DIRECTLY VERIFY consecutive paired-profile row intervals.

Floating calculations propose witnesses only. Every reported covered row is
checked by the independent rational/Arb construction. This driver does not
claim whole-code coverage outside its requested range.
"""

from fractions import Fraction as F
import json

from flint import arb,ctx

from hd_sumrank_profile_pairs import coefficients,exact_region_screen
from hd_sumrank_paired_certificate import receipt


def search(first,last,state=16):
    covered=set(); receipts=[]
    alphas=(.03,.04,.05,.06,.07,.08,.1,.125,.16666666666666666,.25,.5,1.0)
    while any(r not in covered for r in range(first,last+1)):
        gap=next(r for r in range(first,last+1) if r not in covered)
        center=gap
        cutoff=8 if center<=32 else 12 if center<=192 else 16 if center<=512 else 20
        best=None
        for denominator in (480,512,544):
            z=float(F(2.0**(-center/denominator)).limit_denominator(1000000))
            cs=coefficients(z,alphas,state=state,cutoff=cutoff)
            for alpha,coeff in zip(alphas,cs):
                candidate=exact_region_screen(center,z,alpha,coeff)
                if best is None or candidate['log2_upper_screen']<best['log2_upper_screen']: best=candidate
        if best['log2_upper_screen']>=-70:
            print(json.dumps(dict(stopped_at=gap,proposal=best,reason='No accepted numerical witness; proof stays open')),flush=True)
            break
        z=F(best['mark']).limit_denominator(1000000)
        alpha=F(best['alpha']).limit_denominator(10000)
        p=F(1/(1+2**(-best['marker_log_odds']))).limit_denominator(1000000)
        checked=receipt(gap,min(last,center+max(8,center//3)),z,alpha,p,state=state,cutoff=cutoff)
        # Keep the reproducible witness and verified result, not failed scans.
        checked.pop('rejected_rows')
        if not any(a<=gap<=b for a,b in checked['certified_intervals']):
            print(json.dumps(dict(stopped_at=gap,receipt=checked,reason='Proposal did not cover first gap')),flush=True)
            break
        receipts.append(checked)
        for a,b in checked['certified_intervals']: covered.update(range(a,b+1))
        print(json.dumps(checked),flush=True)
    with ctx.workprec(256):
        union=sum((arb(2)**arb(row['union_log2_upper']) for row in receipts),arb(0))
        value=str(union.log()/arb(2).log()) if receipts else None
    missing=[r for r in range(first,last+1) if r not in covered]
    print(json.dumps(dict(requested=[first,last],missing_count=len(missing),first_missing=missing[0] if missing else None,union_log2_upper=value,
                          receipts=len(receipts),directed=True,whole_code_certificate=False)),flush=True)
    return receipts,missing


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first',type=int,default=1)
    parser.add_argument('--last',type=int,default=1024)
    parser.add_argument('--state',type=int,default=16)
    args=parser.parse_args()
    search(args.first,args.last,args.state)
