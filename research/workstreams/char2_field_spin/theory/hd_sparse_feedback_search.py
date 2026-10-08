"""Directed coverage search for the fixed (4,4,4,8) feedback schedule.

Existing witnesses are proposals only. Every receipt recomputes the new
period polynomial. The final claim requires all 8192 active-row counts
and applies the shared-refresh density factor once for all 32768 packets.
"""
from fractions import Fraction as F
import json
from flint import arb,ctx
from hd_bucket_feedback_witnesses import WITNESSES
from hd_randomized_outer_proof import dense_tail_receipt
from hd_sparse_feedback_certificate import receipt
from hd_sparse_feedback_screen import screen


def search(precision=256):
    ctx.prec=precision;covered=set();union=arb(0);accepted=[]
    def accept(witness,cutoff=64):
        nonlocal union
        item=receipt(*witness,precision=precision,cutoff=cutoff)
        if not item['certified_intervals']:return False
        for first,last in item['certified_intervals']:covered.update(range(first,last+1))
        union+=arb(2)**arb(item['union_log2_upper']);accepted.append([witness,cutoff])
        item.pop('rejected_rows',None)
        print(json.dumps(dict(kind='paired',witness=witness,cutoff=cutoff,
                              receipt=item,covered_rows=len(covered))),flush=True)
        return True
    tail=dense_tail_receipt()
    covered.update(range(tail['active_rows_first'],tail['active_rows_last']+1))
    union+=arb(2)**arb(tail['union_log2_upper'])
    print(json.dumps(dict(kind='dense_tail',receipt=tail)),flush=True)
    for kind,witness,cutoff in WITNESSES:
        if kind=='paired':proposal=witness
        elif kind=='global-paired':
            center,*values=witness
            proposal=[center,min(8192,round(center*1.2)),*values]
        else:raise ValueError('Unknown proposal')
        if not all(r in covered for r in range(proposal[0],proposal[1]+1)):accept(proposal,cutoff)
    while len(covered)<8192:
        first=next(r for r in range(1,8193) if r not in covered)
        proposal=screen(first,iterations=130,regional=True)
        print(json.dumps(dict(kind='proposal',proposal=proposal)),flush=True)
        if proposal['log2_upper_screen']>=-85:
            print(json.dumps(dict(kind='incomplete',first_uncovered=first,witnesses=accepted,
                                  whole_code_certificate=False)),flush=True);return
        rational=lambda x:str(F(float(x)).limit_denominator(10000000))
        z=rational(proposal['mark']);a=rational(proposal['alpha'])
        p=rational(1/(1+2**(-proposal['marker_log_odds'])))
        witness=[first,min(8192,max(first+4,round(first*1.2))),z,a,p]
        accept(witness)
        if first not in covered:
            print(json.dumps(dict(kind='incomplete',first_uncovered=first,witnesses=accepted,
                                  whole_code_certificate=False)),flush=True);return
    log2=lambda x:x.log()/arb(2).log();q=arb(65536)
    loss=32768*log2(q/(q-1)*(q**16-1)/(q-1)**16)
    final=log2(union)+loss
    print(json.dumps(dict(kind='complete',construction='global-emission-feedback-4-4-4-8',
                         covered_rows=len(covered),precision=precision,
                         ideal_union_log2_upper=str(log2(union)),density_loss_bits=str(loss),
                         union_log2_upper=str(final),directed=True,
                         whole_code_certificate=bool(final<-40),witnesses=accepted)),flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision',type=int,default=256)
    args=parser.parse_args();search(args.precision)
