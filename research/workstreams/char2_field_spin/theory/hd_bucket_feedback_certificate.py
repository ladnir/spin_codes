"""Directed coverage search for alternating four/eight-check feedback.

FLOAT optimization only proposes witnesses. A whole-code claim requires
complete directed coverage and the final shared-refresh density loss.
"""
from fractions import Fraction as F
import json
from flint import arb,ctx
from hd_bucket_feedback_screen import packet_logs,paired_coefficients
from hd_systematic_mds_screen import optimize
from hd_systematic_mds_witnesses import WITNESSES
from hd_sumrank_paired_certificate import receipt as paired
from hd_sumrank_global_certificate import receipt as global_bound
from hd_randomized_outer_proof import dense_tail_receipt


def rational(x):return F(float(x)).limit_denominator(10000000)


def search(precision=256):
    ctx.prec=precision;covered=set();union=arb(0);accepted=[]
    common=dict(state=16,feedthrough='bucket4',right_feedthrough='bucket8',precision=precision,threshold=-72)
    def accept(item,kind,witness,cutoff=64):
        nonlocal union
        if not item.get('certified_intervals'):return False
        for first,last in item['certified_intervals']:covered.update(range(first,last+1))
        union+=arb(2)**arb(item['union_log2_upper'])
        accepted.append([kind,witness,cutoff]);item.pop('rejected_rows',None)
        print(json.dumps(dict(kind=kind,witness=witness,cutoff=cutoff,receipt=item,covered_rows=len(covered))),flush=True)
        return True
    tail=dense_tail_receipt()
    covered.update(range(tail['active_rows_first'],tail['active_rows_last']+1))
    union+=arb(2)**arb(tail['union_log2_upper'])
    print(json.dumps(dict(kind='dense_tail',receipt=tail)),flush=True)
    # Reusing a witness is not reusing its old guarantee: recompute every
    # operator and coefficient for the new alternating construction.
    for kind,witness in WITNESSES:
        if kind!='paired':continue
        if all(r in covered for r in range(witness[0],witness[1]+1)):continue
        accept(paired(*witness,cutoff=20,**common),'paired',witness,20)
    def propose(r,regional):
        source=lambda z,s,w:packet_logs(z,s,w,4)
        pairs=lambda z,a,s,w:paired_coefficients(z,a,s,w,64,4,8)
        item=optimize(r,16,130,regional,packet_source=source,pair_source=pairs)
        z,a=rational(item['mark']),rational(item['alpha'])
        p=rational(1/(1+2**(-item['marker_log_odds'])))
        return item,z,a,p
    for center in (1536,2048,2560,3072,4096,5500,7000):
        if center in covered:continue
        proposal,z,a,p=propose(center,False)
        witness=[center,str(z),str(a),str(p)]
        accept(global_bound(*witness,paired=True,**common),'global-paired',witness)
    while len(covered)<8192:
        first=next(r for r in range(1,8193) if r not in covered)
        regional=first<1536
        proposal,z,a,p=propose(first,regional)
        print(json.dumps(dict(kind='proposal',proposal=proposal)),flush=True)
        if not regional and proposal['log2_upper_screen']<-85:
            witness=[first,str(z),str(a),str(p)]
            item=global_bound(*witness,paired=True,**common)
            if accept(item,'global-paired',witness) and first in covered:continue
        if not regional:proposal,z,a,p=propose(first,True)
        if proposal['log2_upper_screen']>=-85:
            print(json.dumps(dict(kind='incomplete',first_uncovered=first,proposal=proposal,
                                  witnesses=accepted,whole_code_certificate=False)),flush=True);return
        last=min(8192,max(first+4,round(first*1.2)))
        witness=[first,last,str(z),str(a),str(p)]
        accept(paired(*witness,cutoff=64,**common),'paired',witness)
        if first not in covered:
            print(json.dumps(dict(kind='incomplete',first_uncovered=first,witnesses=accepted,
                                  whole_code_certificate=False)),flush=True);return
    log2=lambda x:x.log()/arb(2).log()
    loss=32768*log2(arb(65536)/65535*(arb(65536)**16-1)/arb(65535)**16)
    result=log2(union)+loss
    print(json.dumps(dict(kind='complete',construction='alternating-block-feedback-4-8',covered_rows=len(covered),
                         precision=precision,ideal_union_log2_upper=str(log2(union)),density_loss_bits=str(loss),
                         union_log2_upper=str(result),directed=True,whole_code_certificate=bool(result<-40),
                         witnesses=accepted)),flush=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--precision',type=int,default=256)
    args=p.parse_args();search(args.precision)
