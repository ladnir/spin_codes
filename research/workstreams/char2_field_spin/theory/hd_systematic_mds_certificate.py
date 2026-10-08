"""Find/replay directed coverage for the systematic RS[64,s] inner.

Only directed arithmetic accepts a row class; floating optimization proposes
rational witnesses. Only a complete directed union can claim a certificate.
"""
from fractions import Fraction as F
import json
from flint import arb,ctx
from hd_systematic_mds_screen import optimize
from hd_sumrank_paired_certificate import receipt as paired
from hd_sumrank_global_certificate import receipt as global_bound
from hd_randomized_outer_proof import dense_tail_receipt


def rational(x):return F(float(x)).limit_denominator(10000000)


def main(state=16,width=64,shared=False):
    ctx.prec=256
    covered=set();union=arb(0);inputs=[]
    def accept(item,kind,witness):
        nonlocal union
        for first,last in item.get('certified_intervals',[]):covered.update(range(first,last+1))
        if item.get('union_log2_upper') is not None:union+=arb(2)**arb(item['union_log2_upper'])
        inputs.append((kind,witness));item.pop('rejected_rows',None)
        print(json.dumps(dict(kind=kind,witness=witness,receipt=item,covered_rows=len(covered))),flush=True)
    tail=dense_tail_receipt()
    covered.update(range(tail['active_rows_first'],tail['active_rows_last']+1))
    union+=arb(2)**arb(tail['union_log2_upper'])
    print(json.dumps(dict(kind='dense_tail',receipt=tail,covered_rows=len(covered))),flush=True)
    def propose(r,regional):
        candidate=optimize(r,state,140,regional,width=width)
        z,a=rational(candidate['mark']),rational(candidate['alpha'])
        p=rational(1/(1+2**(-candidate['marker_log_odds'])))
        return candidate,z,a,p
    for center in (1024,1400,1800,2300,3000,4000,5500,7000):
        if center in covered:continue
        candidate,z,a,p=propose(center,False)
        result=global_bound(center,z,a,p,state=state,threshold=-72,feedthrough='systematic_mds',width=width)
        accept(result,'global',[center,str(z),str(a),str(p)])
    while len(covered)<8192:
        first=next(r for r in range(1,8193) if r not in covered)
        # The marker bound becomes stronger than the paired extraction in
        # dense classes. Choose by directed acceptance, not an assumed split.
        if first>=1000:
            candidate,z,a,p=propose(first,False)
            result=global_bound(first,z,a,p,state=state,threshold=-72,feedthrough='systematic_mds',width=width)
            if result.get('accepted'):
                accept(result,'global',[first,str(z),str(a),str(p)]);continue
        candidate,z,a,p=propose(first,True)
        print(json.dumps(dict(kind='float_proposal',proposal=candidate)),flush=True)
        if candidate['log2_upper_screen']>=-85:
            print(json.dumps(dict(kind='incomplete',first_uncovered=first,reason='screen lacks margin',whole_code_certificate=False)),flush=True);return
        last=min(8192,max(first+4,round(first*1.25)))
        result=paired(first,last,z,a,p,state=state,threshold=-72,feedthrough='systematic_mds',width=width)
        accept(result,'paired',[first,last,str(z),str(a),str(p)])
        if first not in covered:
            print(json.dumps(dict(kind='incomplete',first_uncovered=first,reason='directed witness failed',whole_code_certificate=False)),flush=True);return
    log2=lambda x:x.log()/arb(2).log()
    loss=(2097152//width)*log2((arb(65536)**state-1)/arb(65535)**state)
    if shared:loss+=(2097152//width)*log2(arb(65536)/65535)
    ideal=log2(union);final=ideal+loss
    print(json.dumps(dict(kind='complete',construction='systematic-global-MDS',state=state,packet_width=width,shared_diagonal=shared,
          covered_rows=len(covered),receipts=len(inputs)+1,ideal_union_log2_upper=str(ideal),
          density_loss_bits=str(loss),union_log2_upper=str(final),directed=True,
          whole_code_certificate=bool(final<-40),witnesses=inputs)),flush=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--state',type=int,default=16)
    p.add_argument('--width',type=int,default=64);p.add_argument('--shared',action='store_true')
    args=p.parse_args();main(args.state,args.width,args.shared)
