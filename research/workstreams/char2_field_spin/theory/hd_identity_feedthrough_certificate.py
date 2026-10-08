"""Adaptive directed coverage for the systematic candidate; emits proof inputs.

Only the final record can claim a whole-code certificate. It independently
unions returned intervals and directed failure bounds. The diagonal-Cauchy
density loss is included after the ideal-feedback certificate is summed.
"""
from fractions import Fraction as F
import json
from flint import arb,ctx
from hd_identity_feedthrough_screen import optimize,packet_logs
from hd_identity_feedthrough_pairs import screen
from hd_sumrank_profile_pairs import coefficients,exact_region_screen
from hd_sumrank_paired_certificate import receipt as paired
from hd_sumrank_global_certificate import receipt as global_bound
from hd_randomized_outer_proof import dense_tail_receipt


def rational(x):return F(float(x)).limit_denominator(1000000)


def improve(r,result):
    if result['log2_upper_screen']<-130:return result
    alphas=[min(1,result['alpha']*x) for x in (.9,.95,1,1.05,1.1)]
    for scale in (.96,.98,1,1.02,1.04):
        z=result['mark']**scale
        cs=coefficients(z,alphas,state=16,cutoff=result['cutoff'],packet_source=packet_logs)
        for a,c in zip(alphas,cs):
            candidate=exact_region_screen(r,z,a,c)
            if candidate['log2_upper_screen']<result['log2_upper_screen']:
                candidate['cutoff']=result['cutoff'];result=candidate
    return result


def main():
    ctx.prec=256
    covered=set();union=arb(0);count=0;inputs=[]
    def accept(item,kind,witness):
        nonlocal union,count
        intervals=item.get('certified_intervals',[])
        for first,last in intervals:covered.update(range(first,last+1))
        if item.get('union_log2_upper') is not None:union+=arb(2)**arb(item['union_log2_upper'])
        count+=1;inputs.append((kind,witness))
        item.pop('rejected_rows',None)
        print(json.dumps(dict(kind=kind,witness=witness,receipt=item,covered_rows=len(covered))),flush=True)
    tail=dense_tail_receipt()
    covered.update(range(tail['active_rows_first'],tail['active_rows_last']+1))
    union+=arb(2)**arb(tail['union_log2_upper'])
    print(json.dumps(dict(kind='dense_tail',receipt=tail,covered_rows=len(covered))),flush=True)
    for center in (1536,1800,2048,2500,3072,4096,4500,6144,7400):
        if center in covered:continue
        candidate=optimize(center,16,80)
        z,a=rational(candidate['mark']),rational(candidate['alpha'])
        p=rational(1/(1+2**(-candidate['marker_log_odds'])))
        result=global_bound(center,z,a,p,threshold=-68,feedthrough='identity')
        accept(result,'global',[center,str(z),str(a),str(p)])
    while len(covered)<8192:
        first=next(r for r in range(1,8193) if r not in covered)
        candidate=improve(first,screen(first))
        print(json.dumps(dict(kind='float_proposal',proposal=candidate)),flush=True)
        if candidate['log2_upper_screen']>=-75:
            print(json.dumps(dict(kind='incomplete',first_uncovered=first,reason='screen lacks margin',whole_code_certificate=False)),flush=True)
            return
        z,a=rational(candidate['mark']),rational(candidate['alpha'])
        p=rational(1/(1+2**(-candidate['marker_log_odds'])))
        last=min(8192,max(first+4,round(first*1.25)))
        result=paired(first,last,z,a,p,cutoff=candidate['cutoff'],threshold=-68,feedthrough='identity')
        accept(result,'paired',[first,last,str(z),str(a),str(p),candidate['cutoff']])
        if first not in covered:
            print(json.dumps(dict(kind='incomplete',first_uncovered=first,reason='directed witness failed',whole_code_certificate=False)),flush=True)
            return
    log2=lambda x:x.log()/arb(2).log()
    loss=32768*log2((arb(65536)**16-1)/arb(65535)**16)
    ideal=log2(union);final=ideal+loss
    print(json.dumps(dict(kind='complete',construction='systematic-sumrank-s16',
          covered_rows=len(covered),receipts=count,ideal_union_log2_upper=str(ideal),
          density_loss_bits=str(loss),union_log2_upper=str(final),directed=True,
          whole_code_certificate=bool(final<-40),witnesses=inputs)),flush=True)


if __name__=='__main__':main()
