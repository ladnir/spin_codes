"""Directed GF256 grouped-outer/common-h4 certificate with inverse quiet bound.

This is a new field distribution, not a relabeling of a GF65536 receipt.
The packet entries include the actual diagonal density factors. Search
proposes witnesses; replay accepts only directed arithmetic and full coverage.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb

from flint import arb,ctx

from hd_bucket_feedback import profiles
from hd_common_refresh_packets import local_moments,matrix
from hd_grouped_quiet import quiet_cap
from hd_grouped_subfield_screen import proposal,LENGTH,LIMIT,ROWS
from hd_randomized_outer_proof import Parameters,dense_tail_receipt
from hd_sumrank_paired_certificate import round_up,receipt as regional_receipt
from hd_sumrank_global_certificate import receipt as global_receipt

FIELD=256
RANK=4
THRESHOLD=-80


@lru_cache(2)
def grouped_packets(z,maximum=64):
    z=F(z)
    if not 0<z<1 or not 0<=maximum<=64:raise ValueError('Invalid packet mark or cutoff')
    local=local_moments(z,q=FIELD,h=RANK)
    beta=F(FIELD**16-1,(FIELD-1)**16)
    quiet=quiet_cap(z,q=FIELD,h=RANK)
    groups={}
    for occupancy,count in profiles(16,64,RANK):
        j=sum(occupancy)
        if j>maximum:continue
        k,birth,ret,live=matrix(occupancy,z,local,q=FIELD,h=RANK)
        entries=(k,beta*birth,ret,min(beta*live,quiet) if j==0 else beta*live)
        key=(j,tuple(round_up(x,16) for x in entries))
        groups[key]=groups.get(key,0)+count
    result=tuple((j,count,entries) for (j,entries),count in groups.items())
    for j in range(maximum+1):
        if sum(count for degree,count,_ in result if degree==j)!=comb(64,j):
            raise ArithmeticError('Missing packet support profiles')
    return result


def evaluate(witness,precision=256):
    kind,values=witness
    options=dict(state=16,precision=precision,threshold=THRESHOLD,
        feedthrough='grouped-subfield-common4-inverse-quiet',width=64,paired=False,
        outer_length=256,forbidden_weight=LIMIT,field_order=FIELD,packet_source=grouped_packets)
    if kind=='regional':return regional_receipt(*values,bits=16,**options)
    if kind=='global':return global_receipt(*values,**options)
    raise ValueError('Unknown witness kind')


def initial():
    tail=dense_tail_receipt(Parameters(field=FIELD,rows=ROWS,outer_n=256,outer_k=128,
        width=64,state=16,limit=LIMIT),failure_bits=96)
    return tail,set(range(tail['active_rows_first'],ROWS+1)),arb(2)**arb(tail['union_log2_upper'])


def accept(item,covered,union):
    intervals=item.get('certified_intervals',[])
    if not intervals:raise ArithmeticError('Directed witness accepted no classes')
    for first,last in intervals:covered.update(range(first,last+1))
    return union+arb(2)**arb(item['union_log2_upper'])


def finish(tail,covered,union,witnesses,precision):
    bound=union.log()/arb(2).log()
    missing=[r for r in range(1,ROWS+1) if r not in covered]
    return dict(construction='grouped-subfield-common4-inverse-quiet',field_order=FIELD,
        message_length=LENGTH//2,output_length=LENGTH,forbidden_weight=LIMIT,
        state=16,packet_width=64,local_rank=RANK,precision=precision,
        union_log2_upper=str(bound),uncovered_rows=missing,dense_tail=tail,
        directed=True,whole_code_certificate=not missing and bool(bound<-40),witnesses=witnesses)


def replay(witnesses,precision=256,progress=None):
    with ctx.workprec(precision):
        tail,covered,union=initial()
        for index,witness in enumerate(witnesses):
            item=evaluate(witness,precision);union=accept(item,covered,union)
            if progress:progress(dict(index=index,witness=witness,receipt=item,covered=len(covered)))
        return finish(tail,covered,union,witnesses,precision)


def search(precision=256,progress=None):
    with ctx.workprec(precision):
        tail,covered,union=initial();witnesses=[]
        while len(covered)<ROWS:
            first=next(r for r in range(1,ROWS+1) if r not in covered)
            regional=first<512
            candidate=proposal(first,q=FIELD,h=RANK,quiet_inverse=True,regional=regional,iterations=90)
            values=[str(F(candidate['mark']).limit_denominator(10**7)),
                    str(F(candidate['alpha']).limit_denominator(10**7)),
                    str(F(1/(1+2**(-candidate['marker_log_odds']))).limit_denominator(10**8))]
            witness=['regional',[first,min(ROWS,first+max(4,first//3))]+values] if regional else ['global',[first]+values]
            item=evaluate(witness,precision)
            if not any(a<=first<=b for a,b in item.get('certified_intervals',[])):
                raise ArithmeticError(f'No progress at row {first}: {item}')
            union=accept(item,covered,union);witnesses.append(witness)
            if progress:progress(dict(index=len(witnesses)-1,witness=witness,receipt=item,covered=len(covered)))
        return finish(tail,covered,union,witnesses,precision)


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--precision',type=int,default=256)
    args=p.parse_args()
    print(json.dumps(search(args.precision,lambda item:print(json.dumps(item),flush=True))),flush=True)
