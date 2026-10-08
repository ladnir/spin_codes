"""Directed certificate transfer for shared D2 in M and T; see the proof note."""
from flint import arb,ctx


def transfer(base_bound='-50.06',precision=256):
    with ctx.workprec(precision):
        loss=32768*(arb(65536)/65535).log()/arb(2).log()
        final=arb(base_bound)+loss
        return dict(construction='systematic-RS-shared-diagonal',
                    field_size=65536,state=16,message_length=2**20,code_length=2**21,
                    rate='1/2',forbidden_maximum_weight=838860,
                    source_log2_bound=str(base_bound),additional_loss_bits=str(loss),
                    union_log2_upper=str(final),precision=precision,directed=True,
                    whole_code_certificate=bool(final < -40),
                    underlying_certificate_replayed=False,performance_claim=False)


if __name__=='__main__':
    import argparse,json
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--precision',type=int,default=256)
    p.add_argument('--replay',action='store_true')
    args=p.parse_args();base='-50.06'
    if args.replay:
        from hd_systematic_mds_witnesses import replay
        result=replay(args.precision,lambda i,n,c:print(json.dumps(dict(done=i,total=n,covered=c)),flush=True))
        if not result['whole_code_certificate']:raise ArithmeticError('Source certificate failed')
        base=result['union_log2_upper']
    result=transfer(base,args.precision)
    result['underlying_certificate_replayed']=args.replay
    print(json.dumps(result),flush=True)
