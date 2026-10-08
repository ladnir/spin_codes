"""Directed loss for D2 H D1 feedback; uses the established s16 certificate."""
from flint import arb,ctx


def transfer(base_bound='-59.87'):
    with ctx.workprec(256):
        q,s,packets=65536,16,32768
        density=arb(q**s-1)/arb((q-1)**s)
        loss=packets*density.log()/arb(2).log()
        bound=arb(base_bound)+loss
        return dict(construction='sumrank-s16-diagonal-cauchy-feedback',
                    field_size=q,state=s,message_length=2**20,code_length=2**21,
                    rate='1/2',forbidden_maximum_weight=838860,
                    source_log2_bound=str(base_bound),density_cap=str(density),
                    loss_bits=str(loss),union_log2_upper=str(bound),
                    directed=True,whole_code_certificate=bool(bound < -40),
                    underlying_certificate_replayed=False,performance_claim=False)


if __name__=='__main__':
    import argparse,json
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay',action='store_true')
    args=parser.parse_args()
    base='-59.87'
    if args.replay:
        from hd_sumrank_certificate_witnesses import replay
        result=replay(lambda x:print(json.dumps(x),flush=True))
        if not result['whole_code_certificate']:raise ArithmeticError('Source certificate failed')
        base=result['union_log2_upper']
    result=transfer(base)
    result['underlying_certificate_replayed']=args.replay
    print(json.dumps(result),flush=True)
