"""512-bit conditional lower replay and exact first-moment obstruction ledger."""
import math
from pathlib import Path
from fractions import Fraction as F
from flint import arb,arb_poly,ctx
import bridge as base
import tightened_occupancy as tight


def replay():
    path=base.HERE/'generated/zero_state_q2620_conditional_lower.json';saved=base.read(path)
    assert saved['status']=='CONDITIONAL_ACTUAL_ZERO_STATE_FIRST_MOMENT_LOWER' and saved['configuration']==tight.NAME
    for name,digest in saved['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    assert saved['occupation']==2620 and saved['outer_row_weights']==[38,80] and saved['region_weight_interval']==[64,1536]
    assert saved['input_weight_upper']==209600<=209716
    t,s,spectrum=base.load_map(tight.NAME);assert (t,s)==(128,15);ctx.prec=512
    kernel=tight.kernel_spectrum();poly=arb_poly([arb(kernel.get(j,0)) for j in range(129)])**64
    minimum=min((poly[j]/math.comb(8192,j)).lower() for j in range(64,1537,2));assert minimum>0
    q=2620
    def tail(a,p):return (-q*(a*(a/p).log()+(1-a)*((1-a)/(1-p)).log())).exp()
    bad=256*(tail(arb(63)/q,arb(38)/256)+tail(arb(1537)/q,arb(80)/256))
    parity=(arb(2)**(-255)-bad).lower();assert parity>0
    from audit_bch_q1_full_arb import rational
    factor=rational((math.comb(8192,q)*parity*minimum**256).lower())
    assert factor>=base.decode(saved['factor_lower'])>0
    print('512-bit actual zero-state conditional lower bound replay passed',flush=True)
    tail_path=base.HERE/'generated/joint_tail_lower_80/lower.json';tail_receipt=base.read(tail_path)
    assert tail_receipt['status']=='EXACT_BCH_TAIL_LOWER' and tail_receipt['rational_primal_dual_checks_passed']
    assert tail_receipt['weights']==list(range(38,81,2))
    for name,digest in tail_receipt['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    for name,digest in tail_receipt['outer_sha256'].items():assert base.sha(base.BCH/name)==digest
    lower=base.decode(saved['factor_lower'])*tail_receipt['lower']**q
    assert lower>F(1<<19600)
    print('Actual bad-message first moment at Q2620 exceeds 2^19600; log2 lower',
        math.log2(lower.numerator)-math.log2(lower.denominator),flush=True)
    output=base.HERE/'generated/zero_state_first_moment_obstruction.json'
    result=dict(status='ACTUAL_FIRST_MOMENT_OBSTRUCTION_NOT_FAILURE_PROBABILITY',configuration=tight.NAME,
        occupation=q,first_moment_lower=base.encode(lower),first_moment_exceeds_2_to_19600=True,
        setup_failure_probability_lower_bound=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),path,tail_path]})
    if output.exists():assert result==base.read(output)
    else:base.write_new(output,result)


if __name__=='__main__':replay()
