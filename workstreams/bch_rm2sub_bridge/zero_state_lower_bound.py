"""Outward lower bound for the actual zero-state branch, conditional on a BCH tail count.

Q even independent rows from any nonempty set of even BCH words of weights
38..80 have all 256 column parities even with probability >=2^-255.
Restrict column counts to 64..1536, paying explicit binomial tail bounds.
Uniform region permutations then hit (ker B)^64 with an exact weight ratio.
This yields a FIRST-MOMENT lower bound, not a setup-failure probability.
"""
import math
from pathlib import Path
from fractions import Fraction as F
from flint import arb,arb_poly,ctx
import bridge as base
import tightened_occupancy as tight


def run():
    output=base.HERE/'generated/zero_state_q2620_conditional_lower.json';assert not output.exists()
    t,s,spectrum=base.load_map(tight.NAME);ctx.prec=256
    q=2620;lo=64;hi=1536;assert q%2==0 and 80*q<=209716
    kernel=tight.kernel_spectrum();poly=arb_poly([arb(kernel.get(j,0)) for j in range(129)])**64
    minimum=min((poly[j]/math.comb(8192,j)).lower() for j in range(lo,hi+1,2))
    assert minimum>0
    def tail(a,p):
        divergence=a*(a/p).log()+(1-a)*((1-a)/(1-p)).log()
        return (-q*divergence).exp()
    bad=256*(tail(arb(lo-1)/q,arb(38)/256)+tail(arb(hi+1)/q,arb(80)/256))
    parity=(arb(2)**(-255)-bad).lower();assert parity>0
    from audit_bch_q1_full_arb import rational
    factor=rational((math.comb(8192,q)*parity*minimum**256).lower());assert factor>0
    threshold=-float((arb(factor.numerator)/factor.denominator).log())/(q*math.log(2))
    local=[Path(__file__),Path(tight.__file__),base.HERE/'general_occupancy.py',
        base.HERE/'inputs/t128_s15_selection.json',base.HERE/'inputs/t128_s15_a_spectrum.json',
        base.HERE/'inputs/t128_s15_b_kernel_spectrum.json']
    base.write_new(output,dict(status='CONDITIONAL_ACTUAL_ZERO_STATE_FIRST_MOMENT_LOWER',configuration=tight.NAME,
        occupation=q,outer_row_weights=[38,80],input_weight_upper=80*q,region_weight_interval=[lo,hi],
        factor_lower=base.encode(factor),formula='E[bad messages with Q rows in the tail] >= factor_lower * tail_cardinality^Q',
        required_log2_tail_for_mean_above_one=threshold,
        minimum_region_probability_log2=float(minimum.log())/math.log(2),
        parity_good_probability_log2=float(parity.log())/math.log(2),
        actual_tail_lower_bound_supplied=False,setup_failure_probability_lower_bound=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local}))
    print('Actual zero-state conditional bound: kernel region bits',float(minimum.log())/math.log(2),
        'good parity bits',float(parity.log())/math.log(2),'tail log2 threshold',threshold,flush=True)


if __name__=='__main__':run()
