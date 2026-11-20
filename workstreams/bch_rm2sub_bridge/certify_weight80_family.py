"""An exact-weight actual BCH family for the second-moment investigation.

All Q=2620 occupied rows now have weight exactly 80. This fixes total input
weight and removes variation of row weights without changing the encoder.
"""
import math
from fractions import Fraction as F
from pathlib import Path
from flint import arb,ctx
import bridge as base
import christoffel_caps as christoffel
import screen_exponential_modes as cap_source


def run():
    caps,sources=cap_source.latest_caps()
    tail_path=base.HERE/'generated/joint_tail_lower_80/lower.json';tail=base.read(tail_path)
    for name,digest in tail['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    for name,digest in tail['outer_sha256'].items():assert base.sha(base.BCH/name)==digest
    a80=tail['lower']-sum(caps[w] for w in range(38,80,2));assert 0<a80<=caps[80]
    affine_path=base.HERE/'generated/zero_state_affine_lower.json';affine=base.read(affine_path)
    for name,digest in affine['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    assert affine['occupation']==2620 and affine['configuration']=='t128_s15'
    mean_lower=math.comb(8192,2620)*a80**2620*base.decode(affine['averaged_zero_state_probability_lower'])
    assert mean_lower>F(1<<19000)
    ctx.prec=512;christoffel.build(verify=True)
    from audit_bch_q1_full_arb import rational
    close=F(1,a80);rows=[]
    for c in range(38,81,2):
        i=c//2;j=80-i
        kernel=sum((arb(christoffel.kraw(c,r,i)**2*christoffel.kraw(256-c,s,j)**2)
            /(math.comb(c,r)*math.comb(256-c,s)) for r in range(15) for s in range(15-r)),arb(0)).lower()
        assert kernel>0
        intersection=min(caps[80],math.comb(c,i)*math.comb(256-c,j),math.floor(rational((arb(2)**128/kernel).upper())))
        close+=F(caps[c]*intersection,a80*a80);rows.append(dict(difference_weight=c,intersection_upper=intersection))
    assert close<F(1,128)
    result=dict(status='ACTUAL_EXACT_WEIGHT80_SECOND_MOMENT_FAMILY',configuration='t128_s15',occupation=2620,
        row_weight=80,total_input_weight=209600,shell_size_lower=a80,
        first_moment_lower=base.encode(mean_lower),close_pair_probability_upper=base.encode(close),rows=rows,
        full_second_moment_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),Path(christoffel.__file__),tail_path,affine_path]+sources})
    path=base.HERE/'generated/weight80_family.json'
    if path.exists():assert result==base.read(path)
    else:base.write_new(path,result)
    print('Exact weight80 family: shell lower bits',math.log2(a80),'mean lower bits',math.log2(mean_lower.numerator)-math.log2(mean_lower.denominator),
        'close-pair bits',math.log2(close.denominator)-math.log2(close.numerator),flush=True)


if __name__=='__main__':run()
