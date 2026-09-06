"""All Q>1024 with at most 512 nonconstant rows, if the bound passes.

Every such message has at least 513 all-one rows, hence at least 513 ones
in every region before its permutation. Ordinary rows are bounded together
by one uniform-product counting density, not by 12 group assignments.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
from flint import arb,ctx
import bridge as base
import activation_bridge as q1
import tightened_occupancy as tight
import polynomial_regions as poly
import screen_exponential_modes as caps_module


def run(verify=False):
    output=base.HERE/'generated/constant_heavy_tail_d512_outward.json'
    old=base.read(output) if verify else None
    if not verify:assert not output.exists()
    if old:
        for name,digest in old['local_sha256'].items():assert base.sha(base.HERE/name)==digest
        caps={int(w):v for w,v in old['used_caps'].items()};sources=[]
    else:caps,sources=caps_module.latest_caps()
    gamma=max(F(caps[w]*(1<<256),math.comb(256,w)) for w in base.WEIGHTS if w!=256)
    assert gamma>=1
    t,s,spectrum=base.load_map(tight.NAME);ctx.prec=512 if verify else 256
    lam=(arb(-12)/10).exp();region=poly.regions(t,s,spectrum,tight.kernel_spectrum(),(-lam).exp(),8192)
    envelope=tuple(max(row[k].upper() for row in region[513:]) for k in range(9))
    moment=envelope
    for _ in range(8):moment=q1.positive_mul(moment,moment)
    # The h-count restriction may only decrease 2^(L-d). The moment bound
    # was applied to the eligible h before this count was enlarged.
    counting=sum((math.comb(8192,d)*(1<<(8192-d))*gamma**d for d in range(513)),F(0))
    from audit_bch_q1_full_arb import rational
    bound=counting*rational((sum(moment[:3],arb(0))*(209716*lam).exp()).upper())
    assert bound>0
    if old:
        assert gamma==base.decode(old['ordinary_uniform_density_upper'])
        assert counting==base.decode(old['position_and_row_density_sum'])
        assert bound<=base.decode(old['failure_contribution_upper'])
        print('512-bit constant-heavy tail replay passed',flush=True);return
    margin=math.log2(bound.denominator)-math.log2(bound.numerator)
    local=[Path(__file__),Path(poly.__file__),Path(tight.__file__),Path(caps_module.__file__),
           base.HERE/'general_occupancy.py',base.HERE/'activation_bridge.py',base.HERE/'christoffel_caps.py',
           base.HERE/'generated/christoffel_oa29_caps.json']+sources
    payload=dict(status='OUTWARD_CONSTANT_HEAVY_TAIL_BOUND',configuration=tight.NAME,
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,step_bits=t,state_bits=s,distance_cutoff=209716),
        scope=dict(total_nonzero_rows_minimum=1025,ordinary_nonconstant_rows_maximum=512,all_one_rows_minimum=513),
        used_caps={str(w):caps[w] for w in base.WEIGHTS},ordinary_uniform_density_upper=base.encode(gamma),
        position_and_row_density_sum=base.encode(counting),failure_contribution_upper=base.encode(bound),
        margin_bits_diagnostic=margin,below_2_to_minus_60=bound<F(1,1<<60),all_occupations_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local})
    base.write_new(output,payload)
    print('Constant-heavy full tail margin',margin,'density bits',math.log2(gamma.numerator)-math.log2(gamma.denominator),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
