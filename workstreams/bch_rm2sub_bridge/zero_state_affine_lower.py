"""Sharper actual zero-state first moment using the conserved input weight.

An affine lower bound on log(beta_j) is checked at every eligible even j.
Its negative slope lets the total input-weight cap be used without a
monotonicity or convexity assumption about the sequence beta_j.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
from flint import arb,arb_poly,ctx
import bridge as base
import tightened_occupancy as tight


def run(verify=False):
    path=base.HERE/'generated/zero_state_affine_lower.json';old=base.read(path) if verify else None
    if not verify:assert not path.exists()
    if old:
        for name,digest in old['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    ctx.prec=512 if verify else 256;base.load_map(tight.NAME)
    kernel=tight.kernel_spectrum();poly=arb_poly([arb(kernel.get(j,0)) for j in range(129)])**64
    logs={j:(poly[j]/math.comb(8192,j)).log() for j in range(64,1537,2)}
    slope=base.decode(old['slope']) if old else F.from_float(float((logs[820]-logs[818])/2))
    assert slope<0;b=arb(slope.numerator)/slope.denominator
    intercept=min((value-b*j).lower() for j,value in logs.items())
    q=2620
    def tail(a,p):return (-q*(a*(a/p).log()+(1-a)*((1-a)/(1-p)).log())).exp()
    gamma=arb(2)**-255-256*(tail(arb(63)/q,arb(38)/256)+tail(arb(1537)/q,arb(80)/256))
    assert gamma>arb(2)**-256
    tail_path=base.HERE/'generated/joint_tail_lower_80/lower.json';receipt=base.read(tail_path)
    for name,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    for name,digest in receipt['outer_sha256'].items():assert base.sha(base.BCH/name)==digest
    assert receipt['rational_primal_dual_checks_passed'];a0=receipt['lower']
    from audit_bch_q1_full_arb import rational
    probability=rational((arb(2)**-256*(256*intercept+b*80*q).exp()).lower())
    lower=math.comb(8192,q)*a0**q*probability
    previous=base.read(base.HERE/'generated/zero_state_first_moment_obstruction.json')
    assert lower>base.decode(previous['first_moment_lower'])
    if old:
        assert probability>=base.decode(old['averaged_zero_state_probability_lower'])
        assert lower>=base.decode(old['first_moment_lower'])
        print('512-bit affine first-moment lower replay passed',flush=True);return
    local=[Path(__file__),Path(tight.__file__),base.HERE/'general_occupancy.py',tail_path,
        base.HERE/'inputs/t128_s15_selection.json',base.HERE/'inputs/t128_s15_a_spectrum.json',
        base.HERE/'inputs/t128_s15_b_kernel_spectrum.json']
    base.write_new(path,dict(status='ACTUAL_AFFINE_ZERO_STATE_FIRST_MOMENT_LOWER',configuration=tight.NAME,
        occupation=q,region_even_interval=[64,1536],slope=base.encode(slope),
        intercept_lower=base.encode(rational(intercept)),
        averaged_zero_state_probability_lower=base.encode(probability),first_moment_lower=base.encode(lower),
        setup_failure_probability_lower_bound=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local}))
    print('Affine zero-state lower: probability bits',math.log2(probability.numerator)-math.log2(probability.denominator),
        'first-moment bits',math.log2(lower.numerator)-math.log2(lower.denominator),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
