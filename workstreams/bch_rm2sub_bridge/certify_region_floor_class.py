"""Direct row counting when constant-one rows force a region-weight floor.

For h>=floor, every actual region weight is at least floor, regardless of
ordinary row weights. A uniform entrywise region envelope therefore applies
before averaging any outer rows. Count each ordinary row by 2^128-2.
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


def run(limit,floor,tilt,tag,verify=False):
    assert 0<=limit<=8192 and 1<=floor<=8192 and tag.isidentifier()
    output=base.HERE/'generated'/f'region_floor_{tag}_outward.json'
    old=base.read(output) if verify else None
    if not verify:assert not output.exists()
    if old:
        assert old['scope']==dict(ordinary_rows_maximum=limit,all_one_rows_minimum=floor)
        assert old['tilt_tenth']==tilt
        for name,digest in old['local_sha256'].items():assert base.sha(base.HERE/name)==digest
    t,s,spectrum=base.load_map(tight.NAME);ctx.prec=512 if verify else 256
    lam=(arb(tilt)/10).exp()
    region=poly.regions(t,s,spectrum,tight.kernel_spectrum(),(-lam).exp(),8192)
    maximum=tuple(max(row[k].upper() for row in region[floor:]) for k in range(9))
    value=maximum
    for _ in range(8):value=q1.positive_mul(value,value)
    ordinary=(1<<128)-2
    # Enlarge the number of allowed constant assignments only after applying
    # the uniform moment bound to the eligible assignments with h>=floor.
    count=sum(math.comb(8192,d)*(1<<(8192-d))*ordinary**d for d in range(min(limit,8192-floor)+1))
    from audit_bch_q1_full_arb import rational
    bound=count*rational((sum(value[:3],arb(0))*(209716*lam).exp()).upper())
    assert bound>0
    if old:
        assert count==int(old['row_message_count_upper'])
        assert bound<=base.decode(old['failure_contribution_upper'])
        print('512-bit direct region-floor class replay passed',flush=True);return
    margin=math.log2(bound.denominator)-math.log2(bound.numerator)
    sources=[Path(__file__),Path(poly.__file__),Path(tight.__file__),Path(q1.__file__),base.HERE/'general_occupancy.py',
             base.HERE/'inputs/t128_s15_selection.json',base.HERE/'inputs/t128_s15_a_spectrum.json',base.HERE/'inputs/t128_s15_b_kernel_spectrum.json']
    base.write_new(output,dict(status='OUTWARD_DIRECT_REGION_FLOOR_CLASS',configuration=tight.NAME,
        parameters=dict(message_bits=1<<20,output_bits=1<<21,outer_rows=8192,outer_length=256,step_bits=t,state_bits=s,distance_cutoff=209716),
        scope=dict(ordinary_rows_maximum=limit,all_one_rows_minimum=floor),tilt_tenth=tilt,
        row_message_count_upper=str(count),failure_contribution_upper=base.encode(bound),margin_bits_diagnostic=margin,
        below_2_to_minus_60=bound<F(1,1<<60),all_occupations_certified=False,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in sources}))
    print('Direct region-floor class',limit,floor,'margin',margin,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--ordinary-limit',type=int,required=True);p.add_argument('--floor',type=int,required=True)
    p.add_argument('--tilt',type=int,required=True);p.add_argument('--tag',required=True);p.add_argument('--verify',action='store_true')
    a=p.parse_args();run(a.ordinary_limit,a.floor,a.tilt,a.tag,a.verify)
