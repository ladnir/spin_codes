"""Bounded actual24 event-aligned proposal from a frozen local cache."""
import argparse
from fractions import Fraction
import hashlib
import json
from math import comb,isfinite,log
from pathlib import Path
from time import monotonic

import numpy as np

import cached_fractional24 as base
import event_aligned as event


def run(args):
    theta,alpha=Fraction(args.tilt),Fraction(args.alpha)
    if args.output.exists() or not 1<=args.q<=512 or theta<=0 or not 0<alpha<=1:
        raise ValueError('fresh output and valid occupancy, tilt, fractional power required')
    path=base.find_cache(theta,args.caches)
    local,saved=base.cs.cg.load_local(path)
    if Fraction(saved['tilt'])!=theta:
        raise ArithmeticError('cache weight tilt mismatch')
    pins=event.source_pins()
    pins.update(saved['source_sha256'])
    for source in (Path(__file__).resolve(),Path(base.__file__).resolve(),
                   Path(base.cs.cg.__file__).resolve(),Path(base.cs.cg.variant_gate.__file__).resolve(),path):
        pins[str(source)]=hashlib.sha256(source.read_bytes()).hexdigest()
    if not base.cs.cg.variant_gate.checked_pins(pins):
        raise ArithmeticError('initial source validation failed')
    started=monotonic()
    old=event.marked.potential_operators(local,0.)
    potential,change,algebra=event.birth.rebase_potential(old)
    def progress(stage,elapsed):
        if not base.cs.cg.variant_gate.checked_pins(pins):
            raise ArithmeticError('source changed during event-aligned proposal')
        print(f'actual24 q={args.q}: {stage}; {elapsed:.2f}s',flush=True)
    regional,metadata=event.event_regional(potential,args.q,alpha=float(alpha),progress=progress)
    moment=event.prior.logarithmic.log_power_matrix(regional[args.q],32)
    beta=event.prior.rs_uniform_envelope.UniformInputEnvelope(16,8,8,2).beta
    exponent=(log(comb(512,args.q))+float(alpha)*(
        args.q*(log(beta.numerator)-log(beta.denominator))+float(theta)*13107)+moment)
    margin=-exponent/log(2)
    if not isfinite(margin) or not base.cs.cg.variant_gate.checked_pins(pins):
        raise ArithmeticError('nonfinite result or changed source')
    result=dict(schema='packet8-actual24-cached-event-pair-proposal-1',
        q=args.q,tilt=str(theta),alpha=str(alpha),map_record=saved['map_record'],
        geometry=dict(K=65536,N=131072,groups=512,regions=32,slots_per_region=512,
            physical_steps_per_region=64,state_bits=24,cutoff=13107,
            zero_initial_state=True,continuous_state=True,final_flush=False),
        source_sha256=pins,source_pins_verified_at_finish=True,
        no_new_large_state_census=True,proposal_only=True,whole_code_certificate=False,
        has_outward_endpoints=False,comparison_kernel_not_physical_kernel=True,
        event='nonempty potential occupancy, not nonzero byte labels',
        message_factor='C(512,q)*(beta^q*z^(-cutoff))^alpha',
        subset_union_is_not_powered=True,regional_normalization='one C(512,q), outside alpha',
        margin_bits=margin,fractional_log_moment=moment,
        event_metadata=metadata,algebra_checks=algebra,change_of_basis=change.tolist(),
        potential_operators=potential.tolist(),
        regional_log_operator=[[float(x) if np.isfinite(x) else None for x in row]
                               for row in regional[args.q]],elapsed_seconds=monotonic()-started)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(f'actual24 event q={args.q},theta={theta},alpha={alpha}: {margin:.9f}bits',flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--caches',type=Path,nargs='+',default=[base.cs.cg.HERE/'cache24_v1',base.cs.cg.HERE/'cache24_v2',base.cs.cg.HERE/'cache24_v3'])
    parser.add_argument('--q',type=int,default=16)
    parser.add_argument('--tilt',default='.06')
    parser.add_argument('--alpha',default='.4')
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
