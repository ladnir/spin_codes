"""The single approved actual24 rebased G8 gate, with reusable macro output.

The 64MiB helper setting controls its chunk heuristic. It is not a measured
peak-RSS limit: overlapping temporaries and Python/NumPy allocations remain.
All arithmetic here is a floating proposal, never an outward certificate.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from math import comb,isfinite,log
from pathlib import Path
from time import monotonic

import numpy as np

import cached_fractional24 as base
import long_grouped_gate as long
import potential_birth_gate as birth


def intervals(per_q,threshold):
    result=[]
    for q,value in per_q.items():
        q=int(q)
        if value['margin_bits']>=threshold:
            if result and result[-1][1]+1==q:
                result[-1][1]=q
            else:
                result.append([q,q])
    return result


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    theta,alpha=Fraction(3,50),Fraction(2,5)
    path=base.find_cache(theta,args.caches)
    physical,saved=base.cs.cg.load_local(path)
    if Fraction(saved['tilt'])!=theta:
        raise ArithmeticError('cache weight tilt mismatch')
    pins=long.source_pins()
    pins.update(birth.source_pins())
    pins.update(saved['source_sha256'])
    for source in (Path(__file__).resolve(),Path(base.__file__).resolve(),
                   Path(base.cs.cg.__file__).resolve(),Path(base.cs.cg.variant_gate.__file__).resolve(),path):
        pins[str(source)]=hashlib.sha256(source.read_bytes()).hexdigest()
    if not base.cs.cg.variant_gate.checked_pins(pins):
        raise ArithmeticError('initial source validation failed')
    started=monotonic()
    old=long.marked.potential_operators(physical,0.)
    potential,change,algebra=birth.rebase_potential(old)
    positive=potential[potential>0]
    smallest,largest=float(positive.min()),float(positive.max())
    # This is a conservative range check on the supplied floating matrices,
    # not an exact local endpoint certificate. The generous buffer absorbs
    # insignificant rounding in this diagnostic logarithm.
    product_lower=8*min(0.,log(smallest))
    weighted_lower=float(alpha)*product_lower-log(comb(64,32))
    if min(product_lower,weighted_lower)<log(np.finfo(float).tiny)+100:
        raise ArithmeticError('supplied matrices fail the conservative normal-range gate')
    def progress(done,total):
        if not base.cs.cg.variant_gate.checked_pins(pins):
            raise ArithmeticError('source changed during the G8 macro census')
        print(f'actual24 G8 chunks {done}/{total}; elapsed={monotonic()-started:.2f}s',flush=True)
    print('Building only actual24 rebased G8 theta=.06 alpha=.4',flush=True)
    with np.errstate(under='raise',over='raise',invalid='raise',divide='raise'):
        local,macro=long.long_operators(potential,8,float(alpha),memory_mib=64,progress=progress)
    macro.update(memory_setting_mib=64,memory_setting_is_peak_rss_guarantee=False,
        memory_caveat='nominal chunk heuristic; overlapping temporaries and runtime allocations are additional',
        strict_numpy_underflow_check=True,
        supplied_float_min_positive=smallest,supplied_float_max=largest,
        conservative_positive_product_log_lower=product_lower,
        conservative_weighted_power_log_lower=weighted_lower)
    result=dict(schema='packet8-actual24-cached-rebased-G8-proposal-1',
        tilt=str(theta),alpha=str(alpha),group_steps=8,target_occupancy=16,
        map_record=saved['map_record'],
        geometry=dict(K=65536,N=131072,groups=512,regions=32,slots_per_region=512,
            physical_steps_per_region=64,macro_steps_per_region=8,packets_per_macro=64,
            state_bits=24,cutoff=13107,zero_initial_state=True,continuous_state=True,final_flush=False),
        construction_changed=False,proposal_only=True,whole_code_certificate=False,
        has_outward_endpoints=False,no_new_large_state_census=True,
        source_sha256=pins,source_pins_verified_at_finish=False,
        power_order='potential thinning; full eight-step ordered product; entrywise alpha; tuple averaging',
        message_factor='C(512,q)*beta^(q*alpha)',subset_union_is_not_powered=True,
        birth_basis='one normalized family per potential occupancy',algebra_checks=algebra,
        change_of_basis=change.tolist(),macro_metadata=macro,
        powered_macro_operators=local.tolist(),potential_operators=potential.tolist(),
        occupancy_range=[1,512],per_q={},positive_intervals=[],at_least_40_bit_intervals=[])
    def save():
        if not base.cs.cg.variant_gate.checked_pins(pins):
            raise ArithmeticError('source changed during cached G8 evaluation')
        result['elapsed_seconds']=monotonic()-started
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    save()
    # The frozen helper checked its simultaneous alpha-one tuple sum against
    # this ordinary macro family. Now independently check the whole regional
    # geometry, and then global state propagation, at every occupancy.
    alpha1,_=long.fine.coarse.grouped_operators(potential,8)
    grouped,_=long.prior.placement(alpha1,512,epochs=8,windows=64,force_log=True)
    original,_=long.prior.placement(potential,512,epochs=64,windows=8,force_log=True)
    if not np.array_equal(np.isfinite(grouped),np.isfinite(original)):
        raise ArithmeticError('alpha1 regional structural zeros changed')
    finite=np.isfinite(original)
    regional_error=float(np.max(np.abs(grouped[finite]-original[finite])))
    if regional_error>2e-8:
        raise ArithmeticError('alpha1 regional placement regression failed')
    global_error=0.
    for q in range(513):
        left=long.prior.logarithmic.log_power_matrix(grouped[q],32)
        right=long.prior.logarithmic.log_power_matrix(original[q],32)
        global_error=max(global_error,abs(left-right))
    if global_error>2e-7:
        raise ArithmeticError('alpha1 global propagation regression failed')
    result['alpha1_regression']=dict(occupancy_range=[0,512],
        max_regional_log_error=regional_error,max_global_log_error=global_error,
        all_log_backend=True)
    print(f'alpha1 controls: regional={regional_error:.3g}, global={global_error:.3g}',flush=True)
    save()
    regional,_=long.prior.placement(local,512,epochs=8,windows=64,force_log=True)
    beta=long.prior.rs_uniform_envelope.UniformInputEnvelope(16,8,8,2).beta
    beta_log=log(beta.numerator)-log(beta.denominator)
    for q in range(1,513):
        moment=long.prior.logarithmic.log_power_matrix(regional[q],32)
        exponent=log(comb(512,q))+float(alpha)*(q*beta_log+float(theta)*13107)+moment
        margin=-exponent/log(2)
        if not isfinite(margin):
            raise FloatingPointError('nonfinite all-occupancy proposal')
        result['per_q'][str(q)]=dict(margin_bits=margin,log_moment=moment)
    result['regional_backend']='all-logarithmic floating'
    result['positive_intervals']=intervals(result['per_q'],0.)
    result['at_least_40_bit_intervals']=intervals(result['per_q'],40.)
    result['source_pins_verified_at_finish']=True
    save()
    print(f'actual24 G8 q16={result["per_q"]["16"]["margin_bits"]:.9f}bits; '
          f'positive={result["positive_intervals"]};40bit={result["at_least_40_bit_intervals"]}; '
          f'elapsed={monotonic()-started:.2f}s',flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--caches',type=Path,nargs='+',default=[base.cs.cg.HERE/'cache24_v1',base.cs.cg.HERE/'cache24_v2',base.cs.cg.HERE/'cache24_v3'])
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
