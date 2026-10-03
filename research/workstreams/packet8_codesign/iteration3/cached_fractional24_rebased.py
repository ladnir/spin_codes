"""Replay the frozen actual24 G4 grid in the potential-birth comparison basis."""
import argparse
from fractions import Fraction
import hashlib
import json
from math import comb,log
from pathlib import Path
from time import monotonic

import numpy as np

import cached_fractional24 as base
import potential_birth_gate as birth


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    control=json.loads(args.control.read_text(encoding='utf-8'))
    if (not control['source_pins_verified_at_finish']
            or not base.cs.cg.variant_gate.checked_pins(control['source_sha256'])):
        raise ValueError('control receipt has stale source pins')
    pins=dict(control['source_sha256'])
    pins.update(birth.source_pins())
    for path in (Path(__file__).resolve(),args.control.resolve(),Path(base.__file__).resolve()):
        pins[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    started=monotonic()
    parameter_pairs={(Fraction(t['tilt']),Fraction(t['alpha']))
        for point in control['points'] for t in point['trials']}
    families={}
    algebra={}
    for theta in sorted({theta for theta,_ in parameter_pairs}):
        path=base.find_cache(theta,args.caches)
        local,saved=base.cs.cg.load_local(path)
        if saved['map_record']!=control['map_record'] or Fraction(saved['tilt'])!=theta:
            raise ArithmeticError('cached construction or tilt differs from control')
        old=base.cs.cg.gate.potential_operators(local,0.)
        new,change,checks=birth.rebase_potential(old)
        tuples=base.fine.tuple_products(new,4)
        for candidate_theta,alpha in parameter_pairs:
            if candidate_theta==theta:
                families[(theta,alpha)]=base.fine.fine_operators(tuples,alpha)
        algebra[str(theta)]=dict(checks,change_of_basis=change.tolist())
        pins.update(saved['source_sha256'])
        pins[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    beta=base.cs.cg.gate.prior.rs_uniform_envelope.UniformInputEnvelope(16,8,8,2).beta
    beta_log=log(beta.numerator)-log(beta.denominator)
    result=dict(schema='packet8-actual24-cached-G4-rebased-proposal-1',
        map_record=control['map_record'],geometry=control['geometry'],
        control_receipt=str(args.control.resolve()),construction_changed=False,
        birth_basis='one normalized family per potential occupancy',
        moment_identity='Tnew_j Q=Q Told_j; e0 Q=e0; Q1=1',
        proposal_only=True,whole_code_certificate=False,has_outward_endpoints=False,
        no_new_large_state_census=True,source_sha256=pins,source_pins_verified_at_finish=False,
        message_factor=control['message_factor'],subset_union_is_not_powered=True,
        power_order=control['power_order'],algebra_checks=algebra,points=[])
    for point in control['points']:
        q=point['q']
        choices=[(Fraction(t['tilt']),Fraction(t['alpha'])) for t in point['trials']]
        local=np.asarray([families[key] for key in choices])
        moments,_=base.cs.cd.batch_moments(local,q,epochs=16,windows=32,regions=32)
        if not np.isfinite(moments).all():
            raise FloatingPointError('nonfinite batched fractional moment')
        trials=[]
        for old,(theta,alpha),moment in zip(point['trials'],choices,moments):
            exponent=log(comb(512,q))+float(alpha)*(q*beta_log+float(theta)*13107)+float(moment)
            margin=-exponent/log(2)
            if alpha==1 and abs(margin-old['margin_bits'])>5e-7:
                raise ArithmeticError('alpha1 global comparison changed')
            trials.append(dict(tilt=str(theta),alpha=str(alpha),margin_bits=margin,
                old_margin_bits=old['margin_bits'],gain_bits=margin-old['margin_bits'],
                fractional_log_moment=float(moment)))
        winner=max(range(len(trials)),key=lambda i:trials[i]['margin_bits'])
        replay=base.cs.cd.logarithmic_moment(local[winner],q,epochs=16,windows=32)
        error=float(moments[winner])-replay
        if abs(error)>1e-7:
            raise ArithmeticError('selected winner differs from all-log replay')
        best=dict(trials[winner],all_log_replay=replay,scaled_minus_log_discrepancy=error)
        result['points'].append(dict(q=q,trials=trials,best_choice=best,
            old_best_choice=point['best_choice'],alpha1_global_regression_passed=True))
        if not base.cs.cg.variant_gate.checked_pins(pins):
            raise ArithmeticError('source changed during rebased actual24 gate')
        result['elapsed_seconds']=monotonic()-started
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print(f'actual24 rebasedG4 q={q},theta={best["tilt"]},alpha={best["alpha"]}: '
            f'margin={best["margin_bits"]:.6f},gain={best["gain_bits"]:.6f},logcheck={error:.3g}',flush=True)
    result['source_pins_verified_at_finish']=True
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--control',type=Path,default=base.cs.cg.HERE/'actual24_fractional_g4_v1.json')
    parser.add_argument('--caches',type=Path,nargs='+',default=[base.cs.cg.HERE/'cache24_v1',base.cs.cg.HERE/'cache24_v2'])
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
