"""Fine G4 fractional gate using authenticated actual24 local operators."""
import argparse
from fractions import Fraction
import hashlib
import json
from math import comb,log
from pathlib import Path
from time import monotonic

import numpy as np

import cached_screen as cs
import fine_grouped_gate as fine


def find_cache(tilt,roots):
    paths=[cs.source_path(root,tilt).resolve() for root in roots]
    found=[path for path in paths if path.exists()]
    if len(found)!=1:
        raise ValueError(f'exactly one cache required for tilt{tilt}, found{len(found)}')
    return found[0]


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    points=tuple(map(cs.parse_point,args.points))
    alphas=tuple(map(Fraction,args.alphas))
    if not alphas or min(alphas)<=0 or max(alphas)>1:
        raise ValueError('fractional powers in(0,1] required')
    started=monotonic()
    pins=fine.source_pins()
    for path in (Path(__file__).resolve(),Path(cs.__file__).resolve(),
                 Path(cs.cg.__file__).resolve(),Path(cs.cg.variant_gate.__file__).resolve(),
                 Path(cs.cd.__file__).resolve(),cs.cg.HERE/'test_conditioned_diagnostics.py'):
        pins[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    caches={}
    for tilt in sorted({tilt for _,tilts in points for tilt in tilts}):
        path=find_cache(tilt,args.caches)
        local,saved=cs.cg.load_local(path)
        if Fraction(saved['tilt'])!=tilt:
            raise ArithmeticError('cache weight tilt mismatch')
        potential=cs.cg.gate.potential_operators(local,0.)
        tuples=fine.tuple_products(potential,4)
        families={alpha:fine.fine_operators(tuples,alpha) for alpha in alphas}
        if Fraction(1) not in families:
            families[Fraction(1)]=fine.fine_operators(tuples,1.)
        coarse,_=fine.coarse.grouped_operators(potential,4)
        positive=coarse>0
        if np.any(families[1][~positive]!=0):
            raise ArithmeticError('alpha1 macro has an unsupported positive entry')
        relative=float(np.max(np.abs(families[1][positive]/coarse[positive]-1)))
        if relative>2e-11:
            raise ArithmeticError('alpha1 macro does not match ordinary ordered grouping')
        caches[tilt]=dict(families=families,record=saved['map_record'],
                         alpha1_macro_relative_error=relative)
        pins.update(saved['source_sha256'])
        pins[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    if not cs.cg.variant_gate.checked_pins(pins):
        raise ArithmeticError('initial source validation failed')
    beta=cs.cg.gate.prior.rs_uniform_envelope.UniformInputEnvelope(16,8,8,2).beta
    beta_log=log(beta.numerator)-log(beta.denominator)
    result=dict(schema='packet8-actual24-cached-fine-G4-proposal-1',
        map_record=next(iter(caches.values()))['record'],
        geometry=dict(K=65536,N=131072,groups=512,regions=32,slots_per_region=512,
            physical_steps_per_region=64,state_bits=24,cutoff=13107,
            macro_steps_per_region=16,packets_per_macro=32,group_steps=4,
            zero_initial_state=True,continuous_state=True,final_flush=False),
        message_factor='C(512,q)*beta^(q*alpha)',subset_union_is_not_powered=True,
        power_order='potential thinning; ordered four-step product; entrywise alpha; tuple averaging',
        route_conditioning=False,proposal_only=True,whole_code_certificate=False,
        has_outward_endpoints=False,no_new_large_state_census=True,
        source_sha256=pins,source_pins_verified_at_finish=False,
        alpha1_macro_relative_errors={str(tilt):value['alpha1_macro_relative_error']
            for tilt,value in caches.items()},points=[])
    for q,tilts in points:
        choices=[(tilt,alpha) for tilt in tilts for alpha in alphas]
        families=np.asarray([caches[tilt]['families'][alpha] for tilt,alpha in choices])
        try:
            moments,_=cs.cd.batch_moments(families,q,epochs=16,windows=32,regions=32)
            if not np.isfinite(moments).all():
                raise FloatingPointError('nonfinite batched moment')
            backend='batched scaled float; subnormal underflow allowed for proposals'
        except (ArithmeticError,FloatingPointError) as error:
            moments=np.array([cs.cd.logarithmic_moment(local,q,epochs=16,windows=32)
                              for local in families])
            backend='all-log fallback: '+str(error)
        trials=[]
        for (tilt,alpha),moment in zip(choices,moments):
            exponent=log(comb(512,q))+float(alpha)*(q*beta_log+float(tilt)*13107)+float(moment)
            trials.append(dict(tilt=str(tilt),alpha=str(alpha),fractional_log_moment=float(moment),
                margin_bits=-exponent/log(2)))
        winner=max(range(len(trials)),key=lambda i:trials[i]['margin_bits'])
        replay=cs.cd.logarithmic_moment(families[winner],q,epochs=16,windows=32)
        error=float(moments[winner])-replay
        if abs(error)>1e-7:
            raise ArithmeticError('selected fractional winner differs from all-log replay')
        best=dict(trials[winner],all_log_replay=replay,scaled_minus_log_discrepancy=error)
        result['points'].append(dict(q=q,trials=trials,best_choice=best,search_backend=backend))
        if not cs.cg.variant_gate.checked_pins(pins):
            raise ArithmeticError('source changed during cached fractional screen')
        result['elapsed_seconds']=monotonic()-started
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print(f'actual24G4 q={q},theta={best["tilt"]},alpha={best["alpha"]}: '
              f'margin={best["margin_bits"]:.6f}, logcheck={error:.3g}; '
              f'elapsed={result["elapsed_seconds"]:.2f}s',flush=True)
    result['source_pins_verified_at_finish']=True
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--caches',type=Path,nargs='+',default=[cs.cg.HERE/'cache24_v1',cs.cg.HERE/'cache24_v2'])
    parser.add_argument('--points',nargs='+',default=['64:.2,.3,.4','4:.01,.02,.04','16:.02,.04,.06,.08'])
    parser.add_argument('--alphas',nargs='+',default=['.2','.35','.5','.7','1'])
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
