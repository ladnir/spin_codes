"""Compare route caps at one occupancy using authenticated weight-tilt caches."""
import argparse
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import comb,log
from pathlib import Path
from time import monotonic

import numpy as np

import cached_screen as cs


def candidates(tilts,nus,caps):
    if (not tilts or min(tilts)<=0 or not nus or min(nus)<0
            or not caps or len(set(caps))!=len(caps)
            or any(type(cap) is not int or not 1<=cap<=8 for cap in caps)):
        raise ValueError('positive weight tilts, nonnegative route tilts, distinct caps required')
    return [(tilt,choice) for tilt in tilts for choice in product(nus,repeat=len(caps))]


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    tilts,nus=tuple(map(Fraction,args.tilts)),tuple(map(Fraction,args.nus))
    choices=candidates(tilts,nus,args.caps)
    pins=cs.cg.gate.source_pins()
    for path in (Path(__file__).resolve(),cs.cg.HERE/'test_cap_compare.py',
                 Path(cs.__file__).resolve(),cs.cg.HERE/'test_cached_screen.py',
                 Path(cs.cg.__file__).resolve(),cs.cg.HERE/'test_capped_gate.py',
                 Path(cs.cg.cap_counts.__file__).resolve(),Path(cs.cg.variant_gate.__file__).resolve(),
                 Path(cs.cd.__file__).resolve(),cs.cg.HERE/'test_conditioned_diagnostics.py'):
        pins[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    caches={}
    for tilt in tilts:
        path=cs.source_path(args.cache,tilt).resolve()
        local,saved=cs.cg.load_local(path)
        if Fraction(saved['tilt'])!=tilt:
            raise ArithmeticError('cache weight tilt mismatch')
        caches[tilt]=(local,saved)
        pins.update(saved['source_sha256'])
        pins[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
    witnesses=[cs.cg.search_threshold(args.q,cap) for cap in args.caps]
    bad,bad_record=cs.cg.exact_bad_union(args.q,witnesses)
    beta=cs.cg.gate.prior.rs_uniform_envelope.UniformInputEnvelope(16,8,8,2).beta
    fixed=log(comb(512,args.q))+args.q*(log(beta.numerator)-log(beta.denominator))
    started=monotonic()
    families=[cs.cg.marked(caches[tilt][0],args.caps,tuple(map(float,route_tilts)))
              for tilt,route_tilts in choices]
    moments,backend=cs.evaluate(families,args.q)
    trials=[]
    for (tilt,route_tilts),moment in zip(choices,moments):
        exponent=(fixed+float(tilt)*13107+float(moment)
            -sum(float(nu)*w['h'] for nu,w in zip(route_tilts,witnesses)))
        trials.append(dict(tilt=str(tilt),route_tilts=list(map(str,route_tilts)),
            marked_log_moment=float(moment),good_message_margin_bits=-exponent/log(2),
            combined_margin_bits=-float(np.logaddexp(exponent,log(bad.numerator)-log(bad.denominator)))/log(2)))
    winner=max(range(len(trials)),key=lambda i:trials[i]['good_message_margin_bits'])
    replay=cs.cd.logarithmic_moment(families[winner],args.q)
    error=float(moments[winner])-replay
    if abs(error)>1e-7:
        raise ArithmeticError('selected winner differs from all-log replay')
    if not cs.cg.variant_gate.checked_pins(pins):
        raise ArithmeticError('source changed during cap comparison')
    best=dict(trials[winner],all_log_replay=replay,scaled_minus_log_discrepancy=error)
    result=dict(schema='packet8-iteration3-cap-comparison-proposal-1',
        q=args.q,caps=args.caps,map_record=next(iter(caches.values()))[1]['map_record'],
        route_thresholds=witnesses,bad_route_union=bad_record,search_backend=backend,
        proposal_only=True,whole_code_certificate=False,bad_route_bound_exact=True,
        good_message_bound_outward=False,all_occupancies_checked=False,
        source_sha256=pins,source_pins_verified_at_finish=True,
        search_backend_may_underflow=True,selected_winner_replayed_in_log_space=True,
        no_new_large_state_census=True,trials=trials,best_choice=best,
        elapsed_seconds=monotonic()-started)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(f'q={args.q},caps={args.caps},theta={best["tilt"]},nus={best["route_tilts"]}: '
          f'good={best["good_message_margin_bits"]:.6f},combined={best["combined_margin_bits"]:.6f}, '
          f'logcheck={error:.3g}',flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=cs.cg.HERE/'cache24_v1')
    parser.add_argument('--q',type=int,default=64)
    parser.add_argument('--caps',nargs='+',type=int,default=[2])
    parser.add_argument('--tilts',nargs='+',default=['.2','.3','.4'])
    parser.add_argument('--nus',nargs='+',default=['0','1','2','3','4'])
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
