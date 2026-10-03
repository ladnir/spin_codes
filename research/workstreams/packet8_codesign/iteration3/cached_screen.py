"""Batched H3 proposal search; independently replay each winner in log space."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, log
from pathlib import Path
from time import monotonic

import numpy as np

import capped_gate as cg
import conditioned_diagnostics as cd


def route_threshold(q):
    """Use an exact witness, including deterministic full-slot endpoints."""
    minimum = 32*((q//8)*3+min(q%8,3))
    maximum = 32*min(q,192)
    if minimum == maximum:
        return cg.cap_counts.rational_witness(q,minimum,'1',cap=3,target_bits=60)
    return cg.search_threshold(q,3)


def parse_point(text):
    q,values = text.split(':')
    q,tilts = int(q),tuple(map(Fraction,values.split(',')))
    if not 1<=q<=512 or not tilts or min(tilts)<=0 or len(set(tilts))!=len(tilts):
        raise ValueError('q in1..512 and distinct positive weight tilts required')
    return q,tilts


def source_path(cache,tilt):
    return cache/f'theta_{tilt.numerator}_{tilt.denominator}.json'


def evaluate(local_families,q,*,regions=32):
    """Search aid: use batched scaling, falling back to all-log on failure."""
    try:
        moments,_ = cd.batch_moments(np.asarray(local_families),q,regions=regions)
        if not np.isfinite(moments).all():
            raise FloatingPointError('nonfinite batched moments')
        return moments,'batched scaled float; subnormal underflow allowed for proposals'
    except (FloatingPointError,ArithmeticError) as error:
        values = [cd.logarithmic_moment(local,q,regions=regions) for local in local_families]
        return np.asarray(values),'all-log fallback: '+str(error)


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    points = tuple(map(parse_point,args.points))
    if len({q for q,_ in points})!=len(points):
        raise ValueError('distinct occupancy points required')
    nus = tuple(map(Fraction,args.nus))
    if not nus or min(nus)<0 or len(set(nus))!=len(nus):
        raise ValueError('distinct nonnegative route tilts required')
    pins = cg.gate.source_pins()
    for path in (Path(__file__).resolve(),cg.HERE/'test_cached_screen.py',
                 Path(cg.__file__).resolve(),cg.HERE/'test_capped_gate.py',
                 Path(cg.cap_counts.__file__).resolve(),Path(cg.variant_gate.__file__).resolve(),
                 Path(cd.__file__).resolve(),cg.HERE/'test_conditioned_diagnostics.py'):
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    caches = {}
    for tilt in sorted({tilt for _,tilts in points for tilt in tilts}):
        path = source_path(args.cache,tilt).resolve()
        local,saved = cg.load_local(path)
        if Fraction(saved['tilt'])!=tilt:
            raise ArithmeticError('cache name disagrees with weight tilt')
        caches[tilt]=(local,saved)
        pins.update(saved['source_sha256'])
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    if not cg.variant_gate.checked_pins(pins):
        raise ArithmeticError('initial source validation failed')
    beta = cg.gate.prior.rs_uniform_envelope.UniformInputEnvelope(16,8,8,2).beta
    beta_log = log(beta.numerator)-log(beta.denominator)
    result = dict(schema='packet8-iteration3-cache24-H3-proposal-1',
        map_record=next(iter(caches.values()))[1]['map_record'],
        geometry=dict(K=65536,N=131072,groups=512,regions=32,slots_per_region=512,
            physical_steps_per_region=64,state_bits=24,cutoff=13107,
            zero_initial_state=True,continuous_state=True,final_flush=False),
        source_sha256=pins,source_pins_verified_at_finish=False,
        proposal_only=True,whole_code_certificate=False,bad_route_bound_exact=True,
        good_message_bound_outward=False,all_occupancies_checked=False,
        search_backend_may_underflow=True,every_selected_winner_replayed_in_log_space=True,
        no_new_large_state_census=True,points=[])
    started=monotonic()
    for q,tilts in points:
        witness=route_threshold(q)
        bad=cg.cap_counts.rational_bad_route_bound(q,witness['h'],witness['chernoff_x'],cap=3)
        logbad=log(bad.numerator)-log(bad.denominator) if bad else -float('inf')
        choices=[(tilt,nu) for tilt in tilts for nu in nus]
        families=[cg.marked(caches[tilt][0],[3],[float(nu)]) for tilt,nu in choices]
        moments,backend=evaluate(families,q)
        fixed=log(comb(512,q))+q*beta_log
        trials=[]
        for (tilt,nu),moment in zip(choices,moments):
            exponent=fixed+float(tilt)*13107-float(nu)*witness['h']+float(moment)
            trials.append(dict(tilt=str(tilt),route_tilt=str(nu),
                marked_log_moment=float(moment),good_message_margin_bits=-exponent/log(2),
                combined_margin_bits=-float(np.logaddexp(exponent,logbad))/log(2)))
        winner=max(range(len(trials)),key=lambda i:trials[i]['good_message_margin_bits'])
        checked=cd.logarithmic_moment(families[winner],q)
        discrepancy=float(moments[winner])-checked
        if abs(discrepancy)>1e-7:
            raise ArithmeticError('selected winner disagrees with all-log replay')
        best=dict(trials[winner],all_log_replay=checked,
                  scaled_minus_log_discrepancy=discrepancy)
        result['points'].append(dict(q=q,route_threshold=witness,
            search_backend=backend,trials=trials,best_choice=best))
        if not cg.variant_gate.checked_pins(pins):
            raise ArithmeticError('source changed during cached screen')
        result['elapsed_seconds']=monotonic()-started
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print(f'q={q} h3={witness["h"]} theta={best["tilt"]} nu={best["route_tilt"]}: '
              f'good={best["good_message_margin_bits"]:.6f}, '
              f'combined={best["combined_margin_bits"]:.6f}, logcheck={discrepancy:.3g}; '
              f'elapsed={result["elapsed_seconds"]:.2f}s',flush=True)
    result['source_pins_verified_at_finish']=True
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,default=cg.HERE/'cache24_v1')
    parser.add_argument('--points',nargs='+',required=True)
    parser.add_argument('--nus',nargs='+',default=['3','3.5','4'])
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
