"""H3 gate for the unchanged16-bit-state maps, using fresh local moments."""
import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, log
from pathlib import Path
from time import monotonic

import numpy as np

import capped_gate as cg


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    tilts, nus = tuple(map(Fraction,args.tilts)), tuple(map(Fraction,args.nus))
    if not tilts or min(tilts)<=0 or not nus or min(nus)<0:
        raise ValueError('positive weight tilts and nonnegative route tilts required')
    pins = cg.gate.source_pins()
    for path in (Path(__file__).resolve(),Path(cg.__file__).resolve(),
                 cg.HERE/'test_capped_gate.py',Path(cg.cap_counts.__file__).resolve(),
                 Path(cg.variant_gate.__file__).resolve()):
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    started = monotonic()
    witness = cg.search_threshold(args.q,3)
    bad,bad_record = cg.exact_bad_union(args.q,[witness])
    data,record = cg.gate.maps.prepare('byte_native')
    beta = cg.gate.prior.rs_uniform_envelope.UniformInputEnvelope(16,8,8,2).beta
    fixed = log(comb(512,args.q))+args.q*(log(beta.numerator)-log(beta.denominator))
    result = dict(schema='packet8-iteration3-baseline16-H3-proposal-1',q=args.q,
        map_record=record,geometry=cg.gate.prior.geometry(data),
        unchanged_baseline_construction=True,route_threshold=witness,bad_route_union=bad_record,
        proposal_only=True,whole_code_certificate=False,bad_route_bound_exact=True,
        good_message_bound_outward=False,all_occupancies_checked=False,
        source_sha256=pins,source_pins_verified_at_finish=False,trials=[],
        best_margin_bits=-float('inf'))
    for tilt in tilts:
        weighted,emission,_ = cg.gate.prior.moments(data,np.exp(-float(tilt)))
        local = cg.gate.prior.operators(weighted,emission)
        for nu in nus:
            potential = cg.marked(local,[3],[float(nu)])
            regional,backend = cg.gate.prior.placement(potential,args.q,epochs=64,windows=8)
            moment = cg.gate.prior.logarithmic.log_power_matrix(regional[args.q],32)
            exponent = fixed+float(tilt)*13107-float(nu)*witness['h']+moment
            good = -exponent/log(2)
            combined = -float(np.logaddexp(exponent,log(bad.numerator)-log(bad.denominator)))/log(2)
            trial = dict(tilt=str(tilt),route_tilt=str(nu),backend=backend,
                marked_log_moment=moment,good_message_margin_bits=good,combined_margin_bits=combined)
            result['trials'].append(trial)
            if combined>result['best_margin_bits']:
                result['best_margin_bits'],result['best_choice'] = combined,trial
            if not cg.variant_gate.checked_pins(pins):
                raise ArithmeticError('source changed during baseline control')
            result['elapsed_seconds'] = monotonic()-started
            args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
            print(f'baseline16 q={args.q} theta={tilt} nu={nu}: good={good:.6f}; '
                  f'elapsed={result["elapsed_seconds"]:.2f}s',flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--q',type=int,default=119)
    parser.add_argument('--tilts',nargs='+',default=['.4','.5'])
    parser.add_argument('--nus',nargs='+',default=['3','4','5'])
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
