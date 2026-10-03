"""Bounded floating gate with exact rational capped-route witnesses.

The saved local operators describe the actual 24-bit-state construction.
Only the route-conditioned comparison changes. This script does not certify
the good-message expectation or sum over every outer occupancy.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import comb, exp, floor, isfinite, log
from pathlib import Path
import sys
from time import monotonic

import numpy as np
from scipy.optimize import minimize_scalar

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'route_conditioning'))
import cap_counts
import gate
import variant_gate


def search_threshold(q, cap, *, target_bits=60, denominator_limit=1000):
    """Search a threshold, then verify the rounded rational witness exactly."""
    def continuous_threshold(eta):
        moment, _ = cap_counts.moment_and_mean(q, eta, cap=cap)
        return 1+(-target_bits*log(2)-log(comb(512, q))-32*moment)/eta

    fit = minimize_scalar(lambda eta: -continuous_threshold(eta),
                          bounds=(.1, 5.), method='bounded')
    if not fit.success:
        raise ArithmeticError('route threshold search failed')
    x = Fraction(exp(-fit.x)).limit_denominator(denominator_limit)
    h = floor(continuous_threshold(-log(float(x))))
    witness = cap_counts.rational_witness(q, h, x, cap=cap, target_bits=target_bits)
    while not witness['exact_integer_check_passed']:
        h -= 1
        witness = cap_counts.rational_witness(q, h, x, cap=cap, target_bits=target_bits)
    # One exact check records that the chosen rational tilt does not permit h+1.
    next_witness = cap_counts.rational_witness(q, h+1, x, cap=cap, target_bits=target_bits)
    if next_witness['exact_integer_check_passed']:
        raise ArithmeticError('threshold rounding did not reach the rational maximum')
    witness['floating_search_eta'] = float(fit.x)
    witness['next_integer_threshold_fails_at_same_rational_tilt'] = True
    return witness


def marked(local, caps, nus, *, packet_bits=8):
    """Thin potential slots first, then mark all requested route statistics."""
    if len(caps) != len(nus) or not caps or len(set(caps)) != len(caps):
        raise ValueError('distinct caps and one tilt per cap required')
    if any(type(c) is not int or not 1 <= c < len(local) for c in caps):
        raise ValueError('integer cap in1..physical_packet_count required')
    if any(not isfinite(nu) or nu < 0 for nu in nus):
        raise ValueError('finite nonnegative route tilts required')
    result = gate.potential_operators(local, 0., packet_bits=packet_bits)
    js = np.arange(len(local))
    exponent = sum(float(nu)*np.minimum(js, cap) for cap, nu in zip(caps, nus))
    result *= np.exp(exponent)[:, None, None]
    if not np.isfinite(result).all():
        raise FloatingPointError('marked operator overflow')
    return result


def exact_bad_union(q, witnesses):
    """Union over route conditions; no independence between conditions assumed."""
    bound = sum((cap_counts.rational_bad_route_bound(q, item['h'], item['chernoff_x'],
                 cap=item['cap']) for item in witnesses), Fraction(0))
    if not 0 < bound <= 1:
        raise ArithmeticError('positive route probability upper bound required')
    bits = bound.denominator.bit_length()-bound.numerator.bit_length()
    if bound.numerator << bits > bound.denominator:
        bits -= 1
    identity = hex(bound.numerator)+'/'+hex(bound.denominator)
    return bound, dict(combination='exact rational sum; no independence assumption',
        certified_margin_bits_floor=bits,
        rational_sum_sha256=hashlib.sha256(identity.encode()).hexdigest(),
        estimated_margin_bits=(log(bound.denominator)-log(bound.numerator))/log(2),
        whole_code_certificate=False)


def load_local(path):
    receipt = json.loads(path.read_text(encoding='utf-8'))
    if (not receipt['source_pins_verified_at_finish']
            or not variant_gate.checked_pins(receipt['source_sha256'])):
        raise ValueError('saved local receipt has stale source pins')
    local = np.asarray(receipt['local_operator_matrices'], dtype=float)
    if receipt['map_record']['state_bits'] != 24 or local.shape != (9, 10, 10):
        raise ValueError('actual24-state eight-packet local family required')
    return local, receipt


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    local, saved = load_local(args.receipt)
    caps = [args.cap]
    choices = [tuple(map(Fraction, args.nus))]
    if args.second_cap is not None:
        caps.append(args.second_cap)
        choices.append(tuple(map(Fraction, args.second_nus)))
    if any(not values or min(values) < 0 for values in choices):
        raise ValueError('nonempty nonnegative route tilt grids required')
    pins = gate.source_pins()
    pins.update(saved['source_sha256'])
    for path in (Path(__file__).resolve(), HERE/'test_capped_gate.py',
                 Path(cap_counts.__file__).resolve(), Path(variant_gate.__file__).resolve(),
                 args.receipt.resolve()):
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    started = monotonic()
    witnesses = [search_threshold(args.q, cap) for cap in caps]
    bad, bad_record = exact_bad_union(args.q, witnesses)
    beta = gate.prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2).beta
    theta = float(saved['tilt'])
    fixed = (log(comb(512, args.q))+args.q*(log(beta.numerator)-log(beta.denominator))
             +theta*13107)
    result = dict(schema='packet8-iteration3-capped-route-proposal-1',
        q=args.q, tilt=theta, caps=caps, route_thresholds=witnesses,
        bad_route_union=bad_record, map_record=saved['map_record'],
        saved_local_receipt=str(args.receipt.resolve()), no_new_large_state_census=True,
        geometry=dict(K=65536,N=131072,groups=512,regions=32,slots_per_region=512,
            physical_steps_per_region=64,state_bits=24,cutoff=13107,
            zero_initial_state=True,continuous_state=True,final_flush=False),
        proposal_only=True, whole_code_certificate=False, bad_route_bound_exact=True,
        good_message_bound_outward=False, all_occupancies_checked=False,
        source_sha256=pins, source_pins_verified_at_finish=False,
        trials=[], best_margin_bits=-float('inf'))
    print('exact route thresholds: '+str([(w['cap'], w['h'], w['chernoff_x'])
                                        for w in witnesses]), flush=True)
    for nus in product(*choices):
        potential = marked(local, caps, tuple(map(float, nus)))
        regional, backend = gate.prior.placement(potential, args.q, epochs=64, windows=8)
        moment = gate.prior.logarithmic.log_power_matrix(regional[args.q], 32)
        exponent = fixed+moment-sum(float(nu)*w['h'] for nu, w in zip(nus, witnesses))
        good = -exponent/log(2)
        combined = -float(np.logaddexp(exponent, log(bad.numerator)-log(bad.denominator)))/log(2)
        trial = dict(route_tilts=list(map(str,nus)), backend=backend,
            marked_log_moment=moment, good_message_margin_bits=good,
            combined_margin_bits=combined)
        result['trials'].append(trial)
        if combined > result['best_margin_bits']:
            result['best_margin_bits'], result['best_choice'] = combined, trial
        if not variant_gate.checked_pins(pins):
            raise ArithmeticError('source changed during proposal')
        result['elapsed_seconds'] = monotonic()-started
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print(f'caps={caps} nus={nus}: good={good:.6f}, combined={combined:.6f}; '
              f'elapsed={result["elapsed_seconds"]:.2f}s',flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt',type=Path,default=HERE.parent/'larger_state/trajectory_v1.json')
    parser.add_argument('--q',type=int,default=119)
    parser.add_argument('--cap',type=int,default=3)
    parser.add_argument('--nus',nargs='+',default=['0','1','2','3','4','5'])
    parser.add_argument('--second-cap',type=int)
    parser.add_argument('--second-nus',nargs='+',default=['0'])
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args())
