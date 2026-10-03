"""Last bounded H2 gate using already-authenticated actual24-state operators."""
import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, log
from pathlib import Path
from time import monotonic

import numpy as np

import cap_gate
import variant_gate


def run(args):
    if args.output.exists():
        raise ValueError('fresh output required')
    receipt = json.loads(args.receipt.read_text(encoding='utf-8'))
    if not receipt['source_pins_verified_at_finish'] or not variant_gate.checked_pins(receipt['source_sha256']):
        raise ValueError('stale saved actual24-state receipt')
    local = np.asarray(receipt['local_operator_matrices'], dtype=float)
    if receipt['map_record']['state_bits'] != 24 or local.shape != (9, 10, 10):
        raise ValueError('actual24-state eight-packet local operators required')
    tilt = float(receipt['tilt'])
    nus = tuple(Fraction(value) for value in args.nus)
    if not nus or min(nus) < 0:
        raise ValueError('nonnegative route tilts required')
    pins = cap_gate.sources()
    pins.update(receipt['source_sha256'])
    for path in (Path(__file__).resolve(), Path(variant_gate.__file__).resolve(), args.receipt.resolve()):
        pins[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    if not variant_gate.checked_pins(pins):
        raise ArithmeticError('initial source validation failed')
    witness = cap_gate.cap_counts.rational_witness(119, 2559, '9/35', cap=2, target_bits=60)
    if not witness['exact_integer_check_passed']:
        raise ArithmeticError('exact route threshold failed')
    beta = cap_gate.gate.prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2).beta
    beta_log = log(beta.numerator)-log(beta.denominator)
    result = dict(schema='actual24-state-cap2-conditioned-proposal-1', q=119, h=2559, cap=2,
        tilt=tilt, map_record=receipt['map_record'], source_sha256=pins,
        source_pins_verified_at_finish=False, saved_local_receipt=str(args.receipt.resolve()),
        no_new_large_state_census=True, route_tilts=list(map(str, nus)),
        bad_route_witness=witness, bad_route_bound_exact=True,
        good_message_bound_outward=False, proposal_only=True, whole_code_certificate=False,
        geometry=dict(K=65536, N=131072, groups=512, regions=32, slots_per_region=512,
            physical_steps_per_region=64, state_bits=24, cutoff=13107,
            zero_initial_state=True, continuous_state=True, final_flush=False),
        trials=[], best_margin_bits=-float('inf'))
    started = monotonic()
    for nu in nus:
        marked = cap_gate.marked(local, float(nu))
        regional, backend = cap_gate.gate.prior.placement(marked, 119, epochs=64, windows=8)
        moment = cap_gate.gate.prior.logarithmic.log_power_matrix(regional[119], 32)
        exponent = log(comb(512, 119))+119*beta_log+tilt*13107-float(nu)*2559+moment
        good_margin = -exponent/log(2)
        margin = -float(np.logaddexp(exponent, -witness['certified_margin_bits_floor']*log(2)))/log(2)
        trial = dict(route_tilt=str(nu), backend=backend, marked_log_moment=moment,
            good_message_margin_bits=good_margin, combined_margin_bits=margin)
        result['trials'].append(trial)
        if margin > result['best_margin_bits']:
            result['best_margin_bits'], result['best_choice'] = margin, trial
        if not variant_gate.checked_pins(pins):
            raise ArithmeticError('source changed during saved-operator gate')
        result['elapsed_seconds'] = monotonic()-started
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        print(f'actual24 H2 theta={tilt} nu={nu}: good={good_margin:.6f}, '
              f'best={result["best_margin_bits"]:.6f}; elapsed={result["elapsed_seconds"]:.2f}s', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt', type=Path, default=cap_gate.gate.HERE.parent/'larger_state/trajectory_v1.json')
    parser.add_argument('--nus', nargs='+', default=['3', '3.5', '4'])
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
