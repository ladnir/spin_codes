"""Route-conditioning gates for two already-defined alternative local maps.

The scaled16 mode rebuilds its actual local census. The saved24 mode loads
authenticated local matrices from the existing full24-state receipt; it
does not resample maps or rerun a large Walsh transform.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, log
from pathlib import Path
from time import monotonic

import numpy as np

import gate


def checked_pins(pins):
    return all(Path(path).is_file() and hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
               for path, digest in pins.items())


def run(args):
    thresholds = tuple(gate.threshold(value) for value in args.thresholds)
    nus = tuple(Fraction(value) for value in args.nus)
    tilts = tuple(Fraction(value) for value in args.tilts)
    if args.output.exists() or min(nus) < 0 or min(tilts) <= 0:
        raise ValueError('fresh output and valid positive tilts required')
    started = monotonic()
    pins = gate.source_pins()
    pins[str(Path(__file__).resolve())] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if args.model == 'saved24':
        receipt = json.loads(args.receipt.read_text(encoding='utf-8'))
        if not receipt['source_pins_verified_at_finish'] or not checked_pins(receipt['source_sha256']):
            raise ValueError('saved local operators must retain all authenticated sources')
        if len(tilts) != 1 or float(tilts[0]) != float(receipt['tilt']):
            raise ValueError('saved local matrices support only their original tilt')
        pins.update(receipt['source_sha256'])
        pins[str(args.receipt.resolve())] = hashlib.sha256(args.receipt.read_bytes()).hexdigest()
        record = receipt['map_record']
        saved = np.asarray(receipt['local_operator_matrices'], dtype=float)
        if record['state_bits'] != 24 or saved.shape != (9, 10, 10):
            raise ValueError('actual24-state, eight-packet local receipt required')
    else:
        data, record = gate.maps.prepare('byte_native_A_scaled')
        pins.update(gate.prior.sources('byte_native_A_scaled'))
    if not checked_pins(pins):
        raise ArithmeticError('invalid initial source pins')
    envelope = gate.prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
    beta_log = log(envelope.beta.numerator)-log(envelope.beta.denominator)
    bad = {q: gate.route_counts.bad_route_margin(q, h, float(eta), regions=32) for q, h, eta in thresholds}
    result = dict(schema='packet8-route-conditioned-variant-proposal-1', model=args.model,
        map_record=record, proposal_only=True, whole_code_certificate=False,
        has_outward_endpoints=False, source_sha256=pins, source_pins_verified_at_finish=False,
        saved_local_receipt=str(args.receipt.resolve()) if args.model == 'saved24' else None,
        geometry=dict(K=65536, N=131072, groups=512, regions=32, slots_per_region=512,
            physical_steps_per_region=64, physical_t=64, state_bits=record['state_bits'],
            cutoff=13107, continuous_state=True, zero_initial_state=True, final_flush=False),
        input_tilts=list(map(str, tilts)), route_tilts=list(map(str, nus)),
        thresholds=[dict(q=q, h=h, eta=str(eta), bad_route_margin_bits=bad[q]) for q, h, eta in thresholds],
        best={}, choices={}, trials=[])
    maximum = max(q for q, _, _ in thresholds)
    for tilt in tilts:
        if args.model == 'saved24':
            local = saved
        else:
            weighted, emission, _ = gate.prior.moments(data, np.exp(-float(tilt)))
            local = gate.prior.operators(weighted, emission)
        original, _ = gate.prior.placement(local, maximum, epochs=64, windows=8)
        reference = {q: gate.prior.logarithmic.log_power_matrix(gate.prior.uniform_mixture(original, q), 32)
                     for q, _, _ in thresholds}
        for nu in nus:
            marked = gate.potential_operators(local, float(nu))
            regional, backend = gate.prior.placement(marked, maximum, epochs=64, windows=8)
            trial = dict(tilt=str(tilt), route_tilt=str(nu), backend=backend, witnesses={})
            for q, h, eta in thresholds:
                moment = gate.prior.logarithmic.log_power_matrix(regional[q], 32)
                if nu == 0 and abs(moment-reference[q]) > 2e-8:
                    raise ArithmeticError('zero-route-tilt regression failed')
                exponent = log(comb(512, q))+q*beta_log+float(tilt)*13107-float(nu)*h+moment
                margin = -float(np.logaddexp(exponent, -bad[q]*log(2)))/log(2)
                trial['witnesses'][str(q)] = dict(combined_margin_bits=margin,
                    good_message_margin_bits=-exponent/log(2), bad_route_margin_bits=bad[q],
                    marked_log_moment=moment, old_log_moment=reference[q])
                if str(q) not in result['best'] or margin > result['best'][str(q)]:
                    result['best'][str(q)] = margin
                    result['choices'][str(q)] = dict(tilt=str(tilt), route_tilt=str(nu),
                        good_message_margin_bits=-exponent/log(2), bad_route_margin_bits=bad[q])
            result['trials'].append(trial)
            if not checked_pins(pins):
                raise ArithmeticError('source changed during proposal')
            result['elapsed_seconds'] = monotonic()-started
            args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
            print(f'{args.model} theta={tilt} nu={nu}: best={result["best"]}; '
                  f'elapsed={result["elapsed_seconds"]:.2f}s', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', choices=('scaled16', 'saved24'), required=True)
    parser.add_argument('--receipt', type=Path, default=gate.HERE.parent/'larger_state/trajectory_v1.json')
    parser.add_argument('--tilts', nargs='+', default=['.5'])
    parser.add_argument('--nus', nargs='+', default=['0', '1', '2', '3', '4', '5', '6'])
    parser.add_argument('--thresholds', nargs='+', default=[
        '64:1057:1.5993122190268654', '119:1483:1.8952931340202324', '128:1531:1.9687464380734174'])
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
