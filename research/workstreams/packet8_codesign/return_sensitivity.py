"""Counterfactual return-denominator sensitivity, NOT a larger-state code.

Keep the literal t64/s16 moments, normalized birth shapes, outer measure,
and ordered placement fixed. Only divide nonzero-to-zero entries by the
ratio (2**target_bits-1)/65535. No such transition law is asserted realizable.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, log
from pathlib import Path
from time import monotonic

import numpy as np

import maps
import screen


def change_denominator(local, target_bits):
    if type(target_bits) is not int or not 16 <= target_bits <= 64:
        raise ValueError('counterfactual denominator bit count16..64 required')
    result = local.copy()
    result[:, 1:, 0] *= 65535/((1 << target_bits)-1)
    return result


def sources():
    result = screen.sources('byte_native')
    path = Path(__file__).resolve()
    result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def run(args):
    if args.output.exists():
        raise FileExistsError('fresh output required')
    tilts = tuple(map(Fraction, args.tilts))
    if not tilts or min(tilts) <= 0 or len(set(tilts)) != len(tilts):
        raise ValueError('distinct positive tilts required')
    occupancies = tuple(args.q)
    if not occupancies or any(not 2 <= q <= 512 for q in occupancies):
        raise ValueError('occupancies2..512 required')
    data, record = maps.prepare('byte_native')
    shape = screen.geometry(data)
    beta = screen.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2).beta
    beta_log = log(beta.numerator)-log(beta.denominator)
    start, pins = monotonic(), sources()
    result = dict(schema='packet8-counterfactual-return-sensitivity-1',
        counterfactual_only=True, proposal_only=True, whole_code_certificate=False,
        realizable_larger_state_construction=False,
        actual_state_bits=16, geometry=shape, actual_map=record,
        changed_entries='only local[:,1:,0] multiplied by65535/(2^target_bits-1)',
        unchanged='all emission moments,birth laws,nonzero density coefficients,outer androuting',
        occupancies=list(occupancies), target_bits=args.bits, input_tilts=list(map(str, tilts)),
        source_sha256=pins, source_pins_verified_at_finish=False,
        trials=[], best={str(s): {} for s in args.bits}, choices={str(s): {} for s in args.bits})
    for tilt in tilts:
        weighted, emission, _ = screen.moments(data, np.exp(-float(tilt)))
        baseline = screen.operators(weighted, emission)
        for bits in args.bits:
            local = change_denominator(baseline, bits)
            regional, backend = screen.placement(local, max(occupancies), epochs=64, windows=8)
            witnesses = {}
            for q in occupancies:
                moment = screen.logarithmic.log_power_matrix(screen.uniform_mixture(regional, q), 32)
                value = -(log(comb(512, q))+q*beta_log+float(tilt)*13107+moment)/log(2)
                witnesses[str(q)] = value
                if str(q) not in result['best'][str(bits)] or value > result['best'][str(bits)][str(q)]:
                    result['best'][str(bits)][str(q)] = value
                    result['choices'][str(bits)][str(q)] = str(tilt)
            result['trials'].append(dict(target_bits=bits, tilt=str(tilt), backend=backend,
                                         counterfactual_margin_bits=witnesses))
            if pins != sources():
                raise ArithmeticError('sensitivity source changed during run')
            result['elapsed_seconds'] = monotonic()-start
            args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
            print(f'COUNTERFACTUAL returnbits={bits} tilt={tilt}: '
                  f'best={result["best"][str(bits)]}; elapsed={monotonic()-start:.1f}s', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bits', type=int, nargs='+', default=[20, 24, 28, 32])
    parser.add_argument('--q', type=int, nargs='+', default=[64, 119, 128, 256])
    parser.add_argument('--tilts', nargs='+', default=['.10', '.20', '.30', '.40', '.50', '.60', '.80', '1.2', '1.6'])
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
