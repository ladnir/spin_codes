"""Floating screen for a real 256-to-512-bit outer with routed byte packets.

Eight parallel GF16 RS[16,8] rows form one GF(2^32) MDS[16,8] word.
Each 32-bit symbol receives an independent invertible random map with a
uniform nonzero image for every fixed nonzero input. All 64 byte packets
are shuffled within their group before fresh independent regional routes.
The inner is the unchanged t64/s16 byte-native construction.

These are finite-formula floating proposals, not outward certificates.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, isfinite, log
from pathlib import Path
import sys
from time import monotonic

import numpy as np
from scipy.special import logsumexp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import maps
import screen as local_screen


def geometry(*, symbol_bits=32, message_bits=65536):
    if symbol_bits not in (16, 32) or message_bits % (8*symbol_bits):
        raise ValueError('the compared complete RS16 outer geometries are required')
    group_k, group_n = 8*symbol_bits, 16*symbol_bits
    groups, regions = message_bits//group_k, group_n//8
    if groups % 8:
        raise ValueError('every region must contain whole eight-byte steps')
    return dict(K=message_bits, N=2*message_bits, outer_groups=groups,
        group_dimension=group_k, group_output_bits=group_n,
        outer_base_field_size=16, parallel_RS_rows=symbol_bits//4,
        outer_RS_n=16, outer_RS_k=8, symbol_bits=symbol_bits,
        regions=regions, packet_bits=8, slots_per_region=groups,
        physical_t=64, physical_packet_slots=8,
        physical_steps_per_region=groups//8, physical_steps_total=regions*groups//8,
        state_bits=16, cutoff=(2*message_bits)//10, zero_initial_state=True,
        final_flush=False, continuous_state_across_all_steps_and_regions=True,
        order='y=X+A*a; a_next=M*a+C*X')


def outer(shape):
    envelope = local_screen.rs_uniform_envelope.UniformInputEnvelope(
        16, 8, 8, shape['symbol_bits']//8)
    counts = local_screen.rs_outer.expected_group_support_counts(
        16, 8, 8, shape['symbol_bits']//8)
    if (len(counts) != shape['regions']+1
            or sum(counts) != (1 << shape['group_dimension'])-1
            or envelope.message_bits != shape['group_dimension']
            or envelope.output_bits != shape['group_output_bits']):
        raise ArithmeticError('outer dimensions or exact expected-count mass disagree')
    envelope.verify_shell_domination()
    return envelope, counts


def source_pins():
    result = local_screen.sources('byte_native')
    for path in (Path(__file__).resolve(), HERE/'test_wider.py'):
        result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def run(args):
    tilts = tuple(Fraction(value) for value in args.tilts)
    shape = geometry(symbol_bits=args.symbol_bits)
    if (args.output.exists() or not tilts or min(tilts) <= 0
            or len(set(tilts)) != len(tilts)
            or not 2 <= args.min_q <= args.max_q <= shape['outer_groups']):
        raise ValueError('fresh output, distinct positive tilts, and valid occupancy range required')
    started = monotonic()
    pins = source_pins()
    data, map_record = maps.prepare('byte_native')
    if (data['bits'], data['width'], data['windows'], data['packet_bits']) != (16, 64, 8, 8):
        raise ArithmeticError('unchanged literal byte-native local maps required')
    envelope, counts = outer(shape)
    beta_log = log(envelope.beta.numerator)-log(envelope.beta.denominator)
    count_logs = np.array([log(value.numerator)-log(value.denominator)
                           if value else -np.inf for value in counts])
    best_q1 = np.zeros(shape['regions']+1)
    result = dict(schema='packet8-wider-outer-proposal-1', geometry=shape,
        map_record=map_record, proposal_only=True, whole_code_certificate=False,
        has_outward_endpoints=False, q1_exact_expected_shells=True,
        q2_exact_joint_shells=False, occupancy_range=[args.min_q, args.max_q],
        envelope=envelope.metadata(), local_arithmetic='floating exact finite formulas',
        global_arithmetic='scaled/logarithmic ordered hypergeometric placement',
        source_sha256=pins, source_pins_verified_at_finish=False,
        input_tilts=list(map(str, tilts)), trials=[], best={}, tilt_choices={})
    print(f'symbol{shape["symbol_bits"]}: groups={shape["outer_groups"]}, '
          f'regions={shape["regions"]}, physical steps/region={shape["physical_steps_per_region"]}', flush=True)
    for tilt in tilts:
        weighted, emission, diagnostics = local_screen.moments(data, np.exp(-float(tilt)))
        local = local_screen.operators(weighted, emission)
        regional, backend = local_screen.placement(local, args.max_q,
            epochs=shape['physical_steps_per_region'], windows=8, force_log=args.force_log)
        zero, zero_backend = local_screen.placement(local[:, :1, :1], args.max_q,
            epochs=shape['physical_steps_per_region'], windows=8, force_log=args.force_log)
        trial = dict(tilt=str(tilt), backend=backend, zero_path_backend=zero_backend,
            local_diagnostics=diagnostics, witnesses={})
        for q in range(args.min_q, args.max_q+1):
            moment = local_screen.logarithmic.log_power_matrix(
                local_screen.uniform_mixture(regional, q), shape['regions'])
            zero_moment = shape['regions']*float(local_screen.uniform_mixture(zero, q)[0, 0])
            count = log(comb(shape['outer_groups'], q))+q*beta_log
            margin = -(count+float(tilt)*shape['cutoff']+moment)/log(2)
            if not isfinite(margin):
                raise FloatingPointError('invalid global moment')
            trial['witnesses'][str(q)] = dict(estimated_margin_bits=margin,
                log_moment=moment, zero_state_only_margin_bits=
                -(count+float(tilt)*shape['cutoff']+zero_moment)/log(2))
            if str(q) not in result['best'] or margin > result['best'][str(q)]:
                result['best'][str(q)], result['tilt_choices'][str(q)] = margin, str(tilt)
        best_q1 = np.minimum(best_q1, local_screen.q1_support_logs(regional,
            regions=shape['regions'], tilt=tilt, cutoff=shape['cutoff']))
        q1_terms = log(shape['outer_groups'])+count_logs+best_q1
        result['q1_margin_bits'] = -float(logsumexp(q1_terms))/log(2)
        result['q1_support_log_uppers_floating'] = best_q1.tolist()
        result['q1_dominant_supports'] = np.argsort(q1_terms)[::-1][:5].tolist()
        trial['q1_accumulated_margin_bits'] = result['q1_margin_bits']
        result['trials'].append(trial)
        if pins != source_pins():
            raise ArithmeticError('authenticated source changed during run')
        result['elapsed_seconds'] = monotonic()-started
        worst = min(result['best'], key=result['best'].get)
        result['weakest_checked_q'] = int(worst)
        result['weakest_checked_margin_bits'] = result['best'][worst]
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        print(f'tilt={tilt}: weakest q={worst}, {result["best"][worst]:.6f} bits; '
              f'q1={result["q1_margin_bits"]:.6f}; elapsed={result["elapsed_seconds"]:.2f}s', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--symbol-bits', type=int, choices=(16, 32), default=32)
    parser.add_argument('--min-q', type=int, default=2)
    parser.add_argument('--max-q', type=int, default=256)
    parser.add_argument('--tilts', nargs='+', default=[
        '.00256', '.00512', '.01', '.0256', '.0512', '.10', '.20', '.30', '.40', '.50', '.60', '.80'])
    parser.add_argument('--force-log', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
