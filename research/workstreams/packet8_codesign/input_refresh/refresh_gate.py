"""Floating gates for two actual input-refresh update distributions.

Emission remains y=X+A*s. The shared-refresh order updates s'=M(s+CX).
The independent-refresh order updates s'=a*s+b*CX with independent uniform
nonzero GF(2^16) scalars. Both admit a two-coordinate delta0/uniform closure.
This file does not modify or inherit a distance certificate for the old code.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
from itertools import combinations
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
import screen as prior


def exact_returns(data, maximum=2):
    """Positive integer histograms of wt((I+AC)X), excluding CX=0."""
    if maximum not in (1, 2):
        raise ValueError('only bounded one- and two-packet enumeration is supported')
    b, W, width = data['packet_bits'], data['windows'], data['width']
    labels = np.arange(1, 1 << b, dtype=np.uint64)
    syndromes, output = [], []
    for h in range(W):
        c = data['columns'][b*h:b*(h+1)]
        values = np.array([maps.apply(c, int(x)) for x in labels], dtype=np.uint32)
        syndromes.append(values)
        output.append((labels << np.uint64(b*h)) ^ data['images'][values])
    histograms = np.zeros((maximum+1, width+1), dtype=np.int64)
    for h in range(W):
        histograms[1] += np.bincount(np.bitwise_count(output[h][syndromes[h] != 0]), minlength=width+1)
    if maximum == 2:
        for h, k in combinations(range(W), 2):
            syndrome = syndromes[h][:, None] ^ syndromes[k][None, :]
            words = output[h][:, None] ^ output[k][None, :]
            histograms[2] += np.bincount(np.bitwise_count(words[syndrome != 0]), minlength=width+1)
    denominators = [comb(W, j)*((1 << b)-1)**j for j in range(maximum+1)]
    for j in range(1, maximum+1):
        if j in data['zero_feedback_forbidden'] and histograms[j].sum() != denominators[j]:
            raise ArithmeticError('exact-return census violates structural feedback rank')
    return dict(histograms=histograms, denominators=denominators,
        exact_integer_counts=True, maximum_occupancy=maximum)


def prepare_holder(data, block_bits=16):
    """Four output blocks; only zero-offset moments are needed for returns."""
    width, b, W = data['width'], data['packet_bits'], data['windows']
    if width % block_bits or block_bits > 16:
        raise ValueError('complete output blocks of at most16 bits required')
    columns = tuple((1 << i) ^ int(data['images'][c]) for i, c in enumerate(data['columns']))
    if maps.rank(columns) != width:
        raise ValueError('information-set return bound requires invertible I+AC')
    size, labels = 1 << block_bits, (1 << b)-1
    chars = np.arange(size, dtype=np.uint32)
    coefficients = np.array([[sum(comb(k, ell)*comb(W-k, j-ell)*labels**ell*(-1)**(j-ell)
        for ell in range(max(0, j-W+k), min(j, k)+1))/(comb(W, j)*labels**j)
        for j in range(W+1)] for k in range(W+1)], dtype=float)
    records = []
    for shift in range(0, width, block_bits):
        annihilated = np.zeros(size, dtype=np.uint8)
        for h in range(W):
            active = np.zeros(size, dtype=np.uint8)
            for bit in range(b):
                active |= np.bitwise_count(chars & ((columns[b*h+bit] >> shift) & (size-1))) & 1
            annihilated += active == 0
        records.append(np.bincount(annihilated.astype(np.int64)*(block_bits+1)+np.bitwise_count(chars),
            minlength=(W+1)*(block_bits+1)).reshape(W+1, block_bits+1))
    return dict(block_bits=block_bits, blocks=width//block_bits,
        coefficients=coefficients, character_histograms=records,
        output_rank=maps.rank(columns))


def holder_returns(prepared, z):
    power = float(z)**prepared['blocks']
    b = prepared['block_bits']
    weights = np.arange(b+1)
    fourier = (1+power)**(b-weights)*(1-power)**weights
    moments = np.array([(fourier @ histogram.T @ prepared['coefficients'])/(1 << b)
                        for histogram in prepared['character_histograms']])
    if not np.isfinite(moments).all() or moments.min() < -2e-12:
        raise ArithmeticError('unstable zero-offset Holder moments')
    with np.errstate(divide='ignore'):
        return np.exp(np.log(np.maximum(moments, 0)).mean(axis=0))


def operators(data, weighted, emission, z, *, mode, census=None, holder=None):
    """Entrywise positive bound; never subtract an upper return estimate."""
    W, S, b = data['windows'], 1 << data['bits'], data['packet_bits']
    L = S-1
    total = emission[1:].mean(axis=0)
    source_mass = weighted.sum(axis=0)
    local = np.zeros((W+1, 2, 2))
    local[:, 0, 0] = weighted[0]
    local[:, 0, 1] = source_mass-weighted[0]
    if mode == 'independent':
        lower, upper = np.zeros(W+1), total/L
        upper[0] = 0
        diagnostic = dict(return_bound='T/L; exact return is (T-weighted_CX_zero)/L',
            next_state='a*s+b*CX; a,b independent uniform nonzero GF(2^16) scalars')
    elif mode == 'shared':
        if census is None or holder is None:
            raise ValueError('fresh exact low-occupancy census and Holder data required')
        full_upper = np.minimum(1., ((1+float(z))**b/((1 << b)-1))**np.arange(W+1))
        holder_upper = holder_returns(holder, z)
        full_upper = np.minimum(full_upper, holder_upper)
        lower = np.zeros(W+1)
        upper = np.minimum(total, np.maximum(0., full_upper-weighted[0])/L)
        for j in range(1, census['maximum_occupancy']+1):
            exact = np.dot(census['histograms'][j], float(z)**np.arange(data['width']+1))
            exact /= census['denominators'][j]*L
            lower[j] = upper[j] = exact
        upper[0] = 0
        diagnostic = dict(return_bound='exact j1/j2; min(Holder, information-set,1) minus exact CX0 for higherj',
            next_state='M(s+CX); M fresh independent uniform GL(2,GF256)',
            full_return_moment_upper=full_upper.tolist(), holder_upper=holder_upper.tolist())
    else:
        raise ValueError('unknown realizable refresh order')
    local[:, 1, 0] = upper
    local[:, 1, 1] = total-lower
    if not np.isfinite(local).all() or local.min() < -2e-13:
        raise ArithmeticError('invalid positive comparison operator')
    local = np.maximum(local, 0)
    diagnostic.update(return_lower=lower.tolist(), return_upper=upper.tolist(),
                      uniform_emission=total.tolist())
    return local, diagnostic


def source_pins():
    result = prior.sources('byte_native')
    for path in (Path(__file__).resolve(), HERE/'test_refresh.py'):
        result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def run(args):
    tilts = tuple(Fraction(x) for x in args.tilts)
    if (args.output.exists() or min(tilts) <= 0 or len(set(tilts)) != len(tilts)
            or not 2 <= args.min_q <= args.max_q <= 512):
        raise ValueError('fresh output, distinct positive tilts and valid occupancy range required')
    started, pins = monotonic(), source_pins()
    data, record = maps.prepare('byte_native')
    census = exact_returns(data) if 'shared' in args.modes else None
    holder = prepare_holder(data) if 'shared' in args.modes else None
    shape = prior.geometry(data)
    envelope = prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
    counts = prior.rs_outer.expected_group_support_counts(16, 8, 8, 2)
    count_logs = np.array([log(x.numerator)-log(x.denominator) if x else -np.inf for x in counts])
    beta_log = log(envelope.beta.numerator)-log(envelope.beta.denominator)
    results = {mode: dict(best={}, tilt_choices={}, best_q1=np.zeros(33), trials=[]) for mode in args.modes}
    result = dict(schema='packet8-input-refresh-proposal-1', map_record=record, geometry=shape,
        proposal_only=True, whole_code_certificate=False, has_outward_endpoints=False,
        input_tilts=list(map(str, tilts)), occupancy_range=[args.min_q, args.max_q],
        source_sha256=pins, source_pins_verified_at_finish=False, envelope=envelope.metadata(), modes={})
    print(f'refresh preparation {monotonic()-started:.2f}s', flush=True)
    for tilt in tilts:
        z = np.exp(-float(tilt))
        weighted, emission, prior_diagnostic = prior.moments(data, z)
        for mode, entry in results.items():
            local, diagnostic = operators(data, weighted, emission, z,
                mode=mode, census=census, holder=holder)
            regional, backend = prior.placement(local, args.max_q, epochs=64, windows=8)
            entry['best_q1'] = np.minimum(entry['best_q1'], prior.q1_support_logs(
                regional, regions=32, tilt=tilt, cutoff=13107))
            entry['q1_margin_bits'] = -float(logsumexp(log(512)+count_logs+entry['best_q1']))/log(2)
            witnesses = {}
            for q in range(args.min_q, args.max_q+1):
                moment = prior.logarithmic.log_power_matrix(prior.uniform_mixture(regional, q), 32)
                margin = -(log(comb(512, q))+q*beta_log+float(tilt)*13107+moment)/log(2)
                if not isfinite(margin):
                    raise FloatingPointError('invalid full comparison moment')
                witnesses[str(q)] = dict(estimated_margin_bits=margin, log_moment=moment)
                if str(q) not in entry['best'] or margin > entry['best'][str(q)]:
                    entry['best'][str(q)], entry['tilt_choices'][str(q)] = margin, str(tilt)
            entry['trials'].append(dict(tilt=str(tilt), backend=backend,
                prior_diagnostics=prior_diagnostic, return_diagnostics=diagnostic, witnesses=witnesses))
            worst = min(entry['best'], key=entry['best'].get)
            entry['weakest_checked_q'], entry['weakest_checked_margin_bits'] = int(worst), entry['best'][worst]
            result['modes'][mode] = {key: value.tolist() if isinstance(value, np.ndarray) else value
                                     for key, value in entry.items()}
            print(f'{mode} tilt={tilt}: q1={entry["q1_margin_bits"]:.6f}, '
                  f'weakest q={worst} margin={entry["best"][worst]:.6f}', flush=True)
        if pins != source_pins():
            raise ArithmeticError('source changed during proposal')
        result['elapsed_seconds'] = monotonic()-started
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes', nargs='+', choices=('independent', 'shared'), default=['independent', 'shared'])
    parser.add_argument('--tilts', nargs='+', default=['.00384', '.00768', '.02', '.05', '.10', '.20', '.40', '.60', '.80', '1.0', '1.4', '1.8', '2.2'])
    parser.add_argument('--min-q', type=int, default=2)
    parser.add_argument('--max-q', type=int, default=256)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
