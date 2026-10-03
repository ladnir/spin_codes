"""Bounded-memory fine-count grouping, with a reusable powered macro family.

For group size eight, merge two four-step tuple lists in vectorized chunks.
The two half products are multiplied before entrywise power. Tuple placement
weights stay outside the power. The output stores the powered operator and
all-q floating proposals so the tuple census is never repeated per q.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb, isfinite, log
from pathlib import Path
from time import monotonic

import numpy as np

import fine_grouped_gate as fine

HERE = Path(__file__).resolve().parent
prior, maps, marked = fine.prior, fine.maps, fine.marked


def _half_tuples(potential, steps):
    tuples = fine.tuple_products(potential, steps, cap=min(2, len(potential)-1))
    totals = tuples['totals']
    denominator = np.array([comb(tuples['windows'], int(j)) for j in totals], dtype=float)
    raw = tuples['weights']*denominator
    multiplicities = np.rint(raw)
    if np.any(np.abs(raw-multiplicities) > 1e-6):
        raise ArithmeticError('half-tuple multiplicities are not exact small integers')
    order = np.argsort(totals, kind='stable')
    matrices, totals = tuples['products'][order], totals[order]
    multiplicities = multiplicities[order]
    starts = np.flatnonzero(np.r_[True, totals[1:] != totals[:-1]])
    if not np.array_equal(totals[starts], np.arange(tuples['windows']+1)):
        raise ArithmeticError('missing half occupancy class')
    masses = np.add.reduceat(multiplicities, starts)
    expected = np.array([comb(tuples['windows'], j) for j in range(tuples['windows']+1)])
    if not np.array_equal(masses, expected):
        raise ArithmeticError('half-tuple multiplicities fail Vandermonde normalization')
    return matrices, totals, multiplicities, starts


def long_operators(potential, group_steps, alpha, *, memory_mib=64, progress=None):
    """Fine tuple average of (T_j1...T_jg)^alpha, plus alpha-one regression."""
    potential = np.asarray(potential, dtype=float)
    alpha = float(alpha)
    if (type(group_steps) is not int or group_steps not in (2, 4, 8)
            or not isfinite(alpha) or not 0 < alpha <= 1 or memory_mib < 16):
        raise ValueError('group size2/4/8, alpha in(0,1], and at least16MiB required')
    half = group_steps//2
    matrices, totals, mult, starts = _half_tuples(potential, half)
    W, size = len(potential)-1, potential.shape[1]
    half_windows, full_windows = half*W, group_steps*W
    ends = np.r_[starts[1:], len(totals)]
    local = np.zeros((full_windows+1, size, size))
    regression = np.zeros_like(local)
    # One chunk product dominates memory. Powers and weighted products are
    # in-place. Leave room for the sorted half list, weights, and reductions.
    fixed_bytes = matrices.nbytes+8*len(totals)*size*size+local.nbytes*2
    per_left = len(totals)*(size*size*8+8*3)
    chunk_size = min(8, (int(memory_mib)*2**20-fixed_bytes)//per_left)
    if chunk_size < 1:
        raise ValueError('memory budget cannot hold one vectorized product chunk')
    expected_chunks = sum((int(end-start)+chunk_size-1)//chunk_size for start, end in zip(starts, ends))
    chunks = 0
    for left_total, (start, end) in enumerate(zip(starts, ends)):
        denominator = np.array([comb(full_windows, left_total+int(j)) for j in totals], dtype=float)
        right_weight = mult/denominator
        for first in range(int(start), int(end), chunk_size):
            last = min(first+chunk_size, int(end))
            products = matrices[first:last, None]@matrices[None]
            weights = mult[first:last, None]*right_weight[None]
            # Keep a simultaneous alpha-one control without a second census.
            unpowered_sum = np.einsum('lr,lrab->rab', weights, products, optimize=False)
            regression[left_total:left_total+half_windows+1] += np.add.reduceat(unpowered_sum, starts, axis=0)
            del unpowered_sum
            np.power(products, alpha, out=products)
            products *= weights[:, :, None, None]
            reduced_left = products.sum(axis=0)
            local[left_total:left_total+half_windows+1] += np.add.reduceat(reduced_left, starts, axis=0)
            chunks += 1
            if progress is not None and (chunks % 64 == 0 or chunks == expected_chunks):
                progress(chunks, expected_chunks)
    coarse, backend = fine.coarse.grouped_operators(potential, group_steps)
    positive = coarse > 0
    if (np.any(regression[~positive] != 0)
            or not np.allclose(regression[positive], coarse[positive], rtol=3e-12, atol=1e-300)):
        raise ArithmeticError('simultaneous alpha-one macro regression failed')
    relative_error = float(np.max(np.abs(regression[positive]/coarse[positive]-1)))
    if not np.isfinite(local).all() or np.any(local < 0):
        raise FloatingPointError('invalid powered macro family')
    return local, dict(tuple_count=len(totals)**2, half_tuple_count=len(totals),
        chunk_size=int(chunk_size), chunks=chunks, memory_budget_mib=memory_mib,
        alpha1_max_relative_error=relative_error, alpha1_macro_backend=backend)


def source_pins():
    result = fine.source_pins()
    for path in (Path(__file__).resolve(), HERE/'test_long_grouped_gate.py'):
        result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def run(args):
    theta, alpha = Fraction(args.tilt), Fraction(args.alpha)
    if (args.output.exists() or theta <= 0 or not 0 < alpha <= 1
            or args.group not in (4, 8) or not 2 <= args.min_q <= args.max_q <= 512):
        raise ValueError('fresh output, positive tilt, alpha in(0,1] and feasible q range required')
    start, pins = monotonic(), source_pins()
    data, record = maps.prepare('byte_native')
    shape = prior.geometry(data)
    weighted, emission, diagnostic = prior.moments(data, np.exp(-float(theta)))
    potential = marked.potential_operators(prior.operators(weighted, emission), 0.)
    print(f'Preparing fine g={args.group}, theta={theta}, alpha={alpha}', flush=True)
    def report(done, total):
        if source_pins() != pins:
            raise ArithmeticError('source changed during macro census')
        print(f'macro chunks {done}/{total}; elapsed={monotonic()-start:.2f}s', flush=True)
    local, macro_metadata = long_operators(potential, args.group, alpha,
        memory_mib=args.memory_mib, progress=report)
    envelope = prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
    beta_log = log(envelope.beta.numerator)-log(envelope.beta.denominator)
    result = dict(schema='packet8-long-fine-grouped-proposal-1', geometry=shape,
        map_record=record, construction_changed=False, proposal_only=True,
        whole_code_certificate=False, has_outward_endpoints=False,
        tilt=str(theta), alpha=str(alpha), group_steps=args.group,
        power_order='potential labels, full ordered tuple product, entrywise alpha, tuple average',
        message_factor='C(512,q)*beta^(q*alpha)', subset_union_is_not_powered=True,
        source_sha256=pins, source_pins_verified_at_finish=False,
        envelope=envelope.metadata(), local_diagnostics=diagnostic,
        macro_metadata=macro_metadata, powered_macro_operators=local.tolist(),
        occupancy_range=[args.min_q, args.max_q], per_q={}, positive_intervals=[],
        at_least_50_bit_intervals=[], elapsed_seconds=monotonic()-start)
    def save():
        if source_pins() != pins:
            raise ArithmeticError('source changed during long grouped proposal')
        result['elapsed_seconds'] = monotonic()-start
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    save()
    regional, backend = prior.placement(local, args.max_q,
        epochs=64//args.group, windows=8*args.group)
    result['regional_backend'] = backend
    for q in range(args.min_q, args.max_q+1):
        moment = prior.logarithmic.log_power_matrix(regional[q], 32)
        margin = -(log(comb(512, q))+float(alpha)*(q*beta_log+float(theta)*shape['cutoff'])+moment)/log(2)
        if not isfinite(margin):
            raise FloatingPointError('nonfinite all-q margin')
        result['per_q'][str(q)] = dict(margin_bits=margin, log_moment=moment)
    for threshold, key in ((0., 'positive_intervals'), (50., 'at_least_50_bit_intervals')):
        good = [int(q) for q, value in result['per_q'].items() if value['margin_bits'] >= threshold]
        for q in good:
            if result[key] and q == result[key][-1][1]+1:
                result[key][-1][1] = q
            else:
                result[key].append([q, q])
    result['source_pins_verified_at_finish'] = True
    save()
    selected = result['per_q'].get('119', {})
    print(f'g={args.group}: q119={selected}; positive={result["positive_intervals"]}; '
        f'50bit={result["at_least_50_bit_intervals"]}; elapsed={monotonic()-start:.2f}s', flush=True)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tilt', default='.475')
    parser.add_argument('--alpha', default='.35')
    parser.add_argument('--group', type=int, default=8)
    parser.add_argument('--memory-mib', type=int, default=64)
    parser.add_argument('--min-q', type=int, default=2)
    parser.add_argument('--max-q', type=int, default=512)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
