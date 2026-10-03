"""Pair nonempty potential steps before fractional path expansion.

All local products and dynamic-programming coefficients are logarithmic.
A guarded row/column-scaled product uses ordinary matmul only when every
positive scalar product is provably above the binary64 normal threshold.
These are floating proposals, not outward probability certificates.
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
from scipy.special import logsumexp

import potential_birth_gate as birth

HERE = Path(__file__).resolve().parent
prior, maps, marked = birth.prior, birth.maps, birth.marked
POINTS = {16: ('.06', '.4'), 64: ('.25', '.4')}


def log_product(left, right, *, statistics=None, force_log=False):
    """Ordered batched log product, retaining all positive matrix entries."""
    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    if (left.ndim < 2 or right.ndim < 2 or left.shape[-1] != left.shape[-2]
            or right.shape[-2:] != left.shape[-2:] or np.any(np.isnan(left))
            or np.any(np.isnan(right)) or np.any(np.isposinf(left))
            or np.any(np.isposinf(right))):
        raise ValueError('matching square log matrices required')
    # Row and column scales are restored as logarithms, never exponentiated.
    row = np.max(left, axis=-1, keepdims=True)
    col = np.max(right, axis=-2, keepdims=True)
    row = np.where(np.isfinite(row), row, 0.)
    col = np.where(np.isfinite(col), col, 0.)
    ls, rs = left-row, right-col
    lp, rp = np.isfinite(ls), np.isfinite(rs)
    floor = log(np.finfo(float).tiny)+4.
    lower = float(ls[lp].min(initial=0.))+float(rs[rp].min(initial=0.))
    if force_log or lower < floor:
        if statistics is not None:
            statistics['full_log_products'] += 1
        return prior.logarithmic._multiply(left, right)
    # Every nonzero scalar multiplication is normal; summands cannot cancel.
    with np.errstate(over='raise', invalid='raise', divide='raise', under='raise'):
        linear = np.exp(ls) @ np.exp(rs)
        result = np.full(linear.shape, -np.inf)
        positive = linear > 0
        result[positive] = np.log(linear[positive])
    result += row+col
    if statistics is not None:
        statistics['guarded_linear_products'] += 1
        statistics['smallest_guarded_scalar_log_lower'] = min(
            statistics['smallest_guarded_scalar_log_lower'], lower)
    return result


def event_blocks(potential, epochs, alphas, *, statistics=None, force_log=False):
    """Unnormalized pair A[L,k] and terminal B[L,k], in logarithms."""
    potential = np.asarray(potential, dtype=float)
    alphas = tuple(float(a) for a in alphas)
    if (type(epochs) is not int or epochs < 1 or potential.ndim != 3
            or len(potential) < 2 or potential.shape[1] != potential.shape[2]
            or not np.isfinite(potential).all() or np.any(potential < 0)
            or not alphas or any(not isfinite(a) or not 0 < a <= 1 for a in alphas)):
        raise ValueError('positive family, positive epochs and powers in (0,1] required')
    W, size = len(potential)-1, potential.shape[1]
    logs = prior.logarithmic.array_logs(potential)
    identity = np.full((size, size), -np.inf)
    np.fill_diagonal(identity, 0.)
    zeros = np.empty((epochs+1, size, size))
    zeros[0] = identity
    multiply = lambda a, b: log_product(a, b, statistics=statistics, force_log=force_log)
    for g in range(1, epochs+1):
        zeros[g] = multiply(zeros[g-1], logs[0])
    # left[g,j-1] = T0^g Tj; this shared list has only E*W matrices.
    left = multiply(zeros[:-1, None], logs[None, 1:])
    powers = np.asarray(alphas)[:, None, None, None]
    pair = np.full((len(alphas), epochs+1, 2*W+1, size, size), -np.inf)
    tail = np.full((len(alphas), epochs+1, W+1, size, size), -np.inf)
    tail[:, :, 0] = np.asarray(alphas)[:, None, None, None]*zeros[None]
    choices = np.array([log(comb(W, j)) for j in range(1, W+1)])
    totals = (np.arange(1, W+1)[:, None]+np.arange(1, W+1)[None]).reshape(-1)
    coefficient = (choices[:, None]+choices[None]).reshape(-1)
    masks = [totals == k for k in range(2, 2*W+1)]
    for length in range(1, epochs+1):
        for g0 in range(length):
            # Exactly one event followed by the remaining zero steps.
            block = multiply(left[g0], zeros[length-1-g0])
            value = powers*block[None]+choices[None, :, None, None]
            tail[:, length, 1:] = np.logaddexp(tail[:, length, 1:], value)
        if length >= 2:
            for g0 in range(length-1):
                block = multiply(left[g0, :, None], left[length-2-g0, None])
                block = block.reshape(W*W, size, size)
                value = powers*block[None]+coefficient[None, :, None, None]
                for k, mask in enumerate(masks, start=2):
                    pair[:, length, k] = np.logaddexp(pair[:, length, k],
                        logsumexp(value[:, mask], axis=1))
    return pair, tail


def event_regional(potential, degree, *, epochs=64, alpha=.4,
                   check_alpha1=True, force_log_products=False, progress=None):
    """All regional log operators through degree; no reset at either boundary.

    Return (regional_logs, metadata). The alpha-one control is evaluated in
    the same block census/DP and compared with ordinary exact-q placement.
    """
    potential = np.asarray(potential, dtype=float)
    if (type(degree) is not int or degree < 0 or potential.ndim != 3
            or type(epochs) is not int or epochs < 1 or degree > epochs*(len(potential)-1)):
        raise ValueError('complete family and feasible regional geometry required')
    started = monotonic()
    statistics = dict(full_log_products=0, guarded_linear_products=0,
        smallest_guarded_scalar_log_lower=0.)
    alphas = (float(alpha), 1.) if check_alpha1 and float(alpha) != 1 else (float(alpha),)
    pair, tail = event_blocks(potential, epochs, alphas,
        statistics=statistics, force_log=force_log_products)
    if progress is not None:
        progress('local blocks', monotonic()-started)
    W, size = len(potential)-1, potential.shape[1]
    current = np.full((len(alphas), epochs+1, degree+1, size, size), -np.inf)
    for a in range(len(alphas)):
        np.fill_diagonal(current[a, 0, 0], 0.)
    multiply = lambda a, b: log_product(a, b, statistics=statistics, force_log=force_log_products)
    # D[e,q] enumerates complete pairs only. Every nonempty prefix ends at
    # its second event; no zero-only prefixes are inserted in this table.
    for e in range(2, epochs+1):
        limit = min(degree, W*e)
        for length in range(2, e+1):
            before = e-length
            for k in range(2, min(2*W, limit)+1):
                last = min(limit, k+W*before)
                if last < k:
                    continue
                term = multiply(current[:, before, :last-k+1], pair[:, length, k, None])
                current[:, e, k:last+1] = np.logaddexp(current[:, e, k:last+1], term)
        if progress is not None and (e % 16 == 0 or e == epochs):
            progress(f'pair DP {e}/{epochs}', monotonic()-started)
    regional = np.full((len(alphas), degree+1, size, size), -np.inf)
    for length in range(epochs+1):
        before = epochs-length
        for k in range(min(W, degree)+1):
            if length == 0 and k != 0:
                continue
            last = min(degree, k+W*before)
            if last < k:
                continue
            term = multiply(current[:, before, :last-k+1], tail[:, length, k, None])
            regional[:, k:last+1] = np.logaddexp(regional[:, k:last+1], term)
    regional -= np.array([log(comb(W*epochs, q)) for q in range(degree+1)])[None, :, None, None]
    if np.any(np.isnan(regional)) or np.any(np.isposinf(regional)):
        raise FloatingPointError('invalid event-aligned regional logarithms')
    check = None
    if check_alpha1:
        reference, backend = prior.placement(potential, degree, epochs=epochs, windows=W, force_log=True)
        control = regional[alphas.index(1.)]
        if not np.array_equal(np.isfinite(control), np.isfinite(reference)):
            raise ArithmeticError('alpha1 regional structural zeros changed')
        finite = np.isfinite(reference)
        error = float(np.max(np.abs(control[finite]-reference[finite]), initial=0.))
        if error > 2e-8:
            raise ArithmeticError('alpha1 regional placement regression failed')
        check = dict(max_regional_log_error=error, backend=backend)
    metadata = dict(epochs=epochs, windows=W, degree=degree, alpha=float(alpha),
        pair_tuple_count=W*W*epochs*(epochs-1)//2,
        single_event_tail_tuple_count=W*epochs*(epochs+1)//2,
        local_products_before_fractional_power='logarithmic; structural zeros are -infinity',
        dp_arithmetic='log coefficients with guarded row/column-scaled log products',
        alpha1_regression=check, product_statistics=statistics,
        elapsed_seconds=monotonic()-started)
    return regional[0], metadata


def source_pins():
    result = birth.source_pins()
    for path in (Path(__file__).resolve(), HERE/'test_event_aligned.py'):
        result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def run(args):
    if (args.output.exists() or not args.qs or len(set(args.qs)) != len(args.qs)
            or any(q not in POINTS for q in args.qs)):
        raise ValueError('fresh output and distinct supported occupancies required')
    started, pins = monotonic(), source_pins()
    data, record = maps.prepare('byte_native')
    shape = prior.geometry(data)
    envelope = prior.rs_uniform_envelope.UniformInputEnvelope(16, 8, 8, 2)
    beta_log = log(envelope.beta.numerator)-log(envelope.beta.denominator)
    result = dict(schema='packet8-event-pair-proposal-1', map_record=record,
        geometry=shape, construction_changed=False, proposal_only=True,
        whole_code_certificate=False, has_outward_endpoints=False,
        source_sha256=pins, source_pins_verified_at_finish=False,
        event='nonempty potential occupancy, not nonzero byte labels',
        pairing='within each region; full boundary state matrix retained',
        regional_normalization='one C(512,q), outside the fractional power',
        message_factor='C(512,q)*(beta^q*z^(-cutoff))^alpha',
        envelope=envelope.metadata(), trials=[])
    for q in args.qs:
        theta, alpha = map(Fraction, POINTS[q])
        weighted, emission, diagnostic = prior.moments(data, np.exp(-float(theta)))
        old = marked.potential_operators(prior.operators(weighted, emission), 0.)
        potential, change, algebra = birth.rebase_potential(old)
        def progress(stage, elapsed):
            if source_pins() != pins:
                raise ArithmeticError('source changed during event-pair proposal')
            print(f'q={q}: {stage}; {elapsed:.2f}s', flush=True)
        regional, metadata = event_regional(potential, q, alpha=float(alpha), progress=progress)
        moment = prior.logarithmic.log_power_matrix(regional[q], 32)
        margin = -(log(comb(512, q))+float(alpha)*(q*beta_log+float(theta)*shape['cutoff'])+moment)/log(2)
        if not isfinite(margin):
            raise FloatingPointError('nonfinite event-aligned margin')
        # Alpha-one global agreement follows from full regional agreement;
        # also record the ordinary unpowered global moment for inspection.
        control, _ = prior.placement(potential, q, epochs=64, windows=8, force_log=True)
        control_moment = prior.logarithmic.log_power_matrix(control[q], 32)
        result['trials'].append(dict(q=q, tilt=str(theta), alpha=str(alpha),
            margin_bits=margin, fractional_log_moment=moment,
            alpha1_ordinary_log_moment=control_moment, algebra_checks=algebra,
            event_metadata=metadata, local_diagnostics=diagnostic,
            potential_operators=potential.tolist(),
            regional_log_operator=[[float(x) if np.isfinite(x) else None for x in row]
                for row in regional[q]]))
        if source_pins() != pins:
            raise ArithmeticError('source changed during event-pair proposal')
        result['elapsed_seconds'] = monotonic()-started
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        print(f'q={q} theta={theta} alpha={alpha}: {margin:.9f} bits', flush=True)
    result['source_pins_verified_at_finish'] = True
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qs', nargs='+', type=int, default=[16, 64])
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args())
